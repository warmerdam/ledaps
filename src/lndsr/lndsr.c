#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include "lndsr.h"
#include "keyvalue.h"
#include "const.h"
#include "param.h"
#include "input.h"
#include "prwv_input.h"
#include "lut.h"
#include "output.h"
#include "sr.h"
#include "ar.h"
#include "space.h"
#include "bool.h"
#include "error.h"
#include "clouds.h"

#include "read_grib_tools.h"
#include "sixs_runs.h"

#define AERO_NB_BANDS 3
#define AERO_STATS_NB_BANDS 3
#define SP_INDEX	0
#define WV_INDEX	1
#define ATEMP_INDEX	2
#define OZ_INDEX	0
#define DEBUG_FLAG	0

/* DEM Definition: U_char format, 1 count = 100 meters */
/* 0 = 0 meters */

#define DEMFILE "CMGDEM.hdf"
#define DEM_NBLAT 3600
#define DEM_DLAT 0.05
#define DEM_LATMIN (-90.0)
#define DEM_LATMAX 90.0
#define DEM_NBLON 7200
#define DEM_DLON 0.05
#define DEM_LONMIN (-180.0)
#define DEM_LONMAX 180.0
#define P_DFTVALUE 1013.0

/* Type definitions */

typedef enum {FAILURE = 0, SUCCESS = 1} Status_t;


atmos_t atmos_coef;
#ifdef DEBUG_AR
FILE *fd_ar_diags;
int diags_il_ar;
#endif
#ifdef DEBUG_CLD
FILE *fd_cld_diags;
#endif
/* reading the DEM in hdf */
int32 sds_file_id,sds_id,status;
char sds_name[256];
int32 sds_index;
int32 dim_sizes[2],start[2],stride[2],edges[2];
int32 data_type,n_attrs,rank;
/* Prototypes */
#ifndef  HPUX
#define chand chand_
#define csalbr csalbr_
#endif
void chand(float *phi,float *muv,float *mus,float *tau_ray,float *actual_rho_ray);
void csalbr(float *tau_ray,float *actual_S_r);
int update_atmos_coefs(atmos_t *atmos_coef,Ar_gridcell_t *ar_gridcell, sixs_tables_t *sixs_tables,int ***line_ar,Lut_t *lut,int nband, int bkgd_aerosol);
int update_gridcell_atmos_coefs(int irow,int icol,atmos_t *atmos_coef,Ar_gridcell_t *ar_gridcell, sixs_tables_t *sixs_tables,int **line_ar,Lut_t *lut,int nband, int bkgd_aerosol);
void set_sixs_path_from(const char *path);

#ifdef SAVE_6S_RESULTS
#define SIXS_RESULTS_FILENAME "SIXS_RUN_RESULTS.TXT"
int read_6S_results_from_file(char *filename,sixs_tables_t *sixs_tables);
int write_6S_results_to_file(char *filename,sixs_tables_t *sixs_tables);
#endif
void sun_angles (short jday,float gmt,float flat,float flon,float *ts,float *fs);
/* Functions */

int main (int argc, const char **argv) {
  Param_t *param;
  Input_t *input,*input_b6;
  InputPrwv_t *prwv_input;
  InputOzon_t *ozon_input;
  InputMask_t *input_mask;
  Lut_t *lut;
  Output_t *output;
  int i,j,k, il, is,ib,i_aot,j_aot,ifree;
  int il_start, il_end, il_ar, il_region, is_ar;
  int *line_out[NBAND_SR_MAX];
  int *line_out_buf;
  int ***line_in;
  int **line_in_band_buf;
  int *line_in_buf;
  int ***line_ar;
  int **line_ar_band_buf;
  int *line_ar_buf;
  int ***line_ar_stats;
  int **line_ar_stats_band_buf;
  int *line_ar_stats_buf;
  int** b6_line;
  int* b6_line_buf;
  float *atemp_line;
  int8** qa_line;
  int8* qa_line_buf;
  char** mask_line;
  char* mask_line_buf;
  char **ddv_line;
  char *ddv_line_buf;
  char **rot_cld[3],**ptr_rot_cld[3],**ptr_tmp_cld;
  char **rot_cld_block_buf;
  char *rot_cld_buf;

  Sr_stats_t sr_stats;
  Ar_stats_t ar_stats;
  Ar_gridcell_t ar_gridcell;
  int *prwv_in[NBAND_PRWV_MAX];
  int *prwv_in_buf;
  int *ozon_in;
  
  int nbpts,ipt;
  float faot1;
  int inter_aot[3];

   float lamda[7]={486.,570.,660.,835.,1669.,0.,2207.};

  float scene_gmt;

  Space_t *space;
  Space_def_t space_def;
  char *grid_name = "Grid";
  char *dem_name;
  Img_coord_float_t img;
  Img_coord_int_t loc;
  Geo_coord_t geo;
  Geo_bounds_t bounds;
  int nsds,isds;
  char *sds_names[NBAND_SR_MAX];
  int sds_types[NBAND_SR_MAX];

  t_ncep_ancillary anc_O3,anc_WV,anc_SP,anc_ATEMP;
  double sum_spres_anc,sum_spres_dem;
  int nb_spres_anc,nb_spres_dem;
  float ratio_spres;
  float tmpflt_arr[4],tmpflt;
  double coef;
  int tmpint;
  int nvalues,osize;
  int debug_flag;

	sixs_tables_t sixs_tables;
	int sixs_run;
	char sixs_results_filename[256];

	float center_lat,center_lon;
	float aot550;
	
	char *charptr,tmpfilename[128];
	FILE *fdtmp,*fdtmp2;
	
	short *dem_array;
	int dem_available;

	float actual_rho_ray,actual_T_ray_up,actual_T_ray_down,actual_S_r;
	float mus,muv,phi,tau_ray;
	float rho_ray_P0,T_ray_down_P0,T_ray_up_P0,S_r_P0;
	float tau_ray_sealevel[7]={0.16511,0.08614,0.04716,0.01835,0.00113,0.00037}; /* index=5 => band 7 */

	cld_diags_t cld_diags;
	
float flat,flon,fts,ffs;
double delta_y,delta_x;
float adjust_north;

float sum_value,sumsq_value;
int no_ozone_file;
short jday;

/* Vermote additional variable declaration for the cloud mask May 29 2007 */
int anom,it6;
float t6,t6s_seuil;

float calcuoz(short jday,float flat);
float get_dem_spres(short *dem,float lat,float lon);
void swapbytes(void *val,int nbbytes);



 debug_flag= DEBUG_FLAG;

  no_ozone_file=0;
  
  set_sixs_path_from(argv[0]);

  param = GetParam(argc, argv);
  if (param == (Param_t *)NULL) ERROR("getting runtime parameters", "main");

  /* Open input file */

  input = OpenInput(param->input_file_name);
  if (input == (Input_t *)NULL) ERROR("bad input file", "main");

  if (param->thermal_band) {
    input_b6 = OpenInput(param->temp_file_name);
    if (input_b6 == (Input_t *)NULL) ERROR("bad band 6 input file", "main");
  }

  if ( param->num_prwv_files > 0  && param->num_ncep_files > 0  ) {
     ERROR("both PRWV and PRWV_FIL files specified", "main");
  }
  /* Open prwv input file */

  if ( param->num_prwv_files > 0  ) {
    prwv_input = OpenInputPrwv(param->prwv_file_name);
    if (prwv_input==(InputPrwv_t *)NULL) ERROR("bad input prwv file","main");

    osize= 3 * 
      (prwv_input->size.ntime*prwv_input->size.nlat*prwv_input->size.nlon);

    prwv_in_buf = (int *)calloc((size_t)(osize),sizeof(int));
    if (prwv_in_buf == (int *)NULL) 
      ERROR("allocating input prwv buffer", "main");
    prwv_in[0] = prwv_in_buf;
    for (ib = 1; ib < prwv_input->nband; ib++)
      prwv_in[ib] = prwv_in[ib - 1] + 
        (prwv_input->size.ntime*prwv_input->size.nlat*prwv_input->size.nlon);

    for (ib = 0; ib < prwv_input->nband; ib++)
      {
      nvalues= 
       (prwv_input->size.ntime*prwv_input->size.nlat*prwv_input->size.nlon);

      if (!GetInputPrwv(prwv_input, ib, prwv_in[ib]))
        ERROR("reading input prwv data", "main");
      }
    /**** ozone ***/
   if ( param->num_ozon_files<1 )
     no_ozone_file=1;
   else {
    ozon_input = OpenInputOzon(param->ozon_file_name);
    if (ozon_input==(InputOzon_t *)NULL) ERROR("bad input ozon file", "main");

    osize= 
      (ozon_input->size.ntime*ozon_input->size.nlat*ozon_input->size.nlon);

    ozon_in = (int *)calloc((size_t)(osize),sizeof(int));
    if (ozon_in == (int *)NULL) 
      ERROR("allocating input ozone buffer", "main");

      nvalues= 
       (ozon_input->size.ntime*ozon_input->size.nlat*ozon_input->size.nlon);

      if (!GetInputOzon(ozon_input, 0, ozon_in))
        ERROR("reading input ozone data", "main");
   }
  }

  /* Open input cloudmask */

  if ( param->cloud_flag ) {
    input_mask = OpenInputMask(param->cloud_mask_file);
    if (input == (Input_t *)NULL) ERROR("bad input file", "main");
    }

  /* Get Lookup table */

  lut = GetLut(param->lut_file_name, input->nband, &input->meta, &input->size);
  if (lut == (Lut_t *)NULL) ERROR("bad lut file", "main");
  
  /* Get space definition */
  if (!GetSpaceDefHDF(&space_def, param->input_file_name, grid_name))
    ERROR("getting space metadata from HDF file", "main");
  /* Setup Space */
  space = SetupSpace(&space_def);
  if (space == (Space_t *)NULL) 
    ERROR("setting up space", "main");

  /* compute bounds */
  if (!computeBounds( &bounds, space, input->size.s, input->size.l))
    ERROR("computing bounds", "main");

  /* Create and open output file */

  if (!CreateOutput(param->output_file_name))
    ERROR("creating output file", "main");

  output = OpenOutput(param->output_file_name, input->nband, input->meta.iband, 
                      &input->size);
  if (output == (Output_t *)NULL) ERROR("opening output file", "main");

  /* Open diagnostics files if needed */
#ifdef DEBUG_AR
   strcpy(tmpfilename,param->output_file_name);
	strcat(tmpfilename,".DEBUG_AR");
   fd_ar_diags=fopen(tmpfilename,"w");
   if (fd_ar_diags != NULL) {
		fprintf(fd_ar_diags,"cell_row cell_col total_nb_samples avg_b1 std_b1 avg_b2 std_b2 avg_b3 std_b3 avg_b7 std_b7 szen vzen relaz wv ozone spres fraction_water fraction_clouds fraction_cldshadow fraction_snow spres_ratio tau_ray corrected_T_ray corrected_Sr measured_rho_b1 simulated_b1_01 simulated_b1_02 simulated_b1_03 simulated_b1_04 simulated_b1_05 simulated_b1_06 simulated_b1_07 simulated_b1_08 simulated_b1_09 simulated_b1_10 simulated_b1_11 simulated_b1_12 simulated_b1_13 simulated_b1_14 simulated_b1_15 aot_index coef aot_value new_aot ratio_neg_red\n"); 
	}
#endif
#ifdef DEBUG_CLD
   strcpy(tmpfilename,param->output_file_name);
	strcat(tmpfilename,".DEBUG_CLD");
   fd_cld_diags=fopen(tmpfilename,"w");
   if (fd_cld_diags != NULL) {
		fprintf(fd_cld_diags,"cell_row cell_col nb_samples airtemp_2m avg_t6_clear std_t6_clear avg_b7_clear std_b7_clear\n"); 
	}
#endif
  /* Allocate memory for input lines */

  line_in = (int ***)calloc((size_t)(lut->ar_region_size.l), 
                            sizeof(int **));
  if (line_in == (int ***)NULL) 
    ERROR("allocating input line buffer (a)", "main");

  line_in_band_buf = (int **)calloc((size_t)(lut->ar_region_size.l * 
                                             input->nband), sizeof(int *));
  if (line_in_band_buf == (int **)NULL) 
    ERROR("allocating input line buffer (b)", "main");

  line_in_buf = (int *)calloc((size_t)(input->size.s * lut->ar_region_size.l *
                                       input->nband), sizeof(int));
  if (line_in_buf == (int *)NULL) 
    ERROR("allocating input line buffer (c)", "main");

  for (il = 0; il < lut->ar_region_size.l; il++) {
    line_in[il] = line_in_band_buf;
    line_in_band_buf += input->nband;
    for (ib = 0; ib < input->nband; ib++) {
      line_in[il][ib] = line_in_buf;
      line_in_buf += input->size.s;
    }
  }

  /* Allocate memory for qa line */
    qa_line = (int8**)calloc((size_t)(lut->ar_region_size.l),sizeof(int8 *));
    if (qa_line == (int8**)NULL)ERROR("allocating qa line", "main");
    qa_line_buf = (int8*)calloc((size_t)(input->size.s * lut->ar_region_size.l),sizeof(int8));
    if (qa_line_buf == (char*)NULL)ERROR("allocating qa line buffer", "main");
	for (il = 0; il < lut->ar_region_size.l; il++) {
	 	qa_line[il]=qa_line_buf;
		qa_line_buf += input->size.s;
	 }

  /* Allocate memory for one band 6 line */
  if (param->thermal_band) {
    b6_line = (int**)calloc((size_t)(lut->ar_region_size.l),sizeof(int *));
    if (b6_line == (int**)NULL)ERROR("allocating b6 line", "main");
    b6_line_buf = (int*)calloc((size_t)(input_b6->size.s * lut->ar_region_size.l),sizeof(int));
    if (b6_line_buf == (int*)NULL)ERROR("allocating b6 line buffer", "main");
	 for (il = 0; il < lut->ar_region_size.l; il++) {
	 	b6_line[il]=b6_line_buf;
		b6_line_buf += input_b6->size.s;
	 }
  }
  /* Allocate memory for one air temperature line */
    atemp_line = (float *)calloc((size_t)(input->size.s),sizeof(float));
    if (atemp_line == (float*)NULL)ERROR("allocating atemp line", "main");
	 
  /* Allocate memory for mask line */
    if ( param->cloud_flag ) {
    	mask_line = (char**)calloc((size_t)(lut->ar_region_size.l),sizeof(char *));
    	if (mask_line == (char**)NULL)ERROR("allocating mask line", "main");
    	mask_line_buf = (char*)calloc((size_t)(input_mask->size.s * lut->ar_region_size.l),sizeof(char));
    	if (mask_line_buf == (char*)NULL)ERROR("allocating mask line buffer", "main");
	 	for (il = 0; il < lut->ar_region_size.l; il++) {
	 		mask_line[il]=mask_line_buf;
			mask_line_buf += input->size.s;
	 	}
    }
	 
  /* Allocate memory for ddv line */
    ddv_line = (char**)calloc((size_t)(lut->ar_region_size.l),sizeof(char *));
    if (ddv_line == (char**)NULL)ERROR("allocating ddv line", "main");
    ddv_line_buf = (char*)calloc((size_t)(input->size.s * lut->ar_region_size.l),sizeof(char));
    if (ddv_line_buf == (char*)NULL)ERROR("allocating ddv line buffer", "main");
	 for (il = 0; il < lut->ar_region_size.l; il++) {
	 	ddv_line[il]=ddv_line_buf;
		ddv_line_buf += input->size.s;
	 }

  /* Allocate memory for rotating cloud buffer */

  rot_cld_buf=(char *)calloc((size_t)(input->size.s*lut->ar_region_size.l*3), sizeof(char));
  if (rot_cld_buf == (char *)NULL) 
      ERROR("allocating roatting cloud buffer (a)", "main");
  rot_cld_block_buf=(char **)calloc((size_t)(lut->ar_region_size.l*3), sizeof(char *));
  if (rot_cld_block_buf == (char **)NULL) 
    ERROR("allocating rotating cloud buffer (b)", "main");
  
  for (ib = 0; ib < 3; ib++) {
    rot_cld[ib]=rot_cld_block_buf;
	rot_cld_block_buf += lut->ar_region_size.l;
	for (il = 0; il < lut->ar_region_size.l; il++) {
		rot_cld[ib][il]=rot_cld_buf;
		rot_cld_buf+=input->size.s;
	}
  }


  /* Allocate memory for ar_gridcell */
  ar_gridcell.nbrows=lut->ar_size.l;
  ar_gridcell.nbcols=lut->ar_size.s;
  ar_gridcell.lat=(float *)calloc((size_t)(lut->ar_size.s * lut->ar_size.l),sizeof(float));
  if (ar_gridcell.lat == (float *) NULL)
    ERROR("allocating ar_gridcell.lat", "main");
  ar_gridcell.lon=(float *)calloc((size_t)(lut->ar_size.s * lut->ar_size.l),sizeof(float));
  if (ar_gridcell.lon == (float *)NULL)
    ERROR("allocating ar_gridcell.lon", "main");
  ar_gridcell.sun_zen=(float *)calloc((size_t)(lut->ar_size.s * lut->ar_size.l),sizeof(float));
  if (ar_gridcell.sun_zen == (float *)NULL)
    ERROR("allocating ar_gridcell.sun_zen", "main");
  ar_gridcell.view_zen=(float *)calloc((size_t)(lut->ar_size.s * lut->ar_size.l),sizeof(float));
  if (ar_gridcell.view_zen == (float *)NULL)
    ERROR("allocating ar_gridcell.view_zen", "main");
  ar_gridcell.rel_az=(float *)calloc((size_t)(lut->ar_size.s * lut->ar_size.l),sizeof(float));
  if (ar_gridcell.rel_az == (float *)NULL)
    ERROR("allocating ar_gridcell.rel_az", "main");
  ar_gridcell.wv=(float *)calloc((size_t)(lut->ar_size.s * lut->ar_size.l),sizeof(float));
  if (ar_gridcell.wv == (float *)NULL)
    ERROR("allocating ar_gridcell.wv", "main");
  ar_gridcell.spres=(float *)calloc((size_t)(lut->ar_size.s * lut->ar_size.l),sizeof(float));
  if (ar_gridcell.spres == (float *)NULL)
    ERROR("allocating ar_gridcell.spres", "main");
  ar_gridcell.ozone=(float *)calloc((size_t)(lut->ar_size.s * lut->ar_size.l),sizeof(float));
  if (ar_gridcell.ozone == (float *)NULL)
    ERROR("allocating ar_gridcell.ozone", "main");
  ar_gridcell.spres_dem=(float *)calloc((size_t)(lut->ar_size.s * lut->ar_size.l),sizeof(float));
  if (ar_gridcell.spres_dem == (float *)NULL)
    ERROR("allocating ar_gridcell.spres_dem", "main");


  /* Allocate memory for output lines */

  line_out_buf = (int *)calloc((size_t)(output->size.s * output->nband_tot), 
                               sizeof(int));
  if (line_out_buf == (int *)NULL) 
    ERROR("allocating output line buffer", "main");
  line_out[0] = line_out_buf;
  for (ib = 1; ib < output->nband_tot; ib++)
    line_out[ib] = line_out[ib - 1] + output->size.s;

  /* Allocate memory for the aerosol lines */

  line_ar = (int ***)calloc((size_t)lut->ar_size.l, sizeof(int **));
  if (line_ar == (int ***)NULL) 
    ERROR("allocating aerosol line buffer (a)", "main");

  line_ar_band_buf = (int **)calloc((size_t)(lut->ar_size.l * AERO_NB_BANDS), sizeof(int *));
  if (line_ar_band_buf == (int **)NULL) 
    ERROR("allocating aerosol line buffer (b)", "main");


  line_ar_buf = (int *)calloc((size_t)(lut->ar_size.l * lut->ar_size.s * AERO_NB_BANDS), 
                               sizeof(int));
  if (line_ar_buf == (int *)NULL) 
    ERROR("allocating aerosol line buffer (c)", "main");

  for (il = 0; il < lut->ar_size.l; il++) {
    line_ar[il] = line_ar_band_buf;
    line_ar_band_buf += AERO_NB_BANDS;
    for (ib = 0; ib < AERO_NB_BANDS; ib++) {
      line_ar[il][ib] = line_ar_buf;
      line_ar_buf += lut->ar_size.s;
    }
  }

  line_ar_stats = (int ***)calloc((size_t)lut->ar_size.l, sizeof(int **));
  if (line_ar_stats == (int ***)NULL) 
    ERROR("allocating aerosol stats line buffer (a)", "main");

  line_ar_stats_band_buf = (int **)calloc((size_t)(lut->ar_size.l * AERO_STATS_NB_BANDS), sizeof(int *));
  if (line_ar_stats_band_buf == (int **)NULL) 
    ERROR("allocating aerosol stats line buffer (b)", "main");

  line_ar_stats_buf = (int *)calloc((size_t)(lut->ar_size.l * lut->ar_size.s * AERO_STATS_NB_BANDS), 
                               sizeof(int));
  if (line_ar_stats_buf == (int *)NULL) 
    ERROR("allocating aerosol stats line buffer (c)", "main");

  for (il = 0; il < lut->ar_size.l; il++) {
    line_ar_stats[il] = line_ar_stats_band_buf;
    line_ar_stats_band_buf += AERO_STATS_NB_BANDS;
    for (ib = 0; ib < AERO_STATS_NB_BANDS; ib++) {
      line_ar_stats[il][ib] = line_ar_stats_buf;
      line_ar_stats_buf += lut->ar_size.s;
    }
  }

  /* Initialize the statistics */

  ar_stats.nfill = 0;
  ar_stats.first = true;

  for (ib = 0; ib < output->nband_tot-3; ib++) {
    sr_stats.nfill[ib] = 0;
    sr_stats.nout_range[ib] = 0;
    sr_stats.first[ib] = true;
  }
/****
	Get center lat lon and deviation from true north
****/
	img.l=input->size.l/2.; 
	img.s=input->size.s/2.; 
	img.is_fill=false;
	if (!FromSpace(space, &img, &geo))
    	ERROR("mapping from space (0)", "main");
   center_lat=geo.lat * DEG;
   center_lon=geo.lon * DEG;

  /* Compute scene gmt time */
  if ( (input->meta.acq_date.hour !=0 ) || (input->meta.acq_date.minute != 0 ) || (input->meta.acq_date.second !=0))
   
  scene_gmt=input->meta.acq_date.hour+input->meta.acq_date.minute/60.+input->meta.acq_date.second/3600.;
  else
  scene_gmt=10.5-center_lon/15.;
  if ( scene_gmt < 0.) scene_gmt=scene_gmt+24.;
  
/*  printf("Center: %d %f %f %f %f %f\n",input->meta.acq_date.doy,scene_gmt,center_lat,center_lon,input->meta.sun_zen*DEG,input->meta.sun_az*DEG);
  exit(1); */
   /* Read PRWV Data */
   if ( param->num_prwv_files > 0  ) {

     if (!get_prwv_anc(&anc_SP,prwv_input,prwv_in[SP_INDEX],SP_INDEX))
       ERROR("Can't get PRWV SP data","main");
     if (!get_prwv_anc(&anc_WV,prwv_input,prwv_in[WV_INDEX],WV_INDEX))
       ERROR("Can't get PRWV WV data","main");
     if (!get_prwv_anc(&anc_ATEMP,prwv_input,prwv_in[ATEMP_INDEX],ATEMP_INDEX))
       ERROR("Can't get PRWV ATEMP data","main");
     if (!no_ozone_file)
     	if (!get_ozon_anc(&anc_O3,ozon_input,ozon_in,OZ_INDEX))
       		ERROR("Can't get OZONE data","main");

   } else if ( param->num_ncep_files > 0  ) {

     anc_O3.data[0]=NULL;
     anc_O3.data[1]=NULL;
     anc_O3.data[2]=NULL;
     anc_O3.data[3]=NULL; 
     anc_O3.nblayers=4;
     anc_O3.timeres=6; 
     anc_O3.source=-1;
     strcpy(anc_O3.filename[0],param->ncep_file_name[0]);
     strcpy(anc_O3.filename[1],param->ncep_file_name[1]);
     strcpy(anc_O3.filename[2],param->ncep_file_name[2]);
     strcpy(anc_O3.filename[3],param->ncep_file_name[3]);

     if (read_grib_anc(&anc_O3,TYPE_OZONE_DATA))
          ERROR("Can't read NCEP Ozone data","main");

     anc_WV.data[0]=NULL;
     anc_WV.data[1]=NULL;
     anc_WV.data[2]=NULL;
     anc_WV.data[3]=NULL;
     anc_WV.nblayers=4;
     anc_WV.timeres=6; 
     anc_WV.source=-1;
     strcpy(anc_WV.filename[0],param->ncep_file_name[0]);
     strcpy(anc_WV.filename[1],param->ncep_file_name[1]);
     strcpy(anc_WV.filename[2],param->ncep_file_name[2]);
     strcpy(anc_WV.filename[3],param->ncep_file_name[3]);
     if (read_grib_anc(&anc_WV,TYPE_WV_DATA))
       ERROR("Can't read NCEP WV data","main");

     anc_SP.data[0]=NULL;
     anc_SP.data[1]=NULL;
     anc_SP.data[2]=NULL;
     anc_SP.data[3]=NULL;
     anc_SP.nblayers=4;
     anc_SP.timeres=6; 
     anc_SP.source=-1;
     strcpy(anc_SP.filename[0],param->ncep_file_name[0]);
     strcpy(anc_SP.filename[1],param->ncep_file_name[1]);
     strcpy(anc_SP.filename[2],param->ncep_file_name[2]);
     strcpy(anc_SP.filename[3],param->ncep_file_name[3]);
     if (read_grib_anc(&anc_SP,TYPE_SP_DATA))
       ERROR("Can't read NCEP SP data","main");

     anc_ATEMP.data[0]=NULL;
     anc_ATEMP.data[1]=NULL;
     anc_ATEMP.data[2]=NULL;
     anc_ATEMP.data[3]=NULL;
     anc_ATEMP.nblayers=4;
     anc_ATEMP.timeres=6; 
     anc_ATEMP.source=-1;
     strcpy(anc_ATEMP.filename[0],param->ncep_file_name[0]);
     strcpy(anc_ATEMP.filename[1],param->ncep_file_name[1]);
     strcpy(anc_ATEMP.filename[2],param->ncep_file_name[2]);
     strcpy(anc_ATEMP.filename[3],param->ncep_file_name[3]);
     if (read_grib_anc(&anc_ATEMP,TYPE_ATEMP_DATA))
       ERROR("Can't read NCEP SP data","main");

   } else {
     ERROR("No input NCEP or PRWV data specified","main");
   }

  /* Convert the units */

  for (i=0;i<anc_SP.nblayers;i++)
    for (j=0;j<anc_SP.nbrows*anc_SP.nbcols;j++)
      anc_SP.data[i][j] /= 100.;  /* convert Pascals into millibars */

  for (i=0;i<anc_WV.nblayers;i++)
    for (j=0;j<anc_WV.nbrows*anc_WV.nbcols;j++)
      anc_WV.data[i][j] /= 10.; /* convert original PRWV kg/m2 into g/cm2 */
  if (!no_ozone_file) {
   for (i=0;i<anc_O3.nblayers;i++)
     for (j=0;j<anc_O3.nbrows*anc_O3.nbcols;j++)
      anc_O3.data[i][j] /= 1000.;  /* convert to Dobson */
   }

   /* read DEM file */
   dem_name= (char*)(param->dem_flag ? param->dem_file : DEMFILE );
  /* Open file for SD access */
  sds_file_id = SDstart((char *)dem_name, DFACC_RDONLY);
  if (sds_file_id == HDF_ERROR) {
    ERROR("opening dem_file", "OpenDem");
		   }
  sds_index=0;		   
  sds_id= SDselect(sds_file_id,sds_index);
  status=  SDgetinfo(sds_id, sds_name, &rank, dim_sizes, &data_type,&n_attrs);
  printf("the name of the sds is %s\n",sds_name);
  printf("the rank of the sds is %d\n",rank);
  printf("the dimension of the sds is %d %d\n",dim_sizes[0],dim_sizes[1]);
  start[0]=0;
  start[1]=0;
  edges[0]=3600;
  edges[1]=7200;
  stride[0]=1;
  stride[1]=1;
  dem_array=(short *)malloc(DEM_NBLAT*DEM_NBLON*sizeof(short));
  status=SDreaddata(sds_id,start, stride, edges,dem_array);	   
  printf("the status of the read is %d\n",status);
  if (status != 0 ) {
  printf("Fatal error DEM file not read\n");
  exit(EXIT_FAILURE);
  }
  dem_available=1;


   if ( debug_flag ) {
     print_anc_data(&anc_SP,"SP_DATA");
     print_anc_data(&anc_WV,"WV_DATA");
     print_anc_data(&anc_ATEMP,"ATEMP_DATA");
     if (!no_ozone_file) 
      print_anc_data(&anc_O3,"OZONE_DATA");
   }

      print_anc_data(&anc_O3,"OZONE_DATA");

/****
	Get center lat lon and deviation from true north
****/
	img.l=input->size.l/2.; 
	img.s=input->size.s/2.; 
	img.is_fill=false;
	if (!FromSpace(space, &img, &geo))
    	ERROR("mapping from space (0)", "main");
   center_lat=geo.lat * DEG;
   center_lon=geo.lon * DEG;
	printf ("(y0,x0)=(%d,%d)  (lat0,lon0)=(%f,%f)\n",
		 (int)img.l,(int)img.s,(float)(geo.lat * DEG),(float)(geo.lon * DEG));

	delta_y=img.l;
	delta_x=img.s;

	img.l=input->size.l/2.-100.; 
	img.s=input->size.s/2.; 
	img.is_fill=false;
	if (!FromSpace(space, &img, &geo))
    	ERROR("mapping from space (0)", "main");
/*
	printf ("(y1,x1)=(%d,%d)  (lat1,lon1)=(%f,%f)\n",
		 (int)img.l,(int)img.s,(float)(geo.lat * DEG),(float)(geo.lon * DEG));
*/

	geo.lon=center_lon*RAD;
	geo.is_fill=false;
	if (!ToSpace(space, &geo, &img))
    	ERROR("mapping to space (0)", "main");
/*
	printf ("(y2,x2)=(%d,%d)  (lat2,lon2)=(%f,%f)\n",
		 (int)img.l,(int)img.s,(float)(geo.lat * DEG),(float)(geo.lon * DEG));
*/
	delta_y = delta_y - img.l;
	delta_x = img.s - delta_x;
	adjust_north=(float)(atan(delta_x/delta_y)*DEG);

	printf("True North adjustment = %f\n",adjust_north);


/****
	Recompute solar geometry
****/
/*
   	jday=(short)input->meta.acq_date.doy;    	  
	sun_angles (jday,scene_gmt,center_lat,center_lon,&fts,&ffs);
	input->meta.sun_zen=fts*RAD;
	input->meta.sun_az=ffs*RAD;

	printf("Center: %d %f %f %f %f %f\n",input->meta.acq_date.doy,scene_gmt,center_lat,center_lon,input->meta.sun_zen*DEG,input->meta.sun_az*DEG);
*/	

#ifdef SAVE_6S_RESULTS
	if (read_6S_results_from_file(SIXS_RESULTS_FILENAME,&sixs_tables)) {
#endif
/****
	Run 6S and compute atmcor params
****/
   	interpol_spatial_anc(&anc_WV,center_lat,center_lon,tmpflt_arr);
   	tmpint=(int)(scene_gmt/anc_WV.timeres);
   	if (tmpint>=(anc_WV.nblayers-1))
		tmpint=anc_WV.nblayers-2;
   	coef=(double)(scene_gmt-anc_WV.time[tmpint])/anc_WV.timeres;
   	sixs_tables.uwv=(1.-coef)*tmpflt_arr[tmpint]+coef*tmpflt_arr[tmpint+1];

   	if (!no_ozone_file) {
   		interpol_spatial_anc(&anc_O3,center_lat,center_lon,tmpflt_arr);
   		tmpint=(int)(scene_gmt/anc_O3.timeres);
   		if ( anc_O3.nblayers> 1 ){
      		if (tmpint>=(anc_O3.nblayers-1))tmpint=anc_O3.nblayers-2;
   			coef=(double)(scene_gmt-anc_O3.time[tmpint])/anc_O3.timeres;
    		sixs_tables.uoz=(1.-coef)*tmpflt_arr[tmpint]+coef*tmpflt_arr[tmpint+1];
		} else {
    	  	sixs_tables.uoz=tmpflt_arr[tmpint];
		}
   	} else {
         jday=(short)input->meta.acq_date.doy;    	  
         sixs_tables.uoz=calcuoz(jday,(float)center_lat);
   	}
/***
   interpol_spatial_anc(&anc_SP,center_lat,center_lon,tmpflt_arr);
   tmpint=(int)(scene_gmt/anc_SP.timeres);
   if (tmpint>=(anc_SP.nblayers-1))
		tmpint=anc_SP.nblayers-2;
   coef=(double)(scene_gmt-anc_SP.time[tmpint])/anc_SP.timeres;
   tmpflt=(1.-coef)*tmpflt_arr[tmpint]+coef*tmpflt_arr[tmpint+1];
	sixs_tables.target_alt=(1013.-tmpflt)/-100.; / * target altitude in km (negative=above ground)* /
***/

	sixs_tables.target_alt=0.; /* target altitude in km (sea level) */
	sixs_tables.sza=input->meta.sun_zen*DEG;
	sixs_tables.phi=input->meta.sun_az*DEG;
	sixs_tables.vza=0.;
	sixs_tables.month=9;
	sixs_tables.day=15;
	sixs_tables.srefl=0.14;


/*
	printf ("Center : Lat = %7.2f  Lon = %7.2f \n",center_lat,center_lon);
	printf ("          O3 = %7.2f   SP = %7.2f   WV = %7.2f\n",\
                sixs_tables.uoz,tmpflt,sixs_tables.uwv);

*/
	switch (input->meta.inst) {
		case INST_TM:
			sixs_tables.Inst=SIXS_INST_TM;
			break;
		case INST_ETM:
			sixs_tables.Inst=SIXS_INST_ETM;
			break;
		default:
			ERROR("Unknown Instrument", "main");
	}
	create_6S_tables(&sixs_tables);
#ifdef SAVE_6S_RESULTS
	write_6S_results_to_file(SIXS_RESULTS_FILENAME,&sixs_tables);
	}
#endif
/***
   interpolate ancillary data for AR grid cells
***/
  img.is_fill=false;
   sum_spres_anc=0.;
   sum_spres_dem=0.;
   nb_spres_anc=0;
   nb_spres_dem=0;
  for (il_ar = 0; il_ar < lut->ar_size.l;il_ar++) {
	 img.l=il_ar*lut->ar_region_size.l+lut->ar_region_size.l/2.;
     for (is_ar=0;is_ar < lut->ar_size.s; is_ar++) {
        img.s=is_ar*lut->ar_region_size.s+lut->ar_region_size.s/2.; 
	if (!FromSpace(space, &img, &geo))
    	   ERROR("mapping from space (1)", "main");

        ar_gridcell.lat[il_ar*lut->ar_size.s+is_ar]=geo.lat * DEG;
        ar_gridcell.lon[il_ar*lut->ar_size.s+is_ar]=geo.lon * DEG;
        ar_gridcell.sun_zen[il_ar*lut->ar_size.s+is_ar]=input->meta.sun_zen*DEG;
        ar_gridcell.view_zen[il_ar*lut->ar_size.s+is_ar]=3.5;
        ar_gridcell.rel_az[il_ar*lut->ar_size.s+is_ar]=input->meta.sun_az*DEG;

    	interpol_spatial_anc(&anc_WV,ar_gridcell.lat[il_ar*lut->ar_size.s+is_ar],ar_gridcell.lon[il_ar*lut->ar_size.s+is_ar],tmpflt_arr);
    	tmpint=(int)(scene_gmt/anc_WV.timeres);
    	if (tmpint>=(anc_WV.nblayers-1))
		tmpint=anc_WV.nblayers-2;
    	coef=(double)(scene_gmt-anc_WV.time[tmpint])/anc_WV.timeres;
    	ar_gridcell.wv[il_ar*lut->ar_size.s+is_ar]=(1.-coef)*tmpflt_arr[tmpint]+coef*tmpflt_arr[tmpint+1];

        if (!no_ozone_file) {
           interpol_spatial_anc(&anc_O3,ar_gridcell.lat[il_ar*lut->ar_size.s+is_ar],ar_gridcell.lon[il_ar*lut->ar_size.s+is_ar],tmpflt_arr);
    	   tmpint=(int)(scene_gmt/anc_O3.timeres);
           if ( anc_O3.nblayers> 1 ){
              if (tmpint>=(anc_O3.nblayers-1))
                 tmpint=anc_O3.nblayers-2;
    	      coef=(double)(scene_gmt-anc_O3.time[tmpint])/anc_O3.timeres;
    	      ar_gridcell.ozone[il_ar*lut->ar_size.s+is_ar]=(1.-coef)*tmpflt_arr[tmpint]+coef*tmpflt_arr[tmpint+1];
	   		} else {
    	      ar_gridcell.ozone[il_ar*lut->ar_size.s+is_ar]=tmpflt_arr[tmpint];
	   		}
        } else {
           jday=(short)input->meta.acq_date.doy;
    	   ar_gridcell.ozone[il_ar*lut->ar_size.s+is_ar]=calcuoz(jday,(float)ar_gridcell.lat[il_ar*lut->ar_size.s+is_ar]);
        }

    	interpol_spatial_anc(&anc_SP,ar_gridcell.lat[il_ar*lut->ar_size.s+is_ar],ar_gridcell.lon[il_ar*lut->ar_size.s+is_ar],tmpflt_arr);
    	tmpint=(int)(scene_gmt/anc_SP.timeres);
    	if (tmpint>=(anc_SP.nblayers-1))
			tmpint=anc_SP.nblayers-2;
    	coef=(double)(scene_gmt-anc_SP.time[tmpint])/anc_SP.timeres;
    	ar_gridcell.spres[il_ar*lut->ar_size.s+is_ar]=(1.-coef)*tmpflt_arr[tmpint]+coef*tmpflt_arr[tmpint+1];
        if (ar_gridcell.spres[il_ar*lut->ar_size.s+is_ar] > 0) {
           sum_spres_anc += ar_gridcell.spres[il_ar*lut->ar_size.s+is_ar];
           nb_spres_anc++;
        }
        if (dem_available) {
	   		ar_gridcell.spres_dem[il_ar*lut->ar_size.s+is_ar]=get_dem_spres(dem_array,ar_gridcell.lat[il_ar*lut->ar_size.s+is_ar],ar_gridcell.lon[il_ar*lut->ar_size.s+is_ar]);
           if (ar_gridcell.spres_dem[il_ar*lut->ar_size.s+is_ar] > 0) {
              sum_spres_dem += ar_gridcell.spres_dem[il_ar*lut->ar_size.s+is_ar];
              nb_spres_dem++;
           }
        }

     }
  }
  if (dem_available) {
     ratio_spres=1.;
     if ((nb_spres_anc > 0)&&(nb_spres_dem>0)) 
        ratio_spres=(float)((sum_spres_anc/nb_spres_anc)/(sum_spres_dem/nb_spres_dem));
/*	printf("Ratio pressure %f\n",ratio_spres);*/
     for (il_ar = 0; il_ar < lut->ar_size.l;il_ar++) 
        for (is_ar=0;is_ar < lut->ar_size.s; is_ar++) 
           if ((ar_gridcell.spres[il_ar*lut->ar_size.s+is_ar] > 0)&&(ar_gridcell.spres_dem[il_ar*lut->ar_size.s+is_ar] > 0))
 /*             ar_gridcell.spres[il_ar*lut->ar_size.s+is_ar]=ar_gridcell.spres_dem[il_ar*lut->ar_size.s+is_ar]*ratio_spres; Vermote 11/23/2010*/
              ar_gridcell.spres[il_ar*lut->ar_size.s+is_ar]=ar_gridcell.spres_dem[il_ar*lut->ar_size.s+is_ar]*ar_gridcell.spres[il_ar*lut->ar_size.s+is_ar]/1013.;
  }


/* Compute atmospheric coefs for the whole scene with aot550=0.01 for use in internal cloud screening : NAZMI */
	nbpts=lut->ar_size.l*lut->ar_size.s;
	/***
	Allocate memory for atmos_coeff
	***/
	if (allocate_mem_atmos_coeff(nbpts,&atmos_coef))
        ERROR("Allocating memory for atmos_coef", "main");

    printf("Compute Atmos Params with aot550=0.01\n"); fflush(stdout);
	update_atmos_coefs(&atmos_coef,&ar_gridcell, &sixs_tables,line_ar, lut,input->nband, 1);

  /* Read input first time and compute clear pixels stats for internal cloud screening */

/* allocate memory for cld_diags structure and clear sum and nb of obs */	
	if (allocate_cld_diags(&cld_diags,CLDDIAGS_CELLHEIGHT_5KM, CLDDIAGS_CELLWIDTH_5KM, input->size.l, input->size.s)) {
     		ERROR("couldn't allocate memory from cld_diags","main");
	}

  for (il = 0; il < input->size.l; il++) {
	if (!(il%100)) 
               { printf("Line %d\r",il); fflush(stdout); }
    /* Read each input band */
    for (ib = 0; ib < input->nband; ib++) {
      if (!GetInputLine(input, ib, il, line_in[0][ib]))
        ERROR("reading input data for a line (b)", "main");
    }
    if (!GetInputQALine(input, il, qa_line[0]))
      ERROR("reading input data for qa_line (1)", "main");
    if (param->thermal_band) {
     if (!GetInputLine(input_b6, 0, il, b6_line[0]))
       ERROR("reading input data for b6_line (1)", "main");
    }

    tmpint=(int)(scene_gmt/anc_ATEMP.timeres);
    if (tmpint>=(anc_ATEMP.nblayers-1))
	tmpint=anc_ATEMP.nblayers-2;
    coef=(double)(scene_gmt-anc_ATEMP.time[tmpint])/anc_ATEMP.timeres;
	img.is_fill=false;
/*
	img.l=il_ar*lut->ar_region_size.l+lut->ar_region_size.l/2.;
*/
	img.l=il;
    for (is=0;is < input->size.s; is++) {
	img.s=is;
		if (!FromSpace(space, &img, &geo))
    	   ERROR("mapping from space (2)", "main");
        flat=geo.lat * DEG;
        flon=geo.lon * DEG;

    	interpol_spatial_anc(&anc_ATEMP,flat,flon,tmpflt_arr);
    	atemp_line[is]=(1.-coef)*tmpflt_arr[tmpint]+coef*tmpflt_arr[tmpint+1];
	}
    /* Run Cld Screening Pass1 and compute stats */
	if (param->thermal_band)
		if (!cloud_detection_pass1(lut, input->size.s, il, line_in[0], qa_line[0], b6_line[0], atemp_line,&cld_diags))
      		ERROR("running cloud detection pass 1", "main");

  }
  if (param->thermal_band) {
	for (il=0;il<cld_diags.nbrows;il++) {
    	tmpint=(int)(scene_gmt/anc_ATEMP.timeres);
    	if (tmpint>=(anc_ATEMP.nblayers-1))
			tmpint=anc_ATEMP.nblayers-2;
    	coef=(double)(scene_gmt-anc_ATEMP.time[tmpint])/anc_ATEMP.timeres;
		img.is_fill=false;
		img.l=il*cld_diags.cellheight+cld_diags.cellheight/2.;
		if (img.l >= input->size.l)
			img.l = input->size.l-1;
		for (is=0;is<cld_diags.nbcols;is++) {
			img.s=is*cld_diags.cellwidth+cld_diags.cellwidth/2.;
			if (img.s >= input->size.s)
				img.s = input->size.s-1;
			if (!FromSpace(space, &img, &geo))
    	   		ERROR("mapping from space (3)", "main");
        	flat=geo.lat * DEG;
        	flon=geo.lon * DEG;

    		interpol_spatial_anc(&anc_ATEMP,flat,flon,tmpflt_arr);
			cld_diags.airtemp_2m[il][is]=(1.-coef)*tmpflt_arr[tmpint]+coef*tmpflt_arr[tmpint+1];

			if (cld_diags.nb_t6_clear[il][is] > 0) {
				sum_value=cld_diags.avg_t6_clear[il][is];
				sumsq_value=cld_diags.std_t6_clear[il][is];
				cld_diags.avg_t6_clear[il][is] = sum_value/cld_diags.nb_t6_clear[il][is];
				if (cld_diags.nb_t6_clear[il][is] > 1) {
/**
					cld_diags.std_t6_clear[il][is] = (sumsq_value+cld_diags.nb_t6_clear[il][is]*cld_diags.avg_t6_clear[il][is]*cld_diags.avg_t6_clear[il][is]-2.*cld_diags.avg_t6_clear[il][is]*sum_value)/(cld_diags.nb_t6_clear[il][is]-1);
					cld_diags.std_t6_clear[il][is]=sqrt(fabs(cld_diags.std_t6_clear[il][is]));
**/
					cld_diags.std_t6_clear[il][is] = (sumsq_value-(sum_value*sum_value)/cld_diags.nb_t6_clear[il][is])/(cld_diags.nb_t6_clear[il][is]-1);
					cld_diags.std_t6_clear[il][is]=sqrt(fabs(cld_diags.std_t6_clear[il][is]));

				} else 
					cld_diags.std_t6_clear[il][is] = 0.;
				sum_value=cld_diags.avg_b7_clear[il][is];
				sumsq_value=cld_diags.std_b7_clear[il][is];
				cld_diags.avg_b7_clear[il][is] = sum_value/cld_diags.nb_t6_clear[il][is];
				if (cld_diags.nb_t6_clear[il][is] > 1) {
					cld_diags.std_b7_clear[il][is] = (sumsq_value-(sum_value*sum_value)/cld_diags.nb_t6_clear[il][is])/(cld_diags.nb_t6_clear[il][is]-1);
					cld_diags.std_b7_clear[il][is]=sqrt(fabs(cld_diags.std_t6_clear[il][is]));
				} else
					cld_diags.std_b7_clear[il][is]=0;
			} else {
				cld_diags.avg_t6_clear[il][is]=-9999.;
				cld_diags.avg_b7_clear[il][is]=-9999.;
				cld_diags.std_t6_clear[il][is]=-9999.;
				cld_diags.std_b7_clear[il][is]=-9999.;
			}
		}
	}
	fill_cld_diags(&cld_diags);
#ifdef DEBUG_CLD
	for (il=0;il<cld_diags.nbrows;il++) 
		for (is=0;is<cld_diags.nbcols;is++) 
			if (fd_cld_diags != NULL)
				fprintf(fd_cld_diags,"%d %d %d %f %f %f %f %f\n",il,is,cld_diags.nb_t6_clear[il][is],cld_diags.airtemp_2m[il][is],cld_diags.avg_t6_clear[il][is],cld_diags.std_t6_clear[il][is],cld_diags.avg_b7_clear[il][is],cld_diags.std_b7_clear[il][is]);
	fclose(fd_cld_diags);
#endif
  }
/***
	Create dark target temporary file
***/
	charptr=tempnam(".","temporary_dark_target");
	if (charptr == NULL) ERROR("getting dark target temporary filename", "main");
	strcpy(tmpfilename,charptr);
	if ((fdtmp=fopen(tmpfilename,"w"))==NULL) ERROR("creating dark target temporary file", "main");

  /* Read input second time and create cloud and cloud shadow masks */
  ptr_rot_cld[0]=rot_cld[0];
  ptr_rot_cld[1]=rot_cld[1];
  ptr_rot_cld[2]=rot_cld[2];

  for (il_start = 0, il_ar = 0; 
       il_start < input->size.l; 
       il_start += lut->ar_region_size.l, il_ar++) {

    ar_gridcell.line_lat=&(ar_gridcell.lat[il_ar*lut->ar_size.s]);
    ar_gridcell.line_lon=&(ar_gridcell.lon[il_ar*lut->ar_size.s]);
    ar_gridcell.line_sun_zen=&(ar_gridcell.sun_zen[il_ar*lut->ar_size.s]);
    ar_gridcell.line_view_zen=&(ar_gridcell.view_zen[il_ar*lut->ar_size.s]);
    ar_gridcell.line_rel_az=&(ar_gridcell.rel_az[il_ar*lut->ar_size.s]);
    ar_gridcell.line_wv=&(ar_gridcell.wv[il_ar*lut->ar_size.s]);
    ar_gridcell.line_spres=&(ar_gridcell.spres[il_ar*lut->ar_size.s]);
    ar_gridcell.line_ozone=&(ar_gridcell.ozone[il_ar*lut->ar_size.s]);
    ar_gridcell.line_spres_dem=&(ar_gridcell.spres[il_ar*lut->ar_size.s]);
    
    il_end = il_start + lut->ar_region_size.l - 1;
    if (il_end >= input->size.l) il_end = input->size.l - 1;

    /* Read each input band for each line in region */

    for (il = il_start, il_region = 0; 
         il < (il_end + 1); 
	 il++, il_region++) {

      for (ib = 0; ib < input->nband; ib++) {
        if (!GetInputLine(input, ib, il, line_in[il_region][ib]))
          ERROR("reading input data for a line (a)", "main");
      }
      if (!GetInputQALine(input, il, qa_line[il_region]))
        ERROR("reading input data for qa_line (2)", "main");
	  if (param->thermal_band) {
      	if (!GetInputLine(input_b6, 0, il, b6_line[il_region]))
      	 ERROR("reading input data for b6_line (2)", "main");

        /* Run Cld Screening Pass2 */
	    if (!cloud_detection_pass2(lut, input->size.s, il, line_in[il_region], qa_line[il_region], b6_line[il_region], &cld_diags , ptr_rot_cld[1][il_region]))
      	  ERROR("running cloud detection pass 2", "main");
	  } else {
	    if (!cloud_detection_pass2(lut, input->size.s, il, line_in[il_region], qa_line[il_region], NULL, &cld_diags , ptr_rot_cld[1][il_region]))
      	  ERROR("running cloud detection pass 2", "main");
      }

    }
	if (param->thermal_band) {
		/* Cloud Mask Dilation : 5 pixels */
		if (!dilate_cloud_mask(lut, input->size.s, ptr_rot_cld, 5))
      		ERROR("running cloud mask dilation", "main");

		/* Cloud shadow */
		if (!cast_cloud_shadow(lut, input->size.s, il_start, line_in, b6_line, &cld_diags,ptr_rot_cld,&ar_gridcell,space_def.pixel_size,adjust_north))
      		ERROR("running cloud shadow detection", "main");

		/* Dilate Cloud shadow */
		dilate_shadow_mask(lut, input->size.s, ptr_rot_cld, 5);
	}
/***
	Save cloud and cloud shadow in temporary file
***/
	if (il_ar > 0)
		if (fwrite(ptr_rot_cld[0][0],lut->ar_region_size.l*input->size.s,1,fdtmp)!=1) ERROR("writing dark target to temporary file", "main");
/*	if (fwrite(ptr_rot_cld[1][0],lut->ar_region_size.l*input->size.s,1,fdtmp)!=1) ERROR("writing dark target to temporary file", "main");
*/
/**
	if (fwrite(line_ar[il_ar][0],lut->ar_size.s*2,1,fdtmp2)!=1) ERROR("writing aot1 to temporary file", "main");
**/
	ptr_tmp_cld=ptr_rot_cld[0];
	ptr_rot_cld[0]=ptr_rot_cld[1];
	ptr_rot_cld[1]=ptr_rot_cld[2];
	ptr_rot_cld[2]=ptr_tmp_cld;
	for (i=0;i<lut->ar_region_size.l;i++)
		memset(&ptr_rot_cld[2][i][0],0,input->size.s);
  }
/** Last Block */
  dilate_shadow_mask(lut, input->size.s, ptr_rot_cld, 5);
  if (fwrite(ptr_rot_cld[0][0],lut->ar_region_size.l*input->size.s,1,fdtmp)!=1) ERROR("writing dark target to temporary file", "main");
   fclose(fdtmp);


/***
	Open temporary file for read and write
***/
	if ((fdtmp=fopen(tmpfilename,"r+"))==NULL) ERROR("opening dark target temporary file (r+)", "main");
/*	
	if ((fdtmp2=fopen("Temporary_AOT1_File.dat","w"))==NULL) ERROR("creating temporary aot1 file", "main");
*/
  /* Read input second time and compute the aerosol for each region */

  for (il_start = 0, il_ar = 0; 
       il_start < input->size.l; 
       il_start += lut->ar_region_size.l, il_ar++) {

    ar_gridcell.line_lat=&(ar_gridcell.lat[il_ar*lut->ar_size.s]);
    ar_gridcell.line_lon=&(ar_gridcell.lon[il_ar*lut->ar_size.s]);
    ar_gridcell.line_sun_zen=&(ar_gridcell.sun_zen[il_ar*lut->ar_size.s]);
    ar_gridcell.line_view_zen=&(ar_gridcell.view_zen[il_ar*lut->ar_size.s]);
    ar_gridcell.line_rel_az=&(ar_gridcell.rel_az[il_ar*lut->ar_size.s]);
    ar_gridcell.line_wv=&(ar_gridcell.wv[il_ar*lut->ar_size.s]);
    ar_gridcell.line_spres=&(ar_gridcell.spres[il_ar*lut->ar_size.s]);
    ar_gridcell.line_ozone=&(ar_gridcell.ozone[il_ar*lut->ar_size.s]);
    ar_gridcell.line_spres_dem=&(ar_gridcell.spres[il_ar*lut->ar_size.s]);
    
    il_end = il_start + lut->ar_region_size.l - 1;
    if (il_end >= input->size.l) il_end = input->size.l - 1;
	 
	if (fseek(fdtmp,(long)(il_ar*lut->ar_region_size.l*input->size.s),SEEK_SET)) ERROR("seeking in temporary file (r)", "main");
	if (fread(ddv_line[0],lut->ar_region_size.l*input->size.s,1,fdtmp)!=1) ERROR("reading dark target to temporary file", "main");

    /* Read each input band for each line in region */

    for (il = il_start, il_region = 0; 
         il < (il_end + 1); 
	 il++, il_region++) {
      for (ib = 0; ib < input->nband; ib++) {
        if (!GetInputLine(input, ib, il, line_in[il_region][ib]))
          ERROR("reading input data for a line (a)", "main");
      }

    	if ( param->cloud_flag ) {
/*	printf("reading acca\n");*/
      	if (!GetInputMaskLine(input_mask, il, mask_line[il_region]))
        		ERROR("reading input mask data for a line", "main");
    	}

    }

    /* Compute the aerosol for the regions */

    printf("%d\n",il_ar); fflush(stdout);
/*  printf("%d\r",il_ar); fflush(stdout);    */
#ifdef DEBUG_AR
	diags_il_ar=il_ar;
#endif
    if (!Ar(il_ar,lut, &input->size, line_in, param->cloud_flag, mask_line, ddv_line,line_ar[il_ar], 
	 			line_ar_stats[il_ar], &ar_stats,&ar_gridcell,&sixs_tables))
      ERROR("computing aerosl", "main");
/***
	Save dark target map in temporary file
***/
	if (fseek(fdtmp,il_ar*lut->ar_region_size.l*input->size.s,SEEK_SET)) ERROR("seeking in temporary file (w)", "main");
	if (fwrite(ddv_line[0],lut->ar_region_size.l*input->size.s,1,fdtmp)!=1) ERROR("writing dark target to temporary file", "main");
/**
	if (fwrite(line_ar[il_ar][0],lut->ar_size.s*2,1,fdtmp2)!=1) ERROR("writing aot1 to temporary file", "main");
**/
  }
  printf("\n");
  fclose(fdtmp);
/**
  fclose(fdtmp2);
**/
#ifdef DEBUG_AR
	fclose(fd_ar_diags);
#endif
	/***
	Fill Gaps in the coarse resolution aerosol product for bands 1(0), 2(1) and 3(2)
	**/
    printf("write Fill Gaps ..."); fflush(stdout);
   Fill_Ar_Gaps(lut, line_ar, 0);
/*    printf("WARNING NOT FILLING GAPS IN THE AEROSOL");*/
    printf("Done\n"); fflush(stdout);
/*
  Fill_Ar_Gaps(lut, line_ar, 1);
  Fill_Ar_Gaps(lut, line_ar, 2);
*/
/* Compute atmospheric coefs for the whole scene using retrieved aot : NAZMI */
     nbpts=lut->ar_size.l*lut->ar_size.s;

    printf("Compute Atmos Params\n"); fflush(stdout);
#ifdef NO_AEROSOL_CORRECTION
	update_atmos_coefs(&atmos_coef,&ar_gridcell, &sixs_tables,line_ar, lut,input->nband, 1);
#else
	update_atmos_coefs(&atmos_coef,&ar_gridcell, &sixs_tables,line_ar, lut,input->nband, 0); /*Eric COMMENTED TO PERFORM NO CORRECTION*/
/*        printf("WARNING NO AEROSOL CORRECTION TEST MODE");
	update_atmos_coefs(&atmos_coef,&ar_gridcell, &sixs_tables,line_ar, lut,input->nband, 1); */
#endif

  /* Re-read input and compute surface reflectance */

/***
	Open temporary file for read
***/
	if ((fdtmp=fopen(tmpfilename,"r"))==NULL) ERROR("opening dark target temporary file", "main");

  for (il = 0; il < input->size.l; il++) {
	if (!(il%100)) 
               { printf("Line %d\r",il); fflush(stdout); }
    /* Re-read each input band */

    for (ib = 0; ib < input->nband; ib++) {
      if (!GetInputLine(input, ib, il, line_in[0][ib]))
        ERROR("reading input data for a line (b)", "main");
    }
    
     if (!GetInputLine(input_b6, 0, il, b6_line[0]))
       ERROR("reading input data for b6_line (1)", "main");


    if ( param->cloud_flag ) {
      if (!GetInputMaskLine(input_mask, il, mask_line[0]))
        ERROR("reading input mask data for a line", "main");
    }
    /* Compute the surface reflectance */
	if (param->cloud_flag) {
    	if (!Sr(lut, input->size.s, il, line_in[0], param->cloud_flag, mask_line[0],
            line_out, &sr_stats))
      		ERROR("computing surface reflectance for a line", "main");
	} else {
    	if (!Sr(lut, input->size.s, il, line_in[0], param->cloud_flag, (char *)NULL,
            line_out, &sr_stats))
      		ERROR("computing surface reflectance for a line", "main");
	}

/***
	Read line from dark target temporary file
***/
	if (fread(ddv_line[0],input->size.s,1,fdtmp)!=1) ERROR("reading line from dark target temporary file", "main");

	loc.l=il;
  	 i_aot=il/lut->ar_region_size.l;
	 for (is=0;is<input->size.s;is++) {
	 	loc.s=is;
		j_aot=is/lut->ar_region_size.s;

		if (line_in[0][0][is] != lut->in_fill) {
/* commented by Vermote 12/17/2010 after sucessfull debugging

   		line_out[lut->nband][is]=ar_gridcell.spres[i_aot*lut->ar_size.s+j_aot];

  */
 
			ArInterp(lut, &loc, line_ar, inter_aot); 
			line_out[lut->nband][is] = inter_aot[0];


			if (param->cloud_flag ) {
				line_out[lut->nband+1][is]=0x00;   /* Valid data */
				if (mask_line[0][is] == lut->cloud_cloud)
					line_out[lut->nband+1][is] |= 0x04;  /* set cloud bit */
				if (mask_line[0][is] == lut->cloud_land)
					line_out[lut->nband+1][is] |= 0x20;  /* set land bit */
				if (mask_line[0][is] == lut->cloud_snow)
					line_out[lut->nband+1][is] |= 0x10;  /* set snow bit */
				if (mask_line[0][is] == lut->cloud_shadow)
					line_out[lut->nband+1][is] |= 0x08;  /* set shadow bit */
			} else {
	   		  line_out[lut->nband+1][is]=0x02;  /* QA indicates invalid QA bits */
			}
			if ((ddv_line[0][is]&0x01) == 1)
				line_out[lut->nband+1][is] |= 0x40;  /* set dark target bit */
/**
Set bits for internal cloud mask
bit 8: cloud
bit 9: cloud shadow
bit 10: snow
bit 11: land/water mask
bit 12: adjacent cloud
**/
		if (ddv_line[0][is]&0x10)
			line_out[lut->nband+1][is] |= 0x0800;  /* set internal water mask bit */
		else
			line_out[lut->nband+1][is] &= 0xf7ff;  /* reset internal water mask bit */
		if (ddv_line[0][is]&0x20)
			line_out[lut->nband+1][is] |= 0x0100;  /* set internal cloud mask bit */
		else
			line_out[lut->nband+1][is] &= 0xfeff;  /* reset internal cloud mask bit */
		if (ddv_line[0][is]&0x04)
			line_out[lut->nband+1][is] |= 0x1000;  /* set adjacent cloud mask bit */
		else
			line_out[lut->nband+1][is] &= 0xefff;  /* reset adjacent cloud mask bit */
		if (ddv_line[0][is]&0x40)
			line_out[lut->nband+1][is] |= 0x0200;  /* set internal cloud shadow mask bit */
		else
			line_out[lut->nband+1][is] &= 0xfdff;  /* reset internal cloud shadow mask bit */
		if (ddv_line[0][is]&0x80)
			line_out[lut->nband+1][is] |= 0x0400;  /* set internal snow mask bit */
		else
		line_out[lut->nband+1][is] &= 0xfbff;  /* reset internal snow mask bit */
	   	line_out[lut->nband+2][is]=line_ar_stats[i_aot][0][j_aot];
	  	line_out[lut->nband+3][is]=line_ar_stats[i_aot][1][j_aot];
	   	line_out[lut->nband+4][is]=line_ar_stats[i_aot][2][j_aot];
		} else {
	   	line_out[lut->nband][is]=lut->aerosol_fill;
	   	line_out[lut->nband+1][is]=0x03;  /* QA indicates fill value */
	   	line_out[lut->nband+2][is]=0;
	  	line_out[lut->nband+3][is]=lut->in_fill;
	   	line_out[lut->nband+4][is]=lut->in_fill;
		}
		/* try to redo the cloud mask Vermote May 29 2007 */
/* reset cloud shadow and cloud adjacent bits */
                line_out[lut->nband+1][is] &= 0xfdff;
                line_out[lut->nband+1][is] &= 0xefff;
				
		anom=line_out[0][is]-line_out[2][is]/2.;
		t6=b6_line[0][is]/100.+273.;
		t6s_seuil=280.+(1000.*0.01);
/*		printf ("this is anom %d\n",anom);*/
		if (( ( anom > 300 ) && ( line_out[4][is] > 300) && ( t6 < t6s_seuil) ) || ( (line_out[2][is] > 5000) && ( t6 < t6s_seuil)))
		         {
			line_out[lut->nband+1][is] |= 0x0100;  /* set internal cloud mask bit */
/*		printf ("this is anom >300 and temperature %d %f\n",anom,t6);*/
			}
		else
			line_out[lut->nband+1][is] &= 0xfeff;  /* reset internal cloud mask bit */
		    
	 }
    /* Write each output band */

    for (ib = 0; ib < output->nband_tot-3; ib++) {
      if (!PutOutputLine(output, ib, il, line_out[ib]))
        ERROR("writing output data for a line", "main");
    }
  }
  printf("\n");
	fclose(fdtmp);
	unlink(tmpfilename); 
	
  /* Print the statistics */

  printf(" total pixels %ld\n", ((long)input->size.l * (long)input->size.s));
  printf(" aerosol coarse  nfill %ld  min  %d  max  %d\n", 
         ar_stats.nfill, ar_stats.ar_min, ar_stats.ar_max);

  for (ib = 0; ib < output->nband_tot; ib++) {
    printf(" sr %s  nfill %ld  nout_range %ld  min  %d  max  %d\n", 
            output->sds_sr[ib].name, 
	    sr_stats.nfill[ib], sr_stats.nout_range[ib],
	    sr_stats.sr_min[ib], sr_stats.sr_max[ib]);
  }
  
  /* Write the output metadata */
  
  if (!PutMetadata(output, input->nband, &input->meta, param,lut, &bounds))
    ERROR("writing the metadata", "main");

  /* Close input files */

  if (!CloseInput(input)) ERROR("closing input file", "main");

  if (!CloseOutput(output)) ERROR("closing input file", "main");

  nsds=output->nband_tot-3;
  for (isds = 0; isds < nsds; isds++) {
    sds_names[isds] = output->sds_sr[isds].name;
    sds_types[isds] = output->sds_sr[isds].type;
  }
  if (!PutSpaceDefHDF(&space_def, output->file_name, nsds, 
                      sds_names, sds_types, grid_name))
    ERROR("putting space metadata in HDF file", "main");

  /* Free memory */
  
	free_mem_atmos_coeff(&atmos_coef);

  if (!FreeInput(input)) 
    ERROR("freeing input file stucture", "main");

  if (!FreeInput(input_b6)) 
    ERROR("freeing input file stucture", "main");

  if (!FreeLut(lut)) 
    ERROR("freeing lut file stucture", "main");

  if (!FreeOutput(output)) 
    ERROR("freeing output file stucture", "main");
  if (!FreeSpace(space))
    ERROR("freeing space", "main");


  free(line_out[0]);
  free(line_ar[0][0]);
  free(line_ar[0]);
  free(line_ar);
  free(line_ar_stats[0][0]);
  free(line_ar_stats[0]);
  free(line_ar_stats);
  free(line_in[0][0]);
  free(line_in[0]);
  free(line_in);
  free(qa_line[0]);
  free(qa_line);
  if (param->thermal_band) {
  	free(b6_line[0]);
  	free(b6_line);
  }
  if (param->cloud_flag ) {
    free(mask_line[0]);
    free(mask_line);
  }
  free(ddv_line[0]);
  free(ddv_line);
  free(rot_cld[0][0]);
  free(rot_cld[0]);
  free(ar_gridcell.lat);
  free(ar_gridcell.lon);
  free(ar_gridcell.sun_zen);
  free(ar_gridcell.view_zen);
  free(ar_gridcell.rel_az);
  free(ar_gridcell.wv);
  free(ar_gridcell.spres);
  free(ar_gridcell.ozone);

  for (ifree=0; ifree<(param->num_ncep_files>0?4:1); ifree++) {
     if (anc_O3.data[ifree]!=NULL) free(anc_O3.data[ifree]);
     if (anc_WV.data[ifree]!=NULL) free(anc_WV.data[ifree]);
     if (anc_SP.data[ifree]!=NULL) free(anc_SP.data[ifree]);
  }
  if (dem_available)
     free(dem_array);
  if (!FreeParam(param)) 
    ERROR("freeing parameter stucture", "main");


  /* All done */

  exit(EXIT_SUCCESS);

  return (EXIT_SUCCESS);
}

	  
int allocate_mem_atmos_coeff(int nbpts,atmos_t *atmos_coef) {
	int ib;
	if ((atmos_coef->computed=(int *)malloc(nbpts*sizeof(int)))==(int *)NULL)
		return -1;
	for (ib=0;ib<7;ib++) {
		if ((atmos_coef->tgOG[ib]=(float *)malloc(nbpts*sizeof(float)))==(float *)NULL)
			return -1;
		if ((atmos_coef->tgH2O[ib]=(float *)malloc(nbpts*sizeof(float)))==(float *)NULL)
			return -1;
		if ((atmos_coef->td_ra[ib]=(float *)malloc(nbpts*sizeof(float)))==(float *)NULL)
			return -1;
		if ((atmos_coef->tu_ra[ib]=(float *)malloc(nbpts*sizeof(float)))==(float *)NULL)
			return -1;
		if ((atmos_coef->rho_mol[ib]=(float *)malloc(nbpts*sizeof(float)))==(float *)NULL)
			return -1;
		if ((atmos_coef->rho_ra[ib]=(float *)malloc(nbpts*sizeof(float)))==(float *)NULL)
			return -1;
		if ((atmos_coef->td_da[ib]=(float *)malloc(nbpts*sizeof(float)))==(float *)NULL)
			return -1;
		if ((atmos_coef->tu_da[ib]=(float *)malloc(nbpts*sizeof(float)))==(float *)NULL)
			return -1;
		if ((atmos_coef->S_ra[ib]=(float *)malloc(nbpts*sizeof(float)))==(float *)NULL)
			return -1;
		if ((atmos_coef->td_r[ib]=(float *)malloc(nbpts*sizeof(float)))==(float *)NULL)
			return -1;
		if ((atmos_coef->tu_r[ib]=(float *)malloc(nbpts*sizeof(float)))==(float *)NULL)
			return -1;
		if ((atmos_coef->S_r[ib]=(float *)malloc(nbpts*sizeof(float)))==(float *)NULL)
			return -1;
		if ((atmos_coef->rho_r[ib]=(float *)malloc(nbpts*sizeof(float)))==(float *)NULL)
			return -1;
	}
	return 0;
}

int free_mem_atmos_coeff(atmos_t *atmos_coef) {
	int ib;
	free(atmos_coef->computed);
	for(ib=0;ib<7;ib++) {
		free(atmos_coef->tgOG[ib]);
		free(atmos_coef->tgH2O[ib]);
		free(atmos_coef->td_ra[ib]);
		free(atmos_coef->tu_ra[ib]);
		free(atmos_coef->rho_mol[ib]);
		free(atmos_coef->rho_ra[ib]);
		free(atmos_coef->td_da[ib]);
		free(atmos_coef->tu_da[ib]);
		free(atmos_coef->S_ra[ib]);
		free(atmos_coef->td_r[ib]);
		free(atmos_coef->tu_r[ib]);
		free(atmos_coef->S_r[ib]);
		free(atmos_coef->rho_r[ib]);
	}
	return 0;
}

float calcuoz(short jday,float flat) {
/******************************************************************************
!C
!Routine: calcuoz

!Description:  Gets ozone concentration for a particular day and 
               latitude, interpolating if necessary.

!Revision History:
 Original version:    Nazmi Z El Saleous and Eric Vermote

!Input Parameters:
        jday       Julian day
        lat        latitude (in degrees)
        
!Output Parameters:
        uoz        Ozone concentration (in cm-atm) 

!Return value:
        returns 0 (success)
               -1 (latitude is beyond coverage of data, and so
                  the returned ozone value is an approximation).
                
                
!References and Credits:
      Data is from:
    LONDON J.,BOJKOV R.D.,OLTMANS S. AND KELLEY J.I.,1976,
    ATLAS OF THE GLOBAL DISTRIBUTION OF TOTAL OZONE .JULY 1957-JUNE 1967
    NCAR TECHNICAL NOTE, NCAR/TN/113+STR,PP276
      
!Developers:
      Nazmi Z El Saleous
      Eric Vermote
      University of Maryland / Dept. of Geography
      nazmi.elsaleous@gsfc.nasa.gov

!Design Notes:  
      
!END
*******************************************************************************/
   float t,u,tmpf;
   int i1,i2,j1,j2,Minf,Msup,Latinf,Latsup;
   

/*
CREPARTITION ZONALE PAR BANDE DE 10 DEG DE LATITUDE A PARTIR DE 80 SUD
C            ---   OZONE   ---
*/
   float oz[12][17] = {
   {.315,.320,.315,.305,.300,.280,.260,.240,.240,.240,.250,.280,.320,.350,.375,.380,.380},
   {.280,.300,.300,.300,.280,.270,.260,.240,.240,.240,.260,.300,.340,.380,.400,.420,.420},
   {.280,.280,.280,.280,.280,.260,.250,.240,.250,.250,.270,.300,.340,.400,.420,.440,.440},
   {.280,.280,.280,.280,.280,.260,.250,.250,.250,.260,.280,.300,.340,.380,.420,.430,.430},
   {.280,.290,.300,.300,.280,.270,.260,.250,.250,.260,.270,.300,.320,.360,.380,.400,.400},
   {.280,.300,.300,.305,.300,.280,.260,.250,.250,.260,.260,.280,.310,.330,.360,.370,.370},
   {.290,.300,.315,.320,.305,.280,.260,.250,.240,.240,.260,.270,.290,.310,.320,.320,.320},
   {.300,.310,.320,.325,.320,.300,.270,.260,.240,.240,.250,.260,.280,.290,.300,.300,.290},
   {.300,.320,.325,.335,.320,.300,.280,.260,.240,.240,.240,.260,.270,.280,.280,.280,.280},
   {.320,.340,.350,.345,.330,.300,.280,.260,.240,.240,.240,.260,.260,.280,.280,.280,.280},
   {.360,.360,.360,.340,.320,.300,.280,.260,.240,.240,.240,.260,.280,.300,.310,.310,.300},
   {.340,.350,.340,.320,.310,.280,.260,.250,.240,.240,.240,.260,.300,.320,.330,.340,.330}};
/*     
C      /Begin of interpolation/
C      /Find Neighbours/
C      / First loop for time/
*/
   if (fabs(flat)>=80.) {
    tmpf=.270;
    return -1;
   }
   Minf= (int) ((jday-15.)/30.5);
   if (jday < 15) 
    Minf=Minf-1;
   Latinf=(int) (flat/10.);
   if (flat < 0.) 
    Latinf=Latinf-1;
   t=((jday-15.)-(30.5*Minf))/30.5;
   u=(flat-10.*Latinf)/10.;
   
   Minf=(Minf+12)%12;
   Msup=(Minf+13)%12;
   Latsup=Latinf+1;
   i1=Minf;
   j1=Latinf+8;
   i2=Msup;
   j2=Latsup+8;
   
/*     /Now Calculate Uo3 at the given point Xlat,Jjulien/ */
   tmpf=(1.-t)*(1.-u)*oz[i1][j1]+t*(1.-u)*oz[i2][j1]+t*u*oz[i2][j2]+(1.-t)*u*oz[i1][j2];
   return tmpf;
}

float get_dem_spres(short *dem,float lat,float lon) {
	int i,j,idem,jdem;
	float dem_spres;
		
	idem=(int)((DEM_LATMAX-lat)/DEM_DLAT+0.5);
	if (idem<0)
		idem=0;
	if (idem >= DEM_NBLAT)
		idem=DEM_NBLAT-1;
	jdem=(int)((lon-DEM_LONMIN)/DEM_DLON+0.5);
	if (jdem<0)
		jdem=0;
	if (jdem >= DEM_NBLON)
		jdem=DEM_NBLON-1;
	if (dem[idem*DEM_NBLON+jdem]== -9999)
		dem_spres=1013;
	else
		dem_spres=1013.2*exp(-dem[idem*DEM_NBLON+jdem]/8000.);

	return dem_spres;
}
void swapbytes(void *val,int nbbytes) {
/******************************************************************************
!C
!Routine: swapbytes

!Description:  Does a byte reversal on value "val."

!Revision History:
 Original version:    Nazmi Z Saleous

!Input Parameters:
        val           an array of values (cast to void)
        nbbytes       the length of val (in bytes; cannot be greater than 17 
                      bytes)
        
!Output Parameters:
        val           (updated)

!Return value:
        none
                
!References and Credits:

!Developers:
      Nazmi Z Saleous
      nazmi.saleous@gsfc.nasa.gov

!Design Notes:  
      
!END
*******************************************************************************/
   char *tmpptr1,*tmpptr2,tmpstr[17];
   int i;

   memcpy(tmpstr,val,nbbytes);
   tmpptr1=(char *) val;
   tmpptr2=(char *) tmpstr;
   for (i=0;i<nbbytes;i++)
    tmpptr1[i]=tmpptr2[nbbytes-i-1];
}

int update_atmos_coefs(atmos_t *atmos_coef,Ar_gridcell_t *ar_gridcell, sixs_tables_t *sixs_tables,int ***line_ar,Lut_t *lut,int nband, int bkgd_aerosol) {
	int irow,icol;
	for (irow=0;irow<ar_gridcell->nbrows;irow++)
		for (icol=0;icol<ar_gridcell->nbcols;icol++)
			update_gridcell_atmos_coefs(irow,icol,atmos_coef,ar_gridcell, sixs_tables,line_ar[irow],lut,nband,bkgd_aerosol);
	return 0;
}

int update_gridcell_atmos_coefs(int irow,int icol,atmos_t *atmos_coef,Ar_gridcell_t *ar_gridcell, sixs_tables_t *sixs_tables,int **line_ar,Lut_t *lut,int nband, int bkgd_aerosol) {
	int ib,ipt,k;
	float mus,muv,phi,ratio_spres,tau_ray,aot550;
	double coef;
	float actual_rho_ray,actual_T_ray_up,actual_T_ray_down,actual_S_r;
	float rho_ray_P0,T_ray_up_P0,T_ray_down_P0,S_r_P0;
	float lamda[7]={486.,570.,660.,835.,1669.,0.,2207.};
	float tau_ray_sealevel[7]={0.16511,0.08614,0.04716,0.01835,0.00113,0.00037}; /* index=5 => band 7 */

		ipt=irow*ar_gridcell->nbcols+icol;	
		mus=cos(ar_gridcell->sun_zen[ipt]*RAD);
		muv=cos(ar_gridcell->view_zen[ipt]*RAD);
		phi=ar_gridcell->rel_az[ipt];
		ratio_spres=ar_gridcell->spres[ipt]/1013.;
/* eric is trying something */
/*                if ((irow == 53 )&&(icol == 83)) 
		{
		   printf("let's see what is going on\n"); 
		   }
		   */
		if (bkgd_aerosol) {
			atmos_coef->computed[ipt]=1;
			aot550=0.01;
		} else {
			if (line_ar[0][icol] != lut->aerosol_fill) {
				atmos_coef->computed[ipt]=1;
				aot550=((float)line_ar[0][icol]/1000.)*pow((550./lamda[0]),-1.);
			} else {
/*
				atmos_coef->computed[ipt]=0;
*/
                                
				atmos_coef->computed[ipt]=1;
				aot550=0.01;
			}
		}

		for (k=1;k<SIXS_NB_AOT;k++) {
			if (aot550 < sixs_tables->aot[k])
				break;
		}
		k--;
		if (k>=(SIXS_NB_AOT-1))
			k=SIXS_NB_AOT-2;
		coef=(aot550-sixs_tables->aot[k])/(sixs_tables->aot[k+1]-sixs_tables->aot[k]);


/*                printf("pressure ratio %d %d %f \n",irow,icol,ratio_spres);*/
		for (ib=0;ib < nband; ib++) {
			atmos_coef->tgOG[ib][ipt]=sixs_tables->T_g_og[ib];				
			atmos_coef->tgH2O[ib][ipt]=sixs_tables->T_g_wv[ib];				
			atmos_coef->td_ra[ib][ipt]=(1.-coef)*sixs_tables->T_ra_down[ib][k]+coef*sixs_tables->T_ra_down[ib][k+1];
			atmos_coef->tu_ra[ib][ipt]=(1.-coef)*sixs_tables->T_ra_up[ib][k]+coef*sixs_tables->T_ra_up[ib][k+1];
			atmos_coef->rho_mol[ib][ipt]=sixs_tables->rho_r[ib];				
			atmos_coef->rho_ra[ib][ipt]=(1.-coef)*sixs_tables->rho_ra[ib][k]+coef*sixs_tables->rho_ra[ib][k+1];
			atmos_coef->td_da[ib][ipt]=(1.-coef)*sixs_tables->T_a_down[ib][k]+coef*sixs_tables->T_a_down[ib][k+1];
			atmos_coef->tu_da[ib][ipt]=(1.-coef)*sixs_tables->T_a_up[ib][k]+coef*sixs_tables->T_a_up[ib][k+1];
			atmos_coef->S_ra[ib][ipt]=(1.-coef)*sixs_tables->S_ra[ib][k]+coef*sixs_tables->S_ra[ib][k+1];
/**
			compute DEM-based pressure correction for each grid point
**/
			tau_ray=tau_ray_sealevel[ib]*ratio_spres;
	
			chand(&phi,&muv,&mus,&tau_ray,&actual_rho_ray);

			actual_T_ray_down=((2./3.+mus)+(2./3.-mus)*exp(-tau_ray/mus))/(4./3.+tau_ray); /* downward */
			actual_T_ray_up = ((2./3.+muv)+(2./3.-muv)*exp(-tau_ray/muv))/(4./3.+tau_ray); /* upward */

			csalbr(&tau_ray,&actual_S_r);
						
			rho_ray_P0=sixs_tables->rho_r[ib];
			T_ray_down_P0=sixs_tables->T_r_down[ib];
			T_ray_up_P0=sixs_tables->T_r_up[ib];
			S_r_P0=sixs_tables->S_r[ib];

			atmos_coef->rho_ra[ib][ipt]=actual_rho_ray+(atmos_coef->rho_ra[ib][ipt]-rho_ray_P0); /* will need to correct for uwv/2 */
			atmos_coef->td_ra[ib][ipt] *= (actual_T_ray_down/T_ray_down_P0);
			atmos_coef->tu_ra[ib][ipt] *= (actual_T_ray_up/T_ray_up_P0);
			atmos_coef->S_ra[ib][ipt] = atmos_coef->S_ra[ib][ipt]-S_r_P0+actual_S_r;
			atmos_coef->td_r[ib][ipt] = actual_T_ray_down;
			atmos_coef->tu_r[ib][ipt]=actual_T_ray_up;
			atmos_coef->S_r[ib][ipt]=actual_S_r;
			atmos_coef->rho_r[ib][ipt]=actual_rho_ray;
					
		}
	return 0;
}


void sun_angles (short jday,float gmt,float flat,float flon,float *ts,float *fs) {
      
      double mst,tst,tet,et,ha,delta;
      double dlat,dlon,amuzero,elev,az,caz,azim;
      
      double A1=.000075,A2=.001868,A3=.032077,A4=.014615,A5=.040849;
      double B1=.006918,B2=.399912,B3=.070257,B4=.006758;
      double B5=.000907,B6=.002697,B7=.001480;
      
      dlat=(double)flat*M_PI/180.;
      dlon=(double)flon*M_PI/180.;
		                  
/*      
      SOLAR POSITION (ZENITHAL ANGLE ThetaS,AZIMUTHAL ANGLE PhiS IN DEGREES)
      J IS THE DAY NUMBER IN THE YEAR
 
     MEAN SOLAR TIME (HEURE DECIMALE)
*/
      mst=gmt+(flon)/15.;
      tet=2.*M_PI*(double)jday/365.;

/*    TIME EQUATION (IN MN.DEC) */
      et=A1+A2*cos(tet)-A3*sin(tet)-A4*cos(2.*tet)-A5*sin(2.*tet);
      et=et*12.*60./M_PI;

/*     TRUE SOLAR TIME */

      tst=mst+et/60.;
      tst=(tst-12.);

/*     HOUR ANGLE */

      ha=tst*15.*M_PI/180.;

/*     SOLAR DECLINATION   (IN RADIAN) */

      delta=B1-B2*cos(tet)+B3*sin(tet)-B4*cos(2.*tet)+B5*sin(2.*tet)-
            B6*cos(3.*tet)+B7*sin(3.*tet);

/*    ELEVATION,AZIMUTH */

      amuzero=sin(dlat)*sin(delta)+cos(dlat)*cos(delta)*cos(ha);
      elev=asin(amuzero);
      az=cos(delta)*sin(ha)/cos(elev);
      if (az<-1.)
       az=-1;
      else if (az>1.)
       az=1.;
      caz=(-cos(dlat)*sin(delta)+sin(dlat)*cos(delta)*cos(ha))/cos(elev);
      azim=asin(az);
      if (caz < 0.) 
       azim=M_PI-azim;
      if ((caz > 0.) && (az < 0.)) 
       azim=2*M_PI+azim;
      azim=azim+M_PI;
      if (azim > (2.*M_PI)) 
       azim=azim-2.*M_PI;
      elev=elev*180./M_PI;

/*     CONVERSION IN DEGREES */

      *ts=90.-elev;
      *fs=azim*180./M_PI;
      return;
}

/*
 * make sure the SIXS application exists in the path supplied.
 */
#define MAX_SIXS_PATH_LEN 2048
#define SIXS_APP   "sixsV1.0B"
static char sixs_path[MAX_SIXS_PATH_LEN] = "";
void set_sixs_path_from(const char *path)
{
  struct stat stbuf;
  char *cptr = (char *)NULL;
  
  /*
   * if passed in arg is path to a file, get the directory path
   */
  strncpy(sixs_path, path, MAX_SIXS_PATH_LEN);
  /*
   * if path ends in /, remove
   */
  for (cptr = &(sixs_path[strlen(sixs_path)-1]); *cptr == '/' && cptr >= &sixs_path[0]; --cptr)
  {
    *cptr = (char) 0;
  }
  
  if (stat(sixs_path, &stbuf) != 0)    /* make sure the path exists */
  {
    /*    fprintf(stderr, "find_file: can't stat %s\n", sixs_path);
	  exit(EXIT_FAILURE);*/
  }
  
  /* if it's a file, find the end of the directory */
  cptr = (char *)NULL;
  if (! S_ISDIR(stbuf.st_mode))   /* not a directory */
  {
    cptr = rindex(sixs_path, '/');
    if (cptr == NULL)            /* current directory */
      sixs_path[0] = (char) 0;
  }
  
   
  /*
   * create path to sixs application file
   */
  if (cptr)
  {
    strcpy((cptr+1), SIXS_APP);
  }
  else
  {
    strcat(sixs_path, SIXS_APP);
  }
   
  /*
   * and make sure it exists
   */
  if (stat(sixs_path, &stbuf) != 0)
  {
    /* fprintf(stderr, "find_file: can't stat %s\n", sixs_path);
       exit(EXIT_FAILURE);*/
  }
}

char *get_sixs_path()
{
  if (strlen(sixs_path) < 1)
  {
    printf("%s\n%s\n%s\n",
      "NO SIXS_PATH defined!",
      "  YOU FORGOT TO CALL set_sixs_path_from",
      "  Cannot continue"
    );
    exit(-1);
  }
  
  return &sixs_path[0];
}

