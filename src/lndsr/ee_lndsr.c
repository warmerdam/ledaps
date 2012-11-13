/**
 * lndsr_ee.c - EarthEngine callable "tile" callable version of lndsr.c main.
 *
 * Closely derived from original lndsr.c and related code.
 *
 * Author: Frank Warmerdam <warmerdam@google.com>
 */

#include "ee_lndsr.h"
#include "myproj.h"

atmos_t atmos_coef;

#define AERO_NB_BANDS 3

/************************************************************************/
/*                           ee_lndsr_main()                            */
/************************************************************************/

/**
 * This function is intended to be called from EarthEngine and should
 * only take arguments easily provided by EarthEngine via a simple SWIG
 * wrapped interface.
 * 
 * @param ar_size_s aerosol array width in pixels
 * @param ar_size_l aerosol array height in lines
 * @param aerosol aerosol values (corresponds to line_ar in lndsr)
 * @param anc_SP SP (slope) on aerosol grid x 4 times of day
 * @param anc_WV water vapour on aerosol grid x 4 times of day
 * @param anc_ATEMP air temperature on aerosol grid x 4 times of day
 * @param anc_O3 ozone on aerosol grid 
 * @param anc_dem elevation in meters on aerosol grid, -9999 for nodata.
 * 
 * @return zero on success, non-zero on failure.
 */
int ee_lndsr_main(
    const char *metadata,
    const char *sixs_filename,

    int full_input_size_s, int full_input_size_l,
    int input_offset_s, int input_offset_l,
    int input_size_s, int input_size_l, int nband,
    short *band1, short *band2, short *band3, 
    short *band4, short *band5, short *band7, 
    const short *thermal, const int8* lndcal_QA,
    const float *anc_ATEMP,

    int utm_zone, 
    double ul_x, double ul_y, double pixel_size,

    int ar_size_s, int ar_size_l,
    const short *aerosol,
    const float *anc_SP,
    const float *anc_WV,
    const float *anc_O3,
    const float *anc_dem)

{
    Ar_gridcell_t *ar_gridcell;
    Lut_t *lut;
    Input_meta_t meta;
    Space_def_t space_def;
    Space_t *space = NULL;
    TileDef_t tile_def;
    sixs_tables_t sixs_tables;
    short *landsat_band[6];

    landsat_band[0] = band1;
    landsat_band[1] = band2;
    landsat_band[2] = band3;
    landsat_band[3] = band4;
    landsat_band[4] = band5;
    landsat_band[5] = band7;

    tile_def.full_size.s = full_input_size_s;
    tile_def.full_size.l = full_input_size_l;
    tile_def.offset.s = input_offset_s;
    tile_def.offset.l = input_offset_l;
    tile_def.size.s = input_size_s;
    tile_def.size.l = input_size_l;

    ee_GetInputMeta( &meta, metadata);

    lut = GetLut( "dummy_param_file", 6, &meta, &(tile_def.size) );

    ee_PrepareSpaceDef( &space_def, utm_zone, ul_x, ul_y, pixel_size,
                        &(tile_def.size) );
    space = SetupSpace( &space_def );

/* -------------------------------------------------------------------- */
/*      Compute scene center.                                           */
/* -------------------------------------------------------------------- */
    Img_coord_float_t img;
    Geo_coord_t geo;
    float center_lon, center_lat;
    
    img.l = tile_def.full_size.l / 2.0 - tile_def.offset.l;
    img.s = tile_def.full_size.s / 2.0 - tile_def.offset.s;
    img.is_fill=false;

    if (!FromSpace(space, &img, &geo))
    	RETURN_ERROR("mapping from space (0)", "ee_build_ar_gridcell", 1);

    center_lat=geo.lat * DEG;
    center_lon=geo.lon * DEG;

/* -------------------------------------------------------------------- */
/*      Compute Scene GMT time based on the scene center location.      */
/* -------------------------------------------------------------------- */
    float scene_gmt;

    if ( (meta.acq_date.hour !=0 ) || (meta.acq_date.minute != 0 ) || (meta.acq_date.second !=0))
   
        scene_gmt=meta.acq_date.hour+meta.acq_date.minute/60.+meta.acq_date.second/3600.;
    else
        scene_gmt=10.5-center_lon/15.;

    if ( scene_gmt < 0.) scene_gmt=scene_gmt+24.;

/* -------------------------------------------------------------------- */
/*      For ancillary layers at four times of day, we figure out        */
/*      which two layers the scene's time of day is between and a       */
/*      coefficient indicating how much of each to use.                 */
/* -------------------------------------------------------------------- */
    double time_coef;
    int pixels_ls = tile_def.size.s * tile_def.size.l;
    int time_index = (int)(scene_gmt/ANC_TIMERES);
    if (time_index >= (ANC_NBLAYERS-1))
        time_index=ANC_NBLAYERS-2;
    
    time_coef= (double)(scene_gmt-time_index*ANC_TIMERES)/ANC_TIMERES;
    
/* -------------------------------------------------------------------- */
/*      Build the ar_gridcell structure with ancillary data sampled     */
/*      on the aerosol grid.                                            */
/* -------------------------------------------------------------------- */
    ar_gridcell = ee_build_ar_gridcell( lut, space, &meta, scene_gmt,
                                        anc_SP, anc_WV, anc_O3, anc_dem );
    if (ar_gridcell == NULL) 
        return 1;

/* -------------------------------------------------------------------- */
/*      Compute deviation from true north.                              */
/* -------------------------------------------------------------------- */
    double delta_x, delta_y, adjust_north;

    printf ("(y0,x0)=(%f,%f)  (lat0,lon0)=(%.16f,%.16f)\n",
            img.l,img.s,(float)(geo.lat * DEG),(float)(geo.lon * DEG));

    delta_y=img.l;
    delta_x=img.s;

    img.l = tile_def.full_size.l / 2.0 - tile_def.offset.l - 100.0;
    img.s = tile_def.full_size.s / 2.0 - tile_def.offset.s;
    img.is_fill=false;
    if (!FromSpace(space, &img, &geo))
    	RETURN_ERROR("mapping from space (0)", "main", 1);

    printf ("(y1,x1)=(%f,%f)  (lat1,lon1)=(%.16f,%.16f)\n",
            img.l,img.s,(float)(geo.lat * DEG),(float)(geo.lon * DEG));

    geo.lon=center_lon*RAD;
    geo.is_fill=false;
    if (!ToSpace(space, &geo, &img))
    	RETURN_ERROR("mapping to space (0)", "main", 1);

    printf ("(y2,x2)=(%f,%f)  (lat2,lon2)=(%.16f,%.16f)\n",
            img.l,img.s,(float)(geo.lat * DEG),(float)(geo.lon * DEG));

    delta_y = delta_y - img.l;
    delta_x = img.s - delta_x;

    printf( "delta_y=%g, delta_x=%g\n", delta_y, delta_x );
    adjust_north=(float)(atan(delta_x/delta_y)*DEG);

    printf("True North adjustment = %f\n",adjust_north);

/* -------------------------------------------------------------------- */
/*      Read 6S from data file.  TODO(warmerdam): allow passing         */
/*      these coefficients in as arguments.                             */
/* -------------------------------------------------------------------- */
    if (read_6S_results_from_file(sixs_filename,&sixs_tables) != 0) {
        return 1;
    }

    switch (meta.inst) {
      case INST_TM:
        sixs_tables.Inst=SIXS_INST_TM;
        break;
      case INST_ETM:
        sixs_tables.Inst=SIXS_INST_ETM;
        break;
      default:
        RETURN_ERROR("Unknown Instrument", "main", 1);
    }

/* -------------------------------------------------------------------- */
/*      Prepare line_ar array representation of aerosol values.         */
/* -------------------------------------------------------------------- */
    int ***line_ar = ee_PrepareLineAr( ar_gridcell, aerosol );

/* -------------------------------------------------------------------- */
/*      Compute atmospheric coefs for the whole scene with              */
/*      aot550=0.01 for use in internal cloud screening : NAZMI         */
/* -------------------------------------------------------------------- */
    int nbpts=lut->ar_size.l*lut->ar_size.s;

    if (allocate_mem_atmos_coeff(nbpts,&atmos_coef))
        ERROR("Allocating memory for atmos_coef", "main");

    printf("Compute Atmos Params with aot550=0.01\n"); fflush(stdout);
    update_atmos_coefs(&atmos_coef, ar_gridcell, &sixs_tables, line_ar, lut, nband, 1);
    
    report_timer( "AR Computation Complete" );

/* -------------------------------------------------------------------- */
/*      Prepare working buffers for landsat resolution data.            */
/* -------------------------------------------------------------------- */
    int ***line_in = ee_PrepareLineIn( &tile_def, lut );
    int **b6_line, *b6_line_buf;
    int8 **qa_line, *qa_line_buf;
    int il, is, ib;
    
    b6_line = (int**)calloc((size_t)(lut->ar_region_size.l),sizeof(int *));
    if (b6_line == (int**)NULL)ERROR("allocating b6 line", "main");
    b6_line_buf = (int*)calloc((size_t)(tile_def.size.s * lut->ar_region_size.l),sizeof(int));
    if (b6_line_buf == (int*)NULL)ERROR("allocating b6 line buffer", "main");
    for (il = 0; il < lut->ar_region_size.l; il++) {
        b6_line[il]=b6_line_buf;
        b6_line_buf += tile_def.size.s;
    }
    
    qa_line = (int8**)calloc((size_t)(lut->ar_region_size.l),sizeof(int8 *));
    if (qa_line == (int8**)NULL)ERROR("allocating qa line", "main");
    qa_line_buf = (int8*)calloc((size_t)(tile_def.size.s*lut->ar_region_size.l),sizeof(int8));
    if (qa_line_buf == (char*)NULL)ERROR("allocating qa line buffer", "main");
    for (il = 0; il < lut->ar_region_size.l; il++) {
        qa_line[il]=qa_line_buf;
        qa_line_buf += tile_def.size.s;
    }

/* -------------------------------------------------------------------- */
/*      Read input first time and compute clear pixels stats for        */
/*      internal cloud screening                                        */
/* -------------------------------------------------------------------- */
    cld_diags_t cld_diags;

/* allocate memory for cld_diags structure and clear sum and nb of obs */	
    if (allocate_cld_diags(&cld_diags,CLDDIAGS_CELLHEIGHT_5KM, CLDDIAGS_CELLWIDTH_5KM, 
                           tile_def.size.l,tile_def.size.s)) {
        ERROR("couldn't allocate memory from cld_diags","main");
    }

    for (il = 0; il < tile_def.size.l; il++) {

        // TODO(warmerdam): this implicitly takes time-of-day=0.0 - FIX
        float *atemp_line = (float *) anc_ATEMP + il * tile_def.size.s;

        /* Read each input band */
        for (ib = 0; ib < nband; ib++) {
            for (is = 0; is < tile_def.size.s; is++) {
                line_in[0][ib][is] = landsat_band[ib][il*tile_def.size.s+is];
            }
        }
        
        for (is = 0; is < tile_def.size.s; is++) {
            b6_line[0][is] = thermal[il*tile_def.size.s+is];
            qa_line[0][is] = lndcal_QA[il*tile_def.size.s+is];
        }

        /* Run Cld Screening Pass1 and compute stats */
        if (!cloud_detection_pass1(lut, tile_def.size.s, il, line_in[0], 
                                   qa_line[0], b6_line[0], atemp_line,
                                   &cld_diags)) {
            ERROR("running cloud detection pass 1", "main");
        }
    }

    report_cld_diags(&cld_diags);

    report_timer( "Cloud Detection Pass 1 Complete" );

/* -------------------------------------------------------------------- */
/*      Do some more work with the cloud diag(nostics?) setting air     */
/*      temperature and some other values.                              */
/* -------------------------------------------------------------------- */
    for (il=0;il<cld_diags.nbrows;il++) {
        int is_ls, il_ls;

        il_ls = (int) (il*cld_diags.cellheight+cld_diags.cellheight/2.);
        if (il_ls >= tile_def.size.l )
            il_ls = tile_def.size.l-1;


        for (is=0;is<cld_diags.nbcols;is++) {

            img.s=is*cld_diags.cellwidth+cld_diags.cellwidth/2.;
            is_ls = (int) (is*cld_diags.cellwidth+cld_diags.cellwidth/2.);
            if (is_ls >= tile_def.size.s )
                is_ls = tile_def.size.s-1;

            /*
             * Unlike the original, we will just pull the ATEMP value
             * from the landsat resolution layer. 
             * TODO(warmerdam): Add time-of-day interpolation.
             */

            cld_diags.airtemp_2m[il][is]=
                (1.0-time_coef)*anc_ATEMP[il_ls * tile_def.size.s + is_ls
                                          + time_index * pixels_ls]
                + time_coef*anc_ATEMP[il_ls * tile_def.size.s + is_ls
                                      + (time_index+1) * pixels_ls];

            if (cld_diags.nb_t6_clear[il][is] > 0) {
                float sum_value,sumsq_value;

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
    report_cld_diags(&cld_diags);

    fill_cld_diags(&cld_diags);

    report_cld_diags(&cld_diags);

    report_timer( "Thermal Cloud Diag Complete" );

    return 0;
}

/************************************************************************/
/*                         ee_PrepareSpaceDef()                         */
/************************************************************************/
int ee_PrepareSpaceDef( Space_def_t *space_def, int utm_zone,
                        double ul_x, double ul_y, double pixel_size,
                        Img_coord_int_t *input_size ) {

    memset( space_def, 0, sizeof(Space_def_t));

    space_def->isin_type = SPACE_NOT_ISIN;
    space_def->orientation_angle = 0.0;  /* should get from metadata OrientationAngle? */

    space_def->proj_num = PROJ_UTM;
    space_def->zone = utm_zone;
    space_def->zone_set = true;

    space_def->sphere = 12; /* ? */

    space_def->pixel_size = pixel_size;
    space_def->ul_corner.x = ul_x; /* TODO(warmerdam): center/corner check? */
    space_def->ul_corner.y = ul_y;
    space_def->ul_corner_set = true;

    space_def->img_size.s = input_size->s;
    space_def->img_size.l = input_size->l;

    return 0;
}

/************************************************************************/
/*                          ee_PrepareLineAr()                          */
/*                                                                      */
/*      This produces an "int" line by line representation of the       */
/*      aerosol values which is compatiable with what the other         */
/*      lndsr code expects to receive.                                  */
/************************************************************************/

int ***ee_PrepareLineAr( Ar_gridcell_t *ar_gridcell, const short *aerosol ) {

    int ***line_ar;
    Img_coord_int_t ar_size;
    int il, ib, is;
    int *line_ar_buf;
    int **line_ar_band_buf;

    ar_size.l = ar_gridcell->nbrows;
    ar_size.s = ar_gridcell->nbcols;

    line_ar = (int ***)calloc((size_t)ar_size.l, sizeof(int **));
    if (line_ar == (int ***)NULL) 
        ERROR("allocating aerosol line buffer (a)", "ee_PrepareLineAr");

    line_ar_band_buf = (int **)calloc((size_t)(ar_size.l * AERO_NB_BANDS), sizeof(int *));
    if (line_ar_band_buf == (int **)NULL) 
        ERROR("allocating aerosol line buffer (b)", "ee_PrepareLineAr");

    line_ar_buf = (int *)calloc((size_t)(ar_size.l * ar_size.s * AERO_NB_BANDS), 
                                sizeof(int));
    if (line_ar_buf == (int *)NULL) 
        ERROR("allocating aerosol line buffer (c)", "ee_PrepareLineAr");

    for (il = 0; il < ar_size.l; il++) {
        line_ar[il] = line_ar_band_buf;
        line_ar_band_buf += AERO_NB_BANDS;
        for (ib = 0; ib < AERO_NB_BANDS; ib++) {
            line_ar[il][ib] = line_ar_buf;
            line_ar_buf += ar_size.s;
        }
    }

    return line_ar;
}

/************************************************************************/
/*                          ee_PrepareLineIn()                          */
/*                                                                      */
/*      Prepare the "line_in" buffer which is per-band, per-line        */
/*      pointers to an underlying int buffer but with only enough       */
/*      lines to hold one swath of aerosol data at a time.              */
/************************************************************************/
int ***ee_PrepareLineIn( TileDef_t *tile_def, Lut_t *lut ) {

    int ***line_in;
    int il, ib, is;
    int *line_in_buf;
    int **line_in_band_buf;
    Img_coord_int_t *size = &(tile_def->size);
    const int bands = 5;

    line_in = (int ***)calloc((size_t)size->l, sizeof(int **));
    if (line_in == (int ***)NULL) 
        ERROR("allocating line_in buffer (a)", "ee_PrepareLineIn");

    line_in_band_buf = (int **)calloc((size_t)(size->l * bands), sizeof(int *));
    if (line_in_band_buf == (int **)NULL) 
        ERROR("allocating line_in buffer (b)", "ee_PrepareLineIn");

    line_in_buf = (int *)calloc((size_t)(size->l * size->s * bands), 
                                sizeof(int));
    if (line_in_buf == (int *)NULL) 
        ERROR("allocating line_in buffer (c)", "ee_PrepareLineIn");

    for (il = 0; il < size->l; il++) {
        line_in[il] = line_in_band_buf;
        line_in_band_buf += bands;
        for (ib = 0; ib < bands; ib++) {
            line_in[il][ib] = line_in_buf;
            line_in_buf += size->s;
        }
    }

    return line_in;
}

