/* 
 
 
 Compile with this:
 
 cc -O0 -o bin/SDSreader3.0 src/SDSreader3.0.c \
    -I/home/jim/bin/ERICV/include -I$HDFINC \
    -L$HDFLIB -lmfhdf -ldf -ljpeg -lz -lm


*/
/* Program used to read values from SDSs according to some criteria; 
   Jim Ray, 20-SEP-01. 
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include "mfhdf.h"
/*#include "struct.h"
#include "globals.h"
#include "proj_globals.h"*/
#ifndef M_PI
#define M_PI 3.141592653589793
#endif
#define ROVRDG (M_PI/180.)
#define MAXDIMS 2
#define MAXCRITERIA 10
#define MAXLENGTH 256       /* maximum length of strings */
#define MAXSITES 50
#define MAXBANDS 5120
#define FILL_VALUE -99999
#define EPSLN 1E-3
/*#define DEBUG*/
#define ABS(F) (((F) < 0)? -(F):(F))

/* These are mostly stuff read from the command-line arguments */
short USE_ALL_SDSs = 0;
char SDSlistfile[MAXLENGTH];
char HDFfile[MAXLENGTH];
char SDSnames[MAXBANDS][MAXLENGTH];
int32 SDS_offset[MAXBANDS];                /* for 3-D SDSs */
char SDStofilter[MAXCRITERIA][MAXLENGTH];
short really_boolean[MAXCRITERIA];
double lowend[MAXCRITERIA], highend[MAXCRITERIA];
int n_SDSs ;
short except, proj_ok, n_criteria;
short ULr, ULc, nr, nc, n_sites ;
float overallULR, overallULC, overallNR, overallNC;
short overall_set = 0;
short VERBOSE = 0;
short COLUMNS = 0;
short samples = 0;
char outfile[MAXLENGTH];
char HDFoutfile[MAXLENGTH];
char binfile[MAXLENGTH];
int ref_SDS;
char ref_string[MAXLENGTH];
int aggreg, agg_dim1, agg_dim2, agg_thresh;
float ULR[MAXSITES], ULC[MAXSITES];
float NR[MAXSITES], NC[MAXSITES];

int main(int argc, char **argv)
{
FILE *fdSDSlistfile;
FILE *fd, *fdout;
FILE *fdExtsrc[MAXCRITERIA];
FILE *fddebug;
short ret, i;
int32 which_SDS[MAXBANDS];
short which_filter[MAXCRITERIA];
intn ref;
char *attr;
int32 n_type;
int32 j, ii, jj, file_id, out_id, sds_outid, cr_sds_id, sds_id, n_gattr, n_sets, nn_sets, count, rank;
char label[MAXLENGTH], unit[MAXLENGTH], format[MAXLENGTH];
char coordsys[MAXLENGTH], attrib[MAXLENGTH];
float64 cal[MAXBANDS], cale, off[MAXBANDS], offe;
int32 dims[10], outdims[2], outstart[2], xdims[MAXDIMS], SDS_dims[MAXDIMS], number_type[MAXBANDS], nt, nattr, xsds;
int32 w_number_type[MAXBANDS];
int32 filter_type;
int32 SDS_type[MAXBANDS];
int32 SDS_filtertype[MAXCRITERIA];
int32 index_in_file[MAXBANDS];
long n_rows, n_cols, n_to_skip, save_bin, save_out, save_HDFout;
char name[MAXLENGTH];
char SDSinfile[MAXBANDS][MAXLENGTH];
char sds_table[MAXLENGTH*MAXBANDS];
char resultstring[MAXLENGTH*MAXBANDS];
char debugfile[MAXLENGTH];
int missing, longest, exceeded;
float tmpflt;
long k,li,lj,n_kept;
double **SDS_image_D;
double **SDS_filterimage_D;
double **WHOLESDS_image_D;
double **WHOLESDS_filterimage_D;
short *testimage[MAXBANDS];
float *f_image[MAXBANDS];
double SDS_sum[MAXBANDS], SDS_avg[MAXBANDS], SDS_SD[MAXBANDS];
double tmpdb;
int16 *keepmask;
int warned;
int32 bytes_to_HDFtype(short);
short process_arg(int argc, char **argv, int n_lines);
double *get_images(int32 sd_id, char *name, int32 offset, int32 SDS_index, short SDS_type, FILE *fd, 
                   short, long, short, short, short, short, float subsample_factor);
void sun_angles(short jday, float gmt, double dlat, double dlon, float *ts, float *fs);
void sat_angles(double dlat, double dlon, double dheight, double dlatsub,
                double dlonsub, float *tv, float *fv);
void metareader(int32 sd_id, char *type_of_meta, char *metastring, 
                int32 *count, char *data);
void process_global_attributes(int32 sd_id1, int32 sd_id_out);
void process_local_attributes(int32 sd_id1, int32 sd_id_out);
void process_dim_ids(int32 sd_id1, int32 sd_id_out);
void aggregate(int n_SDSs, int16 *keepmask, double *SDS_image_D[MAXBANDS], short *nr, short *nc);

short jday;
float gmt;
double dlat, dlon, dlatsub, dlonsub, dheight;
float ts, fs, tv, fv;

/* 29-DEC-98: Variables for getting SDS from another HDF as criteria-data */
int32 crit_file_id[MAXCRITERIA];
char otherHDFfile[MAXLENGTH];
float subsample_factor[MAXCRITERIA];
short here, match, skipper;
short take_it;
short go;
int irow, icol;

/* 20-SEP-01 */
int latid, lonid;
char ALLCAPS[70];
char localid[MAXLENGTH] = {'\0'};
int32 n_val;
char tmpstring[MAXLENGTH] = {'\0'};
char ESDT[10];
char JDAY[10];
char GMT[10];
int get_angles;

int16 *i16_data;
int8 *i8_data;
uint16 *ui16_data;
uint8 *ui8_data;
uint32 *ui32_data;
int32 *i32_data;
double *dbl_data;
int n_lines, paramfileline, *newfileindex;
int get_n_lines(char *possibleparamfile, int **newfileindex);
int32 lin, col;
HDF_CHUNK_DEF cdef;
long hcrit_value, lcrit_value;
int truly_boolean, forget_it;
short overall_rows, overall_cols, overall_row0, overall_col0;

/* for now... */
dheight=(705.)*1000.;
jday=250;   /* default */
gmt=12.36;  /* default */
overallULR = overallULC = overallNR = overallNC = -999.0F;




/***************************
 ** Take care of arguments *
 ***************************/ 
/* new, 22-FEB-07: if (argc == 2), we may be using a parameter file;
   Determine ahead of time how many lines it has */
   
/* newfileindex is going to contain either 0 or 1, 1 meaning 
   the HDF input file at this argument has never been seen before,
   0 meaning it has.  newfileindex is meant to limit the reads
   of the input files.  I will set it to NULL in case of 
   memory errors in reading entire files later.  For that reason,
   I have to make sure it is non-NULL at the outset */
   
if (( newfileindex = calloc(1, sizeof(int))) == (int *)NULL) {
    exit(-2);
   }
newfileindex[0] = 1;  /* just in case */

if (argc == 2) n_lines = get_n_lines(argv[1], &newfileindex);
else           n_lines = 1;


if (overall_set == 1) {
   /*printf("XXXXX: oULr = %f, oULc = %f, onr = %f, onc = %f\n", overallULR, overallULC, overallNR, overallNC); */
   overall_row0 = (short)overallULR;
   overall_col0 = (short)overallULC;
   overall_rows = (short)overallNR - overallULR;
   overall_cols = (short)overallNC - overallULC;
   }
   
/*for (paramfileline = 0; paramfileline < n_lines; paramfileline++)  { 
printf("%d\n", newfileindex[paramfileline]);
}
exit(0);*/

WHOLESDS_filterimage_D = (double **)malloc(MAXCRITERIA*sizeof(double *));
SDS_filterimage_D = (double **)malloc(MAXCRITERIA*sizeof(double *));
WHOLESDS_image_D = (double **)malloc(MAXBANDS*sizeof(double *));
SDS_image_D = (double **)malloc(MAXBANDS*sizeof(double *));

for (j=0;j<MAXCRITERIA;j++) {
    WHOLESDS_filterimage_D[j] = (double *)NULL;
    SDS_filterimage_D[j] = (double *)NULL;
    }
for (j=0;j<MAXBANDS;j++)    {
    SDS_image_D[j] = (double *)NULL; 
    WHOLESDS_image_D[j] = (double *)NULL; 
    }       

for (paramfileline = 0; paramfileline < n_lines; paramfileline++)  { /* 22-FEB-07 */ 

ULr = ULc = nr = nc = -999;

proj_ok = warned = exceeded = 0;	
n_SDSs = 0;
except = 0;
n_sites = 0;
n_criteria = aggreg = agg_thresh = 0;
skipper = 1; get_angles = 0;
save_bin = save_out = save_HDFout = ref_SDS = 0;

ret = process_arg(argc, argv, paramfileline);
/*printf("RET=%d\n", ret);*/
if (ret<1) {
     if (ret == 0) printf("\nProgram used to read values from SDSs\n");
     printf("Usage: %s -f <HDF file> -l <SDS namelist> -x -w <criteria> -v -c -o <file> -s -b <prefix> -p <file>\n", argv[0]);
     printf("where -l <SDS namelist> is a file with names of SDSs, or a list of SDSs in HDFfile;\n");
     printf("                        default is to use all SDSs in HDFfile.  Max number: %d\n", MAXBANDS);
     printf("      -x                use all SDSs in file EXCEPT those specifed in <SDS namelist>;\n");
     printf("                        (has no effect if <SDS namelist> is not specified).\n");
     printf("      -v                verbose (diagnostics are printed to screen)\n");
     printf("      -c                results are written in columns instead of extcal-like lines (default)\n");
     printf("      -s                samples are written to stdio\n");
     printf("      -o <file>         output is printed to <file>\n");
     printf("      -h <file>         samples are written to <HDF file>\n");
     printf("      -b <prefix>       samples are written as flat binaries with following \n");
     printf("                          naming convention: <prefix>_<SDS name>_sample_<nn>, \n");
     printf("                        where nn is the site number (up to %d sites can be chosen). \n", 
               MAXSITES);
     printf("      -w <criteria>     is a set of criteria for choosing SDS values. Example:\n");
     printf("                          -w \"50 100 4 4\" \"ndvi LT 1000.4\", \n");
     printf("                        all values within a 4 by 4 box with UL corner at row 50, col 100,\n");
     printf("                        and where NDVI values are below 1000.4.  Criteria are 'and'ed\n");
     printf("                        together.  Each criterion must be enclosed in quotes -- \"...\".  At \n");
     printf("                        least one area criterion (site) must be present, with this format:\n"); 
     printf("                          <UL row> <UL col> <n rows> <n cols> or\n");
     printf("                          <UL lat> <UL lon> <LR lat> <LR lon>;\n");
     printf("                        the first format is the default, the second is used if option -p\n");
     printf("                        (below) is specified.  Up to %d optional boolean criteria can also \n",
               MAXCRITERIA);
     printf("                        be entered, with this format:\n");
     printf("                          <source> <LT|GT|EQ|BT|AND|NOT> <value [2nd val]>\n");
     printf("                        where LT, GT and EQ take one value and BT takes <low> <high> values.\n");
     printf("                        AND and NOT takes <bitmapped> value(s) (2 raised to the target-bit(s)).\n");
     printf("                        \n");
     printf("                        e. g., for AVHRR QAs in which a '1' in bit 4 means sunglint:\n");
     printf("                        to set all glinty pixels to zero, \"<source> NOT 16\"\n");
     printf("                        if a '1' in bit 1 means clouds and '1' in bit 0 means partly cloudy:\n");
     printf("                        to set all glinty and cloudy pixels to zero, \"<source> NOT 19\"\n");
     printf("                        \n");
     printf("                        e. g., for MOD09 State QA in which a '1' in bit 2 means cloud shadow:\n");
     printf("                        to set all cldshd pixels to zero, \"<source> NOT 4\"\n");
     printf("                        to set all BUT the cldshd pixels to zero, \"<source> AND 4\"\n");
     printf("                        \n");
     printf("                        For multibit QAs, supply two arguments; a bitmask and a bitmasked value\n");
     printf("                        e. g., for MOD09 State QA in which bits 6 and 7 encode aerosol quality:\n");
     printf("                        to set all average-aerosol pixels to zero, \"<source> NOT 192 128\"\n");
     printf("                        to set all low-aerosol pixels to zero,     \"<source> NOT 192 64\"\n");
     printf("                        to set all high-aerosol pixels to zero,    \"<source> NOT 192 192\"  (or just \"<source> NOT 192\")\n");
     printf("                        to set all climatological-aerosol pixels to zero, \"<source> NOT 192 0\"\n");
     printf("                        \n");
     printf("                        Argument <source> may be a name of an SDS in the file, (e. g., NDVI), \n");
     printf("                        an external flat binary file name, or a name of an SDS in another\n");
     printf("                        HDF file (HDFfile path:SDSname)\n");
/*     printf("      -p <file>         get projection parameters from <file>.  With -p, user must\n");
     printf("                        enter <UL_lat UL_lon LR_lat LR_lon> instead of row/col values (the\n");
     printf("                        default).\n"); */
     printf("      -r <string>       resolution to process (\"250\", \"500\", \"1km\", etc. Resolution is of first SDS in file by default)\n");
#if defined (DEBUG)     
     printf("                        This version compiled with -DDEBUG, will put out 'keepmask'\n");
#endif
     printf("                        3D SDSs can be processed, but 1) the third dimension cannot be specified in the criteria, \n");
     printf("                        and 2) criteria SDSs (if any) must be 2D.\n");
     printf("      -a <l,c{,t}>      aggregation by the values of <l,c>, where l is for lines, c for columns (required).\n");
     printf("                        optional t is threshold for keep/reject.\n\n");
     printf("New to this version: GOT RID OF THE PROJECTION LIBRARY and PARAMETER FILE HANDLING \n\n");
     printf("\nJim Ray, SSAI, %s, version 3.0.\n\n",__DATE__);
     if (ret == 0) exit(0);
     else          exit(-1);
            }

if (overall_set == 0) {  /* then use just the "-w <whatever>" values as "overall" values */
   printf("XXXXX: oULr = %f, oULc = %f, onr = %f, onc = %f\n", overallULR, overallULC, overallNR, overallNC); 
   overall_row0 = (short)overallULR;
   overall_col0 = (short)overallULC;
   overall_rows = (short)overallNR - overallULR;
   overall_cols = (short)overallNC - overallULC;
   }

	    
/* Echo arguments to user
 **************************/
if (VERBOSE) {
   printf("     HDF file to read = %s\n", HDFfile);
   if (USE_ALL_SDSs) printf("     USING ALL SDSs ");
   if (n_SDSs) {
      if (except) printf("EXCEPT\n");
      else        printf("     Using these SDSs:\n");
        }
   else 
      printf("\n");
   for (i=0;i<n_SDSs;i++) 
      printf("       %s\n", SDSnames[i]);
   if (ref_SDS == -1) printf("     Reference resolution string = %s\n", ref_string);
   printf("     Using these site(s):\n");
   if (proj_ok) 
       for (i=0;i<n_sites;i++) 
           printf("       Upper-left lat, lon; (%f %f); lower-right lat, lon; (%f %f)\n",
                 ULR[i], ULC[i], NR[i], NC[i]);
   else
       for (i=0;i<n_sites;i++) 
           printf("       Upper-left row, column; (%d %d); n rows, columns; (%d %d)\n", 
                 (short)ULR[i], (short)ULC[i], (short)NR[i], (short)NC[i]);
   if (n_criteria) printf("     Using these Boolean criteria:\n");
   for (i=0;i<n_criteria;i++) 
      printf("       %s, low=%f, high=%f\n", SDStofilter[i], lowend[i], highend[i]);
   if (aggreg) printf("     Aggregating %d by %d pixels, %d pixel threshold\n", agg_dim1, agg_dim2, agg_thresh);
}

if (strcmp(binfile,"")) save_bin = 1;
if (strcmp(outfile,"")) {
   save_out = 1;
   if (  (fdout =  (fopen(outfile, "wb"))) == NULL ) {
     printf("Error: cannot open '%s'\n", outfile);
     exit(-7);
                                                     } 
                         }
			 
/* new, 31-OCT-02: set up HDF outputs if necessary */
if (strcmp(HDFoutfile,"")) {
   save_HDFout = 1;
   /*if (  (fdout =  (fopen(HDFoutfile, "rb"))) == NULL ) { */ 
     /* file doesn't exist, create it */
     if (  (out_id = SDstart(HDFoutfile, DFACC_CREATE)) == -1) {
          printf("Error: cannot run SDstart on '%s'; check file\n", HDFoutfile);
          exit(-2);
                                                               }
                                                        /*} */
   /*else {  
     fclose(fdout);
     / * file exists, open for read * /
     if (  (out_id = SDstart(HDFoutfile, DFACC_RDWR)) == -1) {
          printf("Error: cannot run SDstart on '%s'; check file\n", HDFoutfile);
          exit(-2);
                                                             }
        } */

			   }
			 
			 		
/* Make sure HDF file exists and actually is an HDF file
 **********************************************************/
if (  (fd =  (fopen(HDFfile, "rb"))) == NULL ) {
     printf("Error: cannot read '%s'; check file name \n", HDFfile);
     exit(-2);
                                               } 
else fclose(fd);

if (  (file_id = SDstart(HDFfile, DFACC_RDONLY)) == -1) {
     printf("Error: cannot run SDstart on '%s'; check file\n", HDFfile);
     exit(-2);
                                                        }
SDfileinfo(file_id, &n_sets, &n_gattr);
if (n_sets == 0) {
     printf("Error; no SDSs in '%s'\n", HDFfile);
     exit(-9);
                 }
		 
if (save_HDFout == 1) process_global_attributes(file_id, out_id);


/* New, 20-SEP-01: get JDAY, GMT from file's metadata */
n_val = 0;
metareader(file_id, "COREMETADATA", "LOCALGRANULEID", &n_val, localid);
if (n_val == 0) {
   printf("error: cannot find metadatum LOCALGRANULEID in file '%s' -- will use defaults\n", HDFfile);
                }
else {
   /*printf("LOCALGRANULEID is '%s'\n", localid);*/
   for(i=0;i<8;i++)ESDT[i] = localid[i+1];
   ESDT[i] = '\0';
   if (!strstr(ESDT,"CRS")) get_angles = 0;

   strcpy(tmpstring, strchr(localid, '.'));
   for(i=0;i<3;i++)JDAY[i] = tmpstring[i+6];
   JDAY[i] = '\0';
   for(i=0;i<2;i++)GMT[i] = tmpstring[i+10];
   GMT[i] = '\0';
   jday = (short)(atoi(JDAY));
   gmt = atof(GMT);
   for(i=0;i<2;i++)GMT[i] = tmpstring[i+12];
   GMT[i] = '\0';
   gmt += (atof(GMT))/60.0F;
   
   } 


/* Find all the SDSs to be processed in the HDF file, make
 * sure they are there and all have the same ranks
 * and dimensions (for now).
 **********************************************************/
/* What we need is
    n_SDSs -- the total number of SDSs-to-be-processed; and for each 
    SDSnames[] -- its name 
    which_SDS[] -- its index in the HDF file
    SDS_type[] -- its number type.
  */

take_it = 0;
here = 0;
if (USE_ALL_SDSs) {
   for (i=0;i<n_sets;i++)  {
       sds_id = SDselect(file_id, i);
       SDgetinfo(sds_id, name, &rank, dims, &nt, &nattr);
       SDendaccess(sds_id);
       
       /* If except == 1, you are using all SDSs EXCEPT those specified in the list */
       if (except) {
          go = -1;
	  for (j=0;j<n_SDSs;j++) 
	     if (!strcmp(SDSnames[j], name)) 
	        go = j;

          if (go == -1) 
	     go = 1;
	  else 
	     go = 0;
                   }
       else 
          go = 1;

       if ((go == 1) &&
           ((rank == MAXDIMS)||
	    (rank == MAXDIMS+1))) {
	    
          if (!here) {           /* FIX the n_rows/n_cols values to the those of
	                             the first SDS */                 				     
	       /* new to version 2.6 */
	       if ((ref_SDS == -1)&&(!strstr(name, ref_string))) continue;
	   
 	       if (rank == MAXDIMS) {	         
                 n_rows = dims[0];
                 n_cols = dims[1];
	         here=1;
                 which_SDS[take_it] = (int32)i;
                 SDS_type[take_it] = nt;
                 SDS_offset[take_it] = -1;
                 strcpy(SDSinfile[take_it++], name);
		 }
	       else {
                 n_rows = dims[1];
                 n_cols = dims[2];
	         here=1;
		 for(j=0;j<dims[0];j++) {
                   which_SDS[take_it] = (int32)i;
                   SDS_type[take_it] = nt;
                   SDS_offset[take_it] = (int32)j;
                   /*strcpy(SDSinfile[take_it++], name);*/
	           sprintf(SDSinfile[take_it++], "%s_%03d\0", name, j);
		 }
		 }
	              }
	   else {   /* we've already taken the first SDS and fixed the n_rows/n_cols values;
	               make sure the rest have n_rows and n_cols */
	       if (take_it < MAXBANDS) {
	         if (rank == MAXDIMS) {
	            if ((n_rows == dims[0])&&(n_cols == dims[1])) {
                      which_SDS[take_it] = (int32)i;
                      SDS_type[take_it] = nt;
                      SDS_offset[take_it] = -1;
                      strcpy(SDSinfile[take_it++], name);
	                                                          }
		    }
	         else {
	            if ((n_rows == dims[1])&&(n_cols == dims[2])) {
		      for(j=0;j<dims[0];j++) {
                         which_SDS[take_it] = (int32)i;
                         SDS_type[take_it] = nt;
                         SDS_offset[take_it] = (int32)j;
                         /*strcpy(SDSinfile[take_it++], name);*/
			 sprintf(SDSinfile[take_it++], "%s_%03d\0", name, j);
			 }
	                                                          }
		    }
				       }
	       else {
	         if (!warned) {
		      printf("Warning: file '%s' \nhas more SDSs than this program can handle: will only process %d\n",
		          HDFfile, MAXBANDS);
		      warned = 1;
		              }
	            }
	        }
		            } /* if (go)&&(rank == MAXDIMS) */
       /*SDendaccess(sds_id);*/
                           }  /* for (i=0;i<n_sets;i++) */
   n_SDSs = take_it;	   
   if (n_SDSs == 0) {
         printf("Error; apparently, no SDSs in the file are two or three dimensional.\n");
         exit(-10);
                    }
   for (j=0;j<n_SDSs;j++) strcpy(SDSnames[j], SDSinfile[j]);
   
                  }  /* if (USE_ALL_SDSs) */
		  
else {   /* Find all the SDSs-to-be-processed in the HDF file; 
            make sure they are of rank = MAXDIMS and all have the same dimensions */
    take_it = 0;
    for (j=0;j<n_SDSs;j++) {   /* For each SDS specified at the command line ... */
        for (i=0;i<n_sets;i++)  {
            sds_id = SDselect(file_id, i);
            SDgetinfo(sds_id, name, &rank, dims, &nt, &nattr);
	    if (strcmp(name, SDSnames[j]) == 0) {  /* it's here, */
               if ((rank == MAXDIMS)||(rank == MAXDIMS+1)) { 
                  if (!here) {           /* FIX the n_rows/n_cols values to the those of
	                                    the first SDS */
	             /* new to version 2.6 */
	             if ((ref_SDS == -1)&&(!strstr(SDSnames[j], ref_string))) continue;
	   
                     if (rank == MAXDIMS) {
                        n_rows = dims[0];
                        n_cols = dims[1];
	                here=1;
                        which_SDS[take_it] = (int32)i;
                        SDS_type[take_it] = nt;
                        SDS_offset[take_it] = -1;
		        strcpy(SDSnames[take_it++], SDSnames[j]);
			}
                     else {
                        n_rows = dims[1];
                        n_cols = dims[2];
	                here=1;
		        for(jj=0;jj<dims[0];jj++) {
                           which_SDS[take_it] = (int32)i;
                           SDS_type[take_it] = nt;
                           SDS_offset[take_it] = (int32)jj;
		           sprintf(SDSnames[take_it++], "%s_%03d\0", SDSnames[j], jj);
			   }
			}
	                     }
	          else {   /* we've already taken the first SDS and fixed the n_rows/n_cols values;
	                      make sure the rest have n_rows and n_cols */
                     if (rank == MAXDIMS) {
	                if ((n_rows == dims[0])&&(n_cols == dims[1])) {
                           which_SDS[take_it] = (int32)i;
                           SDS_type[take_it] = nt;
                           SDS_offset[take_it] = -1;
		           strcpy(SDSnames[take_it++], SDSnames[j]);
	                                                              }
			  }
                     else if (rank == MAXDIMS+1) {
	                if ((n_rows == dims[1])&&(n_cols == dims[2])) {
		           for(jj=0;jj<dims[0];jj++) {
                               which_SDS[take_it] = (int32)i;
                               SDS_type[take_it] = nt;
                               SDS_offset[take_it] = (int32)jj;
		               sprintf(SDSnames[take_it++], "%s_%03d\0", SDSnames[j], jj);
			       }
	                                                              }
			  }
		       }
			            }
	                                    }
	    SDendaccess(sds_id);
	                        }  /* for (i=0;i<n_sets;i++) */			    
	                   }  /* for (j=0;j<n_SDSs;j++) */
    n_SDSs = take_it;
    /*n_SDSs = 4;*/
    if (n_SDSs == 0) {
         printf("Error, no SDSs in the list can be found, or none in the list are two dimensional.\n");
         exit(-10);
                     }
    }
   
/*
for (j=0;j<n_SDSs;j++) 
    printf("%d SDSname %s, SDS offset %d, index %d, type = %d\n", j, SDSnames[j], SDS_offset[j], which_SDS[j], SDS_type[j]);
 */
longest = strlen(SDSnames[0]);
for (j=1;j<n_SDSs;j++) 
    if (longest < strlen(SDSnames[j])) 
        longest = strlen(SDSnames[j]);
   
/* Make sure all "SDSs" to be filtered upon are 
 * in the HDF file or are external flat binary files 
 * or are in another HDF file...
 *****************************************************/
for (j=0;j<n_criteria;j++) {
    missing = -1;
    for (ii=0;ii<MAXLENGTH;ii++) otherHDFfile[ii] = '\0';
    crit_file_id[j] = 0;
    which_filter[j] = -1; 
    fdExtsrc[j] = 0;
    subsample_factor[j] = 1.0;

    for (i=0;i<n_sets;i++) {
       sds_id = SDselect(file_id, i);
       SDgetinfo(sds_id, name, &rank, dims, &nt, &nattr);
       if (!strcmp(name, SDStofilter[j])) {
          missing = i;
	  SDS_dims[0] = dims[0];
	  SDS_dims[1] = dims[1];
	  filter_type = nt;
	  }
       SDendaccess(sds_id);
                           }   
    
    if (missing != -1) {  /* yes, SDS to be filtered upon is in HDF file */ 
       crit_file_id[j] = file_id;  
       which_filter[j] = missing; 
       SDS_filtertype[j] = filter_type; 
       
       /* make sure it has the right dimensions, n_rows by n_cols */
       if ((n_rows != SDS_dims[0]) || (n_cols != SDS_dims[1])) {
           if (ref_SDS != -1) {
              printf("Error: SDS '%s' has different dimensions from the rest\n", SDStofilter[j]);
              exit(-3);
	      }
           else {
	      /* the SDS-to-filter is a different dimension than the official dimension.  
	         We will need to interpolate or subsample what is read to generate what we need */
		 subsample_factor[j] = (float)SDS_dims[0]/(float)n_rows;
	      }	      
       }
                       }
    else {
      if (strchr(SDStofilter[j], ':')) {  /* New, 29-DEC-98 */  
             /* if SDStofilter[j] has a ':' character, take it as an SDS in another HDF file */
           strcpy(otherHDFfile, SDStofilter[j]);
	   ii=0;
	   while (otherHDFfile[ii] != ':') ii++;
	   otherHDFfile[ii] = '\0';
           if (  (crit_file_id[j] = SDstart(otherHDFfile, DFACC_RDONLY)) == -1) {
              printf("Error: cannot run SDstart on '%s'; check file\n", otherHDFfile);
              exit(-2);
                                                                                }
	   here = ii+1;
	   while (SDStofilter[j][ii++] != '\0') SDStofilter[j][ii-here] = SDStofilter[j][ii];
	   SDStofilter[j][ii-here] = '\0';

           SDfileinfo(crit_file_id[j], &nn_sets, &n_gattr);
           for (ii=0;ii<nn_sets;ii++) {
               cr_sds_id = SDselect(crit_file_id[j], ii);
               SDgetinfo(cr_sds_id, name, &rank, dims, &nt, &nattr);
               if (!strcmp(name, SDStofilter[j])) {
                  which_filter[j] = ii; 
                  SDS_filtertype[j] = nt; 
		  xdims[0] = dims[0];
		  xdims[1] = dims[1];
	       }
	       SDendaccess(cr_sds_id);
	   }
	   if (which_filter[j] == -1) {
              printf("Error: cannot find SDS '%s' in file '%s'; check file\n", SDStofilter[j], otherHDFfile);
              exit(-2);
	   }
           /* make sure it has the right dimensions, n_rows by n_cols */
           if ((n_rows != xdims[0]) || (n_cols != xdims[1])) {
              if (ref_SDS != -1) {
                 printf("Error: SDS '%s' has different dimensions from the rest\n", SDStofilter[j]);
                 exit(-3);
	         }
              else {
	         /* the SDS-to-filter is a different dimension than the official dimension.  
	            We will need to interpolate or subsample what is read to generate what we need */
		    subsample_factor[j] = (float)xdims[0]/(float)n_rows;
	         }	      
           }
       }   /*  if (strchr(SDStofilter[j], ':'))  */  
      else {  /* try an external flat binary file... */
       if (  (fdExtsrc[j] =  (fopen(SDStofilter[j], "rb"))) == NULL ) {
           printf("Error: SDS '%s' is missing, and it's not an external file either\n", SDStofilter[j]);
           exit(-2);
                                                                      }
       /* make sure it has the right dimensions, a multiple of (n_rows*n_cols) */
       fseek(fdExtsrc[j],0L,SEEK_SET);
       fseek(fdExtsrc[j],0L,SEEK_END);
       k = ftell(fdExtsrc[j]);
       /*printf("File is %ld bytes long\n", k);*/
       fseek(fdExtsrc[j],0L,SEEK_SET);
       
       li = (long)n_rows*(long)n_cols;
       tmpflt = (float)k/(float)li;
       if ( (tmpflt - (float)((short)tmpflt)) < EPSLN) {
           SDS_filtertype[j] = bytes_to_HDFtype((short)tmpflt);
           which_filter[j] = -1; 
          /*printf("%d: File is %ld bytes long, %d bytes per pixel\n", j, k, SDS_filtertype[j]);*/
	   }
       else {
           printf("Error: file '%s' has a size that is not an integral ", SDStofilter[j]);
	   printf("multiple of (%d * %d), cannot process\n", n_rows, n_cols);
	   exit(-5);
            }
         }
	  }
                           }  /* for (j=0;j<n_criteria;j++) */
			   

/*for (j=0;j<n_criteria;j++) {
    printf("%p %p %p SDStofilter %s, index %d, type = %d, low = %f, high = %f\n", 
    file_id, crit_file_id[j], fdExtsrc[j], SDStofilter[j], which_filter[j], SDS_filtertype[j],
    lowend[j], highend[j]);
}*/

/* 
 * New, 17-NOV-98: multiple site processing 
 ********************************************/
for (ii=0;ii<n_sites;ii++) {
   if (proj_ok == 1) { 
/*        dlat = ULR[ii]*ROVRDG;
	dlon = ULC[ii]*ROVRDG;
        ProjFOR(&project, &irow, &icol, dlat, dlon);
	ULr = (short)irow;
	ULc = (short)icol;
        dlat = NR[ii]*ROVRDG;
	dlon = NC[ii]*ROVRDG;
        ProjFOR(&project, &irow, &icol, dlat, dlon);
        nr = irow - ULr;
        nc = icol - ULc;
        if ((nc < 0)||(nr < 0)) {
           printf("\tError: there is a coding error in this program. The number of rows/cols to read\n");
           printf("\tfrom the SDS should NOT be (%d,%d).  See developer of the program.\n", nr,nc);
           exit(0);
                                }
        if (nc == 0) {
           printf("\tWarning: UL lon, LR lon difference is under the deltalon for this part of\n");
           printf("\tthe grid, number of columns = 0; will set to 1.\n");
           nc=1;
                     }
        if (nr == 0) {
           printf("\tWarning: UL lat, LR lat difference is under the deltalat for this part of\n");
           printf("\tthe grid, number of rows = 0; will set to 1.\n");
           nr=1;
                     }
/ *
	printf("UL ******** ROW = %d, COL = %d, LAT = %f, LON = %f ******\n", 
	    ULr, ULc, ULR[ii], ULC[ii]); 
	printf("LR ******** NROW = %d, NCOL = %d, LAT = %f, LON = %f ******\n", 
	    nr, nc, NR[ii], NC[ii]); 
 */
                     }
   else {
        ULr = (short)ULR[ii];
	if (ULr < 0) {
            printf("\tWarning: variable ULr rounded to %d, adjusting ULr to zero\n", ULr);
            ULr = 0;
                     }
        ULc = (short)ULC[ii];
	if (ULc < 0) {
            printf("\tWarning: variable ULc rounded to %d, adjusting ULc to zero\n", ULc);
            ULc = 0;
                     }
        nr = (short)NR[ii];
	if (nr <= 0) {
            printf("\tWarning: variable nr rounded to %d, adjusting nr to 1\n", nr);
            nr = 1;
                    }
        nc = (short)NC[ii];
	if (nc <= 0) {
            printf("\tWarning: variable nc rounded to %d, adjusting nc to 1\n", nc);
            nc = 1;
                    }
        }

   /* Adjust 
    *********/
   if (ULr >= n_rows) {
     printf("Error: Requested upper-left row value, %d\n", ULr);
     printf("       Total number of rows in data, %d\n", n_rows-1);
     exit(-5);
                      }
   if (overall_row0 >= n_rows) overall_row0 = 0;
   
   if (ULc >= n_cols) {
     printf("Error: Requested upper-left column value, %d\n", ULc);
     printf("       Total number of pixels in a row of data, %d\n", n_cols-1);
     exit(-5);
                      }
   if (overall_col0 >= n_cols) overall_col0 = 0;
   
   if ((ULc + nc) > n_cols) {
     nc = n_cols - ULc;
     printf("\tWarning: #_cols argument plus UL_col exceeds total number of pixels in a row of data;\n");
     printf("\t         #_cols argument adjusted to %d\n", nc);
                            }
   if ((overall_col0 + overall_cols) > n_cols) overall_cols = n_cols - overall_col0;
   
   if ((ULr + nr) > n_rows) {
     nr = n_rows - ULr;
     printf("\tWarning: #_rows argument plus UL_row exceeds total number of rows of data;\n");
     printf("\t         #_rows argument adjusted to %d\n", nr);
                            }
   if ((overall_row0 + overall_rows) > n_rows) overall_rows = n_rows - overall_row0;

   /* Get enough memory for everything, read everything;
    * assuming, of course, that nc, nr are relatively small...
    *************************************************************/
   k = (long)nr*(long)nc;


   /* New to version 2.5: will TRY to read everything at once,
    * hold all data in memory, and then take the sites from
    * data in memory, instead of read from a file 
    ************************************************************/
   if (newfileindex != (int *)NULL) {
   if (newfileindex[paramfileline] == 1) {
      forget_it = 0;
      for (j=0;j<n_criteria;j++) {
         if (WHOLESDS_filterimage_D[j] != (double *)NULL) {
	    free(WHOLESDS_filterimage_D[j]);
	    WHOLESDS_filterimage_D[j] = (double *)NULL;
	    }
         }
      for (j=0;j<n_SDSs;j++) {
         if (WHOLESDS_image_D[j] != (double *)NULL) {
	    free(WHOLESDS_image_D[j]);
	    WHOLESDS_image_D[j] = (double *)NULL;
	    }
         }
      for (j=0;j<n_criteria;j++) {
          /*WHOLESDS_filterimage_D[j] = get_images(crit_file_id[j], SDStofilter[j], -1, which_filter[j], 
                                         SDS_filtertype[j], fdExtsrc[j], n_cols, 
                                         n_rows*n_cols, 0, 0, n_rows, n_cols, subsample_factor[j]);*/

/*printf("%d %d %d %d %d\n", overall_rows*overall_cols, overall_row0, overall_col0, overall_rows, overall_cols);*/
          WHOLESDS_filterimage_D[j] = get_images(crit_file_id[j], SDStofilter[j], -1, which_filter[j], 
                                         SDS_filtertype[j], fdExtsrc[j], n_cols, 
                                         overall_rows*overall_cols, overall_row0, overall_col0, overall_rows, overall_cols, subsample_factor[j]);
	  if (WHOLESDS_filterimage_D[j] == (double *)NULL) forget_it = 1;			
          SDS_filterimage_D[j] = (double *)malloc(nr*nc*sizeof(double));
	  if (SDS_filterimage_D[j] == (double *)NULL) forget_it = 1;			
                              }  /* for (j=0;j<n_criteria;j++) */

      for (j=0;j<n_SDSs;j++) {
          /*WHOLESDS_image_D[j] = get_images(file_id, SDSnames[j], SDS_offset[j], which_SDS[j], SDS_type[j], NULL, 
                                   n_cols, n_rows*n_cols, 0, 0, n_rows, n_cols, 1.0);*/
/*printf("X %d %d %d %d %d\n", overall_rows*overall_cols, overall_row0, overall_col0, overall_rows, overall_cols);*/
          WHOLESDS_image_D[j] = get_images(file_id, SDSnames[j], SDS_offset[j], which_SDS[j], SDS_type[j], NULL, 
                                   n_cols, overall_rows*overall_cols, overall_row0, overall_col0, overall_rows, overall_cols, 1.0);
/*printf("%p\n", WHOLESDS_image_D[j] );*/

	  if (WHOLESDS_image_D[j] == (double *)NULL) forget_it = 1;			
          SDS_image_D[j] = (double *)malloc(nr*nc*sizeof(double));
	  if (SDS_image_D[j] == (double *)NULL) forget_it = 1;			
                          }   /* for (j=0;j<n_SDSs;j++) */
	
      if (forget_it) {
         printf("Memory errors: will use old method\n");
	 free(newfileindex);
	 newfileindex = (int *)NULL;
         for (j=0;j<n_criteria;j++) {
            if (WHOLESDS_filterimage_D[j] != (double *)NULL) {
	       free(WHOLESDS_filterimage_D[j]);
	       WHOLESDS_filterimage_D[j] = (double *)NULL;
	       }
            }
         for (j=0;j<n_SDSs;j++) {
            if (WHOLESDS_image_D[j] != (double *)NULL) {
	       free(WHOLESDS_image_D[j]);
	       WHOLESDS_image_D[j] = (double *)NULL;
	       }
            }
         for (j=0;j<n_criteria;j++) {
             SDS_filterimage_D[j] = get_images(crit_file_id[j], SDStofilter[j], -1, which_filter[j], 
                                         SDS_filtertype[j], fdExtsrc[j], n_cols, 
                                         k, ULr, ULc, nr, nc, subsample_factor[j]);
				
                                 }  /* for (j=0;j<n_criteria;j++) */

         for (j=0;j<n_SDSs;j++) {
             SDS_image_D[j] = get_images(file_id, SDSnames[j], SDS_offset[j], which_SDS[j], SDS_type[j], NULL, 
                                   n_cols, k, ULr, ULc, nr, nc, 1.0);
                             }   /* for (j=0;j<n_SDSs;j++) */
          }
       else {
	       lj = 0L;
               for (lin=(ULr-overall_row0);lin<(ULr-overall_row0)+nr;lin++) {   /* for each line... */
	          li = lin*overall_cols;
                  for (col=(ULc-overall_col0);col<(ULc-overall_col0)+nc;col++) {   /* for each column... */
	     
                     for (j=0;j<n_criteria;j++) {
                         SDS_filterimage_D[j][lj] = WHOLESDS_filterimage_D[j][li + col];
                                }  /* for (j=0;j<n_criteria;j++) */

                     for (j=0;j<n_SDSs;j++) {
                         SDS_image_D[j][lj] = WHOLESDS_image_D[j][li + col];
                             }   /* for (j=0;j<n_SDSs;j++) */
			     
		     lj++;
		         }
                    }
              for (j=0;j<n_SDSs;j++)free(WHOLESDS_image_D[j]);
              for (j=0;j<n_criteria;j++)free(WHOLESDS_filterimage_D[j]);
			     
	  }		
      }  /*  if (newfileindex[paramfileline] == 1)  */
      else {  /* else, take from the WHOLESDS buffers */
         
         if (forget_it) {
            for (j=0;j<n_criteria;j++) {
                SDS_filterimage_D[j] = get_images(crit_file_id[j], SDStofilter[j], -1, which_filter[j], 
                                         SDS_filtertype[j], fdExtsrc[j], n_cols, 
                                         k, ULr, ULc, nr, nc, subsample_factor[j]);
				
                                 }  /* for (j=0;j<n_criteria;j++) */

            for (j=0;j<n_SDSs;j++) {
                SDS_image_D[j] = get_images(file_id, SDSnames[j], SDS_offset[j], which_SDS[j], SDS_type[j], NULL, 
                                   n_cols, k, ULr, ULc, nr, nc, 1.0);
                             }   /* for (j=0;j<n_SDSs;j++) */
          }
	  else {
	  
	      for (j=0;j<n_criteria;j++) SDS_filterimage_D[j] = (double *)malloc(nr*nc*sizeof(double));
              for (j=0;j<n_SDSs;j++) SDS_image_D[j] = (double *)malloc(nr*nc*sizeof(double));

	      lj = 0L;
              for (lin=(ULr-overall_row0);lin<(ULr-overall_row0)+nr;lin++) {   /* for each line... */
	          li = lin*overall_cols;
                  for (col=(ULc-overall_col0);col<(ULc-overall_col0)+nc;col++) {   /* for each column... */
	     
                     for (j=0;j<n_criteria;j++) {
                         SDS_filterimage_D[j][lj] = WHOLESDS_filterimage_D[j][li + col];
                                }  /* for (j=0;j<n_criteria;j++) */

                     for (j=0;j<n_SDSs;j++) {
                         SDS_image_D[j][lj] = WHOLESDS_image_D[j][li + col];
                             }   /* for (j=0;j<n_SDSs;j++) */
                     lj++;
			     
		         }
                    }
			     
	  }	
      }
      }  /*  if (newfileindex != (int *)NULL)  */
   else {   
      for (j=0;j<n_criteria;j++) {
           SDS_filterimage_D[j] = get_images(crit_file_id[j], SDStofilter[j], -1, which_filter[j], 
                                         SDS_filtertype[j], fdExtsrc[j], n_cols, 
                                         k, ULr, ULc, nr, nc, subsample_factor[j]);
				
                              }  /* for (j=0;j<n_criteria;j++) */

      for (j=0;j<n_SDSs;j++) {
           SDS_image_D[j] = get_images(file_id, SDSnames[j], SDS_offset[j], which_SDS[j], SDS_type[j], NULL, 
                                   n_cols, k, ULr, ULc, nr, nc, 1.0);
                          }   /* for (j=0;j<n_SDSs;j++) */
	}
	
	
	
/* New, 20-SEP-01: if SDStofilter contains lat and lon, then read lats and lons twice each 
   again, into SDS_image_D[n_SDSs+1], SDS_image_D[n_SDSs+2], SDS_image_D[n_SDSs+3], 
   SDS_image_D[n_SDSs+4] -- then OVERWRITE the contents of these SDSs with sun and view
   angles (new for version 1.6). */
   latid = -1;
   lonid = -1;
   for (j=0;j<n_criteria;j++) {
      for (i=0;i<strlen(SDStofilter[j]);i++) ALLCAPS[i] = toupper(SDStofilter[j][i]);
      ALLCAPS[i] = '\0';
      if (strstr(ALLCAPS, "LAT")) latid = j;
      if (strstr(ALLCAPS, "LON")) lonid = j;
   }
   if ((latid != -1)&&(lonid != -1)&&(get_angles)) {
       which_SDS[n_SDSs] = (int32)which_filter[latid];
       SDS_type[n_SDSs] = SDS_filtertype[latid];
       strcpy(SDSnames[n_SDSs], SDStofilter[latid]);
       SDS_image_D[n_SDSs++] = get_images(file_id, SDSnames[n_SDSs], -1, which_SDS[n_SDSs], 
                                          SDS_type[n_SDSs], NULL, n_cols, k, ULr, ULc, nr, nc, subsample_factor[latid]);
       strcpy(SDSnames[n_SDSs-1], "solar zenith angle");
       if (longest < strlen(SDSnames[n_SDSs-1])) longest = strlen(SDSnames[n_SDSs-1]);

       which_SDS[n_SDSs] = (int32)which_filter[latid];
       SDS_type[n_SDSs] = SDS_filtertype[latid];
       strcpy(SDSnames[n_SDSs], SDStofilter[latid]);
       SDS_image_D[n_SDSs++] = get_images(file_id, SDSnames[n_SDSs], -1, which_SDS[n_SDSs], 
                                          SDS_type[n_SDSs], NULL, n_cols, k, ULr, ULc, nr, nc, subsample_factor[latid]);
       strcpy(SDSnames[n_SDSs-1], "solar azimuth angle");
       if (longest < strlen(SDSnames[n_SDSs-1])) longest = strlen(SDSnames[n_SDSs-1]);
					  
       which_SDS[n_SDSs] = (int32)which_filter[lonid];
       SDS_type[n_SDSs] = SDS_filtertype[lonid];
       strcpy(SDSnames[n_SDSs], SDStofilter[lonid]);
       SDS_image_D[n_SDSs++] = get_images(file_id, SDSnames[n_SDSs], -1, which_SDS[n_SDSs], 
                                          SDS_type[n_SDSs], NULL, n_cols, k, ULr, ULc, nr, nc, subsample_factor[lonid]);
       strcpy(SDSnames[n_SDSs-1], "view zenith angle");
       if (longest < strlen(SDSnames[n_SDSs-1])) longest = strlen(SDSnames[n_SDSs-1]);

       which_SDS[n_SDSs] = (int32)which_filter[lonid];
       SDS_type[n_SDSs] = SDS_filtertype[lonid];
       strcpy(SDSnames[n_SDSs], SDStofilter[lonid]);
       SDS_image_D[n_SDSs++] = get_images(file_id, SDSnames[n_SDSs], -1, which_SDS[n_SDSs], 
                                          SDS_type[n_SDSs], NULL, n_cols, k, ULr, ULc, nr, nc, subsample_factor[lonid]);
       strcpy(SDSnames[n_SDSs-1], "view azimuth angle");
       if (longest < strlen(SDSnames[n_SDSs-1])) longest = strlen(SDSnames[n_SDSs-1]);

       /* Assuming we've just read the whole Lat, Lon SDSs, take dlatsub, dlonsub as
          the lat and lon at the middle -- nadir -- of each line */
       latid = n_SDSs-3;
       lonid = n_SDSs-1;
       li = 0L;
       for (irow=0;irow<nr;irow++) { 
          dlatsub=(double)SDS_image_D[latid][(irow*nc)+nc/2]*ROVRDG;
          dlonsub=(double)SDS_image_D[lonid][(irow*nc)+nc/2]*ROVRDG;
          for (icol=0;icol<nc;icol++) { 
             dlat=(double)SDS_image_D[latid][li]*ROVRDG;
             dlon=(double)SDS_image_D[lonid][li]*ROVRDG;
             sun_angles(jday,gmt,dlat,dlon,&ts,&fs);
             sat_angles(dlat,dlon,dheight,dlatsub,dlonsub,&tv,&fv);
	     
	     /* and overwrite */
	     SDS_image_D[n_SDSs-4][li] = ts;
	     SDS_image_D[n_SDSs-3][li] = fs;
	     SDS_image_D[n_SDSs-2][li] = tv;
	     SDS_image_D[n_SDSs-1][li] = fv;
	     li++;
          }
       }

   }


   /* Apply the criteria to a keep/no-keep mask
    *********************************************/
   k = (long)nr*(long)nc;
   if ( (keepmask = (int16 *)malloc(k*sizeof(int16))) == NULL) {
       printf("Error: out of memory, array 'keepmask'\n");
       exit(-7);
   }

   for (j=0;j<n_criteria;j++) {
      lcrit_value = (long)lowend[j];
      hcrit_value = (long)highend[j];
      truly_boolean = really_boolean[j];
      if ((truly_boolean == 1)||(truly_boolean == 2)) {
         if (hcrit_value == -999) hcrit_value = lcrit_value; 
	    } 
     
      for (li=0L;li<k;li++) {   /* for each pixel... */
         if (j == 0) keepmask[li] = 1;      /* Assume you'll keep it */
   
         /* Check the values of the pixel in the criteria arrays against the user-entered
            criteria: */
         if (truly_boolean == 1) { 
            if (((long)SDS_filterimage_D[j][li]&lcrit_value) != hcrit_value) keepmask[li] = 0;
         }
         else if (truly_boolean == 2) {
            if (((long)SDS_filterimage_D[j][li]&lcrit_value) == hcrit_value) keepmask[li] = 0;
        }
         else {
            if ((SDS_filterimage_D[j][li] < lowend[j])||(SDS_filterimage_D[j][li] > highend[j]))  
                keepmask[li] = 0;
              }
           
                              
              } /* for (li=0L;li<k;li++) */
           }

    if (n_criteria == 0) for (li=0L;li<k;li++) keepmask[li] = 1;


   /* Right here -- apply aggregation to 
    * pixels to-be-kept; adjust nr, nc, k accordingly;
    * set keepmask[] to N pixels aggregated (v2.8).
    **************************************************/
   if (aggreg) {
       aggregate(n_SDSs, keepmask, SDS_image_D, &nr, &nc);
       k = (long)nr*(long)nc;
      }

   /* Get mean, SD of all the non-masked pixels 
    * in each of the SDS-to-be-processed 
    *******************************************/
   for (j=0;j<n_SDSs;j++) {
       SDS_sum[j] = 0.0;
       SDS_avg[j] = 0.0;
       SDS_SD[j] = 0.0;
       }
   n_kept = 0L;
   if (save_bin) 
      for (j=0;j<n_SDSs;j++) testimage[j] = (short *)calloc(nc,sizeof(short));  /* one line */
   if (save_HDFout) 
      for (j=0;j<n_SDSs;j++) f_image[j] = (float *)calloc(nc,sizeof(float));

  
   for (j=0;j<n_SDSs;j++) {

      for (li=0L;li<k;li++) {   /* for each pixel... */
          if (keepmask[li] > 0) {
	
              if (j == 0) n_kept++;
              tmpdb = 0.0;
              tmpdb = SDS_image_D[j][li];   
	      SDS_sum[j] += tmpdb;
                                 }  
                            } /* for (li=0L;li<k;li++) */
			 }/* for (j=0;j<n_SDSs;j++) */

   


   if ((save_bin) || (save_HDFout)) {    /* new, v2.3 */

      for (j=0;j<n_SDSs;j++) {
        
	if (save_bin) {
         sprintf(debugfile, "%s_%s_samples_%02d", binfile, SDSnames[j], ii);
         if ((fddebug = fopen(debugfile,"wb") ) == NULL) {
             printf("\tWarning: cannot open file '%s'; will open file named ", debugfile);
             sprintf(debugfile, "%s_samples_%02d", binfile, ii);
	     printf("'%s'\n", debugfile);
	     fddebug = fopen(debugfile,"wb");
	    }
	  }
	  
        if (save_HDFout)  {
           sprintf(debugfile, "%s\0", SDSnames[j]);
           outdims[0] = nr;
	   outdims[1] = nc;
	   outstart[0] = 0;
	   outstart[1] = 0;
           /*printf("%d %s %d %d %d\n", j, SDSnames[j], SDS_type[j], outdims[0], outdims[1]);*/
           sds_outid = SDcreate(out_id, debugfile, SDS_type[j], 2, outdims);
	   
           SDsetdimname( SDgetdimid(sds_outid, 0), "rows\0");	   
           SDsetdimname( SDgetdimid(sds_outid, 1), "columns\0");	   
	   
           cdef.chunk_lengths[0] = 50;
           cdef.chunk_lengths[1] = (int32)outdims[1];
           cdef.comp.comp_type=COMP_CODE_DEFLATE;      /* gzip */
           cdef.comp.cinfo.deflate.level=8;
           SDsetchunk(sds_outid, cdef, HDF_CHUNK|HDF_COMP); 
	        }
	     
         li = -1L;
         for (lin=0;lin<nr;lin++) {   /* for each line... */
           for (col=0;col<nc;col++) {   /* for each column... */
	       li++;
               if (keepmask[li] > 0) {
	
                   if (save_bin) testimage[j][col] = (short)SDS_image_D[j][li];
                   if (save_HDFout) f_image[j][col] = (float)SDS_image_D[j][li];
                                      }
               else {
	
                   if (save_bin) testimage[j][col] = (short)0;
                   if (save_HDFout) f_image[j][col] = (float)0.0;
                    }
				   
                            } /* for (col=0;col<nc;col++) */

           if (save_bin) 
              fwrite(testimage[j], nc*sizeof(short), 1, fddebug);
	      
           if (save_HDFout) {
	   outstart[0] = lin;
           outdims[0] = 1;
	   
	/* 04-DEC-06: we have to convert this data back to the original format! */
	if (SDS_type[j] == DFNT_INT16) {
	     if (lin==0) i16_data = (int16 *)malloc(nc*sizeof(int16));
	     for (col=0;col<nc;col++) i16_data[col] = (int16)f_image[j][col];
             SDwritedata(sds_outid, outstart, NULL, outdims, (VOID *)i16_data);
	     if (lin==(nr-1)) free(i16_data);
	   }
	else if (SDS_type[j] == DFNT_INT8) {
	     if (lin==0) i8_data = (int8 *)malloc(nc);
	     for (col=0;col<nc;col++) i8_data[col] = (int8)f_image[j][col];
             SDwritedata(sds_outid, outstart, NULL, outdims, (VOID *)i8_data);
	     if (lin==(nr-1)) free(i8_data);
	   }
	else if (SDS_type[j] == DFNT_UINT16) {
	     if (lin==0) ui16_data = (uint16 *)malloc(nc*sizeof(uint16));
	     for (col=0;col<nc;col++) ui16_data[col] = (uint16)f_image[j][col];
             SDwritedata(sds_outid, outstart, NULL, outdims, (VOID *)ui16_data);
	     if (lin==(nr-1)) free(ui16_data);
	   }
	else if (SDS_type[j] == DFNT_UINT8) {
	     if (lin==0) ui8_data = (uint8 *)malloc(nc);
	     for (col=0;col<nc;col++) ui8_data[col] = (uint8)f_image[j][col];
             SDwritedata(sds_outid, outstart, NULL, outdims, (VOID *)ui8_data);
	     if (lin==(nr-1)) free(ui8_data);
	   }
	else if (SDS_type[j] == DFNT_UINT32) {
	     if (lin==0) ui32_data = (uint32 *)malloc(nc*sizeof(uint32));
	     for (col=0;col<nc;col++) ui32_data[col] = (uint32)f_image[j][col];
             SDwritedata(sds_outid, outstart, NULL, outdims, (VOID *)ui32_data);
	     if (lin==(nr-1)) free(ui32_data);
	   }
	else if (SDS_type[j] == DFNT_INT32) {
	     if (lin==0) i32_data = (int32 *)malloc(nc*sizeof(int32));
	     for (col=0;col<nc;col++) i32_data[col] = (int32)f_image[j][col];
             SDwritedata(sds_outid, outstart, NULL, outdims, (VOID *)i32_data);
	     if (lin==(nr-1)) free(i32_data);
	   }
	else if ((SDS_type[j] == DFNT_FLOAT64)||(SDS_type[j] == DFNT_DOUBLE)) {
	     if (lin==0) dbl_data = (double *)malloc(nc*sizeof(double));
	     for (col=0;col<nc;col++) dbl_data[col] = (double)f_image[j][col];
             SDwritedata(sds_outid, outstart, NULL, outdims, (VOID *)dbl_data);
	     if (lin==(nr-1)) free(dbl_data);
	   }
	else
           SDwritedata(sds_outid, outstart, NULL, outdims, (VOID *)f_image[j]);
	      }
 

                          } /* for (lin=0;lin<nr;lin++) */
         if (save_bin) {
           fclose(fddebug);
           free(testimage[j]);
             }
	     
         if (save_HDFout)  {
	 	 
	   /* xfer metadata also */
	   xsds = SDselect(file_id, SDnametoindex(file_id, SDSnames[j]));
	   process_local_attributes(xsds, sds_outid);
	   SDendaccess(xsds); 
	
	   SDendaccess(sds_outid);

           free(f_image[j]);
	   }

			 }/* for (j=0;j<n_SDSs;j++) */

         if (save_HDFout)  {
	   /* for now, save keepmask as well */
           outdims[0] = nr;
	   outdims[1] = nc;
	   outstart[0] = 0;
	   outstart[1] = 0;
           sds_outid = SDcreate(out_id, "keepmask\0", DFNT_INT16, 2, outdims);
           cdef.chunk_lengths[0] = 50;
           cdef.chunk_lengths[1] = (int32)outdims[1];
           cdef.comp.comp_type=COMP_CODE_DEFLATE;      /* gzip */
           cdef.comp.cinfo.deflate.level=8;
           SDsetchunk(sds_outid, cdef, HDF_CHUNK|HDF_COMP); 
           SDwritedata(sds_outid, outstart, NULL, outdims, (VOID *)keepmask);
	   SDendaccess(sds_outid);
	   }


       }  /* save */


   if (n_kept > 0)
      for (j=0;j<n_SDSs;j++) {
          SDS_avg[j] = SDS_sum[j]/((double)n_kept);
          SDS_sum[j] = 0.0;
          SDS_SD[j] = 0.0;
          }
   
   if (n_kept > 1) {
      for (j=0;j<n_SDSs;j++) {  
         for (li=0L;li<k;li++) {   /* for each pixel... */
             if (keepmask[li] > 0) {
                 tmpdb = 0.0;
                 tmpdb = (SDS_image_D[j][li]) - SDS_avg[j]; 
	         SDS_sum[j] += tmpdb*tmpdb;
			            }
                               } /* for (li=0L;li<k;li++) */
                             }  /* for (j=0;j<n_SDSs;j++) */
			 
      for (j=0;j<n_SDSs;j++) 
          SDS_SD[j] = sqrt(SDS_sum[j]/((double)(n_kept-1)));
                   }

	  
#if defined (DEBUG)
   fddebug = fopen("keepmask", "wb");
   fwrite(keepmask,  nc*nr*sizeof(int16),  1, fddebug);
   fclose(fddebug);
#endif

   /* Report results to user 
    ***************************/
   if (VERBOSE) {
      sprintf(sds_table,"\n     Site %d\n     Criteria selected %d out of %d pixels %s:\n", 
            ii, n_kept, nr*nc, (aggreg == 1 ? "(aggregated)\0":"\0") );
      for (j=0;j<n_SDSs;j++) {
         if (SDS_offset[j] == -1) 
	    sprintf(sds_table+strlen(sds_table),"       SDS %02d, '%s': ", j, SDSnames[j]);
         else 
	    sprintf(sds_table+strlen(sds_table),"       SDS %02d, '%s' level %d: ", j, SDSnames[j], SDS_offset[j]);
         for (i=0;i<longest-strlen(SDSnames[j]);i++) sprintf(sds_table+strlen(sds_table)," ");
	 sprintf(sds_table+strlen(sds_table),"\tAVG = %f, \tSD = %f\n", SDS_avg[j], SDS_SD[j]);
	     }
      sprintf(sds_table+strlen(sds_table),"\n");
                }

		
   if (COLUMNS) {
       sprintf(resultstring,"Site %d: ",ii);
       if (proj_ok) 
           sprintf(resultstring+strlen(resultstring),"UL lat, lon; (%f %f); LR lat, lon; (%f %f), ",
                 ULR[ii], ULC[ii], NR[ii], NC[ii]);
       else
           sprintf(resultstring+strlen(resultstring),"UL row, column; (%d %d); n rows, columns; (%d %d), ", 
                 (short)ULR[ii], (short)ULC[ii], (short)NR[ii], (short)NC[ii]);
       sprintf(resultstring+strlen(resultstring),"n_samples, %d\n",n_kept);
       sprintf(resultstring+strlen(resultstring),"SDS name");
       for (j=0;j<longest-7;j++) sprintf(resultstring+strlen(resultstring)," "); 
       sprintf(resultstring+strlen(resultstring),"mean      SD\n");

       for (j=0;j<n_SDSs;j++) {
           if (SDS_offset[j] == -1) 
             sprintf(resultstring+strlen(resultstring),"%s ", SDSnames[j]);
	   else
             sprintf(resultstring+strlen(resultstring),"%s level %d ", SDSnames[j], SDS_offset[j]);
           for (i=0;i<longest-strlen(SDSnames[j]);i++) sprintf(resultstring+strlen(resultstring)," ");
           sprintf(resultstring+strlen(resultstring),"%09.3f %09.3f\n", SDS_avg[j], SDS_SD[j]);
	   }
                  
                }
   else {
     if ((ii == 0)||(samples)) {
        sprintf(resultstring,"st_nbs_o/o_N__ ");
        for (j=0;j<n_SDSs;j++) 
           sprintf(resultstring+strlen(resultstring),"SDS_%02d_mn SDS_%02d_sd ",j,j); 
        sprintf(resultstring+strlen(resultstring),"\n"); 
                  }
     else 
       sprintf(resultstring,"\0");
     
     sprintf(resultstring+strlen(resultstring),"%02d %05d %05d ", ii, n_kept, nr*nc);
     for (j=0;j<n_SDSs;j++) 
           sprintf(resultstring+strlen(resultstring),"%09.3f %09.3f ", SDS_avg[j], SDS_SD[j]);
     sprintf(resultstring+strlen(resultstring),"\n");
        }
   
   if (save_out) {
      fprintf(fdout,"%s", resultstring);
      if (VERBOSE) {
          printf("%s", sds_table); 
          printf("%s", resultstring); 
                   }  
                 }               
   else {
      if (VERBOSE) printf("%s", sds_table); 
      printf("%s", resultstring); 
        }  

   /* New ascii dump to stdio */
   if ((samples)&&(n_kept>0)) {
     if (save_out) {
        fprintf(fdout,"\nASCII dump of data samples:\n\n     row  col ");
        if (VERBOSE) printf("\nASCII dump of data samples:\n\n     row  col ");
        for (j=0;j<n_SDSs;j++) {
             fprintf(fdout,"   SDS_%02d",j); 
             if (VERBOSE) printf("   SDS_%02d",j); 
	     }
        fprintf(fdout,"\n"); 
        if (VERBOSE)printf("\n"); 
       
       
        for (li=0L;li<k;li++) {   /* for each pixel... */
           if (keepmask[li] > 0) {
              fprintf(fdout,"ASC %4d %4d ", li/nc, li%nc);
              if (VERBOSE) printf("ASC %4d %4d ", li/nc, li%nc);
              for (j=0;j<n_SDSs;j++) {
		 if ((SDS_type[j] == DFNT_FLOAT32)||(SDS_type[j] == DFNT_FLOAT64)) {
                    fprintf(fdout,"%8.4e ", SDS_image_D[j][li]);
                    if (VERBOSE) printf("%8.4e ", SDS_image_D[j][li]);
		    }
		 else {
                    fprintf(fdout,"%8d ", (int) SDS_image_D[j][li]);
                    if (VERBOSE) printf("%8d ", (int) SDS_image_D[j][li]);
		    }
		 }
              fprintf(fdout,"\n");
              if (VERBOSE) printf("\n");
                                  }
                                   
                              } /* for (li=0L;li<k;li++) */
                   }               
     else {
        printf("\nASCII dump of data samples:\n\n     row  col ");
        for (j=0;j<n_SDSs;j++) {
             printf("   SDS_%02d",j); 
	     }
        printf("\n"); 
       
       
        for (li=0L;li<k;li++) {   /* for each pixel... */
           if (keepmask[li] > 0) {
              printf("ASC %4d %4d ", li/nc, li%nc);
              for (j=0;j<n_SDSs;j++) {
		 if ((SDS_type[j] == DFNT_FLOAT32)||(SDS_type[j] == DFNT_FLOAT64)) {
                    printf("%8.4e ", SDS_image_D[j][li]);
		    }
		 else {
                    printf("%8d ", (int) SDS_image_D[j][li]);
		    }
  		 }
              printf("\n");
                                  }
                                   
                              } /* for (li=0L;li<k;li++) */
          }  
		           } /* if ((samples)&&(n_kept>0)) */


   /* Clean up 
    *************/
   for (j=0;j<n_SDSs;j++)     free(SDS_image_D[j]);
   for (j=0;j<n_criteria;j++) free(SDS_filterimage_D[j]);
   free(keepmask);


                        } /* NEW, 17-NOV-98: for (ii=0;ii<n_sites;ii++) */
if (save_out) fclose(fdout);
for (j=0;j<n_criteria;j++) {
   if (fdExtsrc[j] != 0) fclose(fdExtsrc[j]);   
   if (crit_file_id[j] != file_id) SDend(crit_file_id[j]);
   }
/*SDend(file_id);*/
if (save_HDFout)SDend(out_id);

/* 22-FEB-07 */ }


return(0);
}

/***************************************************************************************/

void aggregate(int n_SDSs, int16 *keepmask, double *SDS_image_D[MAXBANDS], short *nr, short *nc)
{
int i, rowi, coli, agg_rowi, agg_coli;
short newnc, newnr;
double *data, *agg_data;
long l,k,m;
int16 *newkeepmask;

k = (long)(*nr)*(long)(*nc);
data = (double *)malloc(k*sizeof(double));

/*printf("%d %d\n", agg_dim1, agg_dim2);*/

newnr = (*nr)/agg_dim1;
newnc = (*nc)/agg_dim2;
if ( (*nr)%agg_dim1 ) newnr++;
if ( (*nc)%agg_dim2 ) newnc++;

m = (long)newnr*(long)newnc;
agg_data = (double *)malloc(m*sizeof(double));
newkeepmask = (int16 *)malloc(m*sizeof(int16));

for(i=0;i<n_SDSs;i++) {

   for(l=0L;l<m;l++) {
       newkeepmask[l] = 0;
       agg_data[l] = 0.0;
       }
   for(l=0L;l<k;l++) {
      data[l] = SDS_image_D[i][l];
      SDS_image_D[i][l] = 0.0;
      }
   
   l = 0L;
   for(rowi=0;rowi<(*nr);rowi++) {
      agg_rowi = rowi/agg_dim1;
      m = (long)agg_rowi*newnc;
      for(coli=0;coli<(*nc);coli++) {
         agg_coli = coli/agg_dim2;
         if (keepmask[l] > 0) {
	    agg_data[m+agg_coli] += data[l];
	    newkeepmask[m+agg_coli]++;	    
	    }
	 l++;

         }
      }
      
   m = (long)newnr*(long)newnc;
   for(l=0L;l<m;l++) {
      if ((newkeepmask[l] > 1) && (newkeepmask[l] >= agg_thresh)) {
         SDS_image_D[i][l] = agg_data[l];
         SDS_image_D[i][l] /= newkeepmask[l];
         }
      }

   }  /* for(i=0;i<n_SDSs;i++) */
   
m = (long)newnr*(long)newnc;
for(l=0L;l<m;l++) {
   if ((newkeepmask[l] > 1) && (newkeepmask[l] >= agg_thresh)) keepmask[l] = newkeepmask[l];
   else keepmask[l] = 0;
   }
   
*nr = newnr;
*nc = newnc;
free(agg_data);
free(newkeepmask);
free(data);
return;
}

/***************************************************************************************/

int32 bytes_to_HDFtype(short i)
{
/* an attempt to convert between number of bytes (i) and HDF file type */
if (i == 1) {
    printf("bytes_to_HDFtype: number of bytes = %d, will interpret as unsigned 8-bit integer.\n", i);
    return(DFNT_UINT8);   
	    }
else if (i == 2) {
    printf("bytes_to_HDFtype: number of bytes = %d, will interpret as signed short datatype.\n", i);
    return(DFNT_INT16);
                 }
else if (i == 4) {
    printf("bytes_to_HDFtype: number of bytes = %d, will interpret as floating-point datatype.\n", i);
    return(DFNT_FLOAT);
                 }
else if (i == 8) {
    printf("bytes_to_HDFtype: number of bytes = %d, will interpret as double precision datatype.\n", i);
    return(DFNT_DOUBLE);
                 }
else {
   printf("bytes_to_HDFtype, error: number of bytes unrecognized -- %d\n", i);
   exit(-10);
     }
}

/***************************************************************************************/

short SDSread(int32 file_id, int32 offset, int32 SDS_index, void *image, short ULr, short ULc, short nr, short nc)
{  /* this routine assumes SDS_index is the index of the SDS you want, and the SDSs are 2-dim */
int32 sds_id;
int32 start2[MAXDIMS];  
int32 stride2[MAXDIMS] = {1,1};  
int32 edge2[MAXDIMS];
int32 dims2[MAXDIMS];
int32 start3[MAXDIMS+1];  
int32 stride3[MAXDIMS+1] = {1,1,1};  
int32 edge3[MAXDIMS+1];
int32 dims3[MAXDIMS+1];
char name[100];
int32 rank, nt, nattr;

if ((sds_id = SDselect(file_id, SDS_index)) == -1) {
    printf("SDSread Error; finding SDS number %d from file\n", SDS_index);
    SDendaccess(sds_id);
    return(-1);
	      }
	      
if(offset == -1) {     /* two-dimensional */

   start2[0] = ULr;
   start2[1] = ULc;
   edge2[0] = nr;
   edge2[1] = nc;

   if ((SDreaddata(sds_id, start2, NULL, edge2, image)) == -1) {
       printf("SDSread Error; reading SDS number %d from file\n", SDS_index);
       SDendaccess(sds_id);
       return(-1);
	         }
   }	      
else {     /* three-dimensional */

   start3[0] = offset;
   start3[1] = ULr;
   start3[2] = ULc;
   edge3[0] = 1;
   edge3[1] = nr;
   edge3[2] = nc;
   if ((SDreaddata(sds_id, start3, NULL, edge3, image)) == -1) {
       printf("SDSread Error; reading SDS number %d from file (offset = %d)\n", SDS_index, offset);
       SDendaccess(sds_id);
       return(-1);
	         }
   }	      
SDendaccess(sds_id);

return 0;
}

/***************************************************************************************/

void get_site_overall_dims(char *string)
{
short areal_or_boolean(char *);  /* this step populates global float arrays ULR, ULC, NR and NC */

if ( areal_or_boolean(string) ) {
   overall_set = 1;
   /*printf("%s: ULr = %f, ULc = %f, nr = %f, nc = %f\n", string, ULR[n_sites-1], ULC[n_sites-1], NR[n_sites-1], NC[n_sites-1]);  */
   if ((overallULR < 0.0) || (overallULR > ULR[n_sites-1])) overallULR = ULR[n_sites-1];
   if ((overallULC < 0.0) || (overallULC > ULC[n_sites-1])) overallULC = ULC[n_sites-1];
   if ((overallNR < 0.0) || (overallNR < (ULR[n_sites-1]+NR[n_sites-1]))) overallNR = (ULR[n_sites-1]+NR[n_sites-1]);
   if ((overallNC < 0.0) || (overallNC < (ULC[n_sites-1]+NC[n_sites-1]))) overallNC = (ULC[n_sites-1]+NC[n_sites-1]);
   }
return;
}


/***************************************************************************************/

int get_n_lines(char *possibleparamfile, int **newfileindex)
{
FILE *fdparam;
char c;
int n_lines = 0;

char argv[100][150];
char lastfilename[150];
int i, j, list_idx, crit_idx, argc, ai;
int in_list, in_crit;
int in_arg, en_lines, n_quotes, in_quotes;
void get_site_overall_dims(char *string);

if ((fdparam = fopen(possibleparamfile, "r")) == (FILE *)NULL) {
    return(1);  /* silent return */
   }
   
while ( (c = getc(fdparam)) != EOF ) {
   if (c == '\n') n_lines++;   
   }

free(*newfileindex);  /* this has been allocated earlier */

if (( *newfileindex = calloc(n_lines, sizeof(int))) == (int *)NULL) {
    printf("Error trying to allocate memory for %d ints; cannot continue\n", n_lines);
    exit(-2);
   }

rewind(fdparam);

list_idx = crit_idx = 0;
in_list = in_crit = 0;
in_arg = n_quotes = in_quotes = en_lines =  0;
argc = -1;

lastfilename[0] = '\0';

while ( (c = getc(fdparam)) != EOF ) {
   
   
        if (c == '"') {
         n_quotes++;
	 if ((n_quotes%2) == 1){
	     if ((in_quotes == 0)) {
	         argc++;
		 argv[argc][0] = '\0';
		 }
	     in_quotes = 1;
	     }
	 else  in_quotes = 0;
         }
      
      if (in_quotes == 0) {
         if (c != ' ') {
	     if ((in_arg == 0) ) {
	         argc++;
		 argv[argc][0] = '\0';
		 }
	     in_arg = 1;
	     }
         else  in_arg = 0;
	 }
      else  in_arg = 1;
       
      if (c == '\n') {
            
          for(j=0;j<argc;j++) {
	      if ( (argv[j][0] == '-') && (argv[j][1] == 'f')  )   {
	         if (  strcmp(lastfilename, argv[j+1])) {
		     (*newfileindex)[en_lines] = 1;
		     strcpy(lastfilename, argv[j+1]);
		     }
		 }
	      if ( (argv[j][0] == '-') && (argv[j][1] == 'w')  )   {
	         get_site_overall_dims(argv[j+1]);
		 }
	     }
	     
      /*printf("%d\n", en_lines);
      for(ai=0;ai<argc;ai++) printf("* %s *", argv[ai]);
      printf("\n");*/
     
          en_lines++;
	  in_arg = 1;
	  argc = 0;
	  i = 0;
	  }
      else if (in_arg && (c != '"') ) {
          argv[argc][i++] = c;
          argv[argc][i] = '\0';
         }
      else i = 0;


      
      /*printf("c = '%c'; argc = %d, inarg = %d, inquote = %d %d %d\n", c, argc, in_arg, in_quotes, n_lines, en_lines);*/
      }
      
n_sites = 0;
   
fclose(fdparam); 
return(n_lines);  
}

/***************************************************************************************/

short process_arg(int in_argc, char **in_argv, int line_n)
{
FILE *fdparam;
FILE *fdExternSDSlist;
char listargs[MAXLENGTH][150];
char critargs[MAXLENGTH][150];
char argv[100][150];
char c;
int i, list_idx, crit_idx, argc;
int in_list, in_crit;
short ret = 1;
short areal_encountered = 0;
short areal_or_boolean(char *);
char *arg_clean(char *);
/*int initializepar(char *file,t_project *project);*/
int in_arg, n_lines, n_quotes, in_quotes;
int parse_agg_arg(char *arg);
int aggret;

list_idx = crit_idx = 0;
in_list = in_crit = 0;
in_arg = n_lines = n_quotes = in_quotes = 0;

if (in_argc < 2) return(0);

if (in_argc == 2) {   /* new, v2.2 (22-FEB-07) */
   /* let's assume it's a parameter file: open the file and 
      read the commands there */
   if ((fdparam = fopen(in_argv[1], "r")) == (FILE *)NULL) {
       printf("cannot open parameter file '%s'\n", in_argv[1]);
       return(-2);
      }
   argc = -1;
   i = 0;
      
   while ( (c = getc(fdparam)) != EOF ) {
   
      if (c == '"') {
         n_quotes++;
	 if ((n_quotes%2) == 1){
	     if ((in_quotes == 0) && (n_lines == line_n)) {
	         argc++;
		 argv[argc][0] = '\0';
		 }
	     in_quotes = 1;
	     }
	 else  in_quotes = 0;
         }
      
      if (in_quotes == 0) {
         if (c != ' ') {
	     if ((in_arg == 0) && (n_lines == line_n)) {
	         argc++;
		 argv[argc][0] = '\0';
		 }
	     in_arg = 1;
	     }
         else  in_arg = 0;
	 }
      else  in_arg = 1;
       
      if (c == '\n') {
          n_lines++;
	  in_arg = 1;
	  if (n_lines == line_n) argc = 0;
	  i = 0;
	  }
      else if (in_arg && (c != '"') && (n_lines == line_n)) {
          argv[argc][i++] = c;
          argv[argc][i] = '\0';
         }
      else i = 0;
      
      
      /*printf("c = '%c'; argc = %d, inarg = %d, inquote = %d %d %d\n", c, argc, in_arg, in_quotes, n_lines, line_n);   */
      }
   argc++; 
   }  /* end of new material: however, we still
         have to at least copy the incoming arguments
	 to argv */
else {
   argc = in_argc;
   for (i=0;i<argc;i++) strcpy(argv[i], in_argv[i]);
   }

/*printf("line_n = %d, argc = %d\n", line_n, argc);   
for (i=0;i<argc;i++) printf("argv[%d] = '%s'\n", i, argv[i]);*/


for (i=0;i<argc;i++) {
    
    if (argv[i][0] == '-')
       switch(argv[i][1]) {
           case 'f':
               in_list = in_crit = 0;
	       if ((i+1) < argc) strcpy(HDFfile, argv[i+1]);
	       i++;
	       break;
           case 'o':
               in_list = in_crit = 0;
	       if ((i+1) < argc) strcpy(outfile, argv[i+1]);
	       i++;
	       break;
           case 'h':
               in_list = in_crit = 0;
	       if ((i+1) < argc) strcpy(HDFoutfile, argv[i+1]);
	       i++;
	       break;
           case 's':
               in_list = in_crit = 0;
	       samples = 1;
	       break;
           case 'a':
               in_list = in_crit = 0;
	       if ((i+1) < argc) aggret = parse_agg_arg(argv[i+1]);
	       if (aggret == 0) aggreg = 1;
	       i++;
	       break;
           case 'b':
               in_list = in_crit = 0;
	       if ((i+1) < argc) strcpy(binfile, argv[i+1]);
	       i++;
	       break;
          case 'v':
               in_list = in_crit = 0;
	       VERBOSE = 1;
	       break;
          case 'r':
               in_list = in_crit = 0;
	       ref_SDS = -1;  /* this tells rest of program to look into ref_string */
	       if ((i+1) < argc) strcpy(ref_string, argv[i+1]);
	       i++;
	       break;
          case 'c':
               in_list = in_crit = 0;
	       COLUMNS = 1;
	       break;
          case 'p':
/*               in_list = in_crit = 0;
	       if (initializepar(argv[i+1],&project) == -1) {
	          printf("\tWarning: error reading parameter file %s, will not use;\n",argv[i+1]);
	          printf("\twill assume sites are \"UL_Row, UL_Col, #_rows, #_cols\"\n");
	       }
	       else 
	          proj_ok = 1;  / * global */
	       break;
          case 'x':
               in_list = in_crit = 0;
	       except = 1;
	       break;
	   case 'l':
	       in_list = 1;
	       in_crit = 0;
	       break;
	   case 'w':
	       in_list = 0;
	       in_crit = 1;
	       break;
	   default:
               break;
                          }
     if (in_list)
	 strcpy(listargs[list_idx++], argv[i]);
     else if (in_crit)
	 strcpy(critargs[crit_idx++], argv[i]);
                     }
	 
if (!strcmp(HDFfile,"")) {
    printf("Error: HDF file must be specified\n");
    return(-1);
                         }		     
/* Process SDS list arguments 
 * Note: SDS names cannot be verified here since 
 * the HDF file has not yet been opened -- must 
 * do that in main.
 ******************************/
if (list_idx < 2) USE_ALL_SDSs = 1;
else if (list_idx == 2) {
    /* Assume listargs[1] is an external file with a list of SDSs to process; try to open it
       and read it. */
    if ( (fdExternSDSlist = fopen(listargs[1], "rb")) != NULL) {
        printf("Reading %s\n", listargs[1]);

        while ((fscanf(fdExternSDSlist, "%s", SDSnames[n_SDSs]) == 1) && (n_SDSs < MAXBANDS))
            n_SDSs++;

        if ((n_SDSs == MAXBANDS)&&(!feof(fdExternSDSlist))) {
            printf("Warning: list of SDSs to process truncated, only %d can be processed.\n", MAXBANDS);
                                                            }

        if (n_SDSs == 0) {
            printf("Error: something wrong with file '%s'\n", listargs[1]);
            ret = -1;
                         }
	fclose(fdExternSDSlist);
    }  /* Else take it as the name of an SDS in the HDF file */
    else  strcpy(SDSnames[n_SDSs++], arg_clean(listargs[1]));
                        }
else {   /* list_idx > 2, take listarg as a list of SDSs in the HDF file */
    if (list_idx > (MAXBANDS+2)) {
       printf("Warning: list of SDSs to process truncated, only %d can be processed.\n", MAXBANDS);
       list_idx = (MAXBANDS+2);
    }
    
    for (i=1;i<list_idx;i++) strcpy(SDSnames[n_SDSs++], arg_clean(listargs[i]));
     }
if (except) USE_ALL_SDSs = 1;  /* and then subtract out the ones specified in the list -- in main */

/* Process criteria arguments 
 ******************************/
if (crit_idx < 2) {
    printf("Error: no selection criteria entered\n");
    ret = -1;
                  }
else {
    /* Make sure one criterion is areal (four numbers, no strings) */
    for (i=1;i<crit_idx;i++) {
        strcpy(critargs[i], arg_clean(critargs[i]));
        if ( areal_or_boolean(critargs[i]) ) {
	   areal_encountered++;
	   if (overall_set == 0) {
              if ((overallULR < 0.0) || (overallULR > ULR[n_sites-1])) overallULR = ULR[n_sites-1];
              if ((overallULC < 0.0) || (overallULC > ULC[n_sites-1])) overallULC = ULC[n_sites-1];
              if ((overallNR < 0.0) || (overallNR < (ULR[n_sites-1]+NR[n_sites-1]))) overallNR = (ULR[n_sites-1]+NR[n_sites-1]);
              if ((overallNC < 0.0) || (overallNC < (ULC[n_sites-1]+NC[n_sites-1]))) overallNC = (ULC[n_sites-1]+NC[n_sites-1]);
	      }
	   }
	}
    if ( areal_encountered == 0) {
        printf("Error: no area criterion entered\n");
        ret = -1;
        }
     }

return(ret);
}

/***************************************************************************************/

int parse_agg_arg(/*int *agg_dim1, int *agg_dim2, THESE ARE GLOBAL ANYWAY */char *arg)   
{
int i, j, k, l, n, n_comma;
char subarg1[10];
char subarg2[10];
char subarg3[10];

n = strlen(arg);
if (n > 20) return(-1);

n_comma = j = k = l = 0;
for(i=0;i<n;i++) {
   if (arg[i] == ',') {
      n_comma++;
      if (n_comma == 1) subarg1[j] = '\0';
      if (n_comma == 2) subarg2[k] = '\0';
      }
   else {
     if (n_comma == 0) subarg1[j++] = arg[i];
     if (n_comma == 1) subarg2[k++] = arg[i];
     if (n_comma == 2) subarg3[l++] = arg[i];
      }      
   }
  
if ((n_comma == 0)||(n == 0)) return(-1);

agg_dim1 = atoi(subarg1);  /* these are global */
agg_dim2 = atoi(subarg2);
agg_thresh = atoi(subarg3);

return(0);
}

/***************************************************************************************/

short areal_or_boolean(char *string)
{  /* 'string' can be either boolean ('something' LT|GT|BT|AND|EQ value(s)) or areal */
#define MAX_SUBARG 40
char subarg[MAX_SUBARG][MAXLENGTH];
short n_subarg = 0, subargln = 0;
short subarglen[MAX_SUBARG];
short i, j, strlength, arglen;
double hi,lo;
short BOOLEAN = 0;
float ULr, ULc, nr, nc;
double parse_bitstring(char *str);

/* Parse 'string' into subarguments */
for (i=0;i<MAX_SUBARG;i++) {
   strcpy(subarg[i], "\0");
   subarglen[i] = 0;
   }
strlength = strlen(string);
for (i=0;i<strlength;i++) {
     if ((string[i] == ' ')&&(subargln > 0)) { 
        subarg[n_subarg][subargln] = '\0';
        subarglen[n_subarg] = subargln;
        subargln = 0;
        n_subarg++;
	if (n_subarg == MAX_SUBARG) break;
	}
     if (string[i] != ' ') subarg[n_subarg][subargln++] = string[i];
     }
subarg[n_subarg][subargln] = '\0';
     
/* Determine if string is Boolean, and if so, which subargument is the LT|GT|BT|EQ */
for (i=1;i<n_subarg;i++) {   /* the 0th argument CANNOT be boolean */
   if ((!strcmp(subarg[i], "gt")) || (!strcmp(subarg[i], "GT")) ||
       (!strcmp(subarg[i], "lt")) || (!strcmp(subarg[i], "LT")) ||
       (!strcmp(subarg[i], "bt")) || (!strcmp(subarg[i], "BT")) ||
       (!strcmp(subarg[i], "and")) || (!strcmp(subarg[i], "AND")) ||
       (!strcmp(subarg[i], "not")) || (!strcmp(subarg[i], "NOT")) ||
       (!strcmp(subarg[i], "eq")) || (!strcmp(subarg[i], "EQ")))  BOOLEAN = i;
   
}

if (BOOLEAN) {   /* 'string' is Boolean */
  lo = hi = -999.0;
  if (n_criteria < MAXCRITERIA) {
   if ((BOOLEAN+1) <= n_subarg)  lo = atof(subarg[BOOLEAN+1]);
   if ((BOOLEAN+2) <= n_subarg)  hi = atof(subarg[BOOLEAN+2]);

   arglen = 0;
   for (i=0;i<BOOLEAN;i++) {
       for (j=0;j<subarglen[i];j++) 
           SDStofilter[n_criteria][arglen++] = subarg[i][j];	/* May be either SDS in HDF file or external file */
       if (i != (BOOLEAN-1)) SDStofilter[n_criteria][arglen++] = ' ';
   }
   if ((!strcmp(subarg[BOOLEAN], "gt")) || (!strcmp(subarg[BOOLEAN], "GT"))) {
        really_boolean[n_criteria] = 0;
        lowend[n_criteria] = lo+EPSLN;
        highend[n_criteria++] = (double)(-(FILL_VALUE));
      }
   else if ((!strcmp(subarg[BOOLEAN], "lt")) || (!strcmp(subarg[BOOLEAN], "LT"))) {
        really_boolean[n_criteria] = 0;
        lowend[n_criteria] = (double)(FILL_VALUE);
        highend[n_criteria++] = lo-EPSLN;
      }
   else if ((!strcmp(subarg[BOOLEAN], "eq")) || (!strcmp(subarg[BOOLEAN], "EQ"))) {
        really_boolean[n_criteria] = 0;
        lowend[n_criteria] = lo;
        highend[n_criteria++] = lo;
      }
   else if ((!strcmp(subarg[BOOLEAN], "and")) || (!strcmp(subarg[BOOLEAN], "AND"))) {
        really_boolean[n_criteria] = 1;
	if ( ((BOOLEAN+1) <= n_subarg) && subarg[BOOLEAN+1][0] == 'b' ) {
	   lo = parse_bitstring(subarg[BOOLEAN+1]);
	   }
        /*printf("%f %s\n", lo, subarg[BOOLEAN+1]);*/
        lowend[n_criteria] = lo;
	if ( ((BOOLEAN+2) <= n_subarg) && subarg[BOOLEAN+2][0] == 'b' ) {
	   hi = parse_bitstring(subarg[BOOLEAN+2]);
	   }
        highend[n_criteria++] = hi;
      }
   else if ((!strcmp(subarg[BOOLEAN], "not")) || (!strcmp(subarg[BOOLEAN], "NOT"))) {
        really_boolean[n_criteria] = 2;
	if ( ((BOOLEAN+1) <= n_subarg) && subarg[BOOLEAN+1][0] == 'b' ) {
	   lo = parse_bitstring(subarg[BOOLEAN+1]);
	   }
        lowend[n_criteria] = lo;
	if ( ((BOOLEAN+2) <= n_subarg) && subarg[BOOLEAN+2][0] == 'b' ) {
	   hi = parse_bitstring(subarg[BOOLEAN+2]);
	   }
        highend[n_criteria++] = hi;
      }
   else if ((!strcmp(subarg[BOOLEAN], "bt")) || (!strcmp(subarg[BOOLEAN], "BT"))) {
	if (lo < hi) {
           really_boolean[n_criteria] = 0;
           lowend[n_criteria] = lo;
           highend[n_criteria++] = hi;
	}
	else 
	   printf("\tWarning: '%s' has values out of order or a missing high value, ignored\n", string);
      }
    }  /* if (n_criteria < MAXCRITERIA) */
   else {
      printf("\tWarning: criteria buffer exceeded, '%s' ignored\n", string);
        }
   return(0);
        
   }
else {                      /* 'string' is areal (ULr, ULC, nr, nc) */
   ULr = atof(subarg[0]);
   ULc = atof(subarg[1]);
   nr  = atof(subarg[2]);
   nc  = atof(subarg[3]);
   if (proj_ok == 1) {
     /* 29-FEB-00: Note, projection library is screwed up, LATMIN, LATMAX, LONMIN, LONMAX are unreliable.
        Don't do these comparisons.
     printf("%f %f %f %f\n", project.latmax,project.latmin,project.lonmax,project.lonmin);
     
     if ((ULr > project.latmax)||(ULr < project.latmin)) {
          printf("\tWarning: invalid upper-left latitude value, %s, will set to %f\n",subarg[0], project.latmax);
          ULr = project.latmax;
          }
     if ((ULc < project.lonmin)||(ULc > project.lonmax)) {
          printf("\tWarning: invalid upper-left longitude value, %s, will set to %f\n",subarg[1], project.lonmin);
          ULc = project.lonmin;
          }
     if ((nr > project.latmax)||(nr < project.latmin)) {
          printf("\tWarning: invalid lower-right latitude value, %s, will set to %f\n",subarg[2], project.latmin);
          nr = project.latmin;
          }
     if ((nc < project.lonmin)||(nc > project.lonmax)) {
          printf("\tWarning: invalid lower-right longitude value, %s, will set to %f\n",subarg[3], project.lonmax);
          nc = project.lonmax;
          }
     if ((nc < ULc)||(nr > ULr)) {
          printf("\tWarning: invalid lower-right longitude or latitude values, will set to (%f, %f)\n", project.latmin, project.lonmax);
          nr=project.latmin;
	  nc=project.lonmax;
          }*/
                     }
   else {
     if (ULr < 0) {
          printf("\tWarning: invalid upper-left row value, %s, will set to 0\n",subarg[0]);
          ULr = 0;
          }
     if (ULc < 0) {
          printf("\tWarning: invalid upper-left column value, %s, will set to 0\n",subarg[1]);
          ULc = 0;
          }
	
      if (nr < 0) {
          printf("\tWarning: invalid number of rows requested, %s, will set to 0\n",subarg[2]);
          nr = 0;
          }
     if (nr < 0) {
          printf("\tWarning: invalid number of columns requested, %s, will set to 0\n",subarg[3]);
          nc = 0;
          }
	}

   /* Last sanity check */
   if (((short)ULr == 0) && ((short)ULc == 0) && ((short)nr == 0) && ((short)nc == 0)) return(1);
   ULR[n_sites] = ULr;
   ULC[n_sites] = ULc;
   NR[n_sites] = nr;
   NC[n_sites++] = nc;
   return(1);
     }
}

/*****************************************************************************************************************
 **************  start of parse_bitstring() **********************************************************************
 *****************************************************************************************************************/
 
double parse_bitstring(char *str)
{
int i, j, strln, tmp;
double retval;

strln = strlen(str);
retval = 0.0;
j = 0;

for(i=strln-1;i>0;i--,j++) {  /* > 0, not >= 0, because the zeroth character is 'b' */
   tmp = ( 1 << j );
   if ( str[i] == '0' ) {
       /* do nothing */
      }
   else if ( str[i] == '1' ) {
       retval += (double)tmp;
      }
   else {
       /* do nothing */
      }
   /*printf("%c %d %d : %f\n", str[i], tmp, i, retval);*/
   }
return(retval);

}

/***************************************************************************************/

char *arg_clean(char *string)
{
/* Gets rid of parentheses at beginning and end of 'string', and any embedded commas and quotes */
short i, j, strlength;
return(string);  /* THIS IS A SHORT-CIRCUIT; routine is no longer necessary */

strlength = strlen(string);
for (i=0,j=0;i<strlength;i++)
    if ((string[i] != '(')&&(string[i] != ')')&&(string[i] != ',')&&(string[i] != '\'')) 
       string[j++] = string[i];
string[j] = '\0';
return(string);
}

/***************************************************************************************/

void subsampler(double *SDS_image, int32 SDS_type, void *whatever, short nc, short nr, long k, long sub_k, 
                float subsample_factor)
{
long  xl, *SDS_image_l;
short xs, *SDS_image_s;
int   xi, *SDS_image_i;
char  xc, *SDS_image_c;
unsigned long  ul, *SDS_image_ul;
unsigned short us, *SDS_image_us;
unsigned int   ui, *SDS_image_ui;
unsigned char  uc, *SDS_image_uc;
float  f, *SDS_image_f;
double  d, *SDS_image_d;
int factor;
int what_to_do;
long l, ll;
int i, ii, j, jj;

l = ll = 0L;

if (subsample_factor < 1.0) {  /* then we are going from a lower resolution to a higher resolution */
   what_to_do = 1;
   factor = (int)1.0/subsample_factor;
   }
else  {  
   what_to_do = 0;
   factor = (int)subsample_factor;
   }

if ((SDS_type == DFNT_CHAR8)||(SDS_type == DFNT_INT8)) {
    SDS_image_c = (char *)whatever;
    if (what_to_do == 1) {        
       for(i=0;i<nr;i++) {
           for(j=0;j<nc;j++) {
               xc = SDS_image_c[l];
               for(ii=(i*factor);ii<((i*factor)+factor);ii++) {
                  ll = (long)ii*nc*factor;
                  for(jj=(j*factor);jj<((j*factor)+factor);jj++)  SDS_image[ll+jj] = (double)xc;
	          }
                 l++;
               }
           }
        }
    else {        
       for(i=0;i<nr;i=i+factor) {
           for(j=0;j<nc;j=j+factor) {
                 ll = (long)i*nc + j;
                 xc = SDS_image_c[ll];
                 SDS_image[l++] = (double)xc;
	          }
            }
        }       
   }
else if ((SDS_type == DFNT_UCHAR8)||(SDS_type == DFNT_UINT8)) {
    SDS_image_uc = (unsigned char *)whatever;
    if (what_to_do == 1) {        
       for(i=0;i<nr;i++) {
           for(j=0;j<nc;j++) {
               uc = SDS_image_uc[l];
               for(ii=(i*factor);ii<((i*factor)+factor);ii++) {
                  ll = (long)ii*nc*factor;
                  for(jj=(j*factor);jj<((j*factor)+factor);jj++)  SDS_image[ll+jj] = (double)uc;
	          }
                 l++;
               }
           }
        }
    else {        
       for(i=0;i<nr;i=i+factor) {
           for(j=0;j<nc;j=j+factor) {
                 ll = (long)i*nc + j;
                 uc = SDS_image_uc[ll];
                 SDS_image[l++] = (double)uc;
	          }
            }
        }       
  }
else if (SDS_type == DFNT_INT16) {
    SDS_image_s = (short *)whatever;
    if (what_to_do == 1) {        
       for(i=0;i<nr;i++) {
           for(j=0;j<nc;j++) {
               xs = SDS_image_s[l];
               for(ii=(i*factor);ii<((i*factor)+factor);ii++) {
                  ll = (long)ii*nc*factor;
                  for(jj=(j*factor);jj<((j*factor)+factor);jj++)  SDS_image[ll+jj] = (double)xs;
	          }
                 l++;
               }
           }
        }
    else {        
       for(i=0;i<nr;i=i+factor) {
           for(j=0;j<nc;j=j+factor) {
                 ll = (long)i*nc + j;
                 xs = SDS_image_s[ll];
                 SDS_image[l++] = (double)xs;
	          }
            }
        }       
   }
else if (SDS_type == DFNT_UINT16) {
    SDS_image_us = (unsigned short *)whatever;
    if (what_to_do == 1) {        
       for(i=0;i<nr;i++) {
           for(j=0;j<nc;j++) {
               us = SDS_image_us[l];
               for(ii=(i*factor);ii<((i*factor)+factor);ii++) {
                  ll = (long)ii*nc*factor;
                  for(jj=(j*factor);jj<((j*factor)+factor);jj++)  SDS_image[ll+jj] = (double)us;
	          }
                 l++;
               }
           }
        }
    else {        
       for(i=0;i<nr;i=i+factor) {
           for(j=0;j<nc;j=j+factor) {
                 ll = (long)i*nc + j;
                 us = SDS_image_us[ll];
                 SDS_image[l++] = (double)us;
	          }
            }
        }       
   }
else if (SDS_type == DFNT_INT32) {
    SDS_image_i = (int *)whatever;
    if (what_to_do == 1) {        
       for(i=0;i<nr;i++) {
           for(j=0;j<nc;j++) {
               xi = SDS_image_i[l];
               for(ii=(i*factor);ii<((i*factor)+factor);ii++) {
                  ll = (long)ii*nc*factor;
                  for(jj=(j*factor);jj<((j*factor)+factor);jj++)  SDS_image[ll+jj] = (double)xi;
	          }
                 l++;
               }
           }
        }
    else {        
       for(i=0;i<nr;i=i+factor) {
           for(j=0;j<nc;j=j+factor) {
                 ll = (long)i*nc + j;
                 xi = SDS_image_i[ll];
                 SDS_image[l++] = (double)xi;
	          }
            }
        }       
   }
else if (SDS_type == DFNT_UINT32) {
    SDS_image_ui = (unsigned int *)whatever;
    if (what_to_do == 1) {        
       for(i=0;i<nr;i++) {
           for(j=0;j<nc;j++) {
               ui = SDS_image_ui[l];
               for(ii=(i*factor);ii<((i*factor)+factor);ii++) {
                  ll = (long)ii*nc*factor;
                  for(jj=(j*factor);jj<((j*factor)+factor);jj++)  SDS_image[ll+jj] = (double)ui;
	          }
                 l++;
               }
           }
        }
    else {        
       for(i=0;i<nr;i=i+factor) {
           for(j=0;j<nc;j=j+factor) {
                 ll = (long)i*nc + j;
                 ui = SDS_image_ui[ll];
                 SDS_image[l++] = (double)ui;
	          }
            }
        }       
   }
else if (SDS_type == DFNT_INT64) {
    SDS_image_l = (long *)whatever;
    if (what_to_do == 1) {        
       for(i=0;i<nr;i++) {
           for(j=0;j<nc;j++) {
               xl = SDS_image_l[l];
               for(ii=(i*factor);ii<((i*factor)+factor);ii++) {
                  ll = (long)ii*nc*factor;
                  for(jj=(j*factor);jj<((j*factor)+factor);jj++)  SDS_image[ll+jj] = (double)xl;
	          }
                 l++;
               }
           }
        }
    else {        
       for(i=0;i<nr;i=i+factor) {
           for(j=0;j<nc;j=j+factor) {
                 ll = (long)i*nc + j;
                 xl = SDS_image_l[ll];
                 SDS_image[l++] = (double)xl;
	          }
            }
        }       
   }
else if (SDS_type == DFNT_UINT64) {
    SDS_image_ul = (unsigned long *)whatever;
    if (what_to_do == 1) {        
       for(i=0;i<nr;i++) {
           for(j=0;j<nc;j++) {
               ul = SDS_image_ul[l];
               for(ii=(i*factor);ii<((i*factor)+factor);ii++) {
                  ll = (long)ii*nc*factor;
                  for(jj=(j*factor);jj<((j*factor)+factor);jj++)  SDS_image[ll+jj] = (double)ul;
	          }
                 l++;
               }
           }
        }
    else {        
       for(i=0;i<nr;i=i+factor) {
           for(j=0;j<nc;j=j+factor) {
                 ll = (long)i*nc + j;
                 ul = SDS_image_ul[ll];
                 SDS_image[l++] = (double)ul;
	          }
            }
        }       
   }
else if ((SDS_type == DFNT_FLOAT)||(SDS_type == DFNT_FLOAT32)) {
    SDS_image_f = (float *)whatever;
    if (what_to_do == 1) {  /* then we are going from a lower resolution to a higher resolution */
       for(i=0;i<nr;i++) {
           for(j=0;j<nc;j++) {
               f = SDS_image_f[l];
               for(ii=(i*factor);ii<((i*factor)+factor);ii++) {
                  ll = (long)ii*nc*factor;
                  for(jj=(j*factor);jj<((j*factor)+factor);jj++)  SDS_image[ll+jj] = (double)f;
	          }
                 l++;
               }
           }
        }
    else {                         /* then we are going from a higher resolution to a lower resolution */
       for(i=0;i<nr;i=i+factor) {
           for(j=0;j<nc;j=j+factor) {
                 ll = (long)i*nc + j;
                 f = SDS_image_f[ll];
                 SDS_image[l++] = (double)f;
	          }
            }
        }       
   }
else if (SDS_type == DFNT_DOUBLE) {
    SDS_image_d = (double *)whatever;
    if (what_to_do == 1) {        
       for(i=0;i<nr;i++) {
           for(j=0;j<nc;j++) {
               d = SDS_image_d[l];
               for(ii=(i*factor);ii<((i*factor)+factor);ii++) {
                  ll = (long)ii*nc*factor;
                  for(jj=(j*factor);jj<((j*factor)+factor);jj++)  SDS_image[ll+jj] = d;
	          }
                 l++;
               }
           }
        }
    else {        
       for(i=0;i<nr;i=i+factor) {
           for(j=0;j<nc;j=j+factor) {
                 ll = (long)i*nc + j;
                 d = SDS_image_d[ll];
                 SDS_image[l++] = d;
	          }
            }
        }       
   }
return;   
}		

/***************************************************************************************/

double *get_images(int32 file_id, char *SDSnames, int32 offset, int32 which_SDS, short SDS_type, FILE *fd,
                   short n_cols, long k, short ULr, short ULc, short nr, short nc, float subsample_factor)
{
short SDSread(int32 sd_id, int32 offset, int32 SDS_index, void *image, short, short, short, short);
short ext_image(void *image, short size, short ULr, short ULc, 
                short nr, short nc, short n_cols, FILE *fd);
void subsampler(double *SDS_image, int32 SDS_type, void *SDS_image_f, short nc, short nr, long k, long sub_k, 
                float subsample_factor);
long   *SDS_image_l;
short  *SDS_image_s;
int    *SDS_image_i;
char   *SDS_image_c;
unsigned long   *SDS_image_ul;
unsigned short  *SDS_image_us;
unsigned int    *SDS_image_ui;
unsigned char   *SDS_image_uc;
float  *SDS_image_f;
double *SDS_image_d;
double *SDS_image;
long li;
char filename[50];
char extfilename[50];
FILE *fdwr;
int subs = 0;
long subs_k;


#if defined  (DEBUG2)         
sprintf(filename, "SDSreader_debug_%s", SDSnames);
sprintf(extfilename, "SDSreader_debug_ext_%s", SDSnames);
#endif

SDS_image = (double *)NULL;    

if ((SDS_image = (double *)malloc(k*sizeof(double))) == (double *)NULL) {
   printf("get_images, error: out of memory, reading '%s'\n", SDSnames);
   return(SDS_image);
   }

subs_k = k;
if ( fabs(subsample_factor - 1.00000) > 1e-3) {   /* then we are in fact either subsampling or interpolating */
   subs = 1;
   subs_k = (long)((float)k*(subsample_factor*subsample_factor));
   ULr = (short)((float)ULr*subsample_factor);
   ULc = (short)((float)ULc*subsample_factor);
   nr = (short)((float)nr*subsample_factor);
   nc = (short)((float)nc*subsample_factor);
   }
    
SDS_image_s = (short *)NULL;    
    
/* allocate memory for any type; then cast it to double */
if ((SDS_type == DFNT_CHAR8)||(SDS_type == DFNT_INT8)) {
        subs_k *= sizeof(char);
        if ((SDS_image_c = (char *)malloc(subs_k)) == NULL) {
           printf("get_images, error: out of memory, reading '%s'\n", SDSnames);
	   SDS_image = (double *)NULL;
           return(SDS_image);
        }
        if (which_SDS == -1) {    /* an external file */
            if (ext_image(SDS_image_c, SDS_type, ULr, ULc, nr, nc, n_cols, fd) != 1) {
                printf("get_images, error: reading '%s'\n", SDSnames);
                exit(-6);
	    }
                             }
        else {                    /* an SDS */
            SDSread(file_id, offset, which_SDS, (void *)SDS_image_c, ULr, ULc, nr, nc);
#if defined  (DEBUG2)         
            fwrite(SDS_image_c, k, 1, fopen(filename,"wb"));
#endif
             }
        subs_k /= sizeof(char);
	if (subs != 1) {
	   for (li=0L;li<k;li++) 
	      SDS_image[li] = (double)SDS_image_c[li];
	      }
	else {
	    subsampler(SDS_image, SDS_type, (void *)SDS_image_c, nc, nr, k, subs_k, subsample_factor);
	    }
	free(SDS_image_c);
	}
    else if ((SDS_type == DFNT_UCHAR8)||(SDS_type == DFNT_UINT8)) {
        subs_k *= sizeof(unsigned char);
        if ((SDS_image_uc = (unsigned char *)malloc(subs_k)) == NULL) {
           printf("get_images, error: out of memory, reading '%s'\n", SDSnames);
	   SDS_image = (double *)NULL;
           return(SDS_image);
        }
        if (which_SDS == -1) {    /* an external file */
            if (ext_image(SDS_image_uc, SDS_type, ULr, ULc, nr, nc, n_cols, fd) != 1) {
                printf("get_images, error: reading '%s'\n", SDSnames);
                exit(-6);
	    }
                             }
        else {                     /* an SDS */
            SDSread(file_id, offset, which_SDS, (void *)SDS_image_uc, ULr, ULc, nr, nc);
#if defined  (DEBUG2)       
          fwrite(SDS_image_uc, k, 1, fopen(filename,"wb"));
#endif
             }
        subs_k /= sizeof(unsigned char);
	if (subs != 1) {
	   for (li=0L;li<k;li++) 
	      SDS_image[li] = (double)SDS_image_uc[li];
	      }
	else {
	    subsampler(SDS_image, SDS_type, (void *)SDS_image_uc, nc, nr, k, subs_k, subsample_factor);
	    }
	free(SDS_image_uc);
	}
    else if (SDS_type == DFNT_INT16) {
        subs_k *= sizeof(short);
        if ((SDS_image_s = (short *)malloc(subs_k)) == (short *)NULL) {
           printf("get_images, error: out of memory, reading '%s'\n", SDSnames);
	   SDS_image = (double *)NULL;
           return(SDS_image);
        }

        if (which_SDS == -1) {    /* an external file */
            if (ext_image(SDS_image_s, SDS_type, ULr, ULc, nr, nc, n_cols, fd) != 1) {
                printf("get_images, error: reading '%s'\n", SDSnames);
                exit(-6);
	    }
#if defined (DEBUG2)
 	    fdwr = fopen(extfilename,"wb");
	    fwrite(SDS_image_s, k, 1, fdwr);
	    fclose(fdwr);
#endif	       
                             }
        else {                     /* an SDS */
	   /*printf("%p %p %d %d %d %d %d %d\n", SDS_image, SDS_image_s, which_SDS, offset, ULr, ULc, nr, nc);*/
           SDSread(file_id, offset, which_SDS, (void *)SDS_image_s, ULr, ULc, nr, nc);
#if defined  (DEBUG2)         
           fwrite(SDS_image_s, k, 1, fopen(filename,"wb"));
#endif
             }
        subs_k /= sizeof(short);
	if (subs != 1) {
	   for (li=0L;li<k;li++) 
	      SDS_image[li] = (double)(SDS_image_s[li]);
	      }
	else {
	    subsampler(SDS_image, SDS_type, (void *)SDS_image_s, nc, nr, k, subs_k, subsample_factor);
	    }
	free(SDS_image_s);
	SDS_image_s = (short *)NULL;
	}
    else if (SDS_type == DFNT_UINT16) {
        subs_k *= sizeof(unsigned short);
        if ((SDS_image_us = (unsigned short *)malloc(subs_k)) == NULL) {
           printf("get_images, error: out of memory, reading '%s'\n", SDSnames);
	   SDS_image = (double *)NULL;
           return(SDS_image);
        }
        if (which_SDS == -1) {    /* an external file */
            if (ext_image(SDS_image_us, SDS_type, ULr, ULc, nr, nc, n_cols, fd) != 1) {
                printf("get_images, error: reading '%s'\n", SDSnames);
                exit(-6);
	    }
                             }
        else {                     /* an SDS */
            SDSread(file_id, offset, which_SDS, (void *)SDS_image_us, ULr, ULc, nr, nc);
#if defined  (DEBUG2)         
           fwrite(SDS_image_us, k, 1, fopen(filename,"wb"));
#endif
             }
        subs_k /= sizeof(unsigned short);
	if (subs != 1) {
	   for (li=0L;li<k;li++) 
	      SDS_image[li] = (double)SDS_image_us[li];
	      }
	else {
	    subsampler(SDS_image, SDS_type, (void *)SDS_image_us, nc, nr, k, subs_k, subsample_factor);
	    }
	free(SDS_image_us);
	}
    else if (SDS_type == DFNT_INT32) {
        subs_k *= sizeof(int);
        if ((SDS_image_i = (int *)malloc(subs_k)) == NULL) {
           printf("get_images, error: out of memory, reading '%s'\n", SDSnames);
	   SDS_image = (double *)NULL;
           return(SDS_image);
        }
        if (which_SDS == -1) {    /* an external file */
            if (ext_image(SDS_image_i, SDS_type, ULr, ULc, nr, nc, n_cols, fd) != 1) {
                printf("get_images, error: reading '%s'\n", SDSnames);
                exit(-6);
	    }
                             }
        else {                     /* an SDS */
            SDSread(file_id, offset, which_SDS, (void *)SDS_image_i, ULr, ULc, nr, nc);
#if defined  (DEBUG2)         
            fwrite(SDS_image_i, k, 1, fopen(filename,"wb"));
#endif
             }
        subs_k /= sizeof(int);
	if (subs != 1) {
	   for (li=0L;li<k;li++) 
	      SDS_image[li] = (double)SDS_image_i[li];
	      }
	else {
	    subsampler(SDS_image, SDS_type, (void *)SDS_image_i, nc, nr, k, subs_k, subsample_factor);
	    }
	free(SDS_image_i);
	}
    else if (SDS_type == DFNT_UINT32) {
        subs_k *= sizeof(unsigned int);
        if ((SDS_image_ui = (unsigned int *)malloc(subs_k)) == NULL) {
           printf("get_images, error: out of memory, reading '%s'\n", SDSnames);
	   SDS_image = (double *)NULL;
           return(SDS_image);
        }
        if (which_SDS == -1) {    /* an external file */
            if (ext_image(SDS_image_ui, SDS_type, ULr, ULc, nr, nc, n_cols, fd) != 1) {
                printf("get_images, error: reading '%s'\n", SDSnames);
                exit(-6);
	    }
                             }
        else {                     /* an SDS */
            SDSread(file_id, offset, which_SDS, (void *)SDS_image_ui, ULr, ULc, nr, nc);
#if defined  (DEBUG2)         
            fwrite(SDS_image_ui, k, 1, fopen(filename,"wb"));
#endif
             }
        subs_k /= sizeof(unsigned int);
	if (subs != 1) {
	   for (li=0L;li<k;li++) 
	      SDS_image[li] = (double)SDS_image_ui[li];
	      }
	else {
	    subsampler(SDS_image, SDS_type, (void *)SDS_image_ui, nc, nr, k, subs_k, subsample_factor);
	    }
	free(SDS_image_ui);
	}
    else if (SDS_type == DFNT_INT64) {
        subs_k *= sizeof(long);
        if ((SDS_image_l = (long *)malloc(subs_k)) == NULL) {
           printf("get_images, error: out of memory, reading '%s'\n", SDSnames);
	   SDS_image = (double *)NULL;
           return(SDS_image);
        }
        if (which_SDS == -1) {    /* an external file */
            if (ext_image(SDS_image_l, SDS_type, ULr, ULc, nr, nc, n_cols, fd) != 1) {
                printf("get_images, error: reading '%s'\n", SDSnames);
                exit(-6);
	    }
                             }
        else {                     /* an SDS */
            SDSread(file_id, offset, which_SDS, (void *)SDS_image_l, ULr, ULc, nr, nc);
#if defined  (DEBUG2)         
            fwrite(SDS_image_l, k, 1, fopen(filename,"wb"));
#endif
             }
        subs_k /= sizeof(long);
	if (subs != 1) {
	   for (li=0L;li<k;li++) 
	      SDS_image[li] = (double)SDS_image_l[li];
	      }
	else {
	    subsampler(SDS_image, SDS_type, (void *)SDS_image_l, nc, nr, k, subs_k, subsample_factor);
	    }
	free(SDS_image_l);
	}
    else if (SDS_type == DFNT_UINT64) {
        subs_k *= sizeof(unsigned long);
        if ((SDS_image_ul = (unsigned long *)malloc(subs_k)) == NULL) {
           printf("get_images, error: out of memory, reading '%s'\n", SDSnames);
	   SDS_image = (double *)NULL;
           return(SDS_image);
        }
        if (which_SDS == -1) {    /* an external file */
            if (ext_image(SDS_image_ul, SDS_type, ULr, ULc, nr, nc, n_cols, fd) != 1) {
                printf("get_images, error: reading '%s'\n", SDSnames);
                exit(-6);
	    }
                             }
        else {                     /* an SDS */
            SDSread(file_id, offset, which_SDS, (void *)SDS_image_ul, ULr, ULc, nr, nc);
#if defined  (DEBUG2)         
            fwrite(SDS_image_ul, k, 1, fopen(filename,"wb"));
#endif
             }
        subs_k /= sizeof(unsigned long);
	if (subs != 1) {
	   for (li=0L;li<k;li++) 
	      SDS_image[li] = (double)SDS_image_ul[li];
	      }
	else {
	    subsampler(SDS_image, SDS_type, (void *)SDS_image_ul, nc, nr, k, subs_k, subsample_factor);
	    }
	free(SDS_image_ul);
	}
    else if (SDS_type == DFNT_FLOAT) {
        subs_k *= sizeof(float);
        if ((SDS_image_f = (float *)malloc(subs_k)) == NULL) {
           printf("get_images, error: out of memory, reading '%s'\n", SDSnames);
	   SDS_image = (double *)NULL;
           return(SDS_image);
        }
        if (which_SDS == -1) {    /* an external file */
            if (ext_image(SDS_image_f, SDS_type, ULr, ULc, nr, nc, n_cols, fd) != 1) {
                printf("get_images, error: reading '%s'\n", SDSnames);
                exit(-6);
	    }
                             }
        else {                     /* an SDS */
            SDSread(file_id, offset, which_SDS, (void *)SDS_image_f, ULr, ULc, nr, nc);
#if defined  (DEBUG2)      
            fwrite(SDS_image_f, k, 1, fopen(filename,"wb"));
#endif
             }
        subs_k /= sizeof(float);
	if (subs != 1) {
	   for (li=0L;li<k;li++) 
	      SDS_image[li] = (double)SDS_image_f[li];
	      }
	else {
	    subsampler(SDS_image, SDS_type, (void *)SDS_image_f, nc, nr, k, subs_k, subsample_factor);
	    }
		
	free(SDS_image_f);
	}
    else if (SDS_type == DFNT_DOUBLE) {
        subs_k *= sizeof(double);
        if ((SDS_image_d = (double *)malloc(subs_k)) == NULL) {
           printf("get_images, error: out of memory, reading '%s'\n", SDSnames);
	   SDS_image = (double *)NULL;
           return(SDS_image);
        }
        if (which_SDS == -1) {    /* an external file */
            if (ext_image(SDS_image_d, SDS_type, ULr, ULc, nr, nc, n_cols, fd) != 1) {
                printf("get_images,error: reading '%s'\n", SDSnames);
                exit(-6);
	    }
                            }
        else {                     /* an SDS */
            SDSread(file_id, offset, which_SDS, (void *)SDS_image_d, ULr, ULc, nr, nc);
#if defined  (DEBUG2)         
            fwrite(SDS_image, k, 1, fopen(filename,"wb"));
#endif
             }
        subs_k /= sizeof(double);
	if (subs != 1) {
	   for (li=0L;li<k;li++) 
	      SDS_image[li] = SDS_image_d[li];
	      }
	else {
	    subsampler(SDS_image, SDS_type, (void *)SDS_image_d, nc, nr, k, subs_k, subsample_factor);
	    }
		
	free(SDS_image_d);
        } 
return(SDS_image);
}

/***************************************************************************************/

short ext_image(void *image, short size, short ULr, short ULc, short nr, 
                short nc, short n_cols, FILE *fd)
{  /* Extracts image of nr rows by nc cols from fd into image, 
      assuming image is already big enough. */
long k;
short i;
short ret = 1;
char *c;
short *s;
int *n;
long *l;
unsigned char *uc;
unsigned short *us;
unsigned int *un;
unsigned long *ul;
float *f;
double *d;

k = ULr*n_cols+ULc; 
rewind(fd);

if ((size == DFNT_CHAR)||(size == DFNT_INT8)) {    /* are these equivalent? */
   fseek(fd,k,SEEK_CUR);
   c = (char *)image;
   for (i=0;i<nr;i++) {
       if (fread(c, nc, 1, fd) != 1) ret = 0;
       fseek(fd,(n_cols-nc),SEEK_CUR);
       c += nc;
       }
    }
else if ((size == DFNT_UCHAR)||(size == DFNT_UINT8)) {  /* same */
   fseek(fd,k,SEEK_CUR);
   uc = (unsigned char *)image;
   for (i=0;i<nr;i++) {
       if (fread(uc, nc, 1, fd) != 1) ret = 0; 
       fseek(fd,(n_cols-nc),SEEK_CUR);
       uc += nc;
       }
    }
else if (size == DFNT_INT16) {
   fseek(fd,k*sizeof(short),SEEK_CUR);
   s = (short *)image;
   for (i=0;i<nr;i++) { 
       if (fread(s, nc*sizeof(short), 1, fd) != 1) ret = 0; 
       fseek(fd,(n_cols-nc)*sizeof(short),SEEK_CUR);
       s += nc;
       }
    }
else if (size == DFNT_UINT16) {
   fseek(fd,k*sizeof(short),SEEK_CUR);
   us = (unsigned short *)image;
   for (i=0;i<nr;i++) {
       if (fread(us, nc*sizeof(short), 1, fd) != 1) ret = 0; 
       fseek(fd,(n_cols-nc)*sizeof(short),SEEK_CUR);
       us += nc;
       }
    }
else if (size == DFNT_INT32) {
   fseek(fd,k*sizeof(int),SEEK_CUR);
   n = (int *)image;
   for (i=0;i<nr;i++) {
       if (fread(n, nc*sizeof(int), 1, fd) != 1) ret = 0; 
       fseek(fd,(n_cols-nc)*sizeof(int),SEEK_CUR);
       n += nc;
       }
    }
else if (size == DFNT_UINT32) {
   fseek(fd,k*sizeof(int),SEEK_CUR);
   un = (unsigned int *)image;
   for (i=0;i<nr;i++) {
       if (fread(un, nc*sizeof(int), 1, fd) != 1) ret = 0; 
       fseek(fd,(n_cols-nc)*sizeof(int),SEEK_CUR);
       un += nc;
       }
    }
else if (size == DFNT_INT64) {
   fseek(fd,k*sizeof(long),SEEK_CUR);
   l = (long *)image;
   for (i=0;i<nr;i++) {
       if (fread(l, nc*sizeof(long), 1, fd) != 1) ret = 0; 
       fseek(fd,(n_cols-nc)*sizeof(long),SEEK_CUR);
       l += nc;
       }
    }
else if (size == DFNT_UINT64) {
   fseek(fd,k*sizeof(long),SEEK_CUR);
   ul = (unsigned long *)image;
   for (i=0;i<nr;i++) {
       if (fread(ul, nc*sizeof(long), 1, fd) != 1) ret = 0; 
       fseek(fd,(n_cols-nc)*sizeof(long),SEEK_CUR);
       ul += nc;
       }
    }
else if (size == DFNT_FLOAT) {
   fseek(fd,k*sizeof(float),SEEK_CUR);
   f = (float *)image;
   for (i=0;i<nr;i++) {
       if (fread(f, nc*sizeof(float), 1, fd) != 1) ret = 0; 
       fseek(fd,(n_cols-nc)*sizeof(float),SEEK_CUR);
       f += nc;
       }
    }
else if (size == DFNT_DOUBLE) {
   fseek(fd,k*sizeof(double),SEEK_CUR);
   d = (double *)image;
   for (i=0;i<nr;i++) {
       if (fread(d, nc*sizeof(double), 1, fd) != 1) ret = 0; 
       fseek(fd,(n_cols-nc)*sizeof(double),SEEK_CUR);
       d += nc;
       }
    }
return ret;
}

/***** New, 20-SEP-01 ****************************************/

void sun_angles (short jday,float gmt,double dlat,double dlon,
                   float *ts,float *fs) {
      
      double mst,tst,tet,et,ha,delta;
      double amuzero,elev,az,caz,azim;
      
      double A1=.000075,A2=.001868,A3=.032077,A4=.014615,A5=.040849;
      double B1=.006918,B2=.399912,B3=.070257,B4=.006758;
      double B5=.000907,B6=.002697,B7=.001480;
      
                        
/*      
      SOLAR POSITION (ZENITHAL ANGLE ThetaS,AZIMUTHAL ANGLE PhiS IN DEGREES)
      J IS THE DAY NUMBER IN THE YEAR
 
     MEAN SOLAR TIME (HEURE DECIMALE)
*/
      mst=gmt+(dlon/ROVRDG)/15.;
      tet=2.*M_PI*(double)jday/365.;

/*    TIME EQUATION (IN MN.DEC) */
      et=A1+A2*cos(tet)-A3*sin(tet)-A4*cos(2.*tet)-A5*sin(2.*tet);
      et=et*12.*60./M_PI;

/*     TRUE SOLAR TIME */

      tst=mst+et/60.;
      tst=(tst-12.);

/*     HOUR ANGLE */

      ha=tst*15.*ROVRDG;

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
      elev=elev/ROVRDG;

/*     CONVERSION IN DEGREES */

      *ts=90.-elev;
      *fs=azim/ROVRDG;
      return;
}

void sat_angles(double dlat,double dlon,double dheight,double dlatsub,
                   double dlonsub,float *tv,float *fv) {
       
      double psi,cosfv,sinfv;
      
      double EPS=0.0005;  
      double Re= 6367.5*1000.0 ;
          

      psi=acos(sin(dlat)*sin(dlatsub)+cos(dlat)*cos(dlatsub)*cos(dlonsub-dlon));
      if (fabs(psi) > EPS) {
        *tv=fabs(psi+atan(sin(psi)/(dheight/Re+1.-cos(psi))));
        cosfv=(sin(dlatsub)-sin(dlat)*cos(psi))/(sin(psi)*cos(dlat));
        sinfv=(cos(dlatsub)*sin(dlonsub-dlon))/sin(psi);
        *fv=atan2(sinfv,cosfv);
        if (*fv < 0.) 
         *fv=*fv+2.*M_PI;
        *tv=*tv/ROVRDG;
        *fv=*fv/ROVRDG;
      } else {
        *tv=0.;
        *fv=0.;
      }
      return;
}

void metareader(int32 sd_id, char *type_of_meta, char *metastring, 
                int32 *count, char *data)
{
#define MAXATTR 2
int32 i, ii, j, n_attr1, n_sets1, count1, number_type1, n_val, counts[MAXATTR];
int start,within,wasjustobject;
int obj_offset[10],n_obj;
char attrib[MAXATTR][MAXLENGTH], attrib1[MAXLENGTH];
char *charattr, *charatt_array[MAXATTR];
char line[MAXLENGTH], objs[MAXLENGTH];
char lhs[MAXLENGTH], rhs[MAXLENGTH];
void get_a_line(char *text, int *start, char *line);

*count = 0;

/* First of all, look for all metadata with names that
   match type_of_meta.  */
   
SDfileinfo(sd_id, &n_sets1, &n_attr1);
if (n_attr1 == 0) return;


ii = 0;
counts[0] = counts[1] = -1;
for (j=0;j<n_attr1;j++) {
     
     SDattrinfo(sd_id, j, attrib1, &number_type1, &count1);
     n_val = strlen(attrib1);
     for (i=0;i<n_val;i++) attrib1[i] = toupper(attrib1[i]);  
          
     if (strstr(attrib1,type_of_meta)) {
       counts[ii] = count1;  
       attrib[ii][0] = '\0';  
       for (i=0;i<n_val;i++) attrib[ii][i] = attrib1[i];  
       attrib[ii][i] = '\0';  
       if ((charatt_array[ii] = (char *)malloc((count1+10)*sizeof(char))) == NULL) {
	   printf("Out of memory, array 'charatt_array'\n");
           return;
	                                                                   }
        SDreadattr(sd_id, j, charatt_array[ii]);  
	strcat(charatt_array[ii], "\0");
	
	/* Please note: at times, counts[ii]+10 is still less than strlen(charattr_array[ii]).  I don't know 
	   if this is a problem with the HDF library or the kernel.  Reset counts[ii] like this
	   to be on the safe side (at least for now...)  */
	   
	counts[ii] = strlen(charatt_array[ii]); 
	 
	/* If this *isn't* done, on LTDRI3 (but not LTDRI1) we get run-time errors where the C-library
	   gives me this kind of shit (while in get_images()):
	   
	   *** glibc detected *** ericv: free(): invalid next size (fast): 0x0000000015752330 ***
	  
	   presumably because subsequent parsing of the contents of charatt_array[ii] (below)
	   require that counts[ii] has the right value, and if it doesn't, the heap gets corrupted.
	   */
	/*printf ("ii = %d, j = %d, %s: %d %d\n", ii, j, attrib1, count1,strlen(charatt_array[0]));*/
	if (ii < MAXATTR) ii++; 
	}

                        }  /* for (j=0;j<n_attr1;j++) */

/*printf ("ii = %d %d %d\n", ii, counts[0],strlen(charatt_array[0]));*/
			
/* If there are two metadata sets that match (e. g., 
   CoreMetadata and CoreMetadata_1), merge them into one. */
if (ii == 2) {
    if ((counts[0] < 0) || (counts[1] < 0)) return;
    if (VERBOSE) printf("Merging %s (%d) and %s (%d)\n", attrib[0], counts[0], attrib[1], counts[1]);
    if ((charattr = (char *)malloc((counts[0]+counts[1]+10)*sizeof(char))) == NULL) {
	printf("Out of memory, array 'charattr'\n");
        return;
                                                                       }
    charattr[0] = '\0';  
    strcpy(charattr, charatt_array[0]);
    strcat(charattr, charatt_array[1]);
    strcat(charattr, "\0");
    
   }
else if (ii == 1) {
    if (counts[0] < 0) return;
    if ((charattr = (char *)malloc((counts[0]+10)*sizeof(char))) == NULL) {
	printf("Out of memory, array 'charattr'\n");
        return;
                                                                       }
    charattr[0] = '\0';  
    strcpy(charattr, charatt_array[0]);
    strcat(charattr, "\0");
   }
else {
   return;
   }
  

/* Now parse the beast */		       
start = 0;
/*if (strstr(attrib,type_of_meta)) {*/
	 n_obj=0;
	 wasjustobject=0;
	 objs[0] = '\0';
	 do { 
	      line[0] = '\0';
	      get_a_line(charattr, &start, line);
	      
	      /* Get rid of whitespace, if not inside "..." characters. */
	      within=0;
              for (i=0,ii=0;i<strlen(line);i++) {
	          if (line[i] == '"') {
		     if (!within) within=1;
		     else if (within) within=0;
		     }
	          if (within) line[ii++] = line[i]; 
		  else {
		     if (line[i] != ' ') line[ii++] = line[i];  
		     }
		  }
	      line[ii] = '\0';

	      /* Get left-hand-side, rhs of 'equation' */
	      lhs[0] = rhs[0] = '\0';
	      if (strchr(line,'=')) {
                 for (i=0,ii=0;i<(strchr(line,'=')-line);i++) {
		        lhs[ii++]=line[i];
	             }
	         lhs[ii] = '\0';
		 i++;
                 for (ii=0;i<strlen(line)-1;i++,ii++) {
		        rhs[ii]=line[i];
	             }
	         rhs[ii] = '\0';
		 }
		 
	      if (!strcmp(lhs,"OBJECT")) {
	           wasjustobject=1;
	           obj_offset[n_obj++]=strlen(objs);
	           strcat(objs,rhs);/*printf("%s\n", objs);*/
	          }
	      if (!strcmp(lhs,"GROUP")) wasjustobject=0;
	      if ((!strcmp(lhs,"CLASS"))&&(wasjustobject)) {
		   for (i=0;i<strlen(rhs)-1;i++) rhs[i]=rhs[i+1];
	           rhs[i-1] = '\0';
	           strcat(objs,rhs);
	          }
	      if (!strcmp(lhs,"END_OBJECT")) {
	           objs[obj_offset[--n_obj]] = '\0';
	          }
	      if (!strcmp(lhs,"NUM_VAL")) {
	           n_val=atoi(rhs);
	          }
	      
	      if (!strcmp(lhs,"VALUE")) {
		     /*printf("!%s! %d *%s*%s*\n", objs,n_obj,lhs,rhs);*/
	         if (!strcmp(objs, metastring)) {
	             strcpy(data,rhs);
		     *count=n_val;
		     /*printf("!%s! %d *%s*%s*\n", objs,n_obj,lhs,rhs);*/
	                }
	            }
	      if ((!strcmp(lhs,"\0"))&&(wasjustobject)) {
		     /*printf("!%s! %d *%s*%s*\n", objs,n_obj,lhs,rhs);*/
	         if (!strcmp(objs, metastring)) {
	             strcat(data,line);
		     *count+=n_val;
		     /*printf("!%s! %d *%s*%s*\n", objs,n_obj,lhs,rhs);*/
	                }
	            }
	    } while ( line[0] != '\0') ;
	 
	 
	 free(charattr);  
                                     /*}*/


 
}
void get_a_line(char *text, int *start, char *line)
{
int i=0;
int where;
int getout = 0;

where = *start;
line[i] = '\0';
if (text[where] == '\0') return;
else {
    while (getout == 0) {  
	 if ((text[where] == '\0')||(text[where] == '\n')||(i>=MAXLENGTH)) 
	     getout=1;
         line[i++] = text[where++];
      }; 
    *start = where;
    line[i] = '\0';
    return;
     }
}

void process_dim_ids(int32 sd_id1, int32 id_out)
{
/* This routine copies dimension IDs from sd_id1 to id_out. 
 */
int32 id1;
int32 n_attr, n_sets, count, number_type;
char name1[MAXLENGTH];
char dim_name[MAXLENGTH];

id1 = SDgetdimid(sd_id1, 0);
if (SDdiminfo(id1, dim_name, &count, &number_type, &n_attr) == 0) {
  id1 = SDgetdimid(id_out, 0);
  SDsetdimname(id1, dim_name); 
}

id1 = SDgetdimid(sd_id1, 1);
if (SDdiminfo(id1, dim_name, &count, &number_type, &n_attr) == 0) {
  id1 = SDgetdimid(id_out, 1);
  SDsetdimname(id1, dim_name); 
}

}

void process_local_attributes(int32 sd_id1, int32 id_out)
{
/* This routine copies attributes from sd_id1 to id_out. 
 */
      
int32 i, j, jj;
int32 id1;
short equal;
int32 n_attr1, n_sets1, count1, rank1, dims1[5], number_type1;
char name1[MAXLENGTH];
char attrib[MAXLENGTH], attrib1[MAXLENGTH];
char label[MAXLENGTH], tag[MAXLENGTH];
short *contents;
short ii;
char *newstring;
char *charattr;
uchar8 *ucharattr;
int16 *shortattr;
uint16 *ushortattr;
int32 *intattr;
uint32 *uintattr;
float32 *floatattr;
float64 *doubleattr;
char charatt;
uchar8 ucharatt;
int16 shortatt;
uint16 ushortatt;
int32 intatt;
uint32 uintatt;
float32 floatatt;
float64 doubleatt;
void structmetahandler(char *charattr);

id1 = sd_id1;
SDgetinfo(sd_id1, name1, &rank1, dims1, &number_type1, &n_attr1);
     		   		   

for (j=0;j<n_attr1;j++) {
     SDattrinfo(id1, j, attrib1, &number_type1, &count1);

     switch(number_type1) {   
       case DFNT_CHAR8:      
       case DFNT_INT8: 
           if (count1 == 1) {
             SDreadattr(id1, j, &charatt);  
	     SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)&charatt);
	                    }
	   else {
	     if ((charattr = (char *)malloc((count1+1)*sizeof(char))) == NULL) {
	         printf("Out of memory, array 'charattr'\n");
		 return;
	                                                                       }
             SDreadattr(id1, j, charattr);  
	     charattr[count1] = '\0';
	     SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)charattr);
 	     free(charattr);
	        }
           break;
       case DFNT_UCHAR8:    /* treat them like integers */
       case DFNT_UINT8:   
           if (count1 == 1) {
             SDreadattr(id1, j, &ucharatt);  
             SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)&ucharatt);
 	                     }
	   else {
	     if ((ucharattr = 
	       (unsigned char *)malloc(count1*sizeof(unsigned char))) == NULL) {
	         printf("Out of memory, array 'ucharattr'\n");
		 return;
	                                                                       }
             SDreadattr(id1, j, ucharattr);  
             SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)ucharattr);		  
	     free(ucharattr);
	        }
           break;
       case DFNT_INT16: 
           if (count1 == 1) {
             SDreadattr(id1, j, &shortatt);  
             SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)&shortatt);
 	                    }
	   else {
	     if ((shortattr = (short *)malloc(count1*sizeof(short))) == NULL) {
	         printf("Out of memory, array 'shortattr'\n");
		 return;
	                                                                      }
             SDreadattr(id1, j, shortattr);  
             SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)shortattr);
	     free(shortattr);
	        }
           break;
       case DFNT_UINT16: 
           if (count1 == 1) {
             SDreadattr(id1, j, &ushortatt);  
             SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)&ushortatt);
 	                    }
	   else {
	     if ((ushortattr = (unsigned short *)malloc(count1*sizeof(unsigned short))) == NULL) {
	         printf("Out of memory, array 'ushortattr'\n");
		 return;
	                                                                      }
             SDreadattr(id1, j, ushortattr);  
             SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)ushortattr);
	     free(ushortattr);
	        }
           break;
       case DFNT_INT32: 
           if (count1 == 1) {
             SDreadattr(id1, j, &intatt);  
             SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)&intatt);
 	                     }
	   else {
	     if ((intattr = (int32 *)malloc(count1*sizeof(int32))) == NULL) {
	         printf("Out of memory, array 'intattr'\n");
		 return;
	                                                                    }
             SDreadattr(id1, j, intattr);  
             SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)intattr);
	     free(intattr);
	        }
           break;
       case DFNT_UINT32: 
           if (count1 == 1) {
             SDreadattr(id1, j, &uintatt);  
	     SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)&uintatt);
 	                    }
	   else {
	     if ((uintattr = (uint32 *)malloc(count1*sizeof(uint32))) == NULL) {
	         printf("Out of memory, array 'uintattr'\n");
		 return;
	                                                                       }
             SDreadattr(id1, j, uintattr);  
             SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)uintattr);
	     free(uintattr);
	        }
           break;
       case DFNT_FLOAT: 
           if (count1 == 1) {
             SDreadattr(id1, j, &floatatt);  
             SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)&floatatt);
                            } 	                    
	   else {
	     if ((floatattr = (float *)malloc(count1*sizeof(float))) == NULL) {
	         printf("Out of memory, array 'floatattr'\n");
		 return;
	                                                                      }
             SDreadattr(id1, j, floatattr);  
             SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)floatattr);
	     free(floatattr);
	        }
           break;
       case DFNT_DOUBLE: 
           if (count1 == 1) {
             SDreadattr(id1, j, &doubleatt);  
             SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)&doubleatt);
 	                    }
	   else {
	     if ((doubleattr = (double *)malloc(count1*sizeof(double))) == NULL) {
	         printf("Out of memory, array 'doubleattr'\n");
		 return;
	                                                                      }
             SDreadattr(id1, j, doubleattr);  
             SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)doubleattr);
	     free(doubleattr);
	        }
           break;
                   }  /* switch(..) */
                        }  /* for (j=0;j<n_gattr;j++) */


}

void process_global_attributes(int32 sd_id1, int32 id_out)
{
/* This routine copies attributes from sd_id1 to id_out. 
 */
      
int32 i, j, jj;
int32 id1;
short equal;
int start;
int32 n_attr1, n_sets1, count, count1, rank1, dims1[5], number_type1;
char name1[MAXLENGTH];
char attrib[MAXLENGTH], attrib1[MAXLENGTH];
char label[MAXLENGTH], tag[MAXLENGTH];
short *contents;
short ii;
char *newstring;
char *charattr;
uchar8 *ucharattr;
int16 *shortattr;
uint16 *ushortattr;
int32 *intattr;
uint32 *uintattr;
float32 *floatattr;
float64 *doubleattr;
char charatt;
uchar8 ucharatt;
int16 shortatt;
uint16 ushortatt;
int32 intatt;
uint32 uintatt;
float32 floatatt;
float64 doubleatt;
char *newstructmetadata;
char *oldstructmetadata;
int NEWShave[4], newshave, structhave, switch_metadata, n_cols_original, n_lines_original;
float NEWSvals[4], rez;
enum{N,E,W,S};
char projstring[100];
char line[MAXLENGTH], objs[MAXLENGTH];
char lhs[MAXLENGTH], rhs[MAXLENGTH];
void get_a_line(char *text, int *start, char *line);
int get_RHS(char *line);
void put_RHS_str(char *line, char  *x);
void put_RHS_int(char *line, int x);
char MeterString[MAXLENGTH];

id1 = sd_id1;
SDfileinfo(sd_id1, &n_sets1, &n_attr1);

/* This version of SDSreader (v2.1) checks for StructMetadata, and
   if it is found updates its values.  It takes new ULLAT, ULLON, LRLAT
   and LRLON values from the global ULR, LUC, NR, NC values (the zeroth)
   and from the resolution of the original data (since this program doesn't 
   resample in any way).  */
     		   		   
newshave = 1;
structhave = 0;
switch_metadata = 0;
n_cols_original = 0;
n_lines_original = 0;
for (j=0;j<4;j++) NEWShave[j] = 0;

for (j=0;j<n_attr1;j++) {  /* preview the attributes */

     if (structhave == 1) continue;   /* in case of duplicates, take the first */

     SDattrinfo(id1, j, attrib1, &number_type1, &count1);
     for (i=0;i<strlen(attrib1);i++) attrib[i] = toupper(attrib1[i]);  
     attrib[i] = '\0';  
     if (!strcmp(attrib, "WESTBOUNDINGCOORDINATE")) {
         NEWShave[W] = 1;
	 SDreadattr(id1, j, &floatatt);  
	 NEWSvals[W] = floatatt;
        }
     if (!strcmp(attrib, "EASTBOUNDINGCOORDINATE")) {
         NEWShave[E] = 1;
	 SDreadattr(id1, j, &floatatt);  
	 NEWSvals[E] = floatatt;
        }
     if (!strcmp(attrib, "NORTHBOUNDINGCOORDINATE")) {
         NEWShave[N] = 1;
	 SDreadattr(id1, j, &floatatt);  
	 NEWSvals[N] = floatatt;
        }
     if (!strcmp(attrib, "SOUTHBOUNDINGCOORDINATE")) {
         NEWShave[S] = 1;
	 SDreadattr(id1, j, &floatatt);  
	 NEWSvals[S] = floatatt;
        }
     if (strstr(attrib, "STRUCTMETADATA")) {
	 if ((oldstructmetadata = (char *)malloc((count1+100)*sizeof(char))) == NULL) {
	   printf("Out of memory, array 'oldstructmetadata'\n");
           return;
	                                                                            }
         SDreadattr(id1, j, oldstructmetadata);  
	 count = count1 = strlen(oldstructmetadata); 
	 oldstructmetadata[count] = '\0';
         structhave = 1;
        }
     }
for (j=0;j<4;j++) newshave *= NEWShave[j];


/* Check to see if we have the StructMetadata and if
   the projection of the original was LLL. */

if (structhave) {
    start = 0;
    do {
       line[0] = '\0';
       get_a_line(oldstructmetadata, &start, line);
       if (strstr(line, "Projection=GCTP_GEO\0")) switch_metadata = 1;  /* this will only be done for LLL projection */
       if (strstr(line, "XDim=\0")) n_cols_original = get_RHS(line);
       if (strstr(line, "YDim=\0")) n_lines_original = get_RHS(line);
       } while ( line[0] != '\0');
   }
   
switch_metadata *= newshave;
switch_metadata *= (n_cols_original  != 0);
switch_metadata *= (n_lines_original != 0);

if (switch_metadata) {
    /* get memory for new struct metadata */
    if ((newstructmetadata= (char *)malloc((count+100)*sizeof(char))) == NULL) {
        printf("Out of memory, array 'newstructmetadata'\n");
        return;
         }
    newstructmetadata[0] = '\0';
       
    /* Get new bounding coordinates... */
    /* get the rez of the latitudes */
    rez = (float)n_lines_original/(NEWSvals[N] - NEWSvals[S]);
    
    /* get approximate UL lat */
    NEWSvals[N] = NEWSvals[N] - (float)ULR[0]/rez;
    
    /* get approximate LR lat */
    NEWSvals[S] = NEWSvals[N] - (float)(NR[0]-1)/rez;
    
    /* get the rez of the longitudes */
    rez = (float)n_cols_original/(NEWSvals[E] - NEWSvals[W]);
    
    /* get approximate UL lon */
    NEWSvals[W] = NEWSvals[W] + (float)ULC[0]/rez;
    
    /* get approximate LR lon */
    NEWSvals[E] = NEWSvals[W] + (float)(NC[0]-1)/rez;

    
    /* Write a new struct metadata  */
    start = 0;
    do {
       line[0] = '\0';
       get_a_line(oldstructmetadata, &start, line);
       if (strstr(line, "XDim=\0")) put_RHS_int(line, NC[0]);
       if (strstr(line, "YDim=\0")) put_RHS_int(line, NR[0]);
       if (strstr(line, "UpperLeftPointMtrs=\0")) {
          sprintf(MeterString, "(%f,%f)\n", NEWSvals[W]*1000000.0, NEWSvals[N]*1000000.0);  /* true only of LLL */
          put_RHS_str(line, MeterString);
          }
       if (strstr(line, "LowerRightMtrs=\0")) {
          sprintf(MeterString, "(%f,%f)\n", NEWSvals[E]*1000000.0, NEWSvals[S]*1000000.0);
          put_RHS_str(line, MeterString);
	  }
       sprintf(newstructmetadata+strlen(newstructmetadata), "%s", line);
       } while ( line[0] != '\0');
    free(oldstructmetadata);
   }

/* Now process all the attributes as usual, but making substitutions where necessary */

for (j=0;j<n_attr1;j++) {
     SDattrinfo(id1, j, attrib1, &number_type1, &count1);

     switch(number_type1) {   
       case DFNT_CHAR8:      
       case DFNT_INT8: 
           if (count1 == 1) {
             SDreadattr(id1, j, &charatt);  
	     SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)&charatt);
	                    }
	   else {
	     if ((charattr = (char *)malloc((count1+100)*sizeof(char))) == NULL) {
	         printf("Out of memory, array 'charattr'\n");
		 return;
	                                                                       }
             SDreadattr(id1, j, charattr); 
	     charattr[count1] = '\0';
	     for (i=0;i<strlen(attrib1);i++) attrib[i] = toupper(attrib1[i]);  
             attrib[i] = '\0';  
	     if ( strstr(attrib,"STRUCTMETADATA") && (switch_metadata == 1) ) {
	         SDsetattr(id_out, attrib1, number_type1, strlen(newstructmetadata), (VOIDP)newstructmetadata);
 	         free(newstructmetadata);
		 switch_metadata = 0;
		 }
             else {		 
	        SDsetattr(id_out, attrib1, number_type1, strlen(charattr), (VOIDP)charattr);
 		  }
 	     free(charattr);
	        }
           break;
       case DFNT_UCHAR8:    /* treat them like integers */
       case DFNT_UINT8:   
           if (count1 == 1) {
             SDreadattr(id1, j, &ucharatt);  
             SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)&ucharatt);
 	                     }
	   else {
	     if ((ucharattr = 
	       (unsigned char *)malloc(count1*sizeof(unsigned char))) == NULL) {
	         printf("Out of memory, array 'ucharattr'\n");
		 return;
	                                                                       }
             SDreadattr(id1, j, ucharattr);  
             SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)ucharattr);		  
	     free(ucharattr);
	        }
           break;
       case DFNT_INT16: 
           if (count1 == 1) {
             SDreadattr(id1, j, &shortatt);  
             SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)&shortatt);
 	                    }
	   else {
	     if ((shortattr = (short *)malloc(count1*sizeof(short))) == NULL) {
	         printf("Out of memory, array 'shortattr'\n");
		 return;
	                                                                      }
             SDreadattr(id1, j, shortattr);  
             SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)shortattr);
	     free(shortattr);
	        }
           break;
       case DFNT_UINT16: 
           if (count1 == 1) {
             SDreadattr(id1, j, &ushortatt);  
             SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)&ushortatt);
 	                    }
	   else {
	     if ((ushortattr = (unsigned short *)malloc(count1*sizeof(unsigned short))) == NULL) {
	         printf("Out of memory, array 'ushortattr'\n");
		 return;
	                                                                      }
             SDreadattr(id1, j, ushortattr);  
             SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)ushortattr);
	     free(ushortattr);
	        }
           break;
       case DFNT_INT32: 
           if (count1 == 1) {
             SDreadattr(id1, j, &intatt);  
             SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)&intatt);
 	                     }
	   else {
	     if ((intattr = (int32 *)malloc(count1*sizeof(int32))) == NULL) {
	         printf("Out of memory, array 'intattr'\n");
		 return;
	                                                                    }
             SDreadattr(id1, j, intattr);  
             SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)intattr);
	     free(intattr);
	        }
           break;
       case DFNT_UINT32: 
           if (count1 == 1) {
             SDreadattr(id1, j, &uintatt);  
	     SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)&uintatt);
 	                    }
	   else {
	     if ((uintattr = (uint32 *)malloc(count1*sizeof(uint32))) == NULL) {
	         printf("Out of memory, array 'uintattr'\n");
		 return;
	                                                                       }
             SDreadattr(id1, j, uintattr);  
             SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)uintattr);
	     free(uintattr);
	        }
           break;
       case DFNT_FLOAT: 
           if (count1 == 1) {
             SDreadattr(id1, j, &floatatt);  
	     
             if (switch_metadata) {
                for (i=0;i<strlen(attrib1);i++) attrib[i] = toupper(attrib1[i]);  
                attrib[i] = '\0';  
                if (!strcmp(attrib, "WESTBOUNDINGCOORDINATE"))        floatatt = NEWSvals[W];
                else if (!strcmp(attrib, "EASTBOUNDINGCOORDINATE"))   floatatt = NEWSvals[E];
                else if (!strcmp(attrib, "NORTHBOUNDINGCOORDINATE"))  floatatt = NEWSvals[N];
                else if (!strcmp(attrib, "SOUTHBOUNDINGCOORDINATE"))  floatatt = NEWSvals[S];
                 }
		 
             SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)&floatatt);
                            } 	                    
	   else {
	     if ((floatattr = (float *)malloc(count1*sizeof(float))) == NULL) {
	         printf("Out of memory, array 'floatattr'\n");
		 return;
	                                                                      }
             SDreadattr(id1, j, floatattr);  
             SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)floatattr);
	     free(floatattr);
	        }
           break;
       case DFNT_DOUBLE: 
           if (count1 == 1) {
             SDreadattr(id1, j, &doubleatt);  
             SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)&doubleatt);
 	                    }
	   else {
	     if ((doubleattr = (double *)malloc(count1*sizeof(double))) == NULL) {
	         printf("Out of memory, array 'doubleattr'\n");
		 return;
	                                                                      }
             SDreadattr(id1, j, doubleattr);  
             SDsetattr(id_out, attrib1, number_type1, count1, (VOIDP)doubleattr);
	     free(doubleattr);
	        }
           break;
                   }  /* switch(..) */
                        }  /* for (j=0;j<n_gattr;j++) */


}

int get_RHS(char *line)
{
int j, i;
char rhs[10];
int n;

n = strlen(line);

for(j=0;j<n;j++) {
   if (line[j] == '=') break;
  }
i = 0;  
j++;
for(   ;j<n;j++) {
   rhs[i++] = line[j];
  }
rhs[i] = '\0';

return(atoi(rhs));
}


void put_RHS_int(char *line, int x)
{

int j, i;
char rhs[10];
int n;

n = strlen(line);

for(j=0;j<n;j++) {
   if (line[j] == '=') break;
  }
sprintf(rhs, "%d\n\0", x);
i = 0;
j++;
for(   ;j<n;j++) {
   line[j] = rhs[i++];
  }
line[j] = '\0';

return;
}


void put_RHS_str(char *line, char *x)
{

int j, i, n;

n = strlen(line);

for(j=0;j<n;j++) {
   if (line[j] == '=') break;
  }
i = 0;
j++;
for(   ;j<n;j++) {
   line[j] = x[i++];
  }
line[j] = '\0';

return;
}





