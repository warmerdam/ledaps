/**
 * !Description
 * parse metadata from input metadata file (included in the Landsat files) and 
 * create input files for Ledaps modules
 *
 * !Input
 * metadata file from LPGS or NLAPS processing system  
 *
 * !Ouput
 * header file for lndcal module
 *
 * !Developer:
 * Feng Gao (fgao@ltpmail.gsfc.nasa.gov)
 *
 * !Revision:
 *
 *  revision 1.5.8  08/17/2012 Gail Schmidt, USGS EROS
 *  - updated the metadata tags to work with the newly released LPGS metadata
 *    as well as continuing to support the old metadata tags
 *  - cleaned up warning messages from compilation
 *
 *  revision 1.5.7.5  1/18/2012 
 *  - restores the solar zenith angle bug fix from the past for NLAP_W0 format (Greg Ederer)
 *  - fixes a bug when writing the UTM zone (south) into the ENVI hdr file (Greg Ederer)
 *  - added processing for Landsat-4 TM (Feng Gao)
 * 
 *  revision 1.5.7  10/22/2010 Feng Gao
 *  - included global DEM file for lndsr processing
 *  - used actual LMIN, LMAX, QCALMIN and QCALMAX values for ETM+ calibration in LPGS 
 *    if they are provided in metadata, otherwise still use gain_settings (HIGH or LOW)
 *  - cleaned compiling warning messages
 *
 *  revision 1.5.6  9/29/2009 Greg Ederer
 *  - replaced hard coded version strings with macros that can be overridden
 *    at compile time. The build script can substitute current version info
 *    automatically.
 *
 *  revision 1.5.5  9/28/2009 Greg Ederer and Feng Gao 
 *  - replaced fixed ancillary data path with path searched from ANC_PATH 
 *
 *  revision 1.5.4  6/22/2009 Feng Gao
 *  - fixed a bug in computing bais from (lmin, lmax) and (qmin, qmax)  
 *
 *  revision 1.5.3  4/7/2009 Feng Gao
 *  - read GTRasterTypeGeoKey from GeoTiff file and determine the coordinate for the UL of the UL pixel
 *    (note that both LPGS and NLAPS use CENTER of pixel in the provided text metadata)
 *    However, in GeoTiff file, LPGS  uses center (RasterPixelIsArea) 
 *                              NLAPS uses upper left (RasterPixelIsPoint)
 *
 *  revision 1.5.2  4/2/2009 Feng Gao
 *  - added option on re-calibration for public version
 *    in "lndcal*.txt" added:   RE_CAL = "NO" or RE_CAL = "FALSE"
 *    (default: do re-calibration as internal ledaps does)
 *  - included ENVI header for each HDF file so ENVI can load correct geo-reference information
 *
 *  revision 1.5.1 3/19/2009 Feng Gao
 *  - corrected LMIN, LMAX conversion to gain and bias using actual QCALMIN and QCALMAX in the meta file 
 *  - used actual gain and bias values for TM to match recent changes in "lndcal"
 *
 *  revision 1.5  12/16/2008 Feng Gao
 *  - added processing for LPGS TM data (newly EROS released TM and ETM+)    
 *
 *  revision 1.4  8/28/2007  Feng Gao
 *  - added processing for MSS data
 *
 *  revision 1.3  6/22/2007  Feng Gao
 *  - now accept LC-ComPS "standard" metadata from UMD (file name must in *.umd)
 *
 *  revision 1.2  10/21/2006 Feng Gao
 *  - read upper left coordinate from image file if it is a GeoTiff input
 *    note that metadata in GeoCover data are not updated on UL coordinate
 *
 *  revision 1.1  2/24/2006  Feng Gao
 *  - revise to accept metadata from both NLAPS and LPGS Landsat processing system  
 *  - add error handling 
 *  - add log report
 *
 *  original version - 04/05 by Feng Gao (NASA/GSFC through ERT)
 *
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <sys/types.h>
#include <dirent.h>
#include <sys/stat.h>
#include <time.h>
#include <ctype.h>

#include "geotiffio.h"
#include "xtiffio.h"
#include "tiffio.h"

/*
 * the following VERSION definitions are expected to be
 * defined at compile time with the -D option
 */
#ifndef PGE200_VERSION
#define PGE200_VERSION "1.0.9"
#endif

#ifndef PGE201_VERSION
#define PGE201_VERSION "1.0.9"
#endif

#ifndef MOD_CAL_VERSION
#define MOD_CAL_VERSION "1.0.6"
#endif

#ifndef MOD_CSM_VERSION
#define MOD_CSM_VERSION "1.0.4"
#endif

#ifndef MOD_SR_VERSION
#define MOD_SR_VERSION "1.0.4"
#endif

/*
 * define useful constants
 */
#define MAX_STRING_LENGTH 1000
#define DIR_BUF_SIZE (MAXNAMLEN*20)
#define MSS 1
#define TM 2
#define ETM 3
#define SUCCESS 0
#define FAILURE -1

#define MAX_N_BANDS   7
#define N_BANDS_MSS   4

/* define the list of carbon_met PARAMETER IDs */
#define FILE_TYPE            1
#define DATA_PROVIDER        2
#define SATELLITE            3
#define INSTRUMENT           4
#define ACQUISITION_DATE     5
#define PRODUCTION_DATE      6
#define SOLAR_ZENITH         7
#define SOLAR_AZIMUTH        8
#define WRS_SYSTEM           9
#define WRS_PATH            10
#define WRS_ROW             11
#define NBAND               12
#define BANDS               13
#define GAIN                14
#define BIAS                15
#define NSAMPLE             16
#define NLINE               17
#define FILE_NAMES          18
#define NBAND_TH            19
#define BANDS_TH            20
#define GAIN_TH             21
#define BIAS_TH             22
#define NSAMPLE_TH          23
#define NLINE_TH            24
#define FILE_NAMES_TH       25
#define PROJECTION_NUMBER   26
#define PIXEL_SIZE          27
#define UPPER_LEFT_CORNER   28
#define PROJECTION_ZONE     29
#define PROJECTION_SPHERE   30
#define MAX_CM_PARAMS       (PROJECTION_SPHERE + 1)
static char *DELIMIT = ", ";

typedef struct lut {
  int   cm_id;        /* the line number of the cm_keyword in the carbon_met file */
  char *in_name;      /* keyword read from input file */
  int   n_fields;     /* the number of fields to read after the keyword */
} LUT;                /* the last entry in an erray of these MUST BE {0, 0, 0} */

typedef struct metadata {
  char name[MAX_STRING_LENGTH];
  char val[MAX_STRING_LENGTH];
} METADATA;

/* error handling and log file */
void WARN(char *fmt, ...);
void die(char *name, int lineno, char *fmt, ...);
#define DIE(...)         die(__FILE__,__LINE__,__VA_ARGS__)

#define logFile "LogReport"
FILE *LOG_FP;
int  open_log(char name[]);
void LOG(char *fmt, ...);

int conv_date(int *mm, int *dd, int yyyy);
int getMetaFromLPGS(char input[], char scene_name[], char acquisition_date[]);
int getMetaFromNLAPS(char input[], char scene_name[], char acquisition_date[]);
int getMetaFromNLAPS_WO(char input[], char scene_name[], char acquisition_date[]);
int readULfromGeoTiff(char *fname, int *nrows, int *ncols, double ulxy[], double *res);
int getMetaFromLCT(char input[], char scene_name[], char acquisition_date[]);
int writeENVIHeader(char scene_name[], int nrows, int ncols, double ulxy[], double res, int zone);
int create_ENVI_hdr_file(char *scene_name, METADATA *carbon_met);
int find_file(char *path, char *name);
char *cm_name_for_id(int id);
int read_n_fields(char *str, char *key, LUT *lut, char container[MAX_STRING_LENGTH]);
int nvp_to_meta(char *name, char *val, LUT *lut, METADATA *carbon_met);
int validate_carbon_met(METADATA *carbon_met, int param_index);
int create_carbon_met_file(char *file_name, METADATA *carbon_met);

int main(int argc, char *argv[])
{
  int  year, month, day, ret;

  /* file names */
  char  input[MAX_STRING_LENGTH];
  char  scene_name[MAX_STRING_LENGTH];
  char  lndcal_name[MAX_STRING_LENGTH];
  char  lndcsm_name[MAX_STRING_LENGTH];
  char  lndsr_name[MAX_STRING_LENGTH];
  char  dem[MAX_STRING_LENGTH];
  char  ozone[MAX_STRING_LENGTH];
  char  reanalysis[MAX_STRING_LENGTH];
  char  path_buf[DIR_BUF_SIZE];
  char  cal_file[MAX_STRING_LENGTH];
  char  *anc_path;
  char  *tokenptr = NULL;
  char  tmpstr[MAX_STRING_LENGTH] = "\0";
  char  acquisition_date[MAX_STRING_LENGTH] = "\0";

  FILE  *out;

  printf ("\nRunning lndpm ...\n");
  if(open_log("lndpm : Landsat Metadata Parser")==FAILURE) {
    DIE("create log report error");
  }

  if(argc!=2) {
    DIE("Usage: lndpm <input_metadata_file>");
  }
  
  strcpy(input, argv[1]);

  tokenptr = strtok(argv[1], ".");
  sprintf(scene_name, "%s", tokenptr);
  if(strstr(scene_name, "_MTL")) {
    /* remove last 4 chars "_MTL" from scene name */
    strncpy(tmpstr, scene_name, strlen(scene_name)-4);
    strcpy(scene_name, tmpstr);
  }
  
  sprintf(lndcal_name, "lndcal.%s.txt", scene_name);
  sprintf(lndcsm_name, "lndcsm.%s.txt", scene_name);
  sprintf(lndsr_name, "lndsr.%s.txt", scene_name);

  /* generate parameters for LEDAPS module from LPGS metadata file (ETM+) */
  if(strstr(input, "_MTL.") || strstr(input, ".met")) {
    if((ret=getMetaFromLPGS(input, scene_name, acquisition_date))==FAILURE) {
      DIE("reading metadata from LPGS input");
    }
  }

  /* generate parameters for LEDAPS modules from NLAPS metadata file (TM, MSS) */
  else if(strstr(input, ".H1")||strstr(input, ".hdr")) {
    if((ret=getMetaFromNLAPS(input, scene_name, acquisition_date))==FAILURE) {
      DIE("reading metadata from NLAPS input");
    }
  }

  /* generate parameters for LEDAPS modules from NLAPS W0 metadata file (TM, MSS) */
  else if(strstr(input, "_WO") || strstr(input, ".prodReport")) {
    if((ret=getMetaFromNLAPS_WO(input, scene_name, acquisition_date))==FAILURE) {
      DIE("reading metadata from NLAPS input");
    }
  }

  /* generate parameters for LEDAPS modules from UMD metadata file */
  else if(strstr(input, ".umd")||strstr(input, ".UMD")) {
    if((ret=getMetaFromLCT(input, scene_name, acquisition_date))==FAILURE) {
      DIE("reading metadata from LC-ComPS input");
    }
  }

  else {
    DIE("unknown metadata format"); 
  }

  /* write input card for lndcal module */
  if((out=fopen(lndcal_name, "w"))==NULL) {
    DIE("open lndcal input file for write");
  }    

  anc_path=getenv("ANC_PATH");
  if(anc_path == NULL) {
    anc_path = ".";
  }

  fprintf(out, "PARAMETER_FILE\n");
  fprintf(out, "HEADER_FILE = %s.carbon_met.txt\n", scene_name);
  fprintf(out, "REF_FILE = lndcal.%s.hdf\n", scene_name);
  fprintf(out, "THERM_FILE = lndth.%s.hdf\n", scene_name);
  if(ret == TM) {
    /* find, prepare calibration files */
    strcpy(cal_file, "gold.dat");
    strcpy(path_buf, anc_path);
    if (find_file(path_buf, cal_file))
      strcpy(cal_file, path_buf);
    else {
      WARN("could not find gold.dat");
    }
    fprintf(out, "GOLD_FILE = %s\n", cal_file);
    strcpy(cal_file, "gnew.dat");
    strcpy(path_buf, anc_path);
    if (find_file(path_buf, cal_file))
      strcpy(cal_file, path_buf);
    else {
      WARN("could not find gnew.dat");
    }
    fprintf(out, "GNEW_FILE = %s\n", cal_file);
    strcpy(cal_file, "gold_2003.dat");
    strcpy(path_buf, anc_path);
    if (find_file(path_buf, cal_file))
      strcpy(cal_file, path_buf);
    else {
      WARN("could not find gold_2003.dat");
    }
    fprintf(out, "GOLD_2003 = %s\n", cal_file);
  }
  fprintf(out, "PGEVersion = %s\n", PGE200_VERSION);
  fprintf(out, "ProcessVersion = %s\n",MOD_CAL_VERSION);
  fprintf(out, "DNOUT = \"FALSE\"\n");
  fprintf(out, "RE_CAL = \"NO\"\n");
  fprintf(out, "END\n");
  fclose(out);

  /* write input card for lndcsm module */
  if((out=fopen(lndcsm_name, "w"))==NULL) {
    DIE("open lndcsm input file for write");
  }    
  fprintf(out, "PARAMETER_FILE\n");
  fprintf(out, "CSM_FILE = lndcsm.%s.hdf\n", scene_name);
  fprintf(out, "REF_FILE = lndcal.%s.hdf\n", scene_name);
  fprintf(out, "THERM_FILE = lndth.%s.hdf\n", scene_name);
  fprintf(out, "PGEVersion = %s\n", PGE200_VERSION);
  fprintf(out, "ProcessVersion = %s\n", MOD_CSM_VERSION);

  fprintf(out, "END\n");
  fclose(out);

  /* get year, month and day from acquisition date */
  tokenptr = strtok(acquisition_date, "-");
  sscanf(tokenptr, "%d", &year);
  tokenptr = strtok(NULL, "-");
  sscanf(tokenptr, "%d", &month);
  tokenptr = strtok(NULL, "-");
  sscanf(tokenptr, "%d", &day);

  /* convert to day of year */
  conv_date(&month, &day, year);

  /* find, prepare ancillary files */
  strcpy(dem, "CMGDEM.hdf");
  strcpy(path_buf, anc_path);
  if (find_file(path_buf, dem))
  {
    strcpy(dem, path_buf);
    printf("using DEM : %s\n", dem);
  }
  else
  {
    WARN("could not find DEM ancillary data: %s\n  check ANC_PATH environment variable", dem);
  }

  sprintf(ozone, "TOMS_%d%03d.hdf", year, day);
  strcpy(path_buf, anc_path);
  if (find_file(path_buf, ozone))
  {
    strcpy(ozone, path_buf);
    printf("using TOMS : %s\n", ozone);
  }
  else
  {
    WARN("could not find TOMS ancillary data: %s\n  check ANC_PATH environment variable", ozone);
  }
    
  sprintf(reanalysis, "REANALYSIS_%d%03d.hdf", year, day);
  strcpy(path_buf, anc_path);
  if (find_file(path_buf, reanalysis))
  {
    strcpy(reanalysis, path_buf);
    printf("using REANALYSIS : %s\n", reanalysis);
  }
  else
  {
    WARN("could not find REANALYSIS ancillary data: %s\n  check ANC_PATH environment variable", reanalysis);
  }

  /* write input card for lndsr module */
  if((out=fopen(lndsr_name, "w"))==NULL) {
    DIE("open lndsr input file for write");
  }    
  fprintf(out, "PARAMETER_FILE\n");
  fprintf(out, "DEM_FILE = %s\n", dem);
  if(fopen(ozone,"r")!=NULL)  // if not exist then use climatology estimation
    fprintf(out, "OZON_FIL = %s\n", ozone);
  fprintf(out, "PRWV_FIL = %s\n", reanalysis);
  fprintf(out, "CSM_FILE = lndcsm.%s.hdf\n", scene_name);
  fprintf(out, "REF_FILE = lndcal.%s.hdf\n", scene_name);
  fprintf(out, "TEMP_FILE = lndth.%s.hdf\n", scene_name);
  fprintf(out, "SR_FILE = lndsr.%s.hdf\n", scene_name);
  fprintf(out, "PGEVersion = %s\n", PGE201_VERSION);
  fprintf(out, "ProcessVersion = %s\n", MOD_SR_VERSION);
  fprintf(out, "END\n");
  fclose(out);
  fclose(LOG_FP);

  printf ("lndpm complete.\n");
  return SUCCESS;
}


/* retrieve metadata from LPGS format (TM/ETM+) */
/* Modified: Gail Schmidt, USGS EROS
 *   Made changes to the LPGS metadata tags to support the new LPGS metadata
 *   format and tags.  This software will also continue to support the old
 *   metadata tags in case the user wants to use older data in their personal
 *   archive.
 */
int getMetaFromLPGS(char input[], char scene_name[], char acquisition_date[])
{
  int    i, year, pyear, zone, geotiff, instrument;
  int    nrows, ncols, bno;
  float  lmax[8], lmin[8], qmin[8], qmax[8];
  double ulxy[2], res;
  int    chkFlg[34];

 /* vars used in parameter parsing */
  char  buffer[MAX_STRING_LENGTH] = "\0";
  char  tmpstr[MAX_STRING_LENGTH] = "\0";
  char  tmpstr2[MAX_STRING_LENGTH] = "\0";
  char  production_date[MAX_STRING_LENGTH] = "\0";
  char  processing_software[MAX_STRING_LENGTH] = "\0";
  char  sensor[MAX_STRING_LENGTH] = "\0";
  char  met_name[MAX_STRING_LENGTH] = "\0";
  char  *label = NULL;
  char  *tokenptr = NULL;
  char  *seperator = "=\" \t";
  float fnum;

  /* metadata */
  char  meta_list[34][MAX_STRING_LENGTH];
  char  band_fname[8][MAX_STRING_LENGTH];
  char  band_gain[8][MAX_STRING_LENGTH];
  char  ul_xy[2][MAX_STRING_LENGTH];
  char  *meta_lut[34][3] = {
		{"HEADER_FILE",       "", ""},
	    {"FILE_TYPE",         "OUTPUT_FORMAT", ""},
	    {"DATA_PROVIDER",     "STATION_ID", "GROUND_STATION"},
	    {"SATTELITE",         "SPACECRAFT_ID", ""},
	    {"INSTRUMENT",        "SENSOR_ID", ""},
	    {"ACQUISITION_DATE",  "DATE_ACQUIRED", "ACQUISITION_DATE"},
	    {"PRODUCTION_DATE",   "FILE_DATE", "PRODUCT_CREATION_TIME"},
	    {"SOLAR_ZENITH",      "SUN_ELEVATION", ""},
	    {"SOLAR_AziMUTH",     "SUN_AZIMUTH", ""},
	    {"WRS_SYSTEM",        "", ""},
	    {"WRS_PATH",          "WRS_PATH", ""},
	    {"WRS_ROW",           "", ""},
	    {"NBAND",             "", ""},
	    {"BANDS",             "", ""},
	    {"GAIN_SETTINGS",     "", ""},
	    {"NSAMPLE",           "REFLECTIVE_SAMPLES", "PRODUCT_SAMPLES_REF"},
	    {"NLINE",             "REFLECTIVE_LINES", "PRODUCT_LINES_REF"},
	    {"FILE_NAMES",        "", ""},
	    {"NBAND_TH",          "", ""},
	    {"BANDS_TH",          "", ""},
	    {"GAIN_SETTINGS_TH",  "GAIN_BAND_6_VCID_1", "BAND6_GAIN1"},
	    {"NSAMPLE_TH",        "THERMAL_SAMPLES", "PRODUCT_SAMPLES_THM"},
	    {"NLINE_TH",          "THERMAL_LINES", "PRODUCT_LINES_THM"},
	    {"FILE_NAMES_TH",     "", ""},
	    {"PROJECTION_NUMBER", "MAP_PROJECTION", ""},
	    {"PIXEL_SIZE",        "GRID_CELL_SIZE_REFLECTIVE", "GRID_CELL_SIZE_REF"},
	    {"UPPER_LEFT_CORNER", "", ""},
	    {"PROJECTION_ZONE",   "UTM_ZONE", "ZONE_NUMBER"},
	    {"PROJECTION_SPHERE", "", ""},
	    {"GAIN",              "", ""},
	    {"BIAS",              "", ""},
	    {"GAIN_TH",           "", ""},
	    {"BIAS_TH",           "", ""},
	    {"END",               "", ""}};  

  float gain[8], bias[8];
  /*float tm_band_width[8] = {0, 0.066, 0.082, 0.067, 0.128, 0.217, 1.0, 0.252};*/

  FILE *in, *out;
  

  if((in=fopen(input, "r"))==NULL)
    DIE( "Can't open input metadata file %s", input);
  
  sprintf(met_name, "%s.carbon_met.txt", scene_name);
  if((out=fopen(met_name, "w"))==NULL)
    DIE( "Can't open output header file for write %s", met_name);

  for(i=0; i<34; i++) chkFlg[i] = 0;
  for(i=0; i<8; i++)  {qmin[i]=0.0; qmax[i]=255; }

  /* process line by line */
  while(fgets(buffer, MAX_STRING_LENGTH, in) != NULL) {

    /* get string token */
    tokenptr = strtok(buffer, seperator);
    label=tokenptr;
 
    while(tokenptr != NULL) {
 
      tokenptr = strtok(NULL, seperator);

      /* information provided just after label, for some tags we are supporting
         both the old and the new metadata tags */
      for(i=1; i<29; i++)
	if(strcmp(label, meta_lut[i][1]) == 0 || strcmp(label, meta_lut[i][2]) == 0) {
	  chkFlg[i] = 1;
	  switch(i) {      
	    /* some metadata need more work to satisfy input */
	  case 1:
	    if(strcmp(tokenptr,"GEOTIFF") != 0) {
	      sprintf(meta_list[i], "%s = BINARY\n", meta_lut[i][0]);
	      geotiff = 0;
	    }
	    else {
	      sprintf(meta_list[i], "%s = %s\n", meta_lut[i][0], tokenptr);
	      geotiff = 1;
	    }
	    break;
	  case 3:     // SATELLITE
	    if(strcmp(tokenptr,"LANDSAT_7") == 0 || strcmp(tokenptr,"Landsat7") == 0)
	      sprintf(meta_list[i],"%s = LANDSAT_7\n", meta_lut[i][0]);
	    else if(strcmp(tokenptr,"LANDSAT_5") == 0 || strcmp(tokenptr,"Landsat5") == 0)
	      sprintf(meta_list[i],"%s = LANDSAT_5\n", meta_lut[i][0]);
	    else if(strcmp(tokenptr,"LANDSAT_4") == 0 || strcmp(tokenptr,"Landsat4") == 0)
	      sprintf(meta_list[i],"%s = LANDSAT_4\n", meta_lut[i][0]);
	    else
	      chkFlg[i] = 0;
	    break;
	  case 4:    // INSTRUMENT
	    strcpy(sensor, tokenptr);
	    if(strstr(tokenptr,"ETM") != NULL || strstr(tokenptr,"ETM+") != NULL) {
	      sprintf(meta_list[i],"%s = ETM\n", meta_lut[i][0]);
	      instrument = ETM;
	    }
	    else if(strstr(tokenptr,"TM") != NULL) {
	      sprintf(meta_list[i],"%s = TM\n", meta_lut[i][0]);
	      instrument = TM;
	    } 
	    else
	      chkFlg[i] = 0;
	    break;
	  case 6:    // PRODUCTION_DATE
	    // some input give wrong production date and need extra processing   
	    sscanf(tokenptr, "%10s", production_date);
	    break;
	  case 7:    // SOLAR_ZENITH
	    sscanf(tokenptr, "%f", &fnum);
	    sprintf(meta_list[i],"%s = %f\n", meta_lut[i][0], 90-fnum);
	    break;
	  case 20:   // GAIN_SETTINGS_TH
	    if(strcmp(tokenptr,"H") == 0)
	      sprintf(meta_list[i], "%s = HIGH\n", meta_lut[i][0]);
	    else
	      sprintf(meta_list[i], "%s = LOW\n", meta_lut[i][0]);
	    break;
	  case 24:   // projection number
	    if(strstr(tokenptr, "UTM"))
	      sprintf(meta_list[i], "%s = 1\n", meta_lut[i][0]);  
	    else {
	      WARN("Please check projection number");
	      chkFlg[i] = 0;
	    }
	    break; 
	  case 25:   // resolution
	    sscanf(tokenptr, "%lf", &res);
	    sprintf(meta_list[i], "%s = %f\n", meta_lut[i][0], res);
	    break;
	  case 27:
	    sscanf(tokenptr, "%d", &zone);
	    sprintf(meta_list[i], "%s = %d\n", meta_lut[i][0], zone);
	    break;
	  default:  // other parameters can be used directly
	    sscanf(tokenptr, "%s", tmpstr);
	    sprintf(meta_list[i],"%s = %s\n", meta_lut[i][0], tmpstr);
	  }
	  if(i==5) strcpy(acquisition_date, tmpstr);  // save acquisition date for later use	
	}
 
      /* extract coordinates of upleft corner */
      if(strcmp(label, "CORNER_UL_PROJECTION_X_PRODUCT") == 0 || strcmp(label, "PRODUCT_UL_CORNER_MAPX") == 0 || strcmp(label, "SCENE_UL_CORNER_MAPX") == 0){
	chkFlg[26] = 1; 
	sscanf(tokenptr, "%s", ul_xy[0]);
      }
      if(strcmp(label, "CORNER_UL_PROJECTION_Y_PRODUCT") == 0 || strcmp(label, "PRODUCT_UL_CORNER_MAPY") == 0 || strcmp(label, "SCENE_UL_CORNER_MAPY") == 0) 
	sscanf(tokenptr, "%s", ul_xy[1]);

      /* get WRS row */
      if(strcmp(label, "WRS_ROW") == 0 || strcmp(label, "STARTING_ROW") == 0) {
	sscanf(tokenptr, "%s", tmpstr);
	sprintf(meta_list[11],"%s = %s\n", meta_lut[11][0], tmpstr);
	chkFlg[11] = 1;
      }

      if(strcmp(label, "PROCESSING_SOFTWARE_VERSION") == 0 || strcmp(label, "PROCESSING_SOFTWARE") == 0) 
	sscanf(tokenptr, "%s", processing_software);
 
      /* get thermal band filename */
      if(instrument == ETM) {
	if(strcmp(label, "FILE_NAME_BAND_6_VCID_1") == 0 || strcmp(label, "BAND61_FILE_NAME") == 0 || strcmp(label, "BAND6L_FILE_NAME") == 0) {
	  sscanf(tokenptr, "%s", tmpstr);
	  sprintf(meta_list[23],"%s = %s\n", meta_lut[23][0], tmpstr);
	  chkFlg[23] = 1;
	}
      }
      else {
	if(strcmp(label, "FILE_NAME_BAND_6") == 0 || strcmp(label, "BAND6_FILE_NAME") == 0) {
	  sscanf(tokenptr, "%s", tmpstr);
	  sprintf(meta_list[23],"%s = %s\n", meta_lut[23][0], tmpstr);
	  chkFlg[23] = 1;
	}      
      }

      /* extract gain settings and file names for ETM+ band 1-5 and 7 */
      for(i=1; i<=7; i++) {
	if(instrument == ETM) {
	  if(i!=6) {
	    sprintf(tmpstr,"FILE_NAME_BAND_%d",i);
	    sprintf(tmpstr2,"BAND%d_FILE_NAME",i);
	    if(strcmp(label, tmpstr) == 0 || strcmp(label, tmpstr2) == 0) {
	      sscanf(tokenptr, "%s", band_fname[i]);
	      chkFlg[17] = 1;
	    }

	    sprintf(tmpstr,"GAIN_BAND_%d",i);
	    sprintf(tmpstr2,"BAND%d_GAIN",i);
	    if(strcmp(label, tmpstr) == 0 || strcmp(label, tmpstr2) == 0) {
	      if(strcmp(tokenptr,"H") == 0)
		strcpy(band_gain[i], "HIGH");
	      else
		strcpy(band_gain[i], "LOW");
	      chkFlg[14] = 1;
	    }	  
	  }
	}
      }

      /* extract filename and LMIN, LMAX, QCALMIN and QCALMAX if they exist */
      for(i=1; i<=7; i++) {
	if(instrument == ETM && i==6)
    {
	  /* use low gain thermal band for ETM+ */
	  bno = 61;
	  sprintf(tmpstr,"FILE_NAME_BAND_6_VCID_1");
	  sprintf(tmpstr2,"BAND%d_FILE_NAME",bno);
    }
	else
    {
	  bno = i;
	  sprintf(tmpstr,"FILE_NAME_BAND_%d",bno);
	  sprintf(tmpstr2,"BAND%d_FILE_NAME",bno);
    }
	/* get gain/bias info from input (added 12/2008) */
	if(strcmp(label, tmpstr) == 0 || strcmp(label, tmpstr2) == 0) {
	  sscanf(tokenptr, "%s", band_fname[i]);
	  chkFlg[17] = 1;
	}


    if (bno == 61)
	  sprintf(tmpstr,"RADIANCE_MAXIMUM_BAND_6_VCID_1");
    else
	  sprintf(tmpstr,"RADIANCE_MAXIMUM_BAND_%d",bno);
	sprintf(tmpstr2,"LMAX_BAND%d",bno);
	if(strcmp(label, tmpstr) == 0 || strcmp(label, tmpstr2) == 0) {
	  sscanf(tokenptr, "%f", &(lmax[i]));
	  chkFlg[29] = 1; chkFlg[30] = 1;
	  chkFlg[31] = 1; chkFlg[32] = 1;
	  /* disable gain settings */
	  chkFlg[14] = -1; chkFlg[20] = -1;
	}

    if (bno == 61)
	  sprintf(tmpstr,"RADIANCE_MINIMUM_BAND_6_VCID_1");
    else
	  sprintf(tmpstr,"RADIANCE_MINIMUM_BAND_%d",bno);
	sprintf(tmpstr2,"LMIN_BAND%d",bno);
	if(strcmp(label, tmpstr) == 0 || strcmp(label, tmpstr2) == 0)
	  sscanf(tokenptr, "%f", &(lmin[i]));

    if (bno == 61)
	  sprintf(tmpstr,"QUANTIZE_CAL_MIN_BAND_6_VCID_1");
    else
	  sprintf(tmpstr,"QUANTIZE_CAL_MIN_BAND_%d",bno);
	sprintf(tmpstr2,"QCALMIN_BAND%d",bno);
	if(strcmp(label, tmpstr) == 0 || strcmp(label, tmpstr2) == 0)
	  sscanf(tokenptr, "%f", &(qmin[i]));

    if (bno == 61)
	  sprintf(tmpstr,"QUANTIZE_CAL_MAX_BAND_6_VCID_1");
    else
	  sprintf(tmpstr,"QUANTIZE_CAL_MAX_BAND_%d",bno);
	sprintf(tmpstr2,"QCALMAX_BAND%d",bno);
	if(strcmp(label, tmpstr) == 0 || strcmp(label, tmpstr2) == 0)
	  sscanf(tokenptr, "%f", &(qmax[i]));
      }            
      /* in case label (key words) is no the first word in a line */
      label = tokenptr;
    }
  }

  /* start and end */
  sprintf(meta_list[0], "%s\n", meta_lut[0][0]);
  sprintf(meta_list[33], "%s\n", meta_lut[33][0]);
  chkFlg[0] = 1; chkFlg[33] = 1;

  /* some parameters couldn't find from metadata file and are hard-coded for Landsat 7 only */
  sprintf(meta_list[2], "%s = UMD\n", meta_lut[2][0]);   // Only accept UMD/CCRS at this point
  sprintf(meta_list[9], "%s = 2\n", meta_lut[9][0]);     // WRS system
  sprintf(meta_list[12], "%s = 6\n", meta_lut[12][0]);   // NBANDS
  sprintf(meta_list[13], "%s = 1, 2, 3, 4, 5, 7\n", meta_lut[13][0]);  //BANDS
  sprintf(meta_list[18], "%s = 1\n", meta_lut[18][0]);   // NBAND_TH
  sprintf(meta_list[19], "%s = 6\n", meta_lut[19][0]);   // BANDS_TH
  sprintf(meta_list[28], "%s = 12\n", meta_lut[28][0]);  // PROJECTION_SPHERE (WGS84=12 in GCTPC) 
  chkFlg[9] = 1;  chkFlg[12] = 1;
  chkFlg[13] = 1; chkFlg[18] = 1;
  chkFlg[19] = 1; chkFlg[28] = 1;

  /* write band specific parameters */
  /* if defined lmin/lmax/qcalmin/qcalmax then use it instead of gain setting */
  if(chkFlg[29] == 1) {
    /* disable gain settings */
    chkFlg[14] = -1; chkFlg[20] = -1;
    /* use actual defined values */
    for(i=1; i<=7; i++) {
      gain[i] = (lmax[i]-lmin[i])/(qmax[i]-qmin[i]);
      bias[i] = lmin[i]-gain[i]*qmin[i];
    }
    sprintf(meta_list[29],"%s = %f, %f, %f, %f, %f, %f\n", meta_lut[29][0], 
	    gain[1], gain[2], gain[3], gain[4], gain[5], gain[7]);
    sprintf(meta_list[30],"%s = %f, %f, %f, %f, %f, %f\n", meta_lut[30][0], 
	    bias[1], bias[2], bias[3], bias[4], bias[5], bias[7]);
    sprintf(meta_list[31],"%s = %f\n", meta_lut[31][0], gain[6]);
    sprintf(meta_list[32],"%s = %f\n", meta_lut[32][0], bias[6]);
  }
  else {
    sprintf(meta_list[14],"%s = %s, %s, %s, %s, %s, %s\n", meta_lut[14][0], band_gain[1],band_gain[2],
	    band_gain[3],band_gain[4],band_gain[5],band_gain[7]);
    /* use HIGH or LOW gain setting and set gain and bias flag invalid */
    chkFlg[29] = -1; chkFlg[30] = -1;
    chkFlg[31] = -1; chkFlg[32] = -1;
  } 
  sprintf(meta_list[17],"%s = %s, %s, %s, %s, %s, %s\n", meta_lut[17][0], band_fname[1],band_fname[2],
	  band_fname[3],band_fname[4],band_fname[5],band_fname[7]); 

  /* write coordinates of upleft corner */ 
  /*if(strstr(processing_software, "LPGS"))*/
  /* always convert to the UL coordinate system - 4/7/09 */
  sprintf(meta_list[26],"%s = %f, %f\n", 
	  meta_lut[26][0], atof(ul_xy[0])-0.5*res, atof(ul_xy[1])+0.5*res); 
  /*else if(strstr(processing_software, "NLAPS"))
    sprintf(meta_list[26],"%s = %s, %s\n", meta_lut[26][0], ul_xy[0], ul_xy[1]); 
    else {
    WARN("Unknown processing software (neither LPGS nor NLAPS), please check coordinate system");
    } */

  if(geotiff == 1) {
    if(readULfromGeoTiff(band_fname[1], &nrows, &ncols, ulxy, &res) == SUCCESS) {
      sprintf(meta_list[15],"%s = %d\n", meta_lut[15][0], ncols);
      sprintf(meta_list[16],"%s = %d\n", meta_lut[16][0], nrows);
      sprintf(meta_list[25],"%s = %f\n", meta_lut[25][0], res);
      sprintf(meta_list[26],"%s = %f, %f\n", meta_lut[26][0], ulxy[0], ulxy[1]); 
    }
  }
  
  /* some input have production year prior to acquisition date and need special care here */
  /* get year from acquisition date */
  sscanf(acquisition_date, "%d", &year);
  /* get year from production date */
  sscanf(production_date, "%d", &pyear);
  if (pyear<year) strcpy(production_date, acquisition_date);
  sprintf(meta_list[6],"%s = %s\n", meta_lut[6][0], production_date);

  /* write all paramters to input file */
  for(i=0; i<34; i++) {
    if(chkFlg[i]==1) {
      fprintf(out, "%s", meta_list[i]);
      LOG("FOUND: %s", meta_list[i]);
    }
    else if(chkFlg[i]==0) {
      fprintf(out, "%s = \n", meta_lut[i][0]);
      WARN("Metadata %s was not found in metadat. \tPlease enter it in %s manually", meta_lut[i][0], met_name);
    }
  }

  fclose(in);
  fclose(out);

  if(writeENVIHeader(scene_name, nrows, ncols, ulxy, res, zone)==FAILURE)
    return FAILURE;
  else
    return instrument;
}



/* retrieve metadata from NLAPS format (TM & MSS) */
int getMetaFromNLAPS(char input[], char scene_name[], char acquisition_date[])
{
  int   i, year, pyear, month, day, path, row, zone, geotiff;
  int   nrows, ncols;
  int   chkFlg[32];
  double ulxy[2], res;

  /* vars used in parameter parsing */
  char  fileName[MAX_STRING_LENGTH] = "\0";
  char  buffer[MAX_STRING_LENGTH] = "\0";
  char  tmpstr[MAX_STRING_LENGTH] = "\0";
  char  production_date[MAX_STRING_LENGTH] = "\0";
  char  processing_software[MAX_STRING_LENGTH] = "\0";
  char  sensor[MAX_STRING_LENGTH] = "\0";
  char  met_name[MAX_STRING_LENGTH] = "\0";
  char  *label = NULL;
  char  *tokenptr = NULL;
  char  *seperator = ",;=\" \t";
  float fnum;
  float gain[8], bias[8];
  /*float tm_band_width[8] = {0, 0.066, 0.082, 0.067, 0.128, 0.217, 1.0, 0.252};*/

  /* metadata */
  char  meta_list[32][MAX_STRING_LENGTH];
  char  band_fname[8][MAX_STRING_LENGTH];
  char  ul_xy[2][MAX_STRING_LENGTH];
  char  *meta_lut[32][2] = {{"HEADER_FILE",       ""},
			    {"FILE_TYPE",         ""},
			    {"DATA_PROVIDER",     "DATA_SET_TYPE"},
			    {"SATTELITE",         "SATELLITE"},
			    {"INSTRUMENT",        "SATELLITE_INSTRUMENT"},
			    {"ACQUISITION_DATE",  "ACQUISITION_DATE/TIME"},
			    {"PRODUCTION_DATE",   "PROCESSING_DATE/TIME"},
			    {"SOLAR_ZENITH",      "SUN_ELEVATION"},
			    {"SOLAR_AziMUTH",     "SUN_AZIMUTH"},
			    {"WRS_SYSTEM",        ""},
			    {"WRS_PATH",          ""},
			    {"WRS_ROW",           ""},
			    {"NBAND",             ""},
			    {"BANDS",             ""},
			    {"GAIN",              ""},
			    {"BIAS",              ""},
			    {"NSAMPLE",           "PIXELS_PER_LINE"},
			    {"NLINE",             "LINES_PER_DATA_FILE"},
			    {"FILE_NAMES",        ""},
			    {"NBAND_TH",          ""},
			    {"BANDS_TH",          ""},
			    {"GAIN_TH",           "BAND6_RADIOMETRIC_GAINS/BIAS"},
			    {"BIAS_TH",           ""},
			    {"NSAMPLE_TH",        "PIXELS_PER_LINE"},
			    {"NLINE_TH",          "LINES_PER_DATA_FILE"},
			    {"FILE_NAMES_TH",     "BAND6_FILENAME"},
			    {"PROJECTION_NUMBER", "USGS_PROJECTION_NUMBER"},
			    {"PIXEL_SIZE",        "PIXEL_SPACING"},
			    {"UPPER_LEFT_CORNER", "UPPER_LEFT_CORNER"},
			    {"PROJECTION_ZONE",   "USGS_MAP_ZONE"},
			    {"PROJECTION_SPHERE", ""},
			    {"END",               "" }}; 
  FILE *in, *out;
  

  if((in=fopen(input, "r"))==NULL)
    DIE( "Can't open input metadata file %s", input);


  sprintf(met_name, "%s.carbon_met.txt", scene_name);
  if((out=fopen(met_name, "w"))==NULL)
    DIE( "Can't open output header file for write %s", met_name);

  for(i=0; i<32; i++) chkFlg[i] = 0;

  sprintf(meta_list[14], "%s = ", meta_lut[14][0]);  // GAIN
  sprintf(meta_list[15], "%s = ", meta_lut[15][0]);  // BIAS

  /* process line by line */
  while(fgets(buffer, MAX_STRING_LENGTH, in) != NULL) {

    /* get string token */
    tokenptr = strtok(buffer, seperator);
    label=tokenptr;

    while(tokenptr != NULL) {
 
      tokenptr = strtok(NULL, seperator);

      /* information provided just after label */
      for(i=1; i<32; i++)
	if(strcmp(label, meta_lut[i][1]) == 0) {
	  chkFlg[i] = 1;
	  switch(i) {      
	    /* some metadata need more work to satisfy input */
	  case 4:
	    strcpy(sensor, tokenptr);
	    if(strstr(tokenptr,"ETM") != NULL)
	      sprintf(meta_list[i],"%s = ETM\n", meta_lut[i][0]);
	    else
	      if(strstr(tokenptr,"TM") != NULL)
		sprintf(meta_list[i],"%s = TM\n", meta_lut[i][0]);
	      else
		if(strstr(tokenptr,"MSS") != NULL)
		  sprintf(meta_list[i],"%s = MSS\n", meta_lut[i][0]);
	    break;
	  case 5:    // ACQUISITION DATE
	    sscanf(tokenptr, "%10s", acquisition_date);
	    sprintf(meta_list[i], "%s = %s\n", meta_lut[i][0], acquisition_date);
	    break;
	  case 6:    // PRODUCTION_DATE
	    // some input give wrong production date and need extra processing   
	    sscanf(tokenptr, "%10s", production_date);
	    sprintf(meta_list[i], "%s = %s\n", meta_lut[i][0], production_date);
	    break;
	  case 7:    // SOLAR_ZENITH
	    sscanf(tokenptr, "%f", &fnum);
	    sprintf(meta_list[i], "%s = %5.2f\n", meta_lut[i][0], 90-fnum);
	    break;
	  case 28:  // upper left corner
	    tokenptr = strtok(NULL, seperator);
	    tokenptr = strtok(NULL, seperator);
	    sscanf(tokenptr, "%s", ul_xy[0]);
	    tokenptr = strtok(NULL, seperator);
	    sscanf(tokenptr, "%s", ul_xy[1]);
	    break;
	  case 29:
	    sscanf(tokenptr, "%d", &zone);
	    sprintf(meta_list[i], "%s = %d\n", meta_lut[i][0], zone);
	    break;	    
	  default:  // other parameters can be used directly
	    sscanf(tokenptr, "%s", tmpstr);
	    sprintf(meta_list[i],"%s = %s\n", meta_lut[i][0], tmpstr);
	  }
	}
 
      /* get WRS path and row */
      if(strcmp(label, "WRS") == 0) {
	sscanf(tokenptr, "%3d/%3d", &path, &row);
	sprintf(meta_list[10], "%s = %d\n", meta_lut[10][0], path);
	sprintf(meta_list[11], "%s = %d\n", meta_lut[11][0], row);
	chkFlg[10] = 1; chkFlg[11] = 1;
      }

      if(strcmp(label, "PROCESSING_SOFTWARE") == 0) 
	sscanf(tokenptr, "%s", processing_software);
 
      /* extract gain and file names for band 1-7 */
      for(i=1; i<=7; i++) {
	// extract band filename
	sprintf(tmpstr,"BAND%d_FILENAME",i);
	if(strcmp(label, tmpstr) == 0) {
	  chkFlg[18] = 1;	    
	  sscanf(tokenptr, "%s", band_fname[i]);
	  strcpy(fileName, tokenptr);
	}
	sprintf(tmpstr,"BAND%d_RADIOMETRIC_GAINS/BIAS",i);
	if(strcmp(label, tmpstr) == 0) {
	  chkFlg[14] = 1;	    
	  // extract GAIN
	  sscanf(tokenptr, "%f", &(gain[i]));
	  tokenptr = strtok(NULL, seperator);
	  chkFlg[15] = 1;
	  // extract BIAS
	  sscanf(tokenptr, "%f", &(bias[i]));
	  if (i==6) {chkFlg[21] = 1; chkFlg[22] = 1;}
	}
      }
      
      /* in case label (key words) is no the first word in a line */
      label = tokenptr;
    }
  }
  
  /* start and end */
  sprintf(meta_list[0], "%s\n", meta_lut[0][0]);
  sprintf(meta_list[31], "%s\n", meta_lut[31][0]);
  chkFlg[0] = 1;  chkFlg[31] = 1;

  /* some parameters couldn't find from metadata file and are hard-coded for Landsat 7 only */
  /* check file name and see if it is a TIFF file */
  if(strstr(fileName, ".tif") || strstr(fileName, ".TIF")) {
    sprintf(meta_list[1], "%s = GEOTIFF\n", meta_lut[1][0]);   // DATA_FORMAT
    geotiff = 1;
  }
  else {
    sprintf(meta_list[1], "%s = BINARY\n", meta_lut[1][0]);   // DATA_FORMAT
    geotiff = 0;
  }
  sprintf(meta_list[2], "%s = UMD\n", meta_lut[2][0]);   // Only accept UMD/CCRS at this point      
  sprintf(meta_list[19], "%s = 1\n", meta_lut[19][0]);   // NBAND_TH
  sprintf(meta_list[20], "%s = 6\n", meta_lut[20][0]);   // BANDS_TH
  sprintf(meta_list[30], "%s = 12\n", meta_lut[30][0]);  // PROJECTION_SPHERE (WGS84=12 in GCTPC) 
  chkFlg[1] = 1;   chkFlg[9] = 1;
  chkFlg[12] = 1;  chkFlg[13] = 1;
  chkFlg[19] = 1;  chkFlg[20] = 1;
  chkFlg[30] = 1; 

  /* write band specific parameters */
  if(strcmp(sensor, "MSS") == 0) {
    sprintf(meta_list[9], "%s = 1\n", meta_lut[9][0]);     // WRS system
    sprintf(meta_list[12], "%s = 4\n", meta_lut[12][0]);   // NBANDS
    sprintf(meta_list[13], "%s = 1, 2, 3, 4\n", meta_lut[13][0]);  //BANDS
    sprintf(meta_list[14],"%s = %f, %f, %f, %f\n", meta_lut[14][0], gain[1], gain[2], gain[3], gain[4]); 
    sprintf(meta_list[15],"%s = %f, %f, %f, %f\n", meta_lut[15][0], bias[1], bias[2], bias[3], bias[4]); 
    sprintf(meta_list[18],"%s = %s, %s, %s, %s\n", meta_lut[18][0], band_fname[1],band_fname[2], band_fname[3],band_fname[4]); 
  }
  else {
    sprintf(meta_list[9], "%s = 2\n", meta_lut[9][0]);     // WRS system
    sprintf(meta_list[12], "%s = 6\n", meta_lut[12][0]);   // NBANDS
    sprintf(meta_list[13], "%s = 1, 2, 3, 4, 5, 7\n", meta_lut[13][0]);  //BANDS
    sprintf(meta_list[14],"%s = %f, %f, %f, %f, %f, %f\n", meta_lut[14][0], gain[1], gain[2],
	    gain[3], gain[4], gain[5], gain[7]); 
    sprintf(meta_list[15],"%s = %f, %f, %f, %f, %f, %f\n", meta_lut[15][0], bias[1], bias[2],
	    bias[3], bias[4], bias[5], bias[7]); 
    sprintf(meta_list[18],"%s = %s, %s, %s, %s, %s, %s\n", meta_lut[18][0], band_fname[1],band_fname[2],
	    band_fname[3],band_fname[4],band_fname[5],band_fname[7]); 
  } 

  /* thermal band gain/bias */
  sprintf(meta_list[21], "%s = %f\n", meta_lut[21][0], gain[6]);
  sprintf(meta_list[22], "%s = %f\n", meta_lut[22][0], bias[6]);

  /* write coordinates of upleft corner */
  /* printf("processing_software = %s\n", processing_software);*/
  /*if(strstr(processing_software, "LPGS"))*/
  /* always convert to the UL coordinate system - 4/7/09 */
  sprintf(meta_list[28],"%s = %f, %f\n", 
	  meta_lut[28][0], atof(ul_xy[0])-0.5*res, atof(ul_xy[1])+0.5*res); 
  /*else if(strstr(processing_software, "NLAPS"))
    sprintf(meta_list[28],"%s = %s, %s\n", meta_lut[28][0], ul_xy[0], ul_xy[1]); 
    else {
    WARN("Unknown processing software (neither LPGS nor NLAPS), please check coordinate system");
    }*/
  /* printf("%s\n", meta_list[28]);*/

  if(geotiff == 1)  {
    if(readULfromGeoTiff(band_fname[1], &nrows, &ncols, ulxy, &res) == SUCCESS) {
      sprintf(meta_list[16],"%s = %d\n", meta_lut[16][0], ncols);
      sprintf(meta_list[17],"%s = %d\n", meta_lut[17][0], nrows);
      sprintf(meta_list[27],"%s = %f\n", meta_lut[27][0], res);
      sprintf(meta_list[28],"%s = %f, %f\n", meta_lut[28][0], ulxy[0], ulxy[1]); 
    }
    /*  printf("%s\n", meta_list[28]);*/
  }

  /* some input have production year prior to acquisition date and need special care here */
  if(strcmp(sensor, "TM") == 0) {
    if(strstr(acquisition_date, "/")) {
      sscanf(acquisition_date, "%2d%2d%2d", &month, &day, &year);
      sprintf(acquisition_date, "19%2d-%02d-%02d", year, month, day);
    }
    if(strstr(production_date, "/")) {
      sscanf(production_date, "%2d%2d%2d", &month, &day, &year);
      sprintf(production_date, "19%2d-%02d-%02d", year, month, day);
    }
  }

  /* get year from acquisition date */
  sscanf(acquisition_date, "%d", &year);
  /* get year from production date */
  sscanf(production_date, "%d", &pyear);
  if (pyear<year) strcpy(production_date, acquisition_date);
  sprintf(meta_list[5],"%s = %s\n", meta_lut[5][0], acquisition_date);
  sprintf(meta_list[6],"%s = %s\n", meta_lut[6][0], production_date);

  /* write all paramters to input file */
  if(strcmp(sensor, "MSS") == 0) {
    for(i=0; i<32; i++) {
      if(i>=19 && i<=25) continue;
      if(chkFlg[i]) {
	fprintf(out, "%s", meta_list[i]);
	LOG("FOUND: %s", meta_list[i]);
      }
      else {
	fprintf(out, "%s = \n", meta_lut[i][0]);
	WARN("Metadata %s was not found in metadat. \tPlease enter it in %s manually", meta_lut[i][0], met_name);
      }
    }
  }
  else {
    for(i=0; i<32; i++) {
      if(chkFlg[i]) {
	fprintf(out, "%s", meta_list[i]);
	LOG("FOUND: %s", meta_list[i]);
      }
      else {
	fprintf(out, "%s = \n", meta_lut[i][0]);
	WARN("Metadata %s was not found in metadat. \tPlease enter it in %s manually", meta_lut[i][0], met_name);
      }
    }
  }

  fclose(in);
  fclose(out);

  if(writeENVIHeader(scene_name, nrows, ncols, ulxy, res, zone)==FAILURE)
    return FAILURE;

  if(strstr(sensor, "ETM")) return ETM;
  else if(strstr(sensor, "TM")) return TM;
  else if(strstr(sensor, "MSS")) return MSS;
  else return FAILURE;
}

/* given an id for a carbon_met param name, return ptr to the name */
char *cm_name_for_id(int id)
{
  int ii = 0;
  static struct cm_symbol {
    int  id;
    char *name;
  } symbols[] = {
    {FILE_TYPE,          "FILE_TYPE"},
    {DATA_PROVIDER,      "DATA_PROVIDER"},
    {SATELLITE,          "SATTELITE"},          /* yes, this is misspelled */
    {INSTRUMENT,         "INSTRUMENT"},
    {ACQUISITION_DATE,   "ACQUISITION_DATE"},
    {PRODUCTION_DATE,    "PRODUCTION_DATE"},
    {SOLAR_ZENITH,       "SOLAR_ZENITH"},
    {SOLAR_AZIMUTH,      "SOLAR_AziMUTH"},      /* yes, case is correct */
    {WRS_SYSTEM,         "WRS_SYSTEM"},
    {WRS_PATH,           "WRS_PATH"},
    {WRS_ROW,            "WRS_ROW"},
    {NBAND,              "NBAND"},
    {BANDS,              "BANDS"},
    {GAIN,               "GAIN"},
    {BIAS,               "BIAS"},
    {NSAMPLE,            "NSAMPLE"},
    {NLINE,              "NLINE"},
    {FILE_NAMES,         "FILE_NAMES"},
    {NBAND_TH,           "NBAND_TH"},
    {BANDS_TH,           "BANDS_TH"},
    {GAIN_TH,            "GAIN_TH"},
    {BIAS_TH,            "BIAS_TH"},
    {NSAMPLE_TH,         "NSAMPLE_TH"},
    {NLINE_TH,           "NLINE_TH"},
    {FILE_NAMES_TH,      "FILE_NAMES_TH"},
    {PROJECTION_NUMBER,  "PROJECTION_NUMBER"},
    {PIXEL_SIZE,         "PIXEL_SIZE"},
    {UPPER_LEFT_CORNER,  "UPPER_LEFT_CORNER"},
    {PROJECTION_ZONE,    "PROJECTION_ZONE"},
    {PROJECTION_SPHERE,  "PROJECTION_SPHERE"}
  };
  int n_syms = sizeof(symbols) / sizeof(struct cm_symbol);
  
  
  for (ii = 0; ii < n_syms; ++ii)
  {
    if (symbols[ii].id != id) continue;
    return symbols[ii].name;
  }
  
  DIE( "No symbol for id=%d", id);
  return NULL;
}

char *cm_val_for_id(int id, METADATA *metadata)
{
  if (id < 1 || id >= MAX_CM_PARAMS)
    DIE("invalid id=%d", id);
  return metadata[id].val;
}

/* retrieve metadata from NLAPS format (TM & MSS) */
#define IGNORE_STATE (-1)
#define NVP_STATE 1
#define GAIN_BIAS_STATE 2
int getMetaFromNLAPS_WO(char input[], char scene_name[], char acquisition_date[])
{
  struct stat stbuf;
  FILE  *in;
  int   ret = FAILURE;
  int   nrows, ncols;
  double ulxy[2], res;
  int band=0, detector=0, Max_N_Bands = MAX_N_BANDS;
  float gain=0.0, bias=0.0;

  /* vars used in parameter parsing */
  char  buffer[MAX_STRING_LENGTH] = "\0";
  char *str1, *token;
  char *line = (char *)NULL;
  size_t line_len = 0;
  int state = 0;
  int n_read = 0;
  
  /* metadata */
  METADATA *met_list = (METADATA *)NULL;
  char  band_fname[8][MAX_STRING_LENGTH];
  LUT lut[] = {
/*    cm id             W0 param name     n fields */
/*    -------           -------------     -------- */
    {FILE_TYPE,         "Product Format",     1},
    {DATA_PROVIDER,     "DataProvider",       1},
    {SATELLITE,         "Satellite",          1},
    {INSTRUMENT,        "Sensor",             1},
    {ACQUISITION_DATE,  "Scene center date",  3},
    {PRODUCTION_DATE,   "Completion date",    3},
    {SOLAR_ZENITH,      "Sun Elevation",      2},
    {SOLAR_AZIMUTH,     "Sun Azimuth",        2},
    {WRS_SYSTEM,        "WrsSystem",          1},
    {WRS_PATH,          "Path/Strip no.",     1},
    {WRS_ROW,           "Start Row no.",      1},
    {WRS_ROW,           "Row no.",            1},
    {NBAND,             "NBands",             1},
    {BANDS,             "Bands",              1},
    {GAIN,              "BandGain",           0},
    {BIAS,              "BandBias",           0},
    {NSAMPLE,           "N_Cols",             1},
    {NLINE,             "N_Rows",             1},
    {FILE_NAMES,        "BandFileName",       1},
    {NBAND_TH,          "NThermalBands",      1},
    {BANDS_TH,          "ThermalBand",        1},
    {GAIN_TH,           "ThermalGain",        0},
    {BIAS_TH,           "ThermalBias",        0},
    {NSAMPLE_TH,        "N_THCols",           1},
    {NLINE_TH,          "N_THRows",           1},
    {FILE_NAMES_TH,     "THBandFileName",     1},
    {PROJECTION_NUMBER, "Map Projection",     1},
    {PIXEL_SIZE,        "PixelSize",          1},
    {UPPER_LEFT_CORNER, "UpperLeftCorner",    1},
    {PROJECTION_ZONE,   "Zone",               1},
    {PROJECTION_SPHERE, "Earth Ellipsoid",    1},
    {0,                 "Image Orientation",  4},
    {0,                 "Line Spacing",       1},
    {0,                 "Algorithm",          2},
    {0,                 "GainBias",           4},
    {0,                 (char *)NULL,         0}
  };  
  
  /* allocate metadata list */
  met_list = calloc(MAX_CM_PARAMS, sizeof(METADATA));
  if (! met_list)
    DIE("getMetaFromNLAPS_WO : could not allocate %d btytes", MAX_CM_PARAMS * sizeof(METADATA));

  /*
   * Need to know if this is an MSS file!
   * open input _W0 file...
   */
  if((in=fopen(input, "r"))==NULL) {
    DIE("Can't open input metadata file %s", input);
  }
  while(getline(&line, &line_len, in) != -1)
  {
    if (strstr(line, "Sensor:")
    && strstr(line, "MSS"))
    {
      Max_N_Bands = N_BANDS_MSS;
      break;
    }
  }
  fclose(in);        /* done checking MSS */
    

  /*
   * some params are always the same, or are not in the _W0 file,
   * so preload them...
   */
  nvp_to_meta("WrsSystem", "2", lut, met_list);
  nvp_to_meta("DataProvider", "UMD", lut, met_list);
  if (Max_N_Bands > N_BANDS_MSS)  /* not an MSS file */
  {
    nvp_to_meta("NBands", "6", lut, met_list);
    nvp_to_meta("Bands", "1, 2, 3, 4, 5, 7", lut, met_list);
    nvp_to_meta("ThermalBand", "6", lut, met_list);
    nvp_to_meta("NThermalBands", "1", lut, met_list);
  }
  else              /* IS an MSS file */
  {
    nvp_to_meta("NBands", "4", lut, met_list);
    nvp_to_meta("Bands", "1, 2, 3, 4", lut, met_list);
    /* and these fields are not used in carbon.met file... */
    lut[NBAND_TH].cm_id = 0;
    lut[BANDS_TH].cm_id = 0;
    lut[GAIN_TH].cm_id = 0;
    lut[BIAS_TH].cm_id = 0;
    lut[NSAMPLE_TH].cm_id = 0;
    lut[NLINE_TH].cm_id = 0;
    lut[FILE_NAMES_TH].cm_id = 0;
  }
  
  /* band filenames */
  strcpy(buffer, input);
  token = strtok(buffer, "_");
  
  /*
   * unfortunately, filenaming is inconsistent. Files can be
   * Bn.tif, Bn.TIF, Lsomething_Bn.tif, Lsomething_Bn.TIF
   * so try them all before giving up...
   */
  for (band = 1; band <= Max_N_Bands; ++band)
  {
    /* try each combination of the band filename... */
    sprintf(band_fname[band], "%s_B%d.tif", token, band);
    if (stat(band_fname[band], &stbuf) != 0)
    {
      sprintf(band_fname[band], "%s_B%d.TIF", token, band);
      if (stat(band_fname[band], &stbuf) != 0)
      {
        sprintf(band_fname[band], "B%d.tif", band);
        if (stat(band_fname[band], &stbuf) != 0)
        {
	  sprintf(band_fname[band], "B%d.TIF", band);
          if (stat(band_fname[band], &stbuf) != 0)
          {
            DIE("Could not stat file %s", band_fname[band]);
          }
        }
      }
    }

    /* and add it to metadata */
    if (band == 6)       /* thermal band is special */
      nvp_to_meta("THBandFileName", band_fname[band], lut, met_list);
    else
      nvp_to_meta("BandFileName", band_fname[band], lut, met_list);
  }
  
  /*
   * corner points are hard to extract from _W0,
   * so get them from TIF file
   */
  if(readULfromGeoTiff(band_fname[1], &nrows, &ncols, ulxy, &res) == SUCCESS)
  {
    sprintf(buffer,"%d", ncols);
    nvp_to_meta("N_Cols", buffer, lut, met_list);
    sprintf(buffer,"%d", nrows);
    nvp_to_meta("N_Rows", buffer, lut, met_list);
    sprintf(buffer,"%f", res);
    nvp_to_meta("PixelSize", buffer, lut, met_list);
    sprintf(buffer,"%f, %f", ulxy[0], ulxy[1]); 
    nvp_to_meta("UpperLeftCorner", buffer, lut, met_list);
  }
  else
    DIE("Could not read corner points from %s", band_fname[1]);

  /*
   * temperature sometimes has a diffenent size/resolution,
   * so get cols and rows from TIF file
   */
  if (Max_N_Bands > N_BANDS_MSS)  /* not an MSS file */
  {
    if (readULfromGeoTiff(band_fname[6], &nrows, &ncols, ulxy, &res) == SUCCESS)
    {
      sprintf(buffer,"%d", ncols);
      nvp_to_meta("N_THCols", buffer, lut, met_list);
      sprintf(buffer,"%d", nrows);
      nvp_to_meta("N_THRows", buffer, lut, met_list);
      sprintf(buffer,"%f", res);
    }
    else
    {
      DIE("Could not read corner points from %s", band_fname[6]);
    }
  }

  /* process input _W0 file line by line */
  if((in=fopen(input, "r"))==NULL) {
    DIE("Can't open input metadata file %s", input);
  }
  state = NVP_STATE;
  while(getline(&line, &line_len, in) != -1)
  {
       switch (state)
       {
       case NVP_STATE:
         for (str1 = line; ; str1 = NULL)
         {
           /*
	    * keyword
	    */	   
           token = strtok(str1, ":\n");
	   if (token == NULL)
             break;
	   while (*token == ' ') ++token;        /* remove leading whitespace */
	   str1 = (char *)NULL;
	   
	   /*
	    * value(s)...
	    */
	   n_read = read_n_fields((char *)NULL, token, lut, buffer);
	   if (n_read > 0)
             nvp_to_meta(token, buffer, lut, met_list);
	   
	   if ((strcmp(token, "Algorithm") == 0)
	   && ((strcmp(buffer, "NASA") == 0)
	   ||  (strcmp(buffer, "NASA CPF") == 0)
	   ||  (strcmp(buffer, "CCRS") == 0)))
	   {
	     state = GAIN_BIAS_STATE;
	     break;
	   }

         }
        break;
      case GAIN_BIAS_STATE:
        /*
         * value(s)
         */
        n_read = read_n_fields(line, "GainBias", lut, buffer);
        if (n_read < 4)
          break;
	   
        /* handle gain/bias for each band */
	if (! isdigit(buffer[0]))
	  break;
	n_read = sscanf(buffer, "%d %d %f %f", &band, &detector, &gain, &bias);
	if (n_read != 4)
	  DIE("could not parse line=%s", line);
	
	if (band == 6)       /* thermal band */
	{
	  sprintf(buffer, "%f", gain);
          nvp_to_meta("ThermalGain", buffer, lut, met_list);
	  sprintf(buffer, "%f", bias);
          nvp_to_meta("ThermalBias", buffer, lut, met_list);
	}
	else                /* reflectance band */
	{
	  sprintf(buffer, "%f", gain);
          nvp_to_meta("BandGain", buffer, lut, met_list);
	  sprintf(buffer, "%f", bias);
          nvp_to_meta("BandBias", buffer, lut, met_list);
	}
        if (band >= Max_N_Bands) state = IGNORE_STATE;
        break;
      }
       
  }
  free(line);       /* getline mallocs a buffer, so free it */
  fclose(in);        /* done reading metadata */
  
  /*
   * create the carbon_met file
   */
  create_carbon_met_file(scene_name, met_list);
  
  /* set the aquisition date for the rest of the program */
  strcpy(acquisition_date, cm_val_for_id(ACQUISITION_DATE, met_list));
  
  /* write the ENVI header file */
  if (create_ENVI_hdr_file(scene_name, met_list)==FAILURE)
    return FAILURE;
    
  /* set the return value so the rest of the program knows what to do */
  token = cm_val_for_id(INSTRUMENT, met_list);
  if(strstr(token, "ETM")) ret = ETM;
  else if(strstr(token, "TM")) ret = TM;
  else if(strstr(token, "MSS")) ret = MSS;
  
  free(met_list);   /* we allocated this, so free it */
  return ret;
}

int read_n_fields(char *str, char *key, LUT *lut, char container[MAX_STRING_LENGTH])
{
  int ii = 0;
  int n_fields = 0;
  char *value = (char *)NULL;
 
  /* find how many fields to read in */
  for ( ; lut->in_name; lut++ )
  {
    if (strcmp(key, lut->in_name) == 0)
    {
      n_fields = lut->n_fields;
      break;
    }
  }
  
  /* and read in the fields */
  container[0] = 0;                       /* make container empty */
  if (n_fields == 0) return 0;            /* don't read anything if not wanted */
  for (ii = 0; ii < n_fields; ++ii, str=NULL)
  {
    value = strtok(str, " |\t\n");
    if (value == NULL)
      break;
    if (strlen(container) == 0)
      strcpy(container, value);
    else
    {
      strcat(container, " ");
      strcat(container, value);
    }
  }
  
  return ii;
}

int nvp_to_meta(char *name, char *val, LUT *lut, METADATA *carbon_met)
{
  int index = 0;
  METADATA *entry = (METADATA *)NULL;
  
  for ( ; lut->in_name; lut++ )
  {
    if (strcmp(name, lut->in_name) != 0) continue;

    /* found a match -- assign or append the nvp to the corresponding metadata */
    index = lut->cm_id;
    if (index < 1 || index >= MAX_CM_PARAMS)
      return FAILURE;               /* some things are required to be ignored */
    entry = &(carbon_met[index]);
    
    /* name */
    strcpy(entry->name, cm_name_for_id(index));
    
    /* val */
    if (strlen(entry->val) == 0)
      strcpy(entry->val, val);
    else                     /* already some stuff there, so append */
    {
      strcat(entry->val, DELIMIT);
      strcat(entry->val, val);
    }

    /*
     * if this is solar elevation angle, flag it so we can later
     * convert it to solar zenith angle
     */
    if (strcmp(name, "Sun Elevation") == 0)
    {
      strcat(entry->val, " Elevation");
    }
    if (strlen(entry->val) >= MAX_STRING_LENGTH)
      DIE("val too long: '%s'\n", entry->val );
  }
 
  return SUCCESS;
}

int create_carbon_met_file(char *scene_name, METADATA *carbon_met)
{
  int ii = 0;
  FILE *fh = (FILE *)NULL;
  char met_fname[MAX_STRING_LENGTH] = "";
  
  /* build carbon_met filename from scene name */
  sprintf(met_fname, "%s.carbon_met.txt", scene_name);
  
  fh = fopen(met_fname, "w");
  if (! fh)
    DIE( "could not open %s for write", met_fname);
  
  fprintf(fh, "HEADER_FILE\n");
  for (ii = 1; ii < MAX_CM_PARAMS; ++ii)
  {
    if (validate_carbon_met(carbon_met, ii) != SUCCESS)
      DIE( "param %d is invalid", ii);
    LOG("FOUND: %s=%s\n", carbon_met[ii].name, carbon_met[ii].val);
    if (strlen(carbon_met[ii].name) && strlen(carbon_met[ii].val))
    {
      fprintf(fh, "%s = %s\n", carbon_met[ii].name, carbon_met[ii].val);
    }
  }
  fprintf(fh, "END\n");
  
  fclose(fh);
  
  return SUCCESS;
}

int create_ENVI_hdr_file(char *scene_name, METADATA *carbon_met)
{
  int  ii=0, nn=0, zone=0;
  char tmpstr[MAX_STRING_LENGTH] = "\0";
  char *ledaps_file[4]={"lndcal", "lndcsm", "lndth", "lndsr"};
  char *hemisphere = "North";
  FILE *out = (FILE *)NULL;

  nn = sizeof(ledaps_file) / sizeof(char*);
  for(ii=0; ii<nn; ii++) {
    sprintf(tmpstr, "%s.%s.hdf.hdr", ledaps_file[ii], scene_name);
    if((out=fopen(tmpstr, "w"))==NULL) {
      DIE( "Write ENVI header file %s error!", tmpstr);
    }
    fprintf(out,"ENVI\n");
    fprintf(out,"description = {LEDAPS HDF File Imported into ENVI}\n");
    fprintf(out,"samples = %s\n", carbon_met[NSAMPLE].val);
    fprintf(out,"lines   = %s\n", carbon_met[NLINE].val);
    fprintf(out,"bands   = 1\n");
    fprintf(out,"header offset = 0\n");
    fprintf(out,"file type = HDF Scientific Data\n");
    fprintf(out,"data type = 2\n");
    fprintf(out,"interleave = bsq\n");
    fprintf(out,"sensor type = Landsat\n");
    fprintf(out,"byte order = 0\n");
    
    /* zone must be positive */
    sscanf(carbon_met[PROJECTION_ZONE].val, "%d", &zone);
    if (zone < 0) hemisphere = "South";
    fprintf(out,"map info = {UTM, 1.000, 1.000, %s, %s, %s, %d, %s, WGS-84, units=Meters}\n",
      carbon_met[UPPER_LEFT_CORNER].val,
      carbon_met[PIXEL_SIZE].val,
      carbon_met[PIXEL_SIZE].val,
      abs(zone),
      hemisphere
    );
    fclose(out);
  }
  return SUCCESS;
}

/*
 * make sure parameters are defined and formatted correctly
 *
 * NOTE : CURRENTLY THIS HANDLES METADATA INPUT FROM _W0 FILES. IT SHOULD
 *        BE UPDATED TO HANDLE INPUT FROM ANY MET FILE.
 */
int validate_carbon_met(METADATA *carbon_met, int param_index)
{
  char *str = (char *)NULL;
  char *ptrc = (char *)NULL;
  char buffer[MAX_STRING_LENGTH] = "";
  double fval = 0.0;
  int   ival = 0;
  int   year = 0, month=0, day=0;
  
  if (cm_name_for_id(param_index) == NULL)
  {
    WARN("no name for id=%d", param_index);
    return FAILURE;
  }
  
  str = carbon_met[param_index].val;
  if (! str)
  {
    WARN("null val");
    return FAILURE;
  }
  
  switch (param_index)
  {
    case FILE_TYPE:                /* needs to be upper case */
      for ( ; *str; str++) *str = toupper(*str);
      break;
      
    case DATA_PROVIDER:
      break;
    case SATELLITE:                /* use '_' to separate number */
      sscanf(str, "%[a-zA-Z_-]%d", buffer, &ival);
      if ((ptrc=index(buffer, '-')) != NULL) *ptrc = (char)0;
      if ((ptrc=index(buffer, '_')) != NULL) *ptrc = (char)0;
      sprintf(str, "%s_%d", buffer, ival);
      break;
    case INSTRUMENT:
      break;
    case ACQUISITION_DATE:         /* understands 'yyyy-mm-dd' or 'yyyy mm dd' */
      sscanf(str, "%d%[ -]%d%[ -]%d", &year, buffer, &month, buffer, &day);
      sprintf(str, "%04d-%02d-%02d", year, month, day);
      break;
    case PRODUCTION_DATE:          /* understands 'yyyy-mm-dd' or 'yyyy mm dd' */
      sscanf(str, "%d%[ -]%d%[ -]%d", &year, buffer, &month, buffer, &day);
      sprintf(str, "%04d-%02d-%02d", year, month, day);
      break;
    case SOLAR_ZENITH:             /* filter out units string */
      sscanf(str, "%lf", &fval);
      if (strstr(str, "Elevation") != NULL)
        sprintf(str, "%f", 90.0 - fval);/* convert elevation to zenith angle */
      else
        sprintf(str, "%f", fval);
      break;
    case SOLAR_AZIMUTH:            /* filter out units string */
      sscanf(str, "%lf", &fval);
      sprintf(str, "%f", fval);
      break;
    case WRS_SYSTEM:
      break;
    case WRS_PATH:                 /* make sure not float */
      sscanf(str, "%d", &ival);
      sprintf(str, "%d", ival);
      break;
    case WRS_ROW:                  /* make sure not float */
      sscanf(str, "%d", &ival);
      sprintf(str, "%d", ival);
      break;
    case NBAND:
      break;
    case BANDS:
      break;
    case GAIN:
      break;
    case BIAS:
      break;
    case NSAMPLE:
      break;
    case NLINE:
      break;
    case FILE_NAMES:
      break;
    case NBAND_TH:
      break;
    case BANDS_TH:
      break;
    case GAIN_TH:
      break;
    case BIAS_TH:
      break;
    case NSAMPLE_TH:
      break;
    case NLINE_TH:
      break;
    case FILE_NAMES_TH:
      break;
    case PROJECTION_NUMBER:
      if (strcasecmp(str, "UTM") == 0)      /* convert name to GCTP code */
        sprintf(str, "1");
      break;
    case PIXEL_SIZE:
      break;
    case UPPER_LEFT_CORNER:
      break;
    case PROJECTION_ZONE:
      break;
    case PROJECTION_SPHERE:
      if (strcasecmp(str, "WGS84") == 0)    /* convert name to GCTP code */
        sprintf(str, "12");
      break;
  }
  
  return SUCCESS;
}

/* retrieve metadata provided by LC-ComPS from UMD Landsat scene database */
int getMetaFromLCT(char input[], char scene_name[], char acquisition_date[])
{
  int   i, j, year, pyear, zone, geotiff;
  int   nrows, ncols;
  int   chkFlg[34];
  double ulxy[2], res;

  /* vars used in parameter parsing */
  char  buffer[MAX_STRING_LENGTH] = "\0";
  char  tmpstr[MAX_STRING_LENGTH] = "\0";
  char  str[MAX_STRING_LENGTH] = "\0";
  char  production_date[MAX_STRING_LENGTH] = "\0";
  char  sensor[MAX_STRING_LENGTH] = "\0";
  char  hemi[MAX_STRING_LENGTH] = "\0";
  char  met_name[MAX_STRING_LENGTH] = "\0";
  char  *label = NULL;
  char  *tokenptr = NULL;
  char  *seperator = ",;=\" \t";
  float fnum;
  float gain[8], bias[8];
 
  /* metadata */
  char  meta_list[34][MAX_STRING_LENGTH];
  char  band_fname[8][MAX_STRING_LENGTH];
  char  band_gain[8][MAX_STRING_LENGTH];
  char  *meta_lut[34][2] = {{"HEADER_FILE",       ""},
			    {"FILE_TYPE",         "FILE_TYPE"},
			    {"DATA_PROVIDER",     "DATA_PROVIDER"},
			    {"SATTELITE",         "SATELLITE"},
			    {"INSTRUMENT",        "INSTRUMENT"},
			    {"ACQUISITION_DATE",  "ACQUISITION_DATE"},
			    {"PRODUCTION_DATE",   "PROCESSING_DATE"},
			    {"SOLAR_ZENITH",      "SOLAR_ZENITH"},
			    {"SOLAR_AziMUTH",     "SOLAR_AZIMUTH"},
			    {"WRS_SYSTEM",        "WRS_SYSTEM"},
			    {"WRS_PATH",          "WRS_PATH"},
			    {"WRS_ROW",           "WRS_ROW"},
			    {"NBAND",             "NBAND"},
			    {"BANDS",             "BANDS"},
			    {"GAIN",              "GAIN"},
			    {"BIAS",              "BIAS"},
			    {"NSAMPLE",           "NSAMPLE"},
			    {"NLINE",             "NLINE"},
			    {"FILE_NAMES",        "FILE_NAMES"},
			    {"NBAND_TH",          "NBAND_TH"},
			    {"BANDS_TH",          "BANDS_TH"},
			    {"GAIN_TH",           "GAIN_TH"},
			    {"BIAS_TH",           "BIAS_TH"},
			    {"NSAMPLE_TH",        "NSAMPLE_TH"},
			    {"NLINE_TH",          "NLINE_TH"},
			    {"FILE_NAMES_TH",     "FILE_NAMES_TH"},
			    {"PROJECTION_NUMBER", "PROJECTION_NUMBER"},
			    {"PIXEL_SIZE",        "PIXEL_SIZE"},
			    {"UPPER_LEFT_CORNER", "UPPER_LEFT_CORNER"},
			    {"PROJECTION_ZONE",   "PROJECTION_ZONE"},
			    {"PROJECTION_SPHERE", "PROJECTION_SPHERE"},
			    {"GAIN_SETTINGS",     "GAIN_SETTINGS"},
			    {"GAIN_SETTINGS_TH",  "GAIN_SETTINGS_TH"},
			    {"END",               ""}}; 
  FILE *in, *out;
  
  if((in=fopen(input, "r"))==NULL)
    DIE( "Can't open input metadata file %s", input);


  sprintf(met_name, "%s.carbon_met.txt", scene_name);
  if((out=fopen(met_name, "w"))==NULL)
    DIE( "Can't open output header file for write %s", met_name);

  
  for(i=0; i<34; i++) chkFlg[i] = 0;

  /* process line by line */
  while(fgets(tmpstr, MAX_STRING_LENGTH, in) != NULL) {

    strcpy(buffer, tmpstr);
    /* get string token */
    tokenptr = strtok(buffer, seperator);
    label=tokenptr;
    while(tokenptr != NULL) {
 
      tokenptr = strtok(NULL, seperator);

      /* information provided just after label */
      for(i=0; i<34; i++)
	if(strcmp(label, meta_lut[i][1]) == 0) {
	  chkFlg[i] = 1;
	  switch(i) {      
	    /* some metadata need special attentions */
	  case 3:     // satellite, attr name different 
  	    sscanf(tokenptr, "%s", str);
	    sprintf(meta_list[i], "%s = %s\n", meta_lut[i][0], str);
	    break;
	  case 4:     // handle sensor, used to determine BIAS/GAIN or GAIN_SETTINGS
	    sscanf(tokenptr, "%s", sensor);
	    if(strstr(tokenptr,"ETM") != NULL)
	      sprintf(meta_list[i],"%s = ETM\n", meta_lut[i][0]);
	    else
	      sprintf(meta_list[i],"%s = TM\n", meta_lut[i][0]);
	    break;
	  case 5:    // ACQUISITION DATE, need to compare with PRODUCTION_DATE
	    sscanf(tokenptr, "%10s", acquisition_date);
	    sprintf(meta_list[i], "%s = %s\n", meta_lut[i][0], acquisition_date);
	    break;
	  case 6:    // PRODUCTION_DATE
	    // some input give wrong production date and need extra processing   
	    sscanf(tokenptr, "%10s", production_date);
	    sprintf(meta_list[i], "%s = %s\n", meta_lut[i][0], production_date);
	    break;	    
	  case 7:    // SOLAR_ZENITH, need to check if it is solar elevation or zenith 
	    sscanf(tokenptr, "%f", &fnum);
	    sprintf(meta_list[i], "%s = %5.2f\n", meta_lut[i][0], fnum);
	    break;
	  case 14:   // GAIN for TM
	    if(strcmp(sensor, "TM")==0) {
	      for(j=1; j<=6; j++) {
		sscanf(tokenptr, "%f", &gain[j]);
		tokenptr = strtok(NULL, seperator);
	      }
	      gain[7] = gain[6];
	     }
	    break;
	  case 15:   // BIAS for TM
	    if(strcmp(sensor, "TM")==0) {
	      for(j=1; j<=6; j++) {
		sscanf(tokenptr, "%f", &bias[j]);
		tokenptr = strtok(NULL, seperator);
	      }
	      bias[7] = bias[6];
	    }
	    break;
	  case 18:  // FILE_NAMES 
	    for(j=1; j<=6; j++) {
	      sscanf(tokenptr, "%s", str);
	      /* remove .gz extension from filename if there is */
	      if(strstr(str, ".gz"))
		strncpy(band_fname[j], str, strlen(str)-3);
	      else
		strcpy(band_fname[j], str);
	      tokenptr = strtok(NULL, seperator);
	    }
	    strcpy(band_fname[7], band_fname[6]);
	    sprintf(meta_list[i],"%s = %s, %s, %s, %s, %s, %s\n", meta_lut[i][0],band_fname[1],
		    band_fname[2],band_fname[3],band_fname[4],band_fname[5],band_fname[7]); 
	    break;
	  case 21:  // GAIN_TH for TM
	    if(strcmp(sensor, "TM")==0) 
	      sscanf(tokenptr, "%f", &gain[6]);
	    break;
	  case 22:  // BIAS_TH for TM
	    if(strcmp(sensor, "TM")==0) 
	      sscanf(tokenptr, "%f", &bias[6]);
	    break;
	  case 25:  // FILE_NAMES_TH
	    sscanf(tokenptr, "%s", str);
	    /* remove .gz extension from filename if there is */
	    if(strstr(str, ".gz"))
	      strncpy(band_fname[6], str, strlen(str)-3);
	    else
	      strcpy(band_fname[6], str);	
	    sprintf(meta_list[i], "%s = %s\n", meta_lut[i][0], band_fname[6]);	
	    break;
	  case 29:  // ZONE number, different for south and north
	    sscanf(tokenptr, "%d%s", &zone, hemi);
	    // use negative zone number for south and positive for north
	    if(strcmp(hemi, "S")==0 || strcmp(hemi, "s")==0)
	      zone = -zone;
	    sprintf(meta_list[i], "%s = %d\n", meta_lut[i][0], zone);
	    break;	    
	  case 31:  // GAIN_SETTINGS for ETM+
	    if(strcmp(sensor, "TM") != 0) {
	      for(j=1; j<=6; j++) {
		sscanf(tokenptr, "%s", band_gain[j]);
		tokenptr = strtok(NULL, seperator);
	      }
	      strcpy(band_gain[7], band_gain[6]);
	    }
	    break;
	  case 32:  // GAIN_SETTINGS for ETM+
	    if(strcmp(sensor, "TM") != 0)
	      sscanf(tokenptr, "%s", band_gain[6]);
	    break;	    
	  default:  // other parameters can be used directly
	    strcpy(meta_list[i], tmpstr);
	  }
	}
      
      /* in case label (key words) is not the first word in a line */
      if(tokenptr != NULL) 
	label = tokenptr;
    }
  }
  
  /* some parameters couldn't find from metadata file and are hard-coded for Landsat 7 only */

  /* write gain and bias for TM */
  if(strcmp(sensor, "TM") == 0) {
    /* write gain and bias for spectral bands */
    sprintf(meta_list[14],"%s = %f, %f, %f, %f, %f, %f\n", meta_lut[14][0], gain[1], gain[2],
	    gain[3], gain[4], gain[5], gain[7]); 
    sprintf(meta_list[15],"%s = %f, %f, %f, %f, %f, %f\n", meta_lut[15][0], bias[1], bias[2],
	    bias[3], bias[4], bias[5], bias[7]);
    /* write gain and bias for thermal band */ 
    sprintf(meta_list[21],"%s = %f\n", meta_lut[21][0], gain[6]);
    sprintf(meta_list[22],"%s = %f\n", meta_lut[22][0], bias[6]);
  }  
  else {
    sprintf(meta_list[31],"%s = %s, %s, %s, %s, %s, %s\n", meta_lut[31][0], band_gain[1], band_gain[2],
	    band_gain[3], band_gain[4], band_gain[5], band_gain[7]); 
    sprintf(meta_list[32],"%s = %s\n", meta_lut[32][0], band_gain[6]);
  }

  /* check file name and see if it is a TIFF file */
  if(strstr(band_fname[1], ".tif") || strstr(band_fname[1], ".TIF")) {
    sprintf(meta_list[1], "%s = GEOTIFF\n", meta_lut[1][0]);   // DATA_FORMAT
    geotiff = 1;
  }
  else {
    sprintf(meta_list[1], "%s = BINARY\n", meta_lut[1][0]);   // DATA_FORMAT
    geotiff = 0;
  }

  /* use existing coordinate for the UMD metadata - 4/7/09 */
  /* update metadata from GeoTIFF file */
  if(geotiff == 1)  {
    if(readULfromGeoTiff(band_fname[1], &nrows, &ncols, ulxy, &res) == SUCCESS) {
      sprintf(meta_list[16],"%s = %d\n", meta_lut[16][0], ncols);
      sprintf(meta_list[17],"%s = %d\n", meta_lut[17][0], nrows);
      sprintf(meta_list[27],"%s = %f\n", meta_lut[27][0], res);
      sprintf(meta_list[28],"%s = %f, %f\n", meta_lut[28][0], ulxy[0], ulxy[1]); 
    }
  }

  /* check acquisition and production date */
  /* get year from acquisition date */
  sscanf(acquisition_date, "%d", &year);
  /* get year from production date */
  sscanf(production_date, "%d", &pyear);
  /* production date must be later than acquisition date */
  if (pyear<year) {
    strcpy(production_date, acquisition_date); 
    sprintf(meta_list[6],"%s = %s\n", meta_lut[6][0], production_date);
    chkFlg[6] = 2;
  }

  /* start and end */
  sprintf(meta_list[0], "%s\n", meta_lut[0][0]);
  sprintf(meta_list[33], "%s\n", meta_lut[33][0]);
  chkFlg[0] = 1;  chkFlg[33] = 1;

  /* write all paramters to input file */
  for(i=0; i<34; i++) {

    if(strcmp(sensor, "TM")==0) {
      if(i==31 || i==32) continue;
    }
    else  // for ETM+
      if(i==14 || i==15 || i==21 || i==22) continue;

    if(chkFlg[i]==1) {
      fprintf(out, "%s", meta_list[i]);
      LOG("FOUND: %s", meta_list[i]);
    }
    else 
      if(chkFlg[i]==2) {
	fprintf(out, "%s", meta_list[i]);
	LOG("VALUE REPLACED: %s", meta_list[i]);
	WARN("Metadata %s was not found or wrong. \tIt was replaced with %s", meta_lut[i][0], meta_list[i]);
      }
      else {
	fprintf(out, "%s = \n", meta_lut[i][0]);
	WARN("Metadata %s was not found. \tPlease enter it in %s manually", meta_lut[i][0], met_name);
      }
  }

  fclose(in);
  fclose(out);

  if(writeENVIHeader(scene_name, nrows, ncols, ulxy, res, zone)==FAILURE)
    return FAILURE;
  
  if(strstr(sensor, "ETM")) return ETM;
  else return TM;
}



/** read metadata such as nrows, ncols, and upper left corner information from GeoTiff file 
 *  GeoCover header file (*MTL) only contains coordinates from original data source
 */
int readULfromGeoTiff(char *fname, int *nrows, int *ncols, double ulxy[], double *res)
{
  double *tiePoint, *pixelScale;
  uint16 count, coor_sys;

  TIFF *fp_tiff;
  GTIF *gtif;

  if((fp_tiff = XTIFFOpen(fname, "r"))==NULL) {
    DIE("Can't open base GEOTIFF file %s\n", fname);
  } 

  /* get Landsat metadata from tiff file */
  if(TIFFGetField(fp_tiff, TIFFTAG_IMAGEWIDTH, ncols)==0) {
    WARN("Retrieve BASE Landsat file %s error\n", fname);
    return FAILURE;
  }

  if(TIFFGetField(fp_tiff, TIFFTAG_IMAGELENGTH, nrows)==0) {
    WARN("Retrieve BASE Landsat file %s error\n", fname);
    return FAILURE;
  }

  count=6;
  if(TIFFGetField(fp_tiff, TIFFTAG_GEOTIEPOINTS, &count, &tiePoint)==0) {
    WARN("Retrieve BASE Landsat file %s error\n", fname);
    return FAILURE;
  }

  count=3;
  if(TIFFGetField(fp_tiff, TIFFTAG_GEOPIXELSCALE, &count, &pixelScale)==0) {
    WARN("Retrieve BASE Landsat file %s error\n", fname);
    return FAILURE;
  }
  *res = pixelScale[0];

  /* GeoKey 1025 (GTRasterTypeGeoKey) dictates whether the reference
     coordinate is the UL (*RasterPixelIsArea*, code 1) or center
     (*RasterPixelIsPoint*, code 2) of the UL pixel. If this key is missing,
     the default (as defined by the specification) is to be
     *RasterPixelIsArea*, which is the UL of the UL pixel. */
  gtif = GTIFNew(fp_tiff);
  if (GTIFKeyGet(gtif, GTRasterTypeGeoKey, &coor_sys, 0, 1) != 1) {
    WARN("Coordinate system is not defined in %s\n", fname);
    WARN("assume used UL of the UL pixel\n");
  }
  if (coor_sys == RasterPixelIsPoint){
    ulxy[0] = tiePoint[3] - 0.5 * (*res);
    ulxy[1] = tiePoint[4] + 0.5 * (*res);
  }
  else {  /* default use RasterPixelIsArea */
    ulxy[0] = tiePoint[3];
    ulxy[1] = tiePoint[4];
  }
  GTIFFree(gtif);

  return SUCCESS;
}



/**
 * !Description:
 *   converting year-mm-dd to julian day of year. 
 *
 * !References and Credits:
 *   revised from program "jdoy.c" by MODIS LDOPE QA team
 *
 *  Note: *mm and *dd are input/output varialbes
 *  Input julian day number:
 *       input *mm = 0;
 *       input *dd = julian day
 *       output in mm-dd-yyyy format
 *  Input in mm-dd-yyyy format
 *       output *mm = 0;
 *       output *dd = julian day number;
 */
int conv_date(int *mm, int *dd, int yyyy)
{
  int nm, im;
  int st = -1;
  int ndays[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

  if ((yyyy%400 == 0) || ((yyyy%4 == 0) && (yyyy%100 != 0))) ndays[1] = 29; 
  if (*mm == 0)
  {
    if (*dd <= 0) {
      DIE( "Error in input julian date: %d %d", *dd, yyyy);
    }
    else
    {
      for (im=0; ((im<12)&&(*dd>0)); im++)
        *dd -= ndays[im];
      if ((im > 12) || ((im == 12) && (*dd > 0))) {
        DIE( "Error in input julian date: %d %d", *dd, yyyy);
      }
      else
      {
        *mm = im;
        *dd += ndays[*mm - 1];
        st = 1;
      }
    }
  }
  else
  {
    if ((*mm <= 0) || (*dd <= 0)) {
      DIE( "Error in input date: %d %d %d", *mm, *dd, yyyy);
    }
    else
    {
      nm = *mm - 1;
      for (im=0; im<nm; im++)
        *dd += ndays[im];
      *mm = 0;
      st = 1;
    }
  }
  return st;
}


/* write geo-reference information to ENVI header file */
int writeENVIHeader(char scene_name[], int nrows, int ncols, double ulxy[], double res, int zone)
{
  int  i;
  char tmpstr[MAX_STRING_LENGTH] = "\0";
  char *ledaps_file[4]={"lndcal", "lndcsm", "lndth", "lndsr"};
  FILE *out;

  for(i=0; i<4; i++) {
    sprintf(tmpstr, "%s.%s.hdf.hdr", ledaps_file[i], scene_name);
    if((out=fopen(tmpstr, "w"))==NULL) {
      DIE( "Write ENVI header file %s error!", tmpstr);
    }
    fprintf(out,"ENVI\n");
    fprintf(out,"description = {LEDAPS HDF File Imported into ENVI}\n");
    fprintf(out,"samples = %d\n", ncols);
    fprintf(out,"lines   = %d\n", nrows);
    fprintf(out,"bands   = 1\n");
    fprintf(out,"header offset = 0\n");
    fprintf(out,"file type = HDF Scientific Data\n");
    fprintf(out,"data type = 2\n");
    fprintf(out,"interleave = bsq\n");
    fprintf(out,"sensor type = Landsat\n");
    fprintf(out,"byte order = 0\n");
    if(zone > 0) 
      fprintf(out,"map info = {UTM, 1.000, 1.000, %f, %f, %f, %f, %d, North, WGS-84, units=Meters}", 
	      ulxy[0], ulxy[1], res, res, zone);
    else
      fprintf(out,"map info = {UTM, 1.000, 1.000, %f, %f, %f, %f, %d, South, WGS-84, units=Meters}", 
	      ulxy[0], ulxy[1], res, res, -zone);
    fclose(out);
  }
  return SUCCESS;
}


/* open log file for append and add time stamp on it */
int open_log(char name[])
{
  struct tm *currtime;
  time_t t; 
  char str[MAX_STRING_LENGTH];

  if((LOG_FP=fopen(logFile, "a"))==NULL) {
    WARN("cannot open %s", logFile);
    return FAILURE;
  }
  t = time(NULL);
  currtime = (struct tm *)gmtime(&t);
  strftime(str, 100, "%FT%H:%M:%SZ", currtime);
  LOG("\n\n##################################\n");
  LOG("Starting <%s> on %s\n", name, str);
  return SUCCESS;
}

void LOG(char *fmt, ...)
{
  va_list ap;
  
  va_start(ap, fmt);
  vfprintf(LOG_FP, fmt, ap);
  va_end(ap);
}

void WARN(char *fmt, ...)
{
  va_list ap;
  
  fprintf(stderr, "WARNING : ");
  
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
  fprintf(stderr, "\n");
}

void die(char *name, int lineno, char *fmt, ...)
{
  va_list ap;
  
  fprintf(stderr, "FAIL %s:%d : ", name, lineno);
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
  fprintf(stderr, "\n");
  
  exit(FAILURE);
}


/* recursively search a directory tree specified by path, for the
 * first file named name.
 * Returns non-zero if found, and path points to full path,
 * or 0 if not found, and path is unchanged.
 */
 
/* search the specified directory for the specified name */
int _scan_dir(char *path, char *name)
{
  DIR *fd;
  struct dirent *dirent_p;
  char *nbp;
  int found = 0;
  
  nbp = path + strlen(path);
  if (*nbp != '/')
    *nbp++ = '/';       /* add to end of directory name if not already there */
    
  if (nbp+MAXNAMLEN+2 >= path + DIR_BUF_SIZE)
  {
    fprintf(stderr, "path name too long -- cannot search '%s'\n", path);
    return found;       /* path string is too long -- prevent overflow */
  }
  if ((fd = opendir(path)) == NULL)
  {
    fprintf(stderr, "could not read directory '%s'\n", path);
    return found;
  }
  while ((dirent_p = readdir(fd)) != NULL)   /* search directory */
  {
    if (dirent_p->d_ino == 0)                /* slot not in use */
      continue;
    if (strcmp(dirent_p->d_name, ".") == 0   /* ignore current ... */
    || strcmp(dirent_p->d_name, "..") == 0)  /* and parent directory */
      continue;
      
    strcat(path, dirent_p->d_name);          /* check this path */
    if ((found = find_file(path, name)) != 0)
      break;                                /* found it */
    else
      *nbp = '\0';                          /* restore directory name */
  }
  closedir(fd);
  return found;
}

int find_file(char *path, char *name)
{
  struct stat stbuf;
  char pbuf[DIR_BUF_SIZE] = {0};
  int found = 0;
  
  strcpy(pbuf, path);             /* this is the path we are checking */
  
  if (stat(pbuf, &stbuf) != 0)    /* make sure the path exists */
  {
    fprintf(stderr, "find_file: can't stat %s\n", pbuf);
    return found;
  }
  if ((stbuf.st_mode & S_IFMT) == S_IFDIR)   /* directory, so search */
    found = _scan_dir(pbuf, name);
  else                                       /* file, so check */
    found = (strcmp(pbuf + strlen(pbuf) - strlen(name), name) == 0);

  if (found) strcpy(path, pbuf);  /* found file, so remember its location */

  return found;
}
