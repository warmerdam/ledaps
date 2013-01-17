#include "ee_lndsr.h"
#include "keyvalue.h"

static char **ee_ParseMetadataString(const char *metadata)
{
#define MAX_KEYVALS 100
    char **keyvals;
    int  keyval_count = 0;
    const char *start, *end;

    keyvals = (char **) calloc(sizeof(char*),MAX_KEYVALS);

    /* parse the metadata string into a list of keys and values on newlines */

    start = metadata;
    while (*start != '\0') {
        int item_len;

        end = start;
        while (*end != '\0' && *end != '\n') {
            end++;
        }

        item_len = end - start;
        
        keyvals[keyval_count] = (char *) calloc(sizeof(char),item_len+1);
        strncpy(keyvals[keyval_count], start, item_len);
        keyval_count++;

        if( keyval_count == MAX_KEYVALS-2 ) {
            fprintf( stderr, "Too much metadata!  Ignoring extra items.\n");
            break;
        }

        if (*end) 
            start = end+1;
        else
            start = end;
    }

    return keyvals;
}

static char *ee_FindMetadataItem(char **keyvals, const char *key) {
    int i;
    for( i = 0; keyvals[i] != NULL; i += 2 ) {
        if( strcmp(keyvals[i],key) == 0) {
            return keyvals[i+1];
        }
    }
    return NULL;
}

/* these constants are copied from input.c */

#define SDS_PREFIX ("band")
#define SDS_CSM ("cloud_snow_mask")
#define INPUT_FILL (-9999)

#define INPUT_PROVIDER ("DataProvider")
#define INPUT_SAT ("Satellite")
#define INPUT_INST ("Instrument")
#define INPUT_ACQ_DATE ("AcquisitionDate")
#define INPUT_PROD_DATE ("Level1ProductionDate")
#define INPUT_SUN_ZEN ("SolarZenith")
#define INPUT_SUN_AZ ("SolarAzimuth")
#define INPUT_WRS_SYS ("WRS_System")
#define INPUT_WRS_PATH ("WRS_Path")
#define INPUT_WRS_ROW ("WRS_Row")
#define INPUT_NBAND ("NumberOfBands")
#define INPUT_BANDS ("BandNumbers")

#define N_LSAT_WRS1_ROWS  (251)
#define N_LSAT_WRS1_PATHS (233)
#define N_LSAT_WRS2_ROWS  (248)
#define N_LSAT_WRS2_PATHS (233)

bool ee_GetInputMeta(Input_meta_t *meta, const char *metadata) 
{
  double dval[NBAND_REFL_MAX];
  char date[MAX_DATE_LEN + 1];
  int ib, nband, i;
  char **keyvals;
  char *string;
  char *error_string = (char *)NULL;

  keyvals = ee_ParseMetadataString(metadata);
  
  /* Check the parameters */

  meta->fill = INPUT_FILL;

  /* Read the metadata */

  string = ee_FindMetadataItem(keyvals, INPUT_PROVIDER);
  if (string == NULL)
    RETURN_ERROR("reading attribute (data provider)", "ee_GetInputMeta", false);
  meta->provider = (Provider_t)KeyString(string, strlen(string), Provider_string, 
		                         (int)PROVIDER_NULL, (int)PROVIDER_MAX);
  if (meta->provider == PROVIDER_NULL)
    RETURN_ERROR("invalid input metadata (data provider)", 
                 "ee_GetInputMeta", false);

  string = ee_FindMetadataItem(keyvals, INPUT_SAT);
  if (string == NULL)
    RETURN_ERROR("reading attribute (instrument)", "ee_GetInputMeta", false);
  meta->sat = (Sat_t)KeyString(string, strlen(string), Sat_string, 
		               (int)SAT_NULL, (int)SAT_MAX);
  if (meta->sat == SAT_NULL)
    RETURN_ERROR("invalid input metadata (satellite)", "ee_GetInputMeta", false);

  string = ee_FindMetadataItem(keyvals, INPUT_INST);
  if (string == NULL)
    RETURN_ERROR("reading attribute (instrument)", "ee_GetInputMeta", false);
  meta->inst = (Inst_t)KeyString(string, strlen(string), Inst_string, 
		                 (int)INST_NULL, (int)INST_MAX);
  if (meta->inst == INST_NULL)
    RETURN_ERROR("invalid input metadata (instrument)", "ee_GetInputMeta", false);

  string = ee_FindMetadataItem(keyvals, INPUT_ACQ_DATE);
  if (string == NULL)
    RETURN_ERROR("reading attribute (acquisition date)", "ee_GetInputMeta", false);
  if (!DateInit(&meta->acq_date, string, DATE_FORMAT_DATEA_TIME))
    RETURN_ERROR("converting acquisition date", "ee_GetInputMeta", false);

  string = ee_FindMetadataItem(keyvals, INPUT_PROD_DATE);
  if (string == NULL)
    RETURN_ERROR("reading attribute (production date)", "ee_GetInputMeta", false);
  if (!DateInit(&meta->prod_date, string, DATE_FORMAT_DATEA_TIME))
    RETURN_ERROR("converting production date", "ee_GetInputMeta", false);

  string = ee_FindMetadataItem(keyvals, INPUT_SUN_ZEN);
  if (string == NULL)
    RETURN_ERROR("reading attribute (solar zenith)", "ee_GetInputMeta", false);
  dval[0] = atof(string);
  if (dval[0] < -90.0  ||  dval[0] > 90.0)
    RETURN_ERROR("solar zenith angle out of range", "ee_GetInputMeta", false);
  meta->sun_zen = (float)(dval[0] * RAD);

  string = ee_FindMetadataItem(keyvals, INPUT_SUN_AZ);
  if (string == NULL)
    RETURN_ERROR("reading attribute (solar azimuth)", "ee_GetInputMeta", false);
  dval[0] = atof(string);
  if (dval[0] < -360.0  ||  dval[0] > 360.0)
    RETURN_ERROR("solar azimuth angle out of range", "ee_GetInputMeta", false);
  meta->sun_az = (float)(dval[0] * RAD);
  printf( "Sun_az = %.12g (orig = %.12g)\n", meta->sun_az, dval[0]);

  string = ee_FindMetadataItem(keyvals, INPUT_WRS_SYS);
  if (string == NULL)
    RETURN_ERROR("reading attribute (WRS system)", "ee_GetInputMeta", false);
  meta->wrs_sys = (Wrs_t)KeyString(string, 1, Wrs_string, 
		                       (int)WRS_NULL, (int)WRS_MAX);
  if (meta->wrs_sys == WRS_NULL)
    RETURN_ERROR("invalid input metadata (WRS sytem)", "ee_GetInputMeta", false);

  string = ee_FindMetadataItem(keyvals, INPUT_WRS_PATH);
  if (string == NULL)
    RETURN_ERROR("reading attribute (WRS path)", "ee_GetInputMeta", false);
  meta->ipath = atoi(string);
  if (meta->ipath < 1) 
    RETURN_ERROR("WRS path out of range", "ee_GetInputMeta", false);

  string = ee_FindMetadataItem(keyvals, INPUT_WRS_ROW);
  if (string == NULL)
    RETURN_ERROR("reading attribute (WRS row)", "ee_GetInputMeta", false);
  meta->irow = atoi(string);
  if (meta->irow < 1) 
    RETURN_ERROR("WRS path out of range", "ee_GetInputMeta", false);

  string = ee_FindMetadataItem(keyvals, INPUT_NBAND);
  if (string == NULL)
    RETURN_ERROR("reading attribute (number of bands)", "ee_GetInputMeta", false);
  nband = atoi(string);
  if (nband < 1  ||  nband > NBAND_REFL_MAX) 
    RETURN_ERROR("number of bands out of range", "ee_GetInputMeta", false);

  string = ee_FindMetadataItem(keyvals, INPUT_BANDS);
  if (string == NULL)
    RETURN_ERROR("reading attribute (band numbers)", "ee_GetInputMeta", false);

  if (strcmp(string,"1, 2, 3, 4, 5, 7") != 0)
    RETURN_ERROR("unexpected result (band numbers)", 
                 "ee_GetInputMeta", false);

  meta->iband[0] = 1;
  meta->iband[1] = 2;
  meta->iband[2] = 3;
  meta->iband[3] = 4;
  meta->iband[4] = 5;
  meta->iband[5] = 7;

  /* Free parsed results */
  for( i = 0; keyvals[i] != NULL; i++ ) {
    free(keyvals[i]);
  }
  free(keyvals);

  /* Check WRS path/rows */

  error_string = (char *)NULL;

  if (meta->wrs_sys == WRS_1) {
    if (meta->ipath > N_LSAT_WRS1_PATHS)
      error_string = "WRS path number out of range";
    else if (meta->irow > N_LSAT_WRS1_ROWS)
      error_string = "WRS row number out of range";
  } else if (meta->wrs_sys == WRS_2) {
    if (meta->ipath > N_LSAT_WRS2_PATHS)
      error_string = "WRS path number out of range";
    else if (meta->irow > N_LSAT_WRS2_ROWS)
      error_string = "WRS row number out of range";
  } else
    error_string = "invalid WRS system";

  if (error_string != (char *)NULL)
    RETURN_ERROR(error_string, "ee_GetInputMeta", false);

  /* Check satellite/instrument combination */
  
  if (meta->inst == INST_MSS) {
    if (meta->sat != SAT_LANDSAT_1  &&  
        meta->sat != SAT_LANDSAT_2  &&
        meta->sat != SAT_LANDSAT_3  &&  
	  meta->sat != SAT_LANDSAT_4  &&
        meta->sat != SAT_LANDSAT_5)
      error_string = "invalid insturment/satellite combination";
  } else if (meta->inst == INST_TM) {
    if (meta->sat != SAT_LANDSAT_4  &&  
        meta->sat != SAT_LANDSAT_5)
      error_string = "invalid insturment/satellite combination";
  } else if (meta->inst == INST_ETM) {
    if (meta->sat != SAT_LANDSAT_7)
      error_string = "invalid insturment/satellite combination";
  } else
    error_string = "invalid instrument type";

  if (error_string != (char *)NULL)
    RETURN_ERROR(error_string, "ee_GetInputMeta", false);

  /* Check band numbers */

  if (meta->inst == INST_MSS) {
    for (ib = 0; ib < nband; ib++)
      if (meta->iband[ib] > 4) {
        error_string = "band number out of range";
	break;
      }
  } else if (meta->inst == INST_TM  ||
             meta->inst == INST_ETM) {
    for (ib = 0; ib < nband; ib++)
      if ( meta->iband[ib] > 7) {
        error_string = "band number out of range";
	break;
      }
  } else
    error_string = "invalid instrument type";

  if (error_string != (char *)NULL)
    RETURN_ERROR(error_string, "ee_GetInputMeta", false);

  return true;
}
