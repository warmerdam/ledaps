
#include <stdlib.h>
#include "lndsr.h"
#include "lut.h"
#include "input.h"
#include "mystring.h"
#include "error.h"

#define OUTPUT_FILL (-9999)
#define OUTPUT_SATU (20000)
#define MIN_VALID_SR (-2000)
#define MAX_VALID_SR (16000)
#define AEROSOL_FILL (-9999)
#define AEROSOL_REGION_NLINE (40)
#define AEROSOL_REGION_NSAMP (AEROSOL_REGION_NLINE)
#define LONG_NAME_PREFIX ("band %d reflectance")
#define UNITS            ("reflectance")
#define SCALE_FACTOR     (0.0001)
#define ADD_OFFSET       (0.0)
#define SCALE_FACTOR_ERR (0.0)
#define ADD_OFFSET_ERR   (0.0)
#define CALIBRATED_NT    (DFNT_FLOAT32)

/* old definition for ACCA they are incorrect 
#define CLOUD_FILL   (0)
#define CLOUD_LAND   (1)
#define CLOUD_WATER  (2)
#define CLOUD_SHADOW (3)
#define CLOUD_SNOW   (4)
#define CLOUD_CLOUD  (5)
*/
/* correct definition for ACCA */
#define CLOUD_LAND   (0)
#define CLOUD_CLOUD  (1)
#define CLOUD_SHADOW (2)
#define CLOUD_WATER  (4)
#define CLOUD_SNOW   (8)
#define CLOUD_FILL   (255)



Lut_t *GetLut(char *file_name, int nband, Input_meta_t *meta, 
              Img_coord_int_t *input_size) {
  Lut_t *this;

  /* Create the lookup table data structure */

  this = (Lut_t *)malloc(sizeof(Lut_t));
  if (this == (Lut_t *)NULL) 
    RETURN_ERROR("allocating Input data structure", "OpenInput", 
                 (Lut_t *)NULL);

  /* Populate the data structure */

  this->file_name = DupString(file_name);
  if (this->file_name == (char *)NULL) {
    free(this);
    RETURN_ERROR("duplicating file name", "GetLut", (Lut_t *)NULL);
  }

  this->nband = nband;
  this->in_fill = meta->fill;
  this->output_fill = OUTPUT_FILL;
  this->out_satu = OUTPUT_SATU;
  this->aerosol_fill = AEROSOL_FILL;
  this->ar_region_size.l = AEROSOL_REGION_NLINE;
  this->ar_region_size.s = AEROSOL_REGION_NSAMP;
  this->ar_size.l = ((input_size->l - 1) / this->ar_region_size.l) + 1;
  this->ar_size.s = ((input_size->s - 1) / this->ar_region_size.s) + 1;
  this->min_valid_sr = MIN_VALID_SR;
  this->max_valid_sr = MAX_VALID_SR;

  this->cloud_fill=       CLOUD_FILL;       /* cloud_mask_fill_value   */
  this->cloud_land=       CLOUD_LAND;       /* cloud_mask_water_value  */
  this->cloud_water=      CLOUD_WATER;      /* cloud_mask_water_value  */
  this->cloud_shadow=     CLOUD_SHADOW;     /* cloud_mask_shadow_value */
  this->cloud_snow=       CLOUD_SNOW;       /* cloud_mask_snow_value   */
  this->cloud_cloud=      CLOUD_CLOUD;      /* cloud_mask_cloud_value  */
  this->scale_factor=     SCALE_FACTOR;     /* scale factor            */
  this->scale_factor_err= SCALE_FACTOR_ERR; /* scale factor error      */
  this->add_offset=       ADD_OFFSET;       /* add offset              */
  this->add_offset_err=   ADD_OFFSET_ERR;   /* add offset error        */
  this->calibrated_nt=    CALIBRATED_NT;    /* calibrated nt           */

  this->long_name_prefix = DupString(LONG_NAME_PREFIX);
  if (this->long_name_prefix == (char *)NULL) {
    free(this);
    RETURN_ERROR("duplicating long name prefix", "GetLut", (Lut_t *)NULL);
  }

  this->units = DupString(UNITS);
  if (this->units == (char *)NULL) {
    free(this);
    RETURN_ERROR("duplicating ref units", "GetLut", (Lut_t *)NULL);
  }

  if (!InputMetaCopy(meta, nband, &this->meta)) {
    free(this->file_name);
    free(this);
    RETURN_ERROR("copying input metadata", "GetLut", (Lut_t *)NULL);
  }

  return this;
}

bool FreeLut(Lut_t *this) {

  if (this != (Lut_t *)NULL) {
    if (this->file_name != (char *)NULL) free(this->file_name);
    free(this);
  }

  return true;
}
