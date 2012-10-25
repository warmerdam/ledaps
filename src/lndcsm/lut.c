
#include <stdlib.h>
#include "lndcsm.h"
#include "lut.h"
#include "input.h"
#include "mystring.h"
#include "error.h"

#define CLOUD_LAND   (0)
#define CLOUD_CLOUD  (1)    /* bit 0 */
#define CLOUD_SHADOW (2)    /* bit 1 */
#define CLOUD_WATER  (4)    /* bit 2 */
#define CLOUD_SNOW   (8)    /* bit 3 */
#define CLOUD_FILL   (255)

#define old_CLOUD_FILL   (0)
#define old_CLOUD_LAND   (1)
#define old_CLOUD_WATER  (2)
#define old_CLOUD_SHADOW (3)
#define old_CLOUD_SNOW   (4)
#define old_CLOUD_CLOUD  (5)

#define REFL_SCALE_FACTOR (10000.0)
#define THERM_SCALE_FACTOR (100.0)
#define CELC_2_KELVIN (273)
#define TOP_DISTANCE (0.01)
#define	KEEP_RATIO (0.001)
	
#define LONG_NAME        ("cloud snow mask")
#define UNITS            ("flag")
#define SCALE_FACTOR     (0.000)
#define ADD_OFFSET       (0.0)
#define SCALE_FACTOR_ERR (0.0)
#define ADD_OFFSET_ERR   (0.0)
#define CALIBRATED_NT    (DFNT_FLOAT32)
#define MIN_VALID_CSM (0)
#define MAX_VALID_CSM (5)

Lut_t *GetLut(char *file_name, int nband, Input_meta_t *meta) {
  Lut_t *this;
  int ib, iband;

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

  this->nband=      nband;
  this->in_fill=    meta->fill;
  this->csm_fill=   CLOUD_FILL;
  this->csm_land=   CLOUD_LAND;
  this->csm_cloud=  CLOUD_CLOUD;
  this->csm_snow=   CLOUD_SNOW;
  this->csm_shadow= CLOUD_SHADOW;
  this->csm_water=  CLOUD_WATER;
  this->refl_scale_factor=  REFL_SCALE_FACTOR;
  this->therm_scale_factor= THERM_SCALE_FACTOR;
  this->celc2kelvin= CELC_2_KELVIN;
  this->top_distance = TOP_DISTANCE;
  this->keep_ratio = KEEP_RATIO;

  this->long_name = DupString(LONG_NAME);
  if (this->long_name == (char *)NULL) {
    free(this);
    RETURN_ERROR("duplicating long name ref", "GetLut", (Lut_t *)NULL);
  }

  this->units = DupString(UNITS);
  if (this->units == (char *)NULL) {
    free(this);
    RETURN_ERROR("duplicating ref units", "GetLut", (Lut_t *)NULL);
  }

  this->valid_range[0]=     MIN_VALID_CSM;
  this->valid_range[1]=     MAX_VALID_CSM;
  this->scale_factor=       SCALE_FACTOR;
  this->scale_factor_err=   SCALE_FACTOR_ERR;
  this->add_offset=         ADD_OFFSET;
  this->add_offset_err=     ADD_OFFSET_ERR;
  this->calibrated_nt=      CALIBRATED_NT;
  this->output_fill=        CLOUD_FILL;

  if (!InputMetaCopy(meta, nband, &this->meta)) {
    free(this->file_name);
    free(this);
    RETURN_ERROR("copying input metadata", "GetLut", (Lut_t *)NULL);
  }    

  /* Compute the band numbers for TM and ETM */
  
  if (meta->inst == INST_ETM  ||  meta->inst == INST_TM) {

    for (ib = 0; ib < nband; ib++) {
      iband = meta->iband[ib] - 1;
    }

  } else {  /* INST_MSS */

    for (ib = 0; ib < nband; ib++) {
      iband = meta->iband[ib] - 1;
    }

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
