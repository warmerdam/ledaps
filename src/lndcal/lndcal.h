#ifndef LNDCAL_H
#define LNDCAL_H

#define bounded(A,B,C) (A>B?(A<C?A:C):B)
#define min(A,B) (A>B?B:A)
#define max(A,B) (A>B?A:B)

#include "mystring.h"
#include "space.h"

#define NBAND_REFL_MAX (6)
#define NBAND_QA       (1)
#define NBAND_CAL_MAX (NBAND_REFL_MAX + NBAND_QA)
#define QA_BAND_NUM (6)


typedef enum {
  ALG_NASA,
  ALG_NASA_CPF,
  ALG_CCRS 
} Algorithm_t;

/* Data provider type definition */

typedef enum {
  PROVIDER_NULL = -1,
  PROVIDER_UMD = 0, 
  PROVIDER_EDC, 
  PROVIDER_CCRS,
  PROVIDER_MRLC01,
  PROVIDER_MRLC92,
  PROVIDER_USDA,
  PROVIDER_EROS,
  PROVIDER_CI,
  PROVIDER_Scott,
  PROVIDER_MAX
} Provider_t;

extern const Key_string_t Provider_string[PROVIDER_MAX];

/* Satellite type definition */

typedef enum {
  SAT_NULL = -1,
  SAT_LANDSAT_1 = 0, 
  SAT_LANDSAT_2, 
  SAT_LANDSAT_3, 
  SAT_LANDSAT_4, 
  SAT_LANDSAT_5, 
  SAT_LANDSAT_7, 
  SAT_MAX
} Sat_t;

extern const Key_string_t Sat_string[SAT_MAX];

/* Instrument type definition */

typedef enum {
  INST_NULL = -1,
  INST_MSS = 0, 
  INST_TM,
  INST_ETM, 
  INST_MAX
} Inst_t;

extern const Key_string_t Inst_string[INST_MAX];

/* World Reference System (WRS) type definition */

typedef enum {
  WRS_NULL = -1,
  WRS_1 = 0, 
  WRS_2,
  WRS_MAX
} Wrs_t;

extern const Key_string_t Wrs_string[WRS_MAX];

/* Band gain settings (ETM+ only) */

typedef enum {
  GAIN_NULL = -1,
  GAIN_HIGH = 0, 
  GAIN_LOW, 
  GAIN_MAX
} Gain_t;

extern const Key_string_t Gain_string[GAIN_MAX];

#ifndef IMG_COORD_INT_TYPE_DEFINED

#define IMG_COORD_INT_TYPE_DEFINED

/* Integer image coordinates data structure */

typedef struct {
  int l;                /* line number */
  int s;                /* sample number */
} Img_coord_int_t;

#endif

typedef struct {
  float DN_to_Radiance_gain[7];
  float DN_to_Radiance_offset[7];
  float final_gain[7];
  float final_offset[7];
  float reflective_forward_gains[6][16];
  float reflective_forward_offsets[6][16];
  float reflective_reverse_gains[6][16];
  float reflective_reverse_offsets[6][16];
  float thermal_forward_gains[4];
  float thermal_forward_offsets[4];
  float thermal_reverse_gains[4];
  float thermal_reverse_offsets[4];
  int reference_detector[7];
  bool nasa_flag;
  Algorithm_t algorithm;
  int completion_year;
  int completion_month;
  int completion_day;
} Param_wo_t;

typedef struct {
  float gains[7];
  bool valid_flag;
} Gains_t;

#endif
