#ifndef LNDCSM_H
#define LNDCSM_H

#include "mystring.h"

enum {BRIGHT = 0, GREEN, WET, HIGH, LOW, NBAND_KT};
enum {CLMASK = 0, NBAND_CSM};

#define NBAND_REFL_MAX (6)

/* Data provider type definition */

typedef enum {
  PROVIDER_NULL = -1,
  PROVIDER_UMD = 0, 
  PROVIDER_EDC,
  PROVIDER_CCRS,
  PROVIDER_MRLC01,
  PROVIDER_MRLC92,
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

/* Integer image coordinates data structure */

#ifndef IMG_COORD_INT_TYPE_DEFINED
#define IMG_COORD_INT_TYPE_DEFINED

typedef struct {
  int l;                /* line number */
  int s;                /* sample number */
} Img_coord_int_t;

#endif

#endif
