#ifndef LNDSR_H
#define LNDSR_H

#define bounded(A,B,C) (A>B?(A<C?A:C):B)
#define min(A,B) (A>B?B:A)
#define max(A,B) (A>B?A:B)
#define nint(A)(A<0?(int)(A-0.5):(int)(A+0.5))

#include "mystring.h"
#include "space.h"

#define NBAND_SR_EXTRA (5)
#define NBAND_REFL_MAX (6)
#define NBAND_PRWV_MAX (3)
#define NBAND_SR_MAX (NBAND_REFL_MAX + NBAND_SR_EXTRA)

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

/* Ozone source */

typedef enum {
  OZSRC_NULL = -1,
  OZSRC_NIMBUS7 = 0, 
  OZSRC_METEOR3, 
  OZSRC_EARTHPROBE, 
  OZSRC_NIMBUS7_FILL, 
  OZSRC_METEOR3_FILL, 
  OZSRC_EARTHPROBE_FILL, 
  OZSRC_FILL, 
  OZSRC_MAX
} Ozsrc_t;

extern const Key_string_t Ozsrc_string[OZSRC_MAX];

/* Integer image coordinates data structure */
#ifndef IMG_COORD_INT_TYPE_DEFINED
#define IMG_COORD_INT_TYPE_DEFINED
typedef struct {
  int l;                /* line number */
  int s;                /* sample number */
} Img_coord_int_t;
#endif
typedef struct {
  int nbrows,nbcols;
  float *lat,*lon;
  float *sun_zen,*view_zen,*rel_az;
  float *wv,*spres,*ozone,*spres_dem;
  float *line_lat,*line_lon,*line_sun_zen,*line_view_zen,*line_rel_az;
  float *line_wv,*line_spres,*line_ozone,*line_spres_dem;
} Ar_gridcell_t;

typedef struct {
  double min_lon;  /* Geodetic longitude coordinate (degrees) */ 
  double min_lat;  /* Geodetic latitude coordinate (degrees) */ 
  double max_lon;  /* Geodetic longitude coordinate (degrees) */ 
  double max_lat;  /* Geodetic latitude coordinate (degrees) */ 
  bool is_fill;    /* Flag to indicate whether the point is a fill value; */                            
} Geo_bounds_t;

typedef struct {
int *computed;
float *tgOG[7],*tgH2O[7],*td_ra[7],*tu_ra[7],*rho_mol[7],*rho_ra[7],*td_da[7],*tu_da[7],*S_ra[7];
float *td_r[7],*tu_r[7],*S_r[7],*rho_r[7];
} atmos_t;

int allocate_mem_atmos_coeff(int nbpts,atmos_t *atmos_coef);
int free_mem_atmos_coeff(atmos_t *atmos_coef);

#endif
