#ifndef LNDSR_EE_H
#define LNDSR_EE_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <time.h>

#include "lndsr.h"
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

#include "ee_lndsr_main.h"
#include "read_grib_tools.h"
#include "sixs_runs.h"

/* We hard code things to 4 layers of daily ancillary data, 6 hours apart. */

#define ANC_TIMERES 6.0 /* 6.0 hours time resolution for most ancillary layers*/
#define ANC_NBLAYERS 4  /* 4 "time layers" for each ancillary layer */


Ar_gridcell_t *ee_build_ar_gridcell(Lut_t *lut,
                                    Space_t *space,
                                    Input_meta_t *meta,
                                    double scene_gmt,
                                    const float *anc_SP,
                                    const float *anc_WV,
                                    const float *anc_O3,
                                    const float *anc_dem);


bool ee_GetInputMeta(Input_meta_t *meta, const char *metadata);

int ee_PrepareSpaceDef( Space_def_t *space_def, int utm_zone,
                        double ul_x, double ul_y, double pixel_size,
                        Img_coord_int_t *input_size );
int ***ee_PrepareLineAr( Ar_gridcell_t *ar_gridcell, const short *aerosol );
int ***ee_PrepareLineIn( TileDef_t *tile_def, Lut_t *lut );

#endif /* ndef LNDSR_EE_H */
