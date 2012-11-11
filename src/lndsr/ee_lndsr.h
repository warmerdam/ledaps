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

#include "read_grib_tools.h"
#include "sixs_runs.h"

Ar_gridcell_t *ee_build_ar_gridcell(Lut_t *lut,
                                    Space_t *space,
                                    Input_meta_t *meta,
                                    double scene_gmt,
                                    const float *anc_SP,
                                    const float *anc_WV,
                                    const float *anc_O3,
                                    const float *anc_dem);

int ee_lndsr_main(
    const char *metadata,
    const char *sixs_filename,

    int full_input_size_s, int full_input_size_l,
    int input_offset_s, int input_offset_l,
    int input_size_s, int input_size_l, int nband,
    short *band1, short *band2, short *band3, 
    short *band4, short *band5, short *band7,

    int utm_zone, 
    double ul_x, double ul_y, double pixel_size,

    int ar_size_s, int ar_size_l,
    const short *aerosol,
    const float *anc_SP,
    const float *anc_WV,
    const float *anc_ATEMP,
    const float *anc_O3,
    const float *anc_dem);

bool ee_GetInputMeta(Input_meta_t *meta, const char *metadata);

int ee_PrepareSpaceDef( Space_def_t *space_def, int utm_zone,
                        double ul_x, double ul_y, double pixel_size,
                        Img_coord_int_t *input_size );

#endif /* ndef LNDSR_EE_H */

