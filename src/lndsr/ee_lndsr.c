/**
 * lndsr_ee.c - EarthEngine callable "tile" callable version of lndsr.c main.
 *
 * Closely derived from original lndsr.c and related code.
 *
 * Author: Frank Warmerdam <warmerdam@google.com>
 */

#include "ee_lndsr.h"
#include "myproj.h"

atmos_t atmos_coef;

/* 
   dem_array - only needed to compute spres in ar_gridcell. 
   
 */

/**
 * This function is intended to be called from EarthEngine and should
 * only take arguments easily provided by EarthEngine via a simple SWIG
 * wrapped interface.
 * 
 * @param ar_size_s aerosol array width in pixels
 * @param ar_size_l aerosol array height in lines
 * @param aerosol aerosol values (corresponds to line_ar in lndsr)
 * @param anc_SP SP (slope) on aerosol grid x 4 times of day
 * @param anc_WV water vapour on aerosol grid x 4 times of day
 * @param anc_ATEMP air temperature on aerosol grid x 4 times of day
 * @param anc_O3 ozone on aerosol grid 
 * @param anc_dem elevation in meters on aerosol grid, -9999 for nodata.
 * 
 * @return zero on success, non-zero on failure.
 */
int ee_lndsr_main(
    const char *metadata,

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
    const float *anc_dem)

{
    Ar_gridcell_t *ar_gridcell;
    Lut_t *lut;
    Img_coord_int_t input_size;
    Input_meta_t meta;
    Space_def_t space_def;
    Space_t *space = NULL;

    input_size.l = input_size_l;
    input_size.s = input_size_s;

    ee_GetInputMeta( &meta, metadata);

    lut = GetLut( "dummy_param_file", 6, &meta, &input_size );

    ee_PrepareSpaceDef( &space_def, utm_zone, ul_x, ul_y, pixel_size,
                        &input_size );
    space = SetupSpace( &space_def );

    ar_gridcell = ee_build_ar_gridcell( lut, space, &meta, 
                                        anc_SP, anc_WV, anc_O3, anc_dem );

    return 0;
}

int ee_PrepareSpaceDef( Space_def_t *space_def, int utm_zone,
                        double ul_x, double ul_y, double pixel_size,
                        Img_coord_int_t *input_size ) {

    memset( space_def, 0, sizeof(Space_def_t));

    space_def->isin_type = SPACE_NOT_ISIN;
    space_def->orientation_angle = 0.0;  /* should get from metadata OrientationAngle? */

    space_def->proj_num = PROJ_UTM;
    space_def->zone = utm_zone;
    space_def->zone_set = true;

    space_def->sphere = 12; /* ? */

    space_def->pixel_size = pixel_size;
    space_def->ul_corner.x = ul_x; /* TODO(warmerdam): center/corner check? */
    space_def->ul_corner.y = ul_y;
    space_def->ul_corner_set = true;

    space_def->img_size.s = input_size->s;
    space_def->img_size.l = input_size->l;

    return 0;
}
