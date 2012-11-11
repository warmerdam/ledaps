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

/************************************************************************/
/*                           ee_lndsr_main()                            */
/************************************************************************/

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
    const float *anc_dem)

{
    Ar_gridcell_t *ar_gridcell;
    Lut_t *lut;
    Input_meta_t meta;
    Space_def_t space_def;
    Space_t *space = NULL;
    TileDef_t tile_def;
    sixs_tables_t sixs_tables;

    tile_def.full_size.s = full_input_size_s;
    tile_def.full_size.l = full_input_size_l;
    tile_def.offset.s = input_offset_s;
    tile_def.offset.l = input_offset_l;
    tile_def.size.s = input_size_s;
    tile_def.size.l = input_size_l;

    ee_GetInputMeta( &meta, metadata);

    lut = GetLut( "dummy_param_file", 6, &meta, &(tile_def.size) );

    ee_PrepareSpaceDef( &space_def, utm_zone, ul_x, ul_y, pixel_size,
                        &(tile_def.size) );
    space = SetupSpace( &space_def );

/* -------------------------------------------------------------------- */
/*      Compute scene center.                                           */
/* -------------------------------------------------------------------- */
    Img_coord_float_t img;
    Geo_coord_t geo;
    float center_lon, center_lat;
    
    img.l = tile_def.full_size.l / 2.0 - tile_def.offset.l;
    img.s = tile_def.full_size.s / 2.0 - tile_def.offset.s;
    img.is_fill=false;

    if (!FromSpace(space, &img, &geo))
    	RETURN_ERROR("mapping from space (0)", "ee_build_ar_gridcell", 1);

    center_lat=geo.lat * DEG;
    center_lon=geo.lon * DEG;

/* -------------------------------------------------------------------- */
/*      Compute Scene GMT time based on the scene center location.      */
/* -------------------------------------------------------------------- */
    float scene_gmt;

    if ( (meta.acq_date.hour !=0 ) || (meta.acq_date.minute != 0 ) || (meta.acq_date.second !=0))
   
        scene_gmt=meta.acq_date.hour+meta.acq_date.minute/60.+meta.acq_date.second/3600.;
    else
        scene_gmt=10.5-center_lon/15.;

    if ( scene_gmt < 0.) scene_gmt=scene_gmt+24.;

/* -------------------------------------------------------------------- */
/*      Build the ar_gridcell structure with ancillary data sampled     */
/*      on the aerosol grid.                                            */
/* -------------------------------------------------------------------- */
    ar_gridcell = ee_build_ar_gridcell( lut, space, &meta, scene_gmt,
                                        anc_SP, anc_WV, anc_O3, anc_dem );
    if (ar_gridcell == NULL) 
        return 1;

/* -------------------------------------------------------------------- */
/*      Compute deviation from true north.                              */
/* -------------------------------------------------------------------- */
    double delta_x, delta_y, adjust_north;

    printf ("(y0,x0)=(%f,%f)  (lat0,lon0)=(%.16f,%.16f)\n",
            img.l,img.s,(float)(geo.lat * DEG),(float)(geo.lon * DEG));

    delta_y=img.l;
    delta_x=img.s;

    img.l = tile_def.full_size.l / 2.0 - tile_def.offset.l - 100.0;
    img.s = tile_def.full_size.s / 2.0 - tile_def.offset.s;
    img.is_fill=false;
    if (!FromSpace(space, &img, &geo))
    	RETURN_ERROR("mapping from space (0)", "main", 1);

    printf ("(y1,x1)=(%f,%f)  (lat1,lon1)=(%.16f,%.16f)\n",
            img.l,img.s,(float)(geo.lat * DEG),(float)(geo.lon * DEG));

    geo.lon=center_lon*RAD;
    geo.is_fill=false;
    if (!ToSpace(space, &geo, &img))
    	RETURN_ERROR("mapping to space (0)", "main", 1);

    printf ("(y2,x2)=(%f,%f)  (lat2,lon2)=(%.16f,%.16f)\n",
            img.l,img.s,(float)(geo.lat * DEG),(float)(geo.lon * DEG));

    delta_y = delta_y - img.l;
    delta_x = img.s - delta_x;

    printf( "delta_y=%g, delta_x=%g\n", delta_y, delta_x );
    adjust_north=(float)(atan(delta_x/delta_y)*DEG);

    printf("True North adjustment = %f\n",adjust_north);

/* -------------------------------------------------------------------- */
/*      Read 6S from data file.  TODO(warmerdam): allow passing         */
/*      these coefficients in as arguments.                             */
/* -------------------------------------------------------------------- */
    if (read_6S_results_from_file(sixs_filename,&sixs_tables) != 0) {
        return 1;
    }

    switch (meta.inst) {
      case INST_TM:
        sixs_tables.Inst=SIXS_INST_TM;
        break;
      case INST_ETM:
        sixs_tables.Inst=SIXS_INST_ETM;
        break;
      default:
        RETURN_ERROR("Unknown Instrument", "main", 1);
    }

/* -------------------------------------------------------------------- */
/*      Compute atmospheric coefs for the whole scene with              */
/*      aot550=0.01 for use in internal cloud screening : NAZMI         */
/* -------------------------------------------------------------------- */
    int nbpts=lut->ar_size.l*lut->ar_size.s;

    if (allocate_mem_atmos_coeff(nbpts,&atmos_coef))
        ERROR("Allocating memory for atmos_coef", "main");

//    printf("Compute Atmos Params with aot550=0.01\n"); fflush(stdout);
//    update_atmos_coefs(&atmos_coef,ar_gridcell, &sixs_tables,line_ar, lut,input->nband, 1);
    
    return 0;
}

/************************************************************************/
/*                         ee_PrepareSpaceDef()                         */
/************************************************************************/
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
