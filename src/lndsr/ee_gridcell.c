/**
 * ee_gridcell - EarthEngine support code to prepare the Ar_gridcell_t 
 *               data structure.
 *
 * Closely derived from original lndsr.c and related code.
 *
 * Author: Frank Warmerdam <warmerdam@google.com>
 */

#include "ee_lndsr.h"

#define ANC_TIMERES 6.0 /* 6.0 hours time resolution for most ancillary layers*/
#define ANC_NBLAYERS 4  /* 4 "time layers" for each ancillary layer */


Ar_gridcell_t *ee_build_ar_gridcell(Lut_t *lut,
                                    Space_t *space,
                                    Input_meta_t *meta,
                                    double scene_gmt,
                                    const float *anc_SP,
                                    const float *anc_WV,
                                    const float *anc_O3,
                                    const float *anc_dem)

{
    Ar_gridcell_t *ar_gridcell;

/* -------------------------------------------------------------------- */
/*      Allocate structure and attached buffers.                        */
/* -------------------------------------------------------------------- */
    ar_gridcell = (Ar_gridcell_t *) calloc(sizeof(Ar_gridcell_t),1);

    ar_gridcell->nbrows=lut->ar_size.l;
    ar_gridcell->nbcols=lut->ar_size.s;
    ar_gridcell->lat=(float *)calloc((size_t)(lut->ar_size.s * lut->ar_size.l),sizeof(float));
    if (ar_gridcell->lat == (float *) NULL)
        ERROR("allocating ar_gridcell->lat", "main");
    ar_gridcell->lon=(float *)calloc((size_t)(lut->ar_size.s * lut->ar_size.l),sizeof(float));
    if (ar_gridcell->lon == (float *)NULL)
        ERROR("allocating ar_gridcell->lon", "main");
    ar_gridcell->sun_zen=(float *)calloc((size_t)(lut->ar_size.s * lut->ar_size.l),sizeof(float));
    if (ar_gridcell->sun_zen == (float *)NULL)
        ERROR("allocating ar_gridcell->sun_zen", "main");
    ar_gridcell->view_zen=(float *)calloc((size_t)(lut->ar_size.s * lut->ar_size.l),sizeof(float));
    if (ar_gridcell->view_zen == (float *)NULL)
        ERROR("allocating ar_gridcell->view_zen", "main");
    ar_gridcell->rel_az=(float *)calloc((size_t)(lut->ar_size.s * lut->ar_size.l),sizeof(float));
    if (ar_gridcell->rel_az == (float *)NULL)
        ERROR("allocating ar_gridcell->rel_az", "main");
    ar_gridcell->wv=(float *)calloc((size_t)(lut->ar_size.s * lut->ar_size.l),sizeof(float));
    if (ar_gridcell->wv == (float *)NULL)
        ERROR("allocating ar_gridcell->wv", "main");
    ar_gridcell->spres=(float *)calloc((size_t)(lut->ar_size.s * lut->ar_size.l),sizeof(float));
    if (ar_gridcell->spres == (float *)NULL)
        ERROR("allocating ar_gridcell->spres", "main");
    ar_gridcell->ozone=(float *)calloc((size_t)(lut->ar_size.s * lut->ar_size.l),sizeof(float));
    if (ar_gridcell->ozone == (float *)NULL)
        ERROR("allocating ar_gridcell->ozone", "main");
    ar_gridcell->spres_dem=(float *)calloc((size_t)(lut->ar_size.s * lut->ar_size.l),sizeof(float));
    if (ar_gridcell->spres_dem == (float *)NULL)
        ERROR("allocating ar_gridcell->spres_dem", "main");

/* -------------------------------------------------------------------- */
/*      interpolate ancillary data for AR grid cells (in time only)     */
/* -------------------------------------------------------------------- */
    Img_coord_int_t loc;
    Img_coord_float_t img;
    Geo_coord_t geo;
    int il_ar, is_ar;
    double sum_spres_anc,sum_spres_dem;
    int nb_spres_anc,nb_spres_dem;

    img.is_fill=false;
    sum_spres_anc=0.;
    sum_spres_dem=0.;
    nb_spres_anc=0;
    nb_spres_dem=0;
    for (il_ar = 0; il_ar < lut->ar_size.l;il_ar++) {
        img.l=il_ar*lut->ar_region_size.l+lut->ar_region_size.l/2.;
        for (is_ar=0;is_ar < lut->ar_size.s; is_ar++) {
            float tmpflt_arr[4],tmpflt;
            double coef;
            int tmpint;
            int jday;
            int ipt = is_ar + il_ar * lut->ar_size.s;
            int layer_size = lut->ar_size.s * lut->ar_size.l;

            img.s=is_ar*lut->ar_region_size.s+lut->ar_region_size.s/2.; 
            if (!FromSpace(space, &img, &geo))
                ERROR("mapping from space (1)", "main");

            ar_gridcell->lat[ipt]=geo.lat * DEG;
            ar_gridcell->lon[ipt]=geo.lon * DEG;
            ar_gridcell->sun_zen[ipt]=meta->sun_zen*DEG;
            ar_gridcell->view_zen[ipt]=3.5;
            ar_gridcell->rel_az[ipt]=meta->sun_az*DEG;

            int time_index = (int)(scene_gmt/ANC_TIMERES);
            if (time_index >= tmpint>=(ANC_NBLAYERS-1))
		time_index=ANC_NBLAYERS-2;

            coef=(double)(scene_gmt-time_index*ANC_TIMERES)/ANC_TIMERES;

            ar_gridcell->wv[ipt]=
                (1.-coef)*anc_WV[ipt + time_index*layer_size]
                + coef*anc_WV[ipt + (time_index+1)*layer_size];
            ar_gridcell->wv[ipt] /= 10.0; /* convert original PRWV kg/m2 into g/cm2 */

            if (anc_O3 != NULL) {
                ar_gridcell->ozone[ipt] = anc_O3[ipt];
                ar_gridcell->ozone[ipt] /= 1000.0; /* convert to Dobson */
            } else {
                jday=(short)meta->acq_date.doy;
                ar_gridcell->ozone[il_ar*lut->ar_size.s+is_ar]=calcuoz(jday,(float)ar_gridcell->lat[il_ar*lut->ar_size.s+is_ar]);
            }

            ar_gridcell->spres[ipt]=
                (1.-coef)*anc_SP[ipt + time_index*layer_size]
                + coef*anc_SP[ipt + (time_index+1)*layer_size];
            ar_gridcell->spres[ipt] /= 100.0; /* convert Pascals into millibars */

            if (ar_gridcell->spres[ipt] > 0) {
                sum_spres_anc += ar_gridcell->spres[ipt];
                nb_spres_anc++;
            }

            if (anc_dem[ipt] == -9999)
                ar_gridcell->spres_dem[ipt]=1013;
            else
                ar_gridcell->spres_dem[ipt]=
                    1013.2*exp(-anc_dem[ipt]/8000.);

            /* TODO(warmerdam): debugging */
            if ( is_ar_check_pixel(is_ar, il_ar) ) {
                printf( "  DEBUG: ar cell %d [%d,%d] (%.15g,%.15g)\n", 
                        ipt, is_ar, il_ar,
                        ar_gridcell->lat[ipt], ar_gridcell->lon[ipt] );
                printf( "         sun_zen = %.15g\n", ar_gridcell->sun_zen[ipt] );
                printf( "         view_zen = %.15g\n", ar_gridcell->view_zen[ipt] );
                printf( "         rel_az = %.15g\n", ar_gridcell->rel_az[ipt] );
                printf( "         wv = %.15g\n", ar_gridcell->wv[ipt] );
                printf( "         ozone = %.15g\n", ar_gridcell->ozone[ipt] );
                printf( "         spres = %.15g\n", ar_gridcell->spres[ipt] );
                printf( "         spres_dem = %.15g\n", ar_gridcell->spres_dem[ipt] );
            }
        } /* is_ar */
    } /* il_ar */

    return ar_gridcell;
}
