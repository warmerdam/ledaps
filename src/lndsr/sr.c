
#include "sr.h"
#include "ar.h"
#include "const.h"
#include "sixs_runs.h"


extern atmos_t atmos_coef;
int SrInterpAtmCoef(Lut_t *lut, Img_coord_int_t *input_loc, atmos_t *atmos_coef,atmos_t *interpol_atmos_coef); 

bool Sr(Lut_t *lut, int nsamp, int il, int **line_in,  bool mask_flag, 
        char* mask_line, int **line_out, Sr_stats_t *sr_stats) 
{
  int is;
  bool is_fill, not_fill;
  int ib;
  int iband_ar;
  Img_coord_int_t loc;
  float rho,tmpflt;
  int i_aot,j_aot,ipt;
	atmos_t interpol_atmos_coef;


	allocate_mem_atmos_coeff(1,&interpol_atmos_coef);

  iband_ar = lut->nband;
  loc.l = il;

  i_aot=il/lut->ar_region_size.l;

  for (is = 0; is < nsamp; is++) {
    loc.s = is;
    
    j_aot=is/lut->ar_region_size.s;
    ipt=i_aot*lut->ar_size.s+j_aot;

    not_fill = true;
    for (ib = 0; ib < lut->nband; ib++) {
      if (line_in[ib][is] == lut->in_fill) {
        not_fill = false;
	break;
      }
    }
/*
    if (not_fill)
     / * line_out[iband_ar][is] = ArInterp(lut, &loc, line_ar); * /
      line_out[iband_ar][is] = line_ar[i_aot][j_aot];
    else
      line_out[iband_ar][is] = lut->aerosol_fill;
*/
/*
    if (not_fill) {
      for (ib = 0; ib < lut->nband; ib++) { 
        if (line_in[ib][is] == lut->in_fill) {
	  is_fill = true;
          line_out[ib][is] = lut->output_fill;
          sr_stats->nfill[ib]++;
	} else  
	  line_out[ib][is] = line_in[ib][is];
      }
      sr_stats->nfill[iband_ar]++;
      continue;
    } else {
      if (sr_stats->first[iband_ar]) {

        sr_stats->sr_min[iband_ar] = sr_stats->sr_max[iband_ar] = 
	  line_out[iband_ar][is];
        sr_stats->first[iband_ar] = false;

      } else {

        if (line_out[iband_ar][is] < sr_stats->sr_min[iband_ar])
          sr_stats->sr_min[iband_ar] = line_out[iband_ar][is];

        if (line_out[iband_ar][is] > sr_stats->sr_max[iband_ar])
          sr_stats->sr_max[iband_ar] = line_out[iband_ar][is];
      } 
    }
*/
    is_fill = false;

/*
NAZMI 6/2/04 : correct even cloudy pixels

    if ( mask_flag && (
         mask_line[is]==lut->cloud_snow || 
         mask_line[is]==lut->cloud_cloud ) )is_fill= true; 

*/
	 SrInterpAtmCoef(lut, &loc, &atmos_coef,&interpol_atmos_coef);

    for (ib = 0; ib < lut->nband; ib++) {
      if (line_in[ib][is] == lut->in_fill) {
        is_fill = true;
        line_out[ib][is] = lut->output_fill;
        sr_stats->nfill[ib]++;
      } else {
/* BEGIN introduce correction to handle saturated data Fri-Sep-11-17:49:35-EDT-2009 EV */

      if (line_in[ib][is] != lut->max_valid_sr) {
      rho=(float)line_in[ib][is]/10000.;
      rho=(rho/interpol_atmos_coef.tgOG[ib][0]-interpol_atmos_coef.rho_ra[ib][0]);
		tmpflt=(interpol_atmos_coef.tgH2O[ib][0]*interpol_atmos_coef.td_ra[ib][0]*interpol_atmos_coef.tu_ra[ib][0]);
		rho /= tmpflt;
		rho /= (1.+interpol_atmos_coef.S_ra[ib][0]*rho);


        line_out[ib][is] = (short)(rho*10000.);

        } else
	{
       line_out[ib][is] = lut->max_valid_sr;
	}
	
/* END introduce correction to handle saturated data Fri-Sep-11-17:49:35-EDT-2009 EV */

	if (line_out[ib][is] < lut->min_valid_sr) {
	  sr_stats->nout_range[ib]++;
	  line_out[ib][is] = lut->min_valid_sr;
	}
	if (line_out[ib][is] > lut->max_valid_sr) {
	  sr_stats->nout_range[ib]++;
	  line_out[ib][is] = lut->max_valid_sr;
	}
      }

      if (is_fill) continue;

      if (sr_stats->first[ib]) {

        sr_stats->sr_min[ib] = sr_stats->sr_max[ib] = line_out[ib][is];
        sr_stats->first[ib] = false;

      } else {

        if (line_out[ib][is] < sr_stats->sr_min[ib])
          sr_stats->sr_min[ib] = line_out[ib][is];

        if (line_out[ib][is] > sr_stats->sr_max[ib])
          sr_stats->sr_max[ib] = line_out[ib][is];
      } 
    }

  }
	free_mem_atmos_coeff(&interpol_atmos_coef);

  return true;
}


int SrInterpAtmCoef(Lut_t *lut, Img_coord_int_t *input_loc, atmos_t *atmos_coef,atmos_t *interpol_atmos_coef) 
/* 
  Point order:

    0 ---- 1    +--> sample
    |      |    |
    |      |    v
    2 ---- 3   line

 */

{
  Img_coord_int_t p[4];
  int iar;
  int i, n,ipt, ib;
  double dl, ds, w;
  double sum[7][13], sum_w;
  Img_coord_int_t ar_region_half;

  iar = lut->aerosol_fill;

  ar_region_half.l = (lut->ar_region_size.l + 1) / 2;
  ar_region_half.s = (lut->ar_region_size.s + 1) / 2;

  p[0].l = (input_loc->l - ar_region_half.l) / lut->ar_region_size.l;

  p[2].l = p[0].l + 1;
  if (p[2].l >= lut->ar_size.l) {
    p[2].l = lut->ar_size.l - 1;
    if (p[0].l > 0) p[0].l--;
  }    
      
  p[1].l = p[0].l;
  p[3].l = p[2].l;

  p[0].s = (input_loc->s - ar_region_half.s) / lut->ar_region_size.s;
  p[1].s = p[0].s + 1;

  if (p[1].s >= lut->ar_size.s) {
    p[1].s = lut->ar_size.s - 1;
    if (p[0].s > 0) p[0].s--;
  }    

  p[2].s = p[0].s;
  p[3].s = p[1].s;

  n = 0;
  sum_w = 0.0;
  for (ipt=0;ipt<9;ipt++)
	 for (ib=0;ib<7;ib++)
		sum[ib][ipt]=0.;

  for (i = 0; i < 4; i++) {
    if (p[i].l != -1  &&  p[i].s != -1) {
		ipt=p[i].l * lut->ar_size.s + p[i].s;
		if (!(atmos_coef->computed[ipt])) continue; 

      dl = (input_loc->l - ar_region_half.l) - (p[i].l * lut->ar_region_size.l);
      dl = fabs(dl) / lut->ar_region_size.l;
      ds = (input_loc->s - ar_region_half.s) - (p[i].s * lut->ar_region_size.s);
      ds = fabs(ds) / lut->ar_region_size.s;
      w = (1.0 - dl) * (1.0 - ds);

      n++;
      sum_w += w;

		for (ib=0;ib<6;ib++) {
      	sum[ib][0] += (atmos_coef->tgOG[ib][ipt] * w);
      	sum[ib][1] += (atmos_coef->tgH2O[ib][ipt] * w);
      	sum[ib][2] += (atmos_coef->td_ra[ib][ipt] * w);
      	sum[ib][3] += (atmos_coef->tu_ra[ib][ipt] * w);
       	sum[ib][4] += (atmos_coef->rho_mol[ib][ipt] * w);
      	sum[ib][5] += (atmos_coef->rho_ra[ib][ipt] * w);
      	sum[ib][6] += (atmos_coef->td_da[ib][ipt] * w);
      	sum[ib][7] += (atmos_coef->tu_da[ib][ipt] * w);
     		sum[ib][8] += (atmos_coef->S_ra[ib][ipt] * w);
     		sum[ib][9] += (atmos_coef->td_r[ib][ipt] * w);
     		sum[ib][10] += (atmos_coef->tu_r[ib][ipt] * w);
     		sum[ib][11] += (atmos_coef->S_r[ib][ipt] * w);
     		sum[ib][12] += (atmos_coef->rho_r[ib][ipt] * w);
		}
    }
  }

  if (n > 0) {
	for (ib=0;ib<6;ib++) {
		interpol_atmos_coef->tgOG[ib][0]=sum[ib][0]/sum_w;
		interpol_atmos_coef->tgH2O[ib][0]=sum[ib][1]/sum_w;
		interpol_atmos_coef->td_ra[ib][0]=sum[ib][2]/sum_w;
		interpol_atmos_coef->tu_ra[ib][0]=sum[ib][3]/sum_w;
		interpol_atmos_coef->rho_mol[ib][0]=sum[ib][4]/sum_w;
		interpol_atmos_coef->rho_ra[ib][0]=sum[ib][5]/sum_w;
		interpol_atmos_coef->td_da[ib][0]=sum[ib][6]/sum_w;
		interpol_atmos_coef->tu_da[ib][0]=sum[ib][7]/sum_w;
		interpol_atmos_coef->S_ra[ib][0]=sum[ib][8]/sum_w;
		interpol_atmos_coef->td_r[ib][0]=sum[ib][9]/sum_w;
		interpol_atmos_coef->tu_r[ib][0]=sum[ib][10]/sum_w;
		interpol_atmos_coef->S_r[ib][0]=sum[ib][11]/sum_w;
		interpol_atmos_coef->rho_r[ib][0]=sum[ib][12]/sum_w;

	}
  }

  return 0;
}

/* 
!C******************************************************************************


!Description: 'computeBounds' computes the boundry corners of the output iamge

!Prototype : 

  bool computeBounds(Geo_bounds_t *bounds, Space_t *space, int nps, int nls)
 
!Input Parameters:
 space            space (transformation) definition structure that defines
                  the transformation between geo to map coordinates
                  (Space_t*)

 nps              image size number of samples (int)
 nls              image size number of lines   (int)


!Output Parameters:
 bounds           output boundry structure (Geo_bounds_t)

!Team Unique Header:

 ! Design Notes:
   1. An error status is returned when:
       a. the transformaion can not proceed

   2. Error messages are handled with the 'RETURN_ERROR' macro.

!END****************************************************************************
*/
bool computeBounds(Geo_bounds_t *bounds, Space_t *space, int nps, int nls)
{
  const float pixcorn_x[4]={-0.5,-0.5, 0.5, 0.5}; /* 4 corners of a pixel */
  const float pixcorn_y[4]={-0.5, 0.5,-0.5, 0.5}; /* 4 corners of a pixel */
  Img_coord_float_t img;
  Geo_coord_t geo;
  int i,ix,iy,ic;
  bounds->min_lat=  99999.9999;
  bounds->min_lon=  99999.9999;
  bounds->max_lat= -99999.9999;
  bounds->max_lon= -99999.9999;

  for ( i=0; i<(nps*2+nls*2); i++ )
    {
    if ( i<nps ) /* top edge */
      {
      ix= i;
      iy= 0;
      }
    else if ( i<(nps*2) )
      {
      ix= i-nps;
      iy= (nls-1);
      }
    else if ( i <(nps*2+nls) )
      {
      ix= 0;
      iy= i-2*nps;
      }
    else
      {
      ix= nps-1;
      iy= i-(2*nps+nls);
      }
    for ( ic=0; ic<4; ic++ )
      {
      img.l = (double)iy + pixcorn_y[ic];
      img.s = (double)ix + pixcorn_x[ic];
      img.is_fill = false;
      if (!FromSpace(space, &img, &geo))
        RETURN_ERROR("mapping from sapce", "computeBounds", false);
      bounds->max_lat= max(bounds->max_lat,geo.lat);
      bounds->min_lat= min(bounds->min_lat,geo.lat);
      bounds->max_lon= max(bounds->max_lon,geo.lon);
      bounds->min_lon= min(bounds->min_lon,geo.lon);
      }
    }
  return true;
 }
