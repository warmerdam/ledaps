#include "ar.h"
#include "const.h"
#include "error.h"
#include "sixs_runs.h"
#include "clouds.h"

/* #define VRA_THRESHOLD 0.1 */
#define VRA_THRESHOLD 0.08

extern atmos_t atmos_coef;

int allocate_mem_atmos_coeff(int nbpts,atmos_t *atmos_coef);
int free_mem_atmos_coeff(atmos_t *atmos_coef);
int SrInterpAtmCoef(Lut_t *lut, Img_coord_int_t *input_loc, atmos_t *atmos_coef,atmos_t *interpol_atmos_coef);

bool cloud_detection_pass1(Lut_t *lut, int nsamp, int il, int **line_in, int8 *qa_line, int *b6_line,float *atemp_line,
		  cld_diags_t *cld_diags) {

	int is;
	bool is_fill, not_fill;
	int ib;
	Img_coord_int_t loc;
	float tmpflt;
	float rho1,rho3,rho4,rho5,rho7,t6;
	atmos_t interpol_atmos_coef;
   	int C1,C1p,C2,C2p,C3,C3p,C4,C5,water;
	int cld_row,cld_col;
   float vra,ndvi;


	allocate_mem_atmos_coeff(1,&interpol_atmos_coef);
  	loc.l = il;
	cld_row=il/cld_diags->cellheight;

	for (is = 0; is < nsamp; is++) {
    	loc.s = is;
		cld_col=is/cld_diags->cellwidth;

/*
		if ((line_in[0][is]== lut->in_fill) || (line_in[1][is]== lut->in_fill) || (line_in[2][is]== lut->in_fill) ||
			(line_in[3][is]== lut->in_fill) || (line_in[4][is]== lut->in_fill) ||
			(line_in[5][is]== lut->in_fill) || (b6_line[is]== lut->in_fill))
*/
		if ((qa_line[is]&0x01)==0x01)
			is_fill=true;
		else
			is_fill=false;
		if (! is_fill) {
		  if (((qa_line[is]&0x08)==0x00)||((lut->meta.inst==INST_TM)&&(line_in[2][is]<5000))) { /* no saturation in band 3 */
	 		SrInterpAtmCoef(lut, &loc, &atmos_coef, &interpol_atmos_coef);
			
	 		rho1=line_in[0][is]/10000.;
      		rho1=(rho1/interpol_atmos_coef.tgOG[0][0]-interpol_atmos_coef.rho_ra[0][0]);
			tmpflt=(interpol_atmos_coef.tgH2O[0][0]*interpol_atmos_coef.td_ra[0][0]*interpol_atmos_coef.tu_ra[0][0]);
			rho1 /= tmpflt;
			rho1 /= (1.+interpol_atmos_coef.S_ra[0][0]*rho1);
	 		rho3=line_in[2][is]/10000.;
      		rho3=(rho3/interpol_atmos_coef.tgOG[2][0]-interpol_atmos_coef.rho_ra[2][0]);
			tmpflt=(interpol_atmos_coef.tgH2O[2][0]*interpol_atmos_coef.td_ra[2][0]*interpol_atmos_coef.tu_ra[2][0]);
			rho3 /= tmpflt;
			rho3 /= (1.+interpol_atmos_coef.S_ra[2][0]*rho3);
	 		rho4=line_in[3][is]/10000.;
      		rho4=(rho4/interpol_atmos_coef.tgOG[3][0]-interpol_atmos_coef.rho_ra[3][0]);
			tmpflt=(interpol_atmos_coef.tgH2O[3][0]*interpol_atmos_coef.td_ra[3][0]*interpol_atmos_coef.tu_ra[3][0]);
			rho4 /= tmpflt;
			rho4 /= (1.+interpol_atmos_coef.S_ra[3][0]*rho4);
	 		rho5=line_in[4][is]/10000.;
      		rho5=(rho5/interpol_atmos_coef.tgOG[4][0]-interpol_atmos_coef.rho_ra[4][0]);
			tmpflt=(interpol_atmos_coef.tgH2O[4][0]*interpol_atmos_coef.td_ra[4][0]*interpol_atmos_coef.tu_ra[4][0]);
			rho5 /= tmpflt;
			rho5 /= (1.+interpol_atmos_coef.S_ra[4][0]*rho5);
	 		rho7=line_in[5][is]/10000.;
      		rho7=(rho7/interpol_atmos_coef.tgOG[5][0]-interpol_atmos_coef.rho_ra[5][0]);
			tmpflt=(interpol_atmos_coef.tgH2O[5][0]*interpol_atmos_coef.td_ra[5][0]*interpol_atmos_coef.tu_ra[5][0]);
			rho7 /= tmpflt;
			rho7 /= (1.+interpol_atmos_coef.S_ra[5][0]*rho7);
			t6=b6_line[is]/100.+273.;

			vra=rho1-rho3/2;
					
			C1=(int)(vra>VRA_THRESHOLD);
			C1p=!C1;
/*					
			C2=(t6 < (atemp_line[is]-10.));  
			C2p=(t6 > (atemp_line[is]-10.));
*/
			C2=(t6 < (atemp_line[is]-7.));  
			C2p=(t6 > (atemp_line[is]-7.));
   
			tmpflt=rho4/rho3;
			C3=((tmpflt>=0.9)&&(tmpflt<=1.3));
			C3p=!C3;

			C4=(rho7 > 0.03);

			C5=((rho3 > 0.6)||(rho4 > 0.6));
					
/**
			Water test :
			ndvi < 0 => water
			((0<ndvi<0.1) or (b4<5%)) and b5 < 0.01 => turbid water
**/
			if ((rho4+rho3) != 0)
				ndvi=(rho4-rho3)/(rho4+rho3);
			else
				ndvi=0.01;
			water=(ndvi < 0)||((((ndvi>0)&&(ndvi<0.1))||(rho4<0.05))&&(rho5<0.02));

			if (!water) { /* if not water */
				if ((t6 > (atemp_line[is]-20.)) && (!C5)) { 
					if (!((C1||C3)&&C2&&C4)) { /* clear */
						cld_diags->avg_t6_clear[cld_row][cld_col] += t6;
						cld_diags->std_t6_clear[cld_row][cld_col] += (t6*t6);
						cld_diags->avg_b7_clear[cld_row][cld_col] += rho7;
						cld_diags->std_b7_clear[cld_row][cld_col] += (rho7*rho7);
						cld_diags->nb_t6_clear[cld_row][cld_col] ++;
					}
				} 
			}
		  } 
		} 
	}
	free_mem_atmos_coeff(&interpol_atmos_coef);
	return true;
}

bool cloud_detection_pass2(Lut_t *lut, int nsamp, int il, int **line_in, int8 *qa_line, int *b6_line,
		  cld_diags_t *cld_diags,char *ddv_line) {

/**
use ddv_line to store internal cloud screening info
bit 2 = adjacent cloud 1=yes 0=no
bit 3 = fill value 1=fill 0=valid
bit 4 = land/water mask 1=land 0=water
bit 5 = cloud 0=clear 1=cloudy
bit 6 = cloud shadow 
bit 7 = snow
**/
	int is;
	int il_ar,is_ar;
	bool is_fill, not_fill,thermal_band;
	int ib;
	Img_coord_int_t loc;
	float rho1,rho2,rho3,rho4,rho5,rho7,t6;
	atmos_t interpol_atmos_coef;
   	int C1,C1p,C2,C2p,C3,C3p,C4,C5,water;
	int cld_row,cld_col;
	float vra,ndvi,ndsi,temp_snow_thshld;
	float temp_b6_clear,temp_thshld1,temp_thshld2,avg_b7_clear,atemp_ancillary;
	float tmpflt,tmpflt_arr[10];



	thermal_band=true;
	if (b6_line == NULL) 
		thermal_band=false;
	temp_snow_thshld=380.; /* now flag snow and possibly salt pan */
	allocate_mem_atmos_coeff(1,&interpol_atmos_coef);
  	loc.l = il;
	cld_row=il/cld_diags->cellheight;
	il_ar=il/lut->ar_region_size.l;
	if (il_ar >= lut->ar_size.l)
		il_ar=lut->ar_size.l-1;
	for (is = 0; is < nsamp; is++) {
    	loc.s = is;
		cld_col=is/cld_diags->cellwidth;
		is_ar=is/lut->ar_region_size.s;
		if (is_ar >= lut->ar_size.s)
			is_ar=lut->ar_size.s-1;

		is_fill=false;
		if (thermal_band) {
			if (b6_line[is]== lut->in_fill) {
				is_fill=true;
				ddv_line[is] = 0x08;
			}
		}
/*
		if ((line_in[0][is]== lut->in_fill) || (line_in[1][is]== lut->in_fill) || (line_in[2][is]== lut->in_fill) ||
			(line_in[3][is]== lut->in_fill) || (line_in[4][is]== lut->in_fill) ||
			(line_in[5][is]== lut->in_fill)) {
*/
		if ((qa_line[is]&0x01)==0x01) {
			is_fill=true;
			ddv_line[is] = 0x08;
		}
		if (! is_fill) {
			ddv_line[is] &= 0x44; /* reset all bits except cloud shadow and adjacent cloud */ 

		  if (((qa_line[is]&0x08)==0x08)||((lut->meta.inst==INST_TM)&&(line_in[2][is]>=5000))) {  /* saturated band 3 */
			if (thermal_band) {
				t6=b6_line[is]/100.+273.;

				interpol_clddiags_1pixel(cld_diags, il,is,tmpflt_arr);
				temp_b6_clear=tmpflt_arr[0];
				avg_b7_clear=tmpflt_arr[1];
				atemp_ancillary=tmpflt_arr[2];
				if (temp_b6_clear < 0.) {
					temp_thshld1=atemp_ancillary-20.;
					temp_thshld2=atemp_ancillary-20.;
				} else {
					if (cld_diags->std_t6_clear[cld_row][cld_col] > 0.) {
						temp_thshld1=temp_b6_clear-(cld_diags->std_t6_clear[cld_row][cld_col]+4.);
						temp_thshld2=temp_b6_clear-(cld_diags->std_t6_clear[cld_row][cld_col]);
					} else {
						temp_thshld1=temp_b6_clear-4.;
						temp_thshld2=temp_b6_clear-2.;
					}
				}
		    	if ((((qa_line[is]&0x20)==0x20)||((lut->meta.inst==INST_TM)&&(line_in[4][is]>=5000)))&&(t6 < temp_thshld1)) {  /* saturated band 5 and t6 < threshold => cloudy */
					ddv_line[is] &= 0xbf; /* reset shadow bit */
					ddv_line[is] &= 0xfb; /* reset adjacent cloud bit */
					ddv_line[is] |= 0x20; /* set cloud bit */
				} else if ((line_in[4][is]<2000)&&(t6 < temp_snow_thshld)) { /* snow */
					ddv_line[is] |= 0x80;
				} else { /* asuume cloudy */
					ddv_line[is] &= 0xbf; /* reset shadow bit */
					ddv_line[is] &= 0xfb; /* reset adjacent cloud bit */
					ddv_line[is] |= 0x20; /* set cloud bit */
				}
			}
		  } else {
	 		SrInterpAtmCoef(lut, &loc, &atmos_coef, &interpol_atmos_coef);
			
	 		rho1=line_in[0][is]/10000.;
      		rho1=(rho1/interpol_atmos_coef.tgOG[0][0]-interpol_atmos_coef.rho_ra[0][0]);
			tmpflt=(interpol_atmos_coef.tgH2O[0][0]*interpol_atmos_coef.td_ra[0][0]*interpol_atmos_coef.tu_ra[0][0]);
			rho1 /= tmpflt;
			rho1 /= (1.+interpol_atmos_coef.S_ra[0][0]*rho1);
	 		rho2=line_in[1][is]/10000.;
      		rho2=(rho2/interpol_atmos_coef.tgOG[1][0]-interpol_atmos_coef.rho_ra[1][0]);
			tmpflt=(interpol_atmos_coef.tgH2O[1][0]*interpol_atmos_coef.td_ra[1][0]*interpol_atmos_coef.tu_ra[1][0]);
			rho2 /= tmpflt;
			rho2 /= (1.+interpol_atmos_coef.S_ra[1][0]*rho2);
	 		rho3=line_in[2][is]/10000.;
      		rho3=(rho3/interpol_atmos_coef.tgOG[2][0]-interpol_atmos_coef.rho_ra[2][0]);
			tmpflt=(interpol_atmos_coef.tgH2O[2][0]*interpol_atmos_coef.td_ra[2][0]*interpol_atmos_coef.tu_ra[2][0]);
			rho3 /= tmpflt;
			rho3 /= (1.+interpol_atmos_coef.S_ra[2][0]*rho3);
	 		rho4=line_in[3][is]/10000.;
      		rho4=(rho4/interpol_atmos_coef.tgOG[3][0]-interpol_atmos_coef.rho_ra[3][0]);
			tmpflt=(interpol_atmos_coef.tgH2O[3][0]*interpol_atmos_coef.td_ra[3][0]*interpol_atmos_coef.tu_ra[3][0]);
			rho4 /= tmpflt;
			rho4 /= (1.+interpol_atmos_coef.S_ra[3][0]*rho4);
	 		rho5=line_in[4][is]/10000.;
      		rho5=(rho5/interpol_atmos_coef.tgOG[4][0]-interpol_atmos_coef.rho_ra[4][0]);
			tmpflt=(interpol_atmos_coef.tgH2O[4][0]*interpol_atmos_coef.td_ra[4][0]*interpol_atmos_coef.tu_ra[4][0]);
			rho5 /= tmpflt;
			rho5 /= (1.+interpol_atmos_coef.S_ra[4][0]*rho5);
	 		rho7=line_in[5][is]/10000.;
      		rho7=(rho7/interpol_atmos_coef.tgOG[5][0]-interpol_atmos_coef.rho_ra[5][0]);
			tmpflt=(interpol_atmos_coef.tgH2O[5][0]*interpol_atmos_coef.td_ra[5][0]*interpol_atmos_coef.tu_ra[5][0]);
			rho7 /= tmpflt;
			rho7 /= (1.+interpol_atmos_coef.S_ra[5][0]*rho7);
			if (thermal_band)
				t6=b6_line[is]/100.+273.;

			interpol_clddiags_1pixel(cld_diags, il,is,tmpflt_arr);
			temp_b6_clear=tmpflt_arr[0];
/*			printf("temp_b6_clear %f\n",temp_b6_clear);*/
			avg_b7_clear=tmpflt_arr[1];
			atemp_ancillary=tmpflt_arr[2];
			if (temp_b6_clear < 0.) {
				temp_thshld1=atemp_ancillary-20.;
				temp_thshld2=atemp_ancillary-20.;
			} else {
/**
				temp_thshld1=temp_b6_clear-5.;
				temp_thshld2=temp_b6_clear-5.;
**/
				if (cld_diags->std_t6_clear[cld_row][cld_col] > 0.) {
					temp_thshld1=temp_b6_clear-(cld_diags->std_t6_clear[cld_row][cld_col]+4.);
					temp_thshld2=temp_b6_clear-(cld_diags->std_t6_clear[cld_row][cld_col]);
/*EV add the two following lines on 5-21/-07 and commented them on 5-22-07*/	
/*				
					temp_thshld1=temp_b6_clear+(cld_diags->std_t6_clear[cld_row][cld_col]+4.);
					temp_thshld2=temp_b6_clear+(cld_diags->std_t6_clear[cld_row][cld_col]);
					*/
				} else {
					temp_thshld1=temp_b6_clear-4.;
					temp_thshld2=temp_b6_clear-2.;
				}
			}

			if (thermal_band) {
				vra=rho1-rho3/2.;
					
				C1=(int)(vra>VRA_THRESHOLD);
				C1p=!C1;
					
           		C2=(t6 < temp_thshld1);  
            	C2p=(t6 > temp_thshld1);
   
				tmpflt=rho4/rho3;
				C3=((tmpflt>=0.9)&&(tmpflt<=1.2));
				C3p=!C3;

				C4=(rho7 > 0.03);

            	C5=(t6 < temp_thshld2)&&C1;
			}
/*
			printf ("PASS2: T6=%f THRESH=%f R7=%f  C2=%d  C4=%d\n",t6,temp_thshld1,rho7,C2,C4);
*/					

					
/**
			Water test :
			ndvi < 0 => water
			((0<ndvi<0.1) or (b4<5%)) and b5 < 0.01 => turbid water
**/
			if ((rho4+rho3) != 0)
				ndvi=(rho4-rho3)/(rho4+rho3);
			else
				ndvi=0.01;
			water=(ndvi < 0)||((((ndvi>0)&&(ndvi<0.1))||(rho4<0.05))&&(rho5<0.02));
			if (thermal_band) {
			if (!water) { /* if not water */
				ddv_line[is] |= 0x10;
				if ((C2||C5)&&C4) { /* cloudy */
					ddv_line[is] &= 0xbf; /* reset shadow bit */
					ddv_line[is] &= 0xfb; /* reset adjacent cloud bit */
					ddv_line[is] |= 0x20; /* set cloud bit */
				} else { /* clear */
					ddv_line[is] &= 0xdf;
					ndsi=(rho2-rho5)/(rho2+rho5);
					if ((ndsi > 0.3)&&(t6 < temp_snow_thshld)&&(rho4 > 0.2))
						ddv_line[is] |= 0x80;
				}
			} else 
				ddv_line[is] &= 0xef; 
			} else { /* no thermal band - cannot run cloud mask */
				ddv_line[is] &= 0xdf; /* assume clear */
				if (!water) { /* if not water */
					ddv_line[is] |= 0x10;
				} else {
					ddv_line[is] &= 0xef; 
				}
			}
		  }
		} 
	}
	free_mem_atmos_coeff(&interpol_atmos_coef);
	return true;
}


bool dilate_cloud_mask(Lut_t *lut, int nsamp, char ***cloud_buf, int dilate_dist) {
/**
use ddv_line to store internal cloud screening info
bit 2 = adjacent cloud 1=yes 0=no
bit 3 = fill value 1=fill 0=valid
bit 4 = land/water mask 1=land 0=water
bit 5 = cloud 0=clear 1=cloudy
bit 6 = cloud shadow 
bit 7 = snow
**/

	int il,is,il_adj,is_adj,buf_ind;
	int k,l;

	for (il=0;il<lut->ar_region_size.l;il++) {
		for (is=0;is<nsamp;is++) {

			if (cloud_buf[1][il][is] & 0x20) { /* if cloudy dilate */

				for (k=(il-dilate_dist);k<(il+dilate_dist);k++) {
					il_adj=k;
					buf_ind=1;
					if (k < 0) {
						buf_ind--;
						il_adj += lut->ar_region_size.l; 
					}
					if (k >= lut->ar_region_size.l) {
						buf_ind++;
						il_adj -= lut->ar_region_size.l; 
					}
					if ((il_adj>=0)&&(il_adj<lut->ar_region_size.l)) {
						for (is_adj=(is-dilate_dist);is_adj<(is+dilate_dist);is_adj++) {
							if ((is_adj>=0)&&(is_adj<nsamp)) {
								if (!(cloud_buf[buf_ind][il_adj][is_adj] & 0x20)) { /* if not cloudy */
									cloud_buf[buf_ind][il_adj][is_adj] &= 0xfb; /* reset adjacent cloud bit */
									cloud_buf[buf_ind][il_adj][is_adj] &= 0xbf; /* reset shadow bit */
									cloud_buf[buf_ind][il_adj][is_adj] |= 0x04; /* set adjacent cloud bit */
								}
							}
						}
					}
				}
			} /* if cloudy */
		}
	}
	return true;
}


bool cast_cloud_shadow(Lut_t *lut, int nsamp, int il_start, int ***line_in, int **b6_line,
		  cld_diags_t *cld_diags,char ***cloud_buf, Ar_gridcell_t *ar_gridcell,float pixel_size,float adjust_north) {

	int il,is,il_ar,is_ar,il_clddiags,is_clddiags,shd_buf_ind;
	float t6,temp_b6_clear,atemp_ancillary,tmpflt_arr[10];
	float conv_factor,cld_height,ts,tv,fs,fv,dx,dy;
	int shd_x,shd_y;

/***
			Cloud Shadow
***/
	il_ar=il_start/lut->ar_region_size.l;
	if (il_ar >= lut->ar_size.l)
		il_ar=lut->ar_size.l-1;
	for (il=0;il<lut->ar_region_size.l;il++) {
		for (is=0;is<nsamp;is++) {
			is_ar=is/lut->ar_region_size.s;
			if (is_ar >= lut->ar_size.s)
				is_ar = lut->ar_size.s - 1;

			t6=b6_line[il][is]/100.+273.;
			interpol_clddiags_1pixel(cld_diags, il+il_start,is,tmpflt_arr);
			temp_b6_clear=tmpflt_arr[0];
			atemp_ancillary=tmpflt_arr[2];

			if (cloud_buf[1][il][is] & 0x20) { /* if cloudy cast shadow */
				conv_factor=6.;
				while (conv_factor <= 6.) {
					if (temp_b6_clear>0)
      					cld_height=(temp_b6_clear-t6)/conv_factor;
					else
      					cld_height=(atemp_ancillary-t6)/conv_factor;
      				if (cld_height > 0. ) {
      					ts=ar_gridcell->sun_zen[il_ar*lut->ar_size.s+is_ar]/DEG;
      					fs=(ar_gridcell->rel_az[il_ar*lut->ar_size.s+is_ar]-adjust_north)/DEG;
      					tv=ar_gridcell->view_zen[il_ar*lut->ar_size.s+is_ar]/DEG;
      					fv=0.;
/*
      					dy=sin(fv-2.*M_PI)*tan(tv)*cld_height;
      					dx=cos(fv-2.*M_PI)*tan(tv)*cld_height;
*/
      					dy=sin(2.*M_PI-fv)*tan(tv)*cld_height;
      					dx=cos(2.*M_PI-fv)*tan(tv)*cld_height;
	
						dy=0;
						dx=0;

      					shd_x=is-dx*1000./pixel_size;
      					shd_y=il+dy*1000./pixel_size;
/*
      					dy=sin(fs-2.*M_PI)*tan(ts)*cld_height;
      					dx=cos(fs-2.*M_PI)*tan(ts)*cld_height;
*/
      					dy=cos(fs)*tan(ts)*cld_height;
      					dx=sin(fs)*tan(ts)*cld_height;

      					shd_x=shd_x-dx*1000./pixel_size;
      					shd_y=shd_y+dy*1000./pixel_size;

      					if ((shd_x>=0)&&(shd_x<nsamp)) {
							shd_buf_ind=1;
							if (shd_y<0) {
								shd_buf_ind--;
								shd_y += lut->ar_region_size.l;
							}
							if (shd_y >= lut->ar_region_size.l) {
								shd_buf_ind++;
								shd_y -= lut->ar_region_size.l;
							}
       						if ((shd_y >= 0)&&(shd_y < lut->ar_region_size.l)) {
         						if ((cloud_buf[shd_buf_ind][shd_y][shd_x]&0x20)==0x00) /* if clear observation */
         							cloud_buf[shd_buf_ind][shd_y][shd_x] |= 0x40;
       							}
						}
      				} /*end if cld_height >0*/
					conv_factor += 1.;
				} /* while (conv_fact <= 6.) */
			}
		}
	}
	return true;
}

bool dilate_shadow_mask(Lut_t *lut, int nsamp, char ***cloud_buf, int dilate_dist) {
/**
use ddv_line to store internal cloud screening info
bit 2 = adjacent cloud 1=yes 0=no
bit 3 = fill value 1=fill 0=valid
bit 4 = land/water mask 1=land 0=water
bit 5 = cloud 0=clear 1=cloudy
bit 6 = cloud shadow 
bit 7 = snow
**/

	int il,is,il_adj,is_adj,buf_ind;
	int k,l;
	char *fill_mask;

	if ((fill_mask=(char *)calloc(lut->ar_region_size.l*nsamp,sizeof(char)))==NULL)
		return false;
	for (il=0;il<lut->ar_region_size.l;il++) {
		for (is=0;is<nsamp;is++) {
			if ((cloud_buf[0][il][is] & 0x40)&&(!fill_mask[il*nsamp+is])) { /* if cloud shadow dilate */
				for (k=(il-dilate_dist);k<=(il+dilate_dist);k++) {
					il_adj=k;
					buf_ind=0;
					if (k >= lut->ar_region_size.l) {
						buf_ind++;
						il_adj -= lut->ar_region_size.l; 
					}
					if (k>=0) {
						if ((il_adj>=0)&&(il_adj<lut->ar_region_size.l)) {
							for (is_adj=(is-dilate_dist);is_adj<=(is+dilate_dist);is_adj++) {
								if ((is_adj>=0)&&(is_adj<nsamp)) {
									if (!((cloud_buf[buf_ind][il_adj][is_adj] & 0x20)||
									  	(cloud_buf[buf_ind][il_adj][is_adj] & 0x04)||
									  	(cloud_buf[buf_ind][il_adj][is_adj] & 0x40))) { /* if not cloud, adjacent cloud or cloud shadow */
										cloud_buf[buf_ind][il_adj][is_adj] |= 0x40; /* set adjacent cloud shadow bit */
										fill_mask[il_adj*nsamp+is_adj]=1;
									}
								}
							}
						}
					}
				}
			} /* if cloud shadow */
		}
	}
	free(fill_mask);


	return true;
}

int allocate_cld_diags(cld_diags_t *cld_diags,int cell_height, int cell_width, int scene_height, int scene_width) {
	
	int i,j;
	
	cld_diags->cellheight=cell_height;
	cld_diags->cellwidth=cell_width;
	cld_diags->nbrows=(scene_height-1)/cell_height+1;
	cld_diags->nbcols=(scene_width-1)/cell_width+1;
	if ((cld_diags->avg_t6_clear=(float **)malloc(cld_diags->nbrows*sizeof(float *)))==NULL)
		return -1;
	for (i=0;i<cld_diags->nbrows;i++)
		if ((cld_diags->avg_t6_clear[i]=(float *)calloc(cld_diags->nbcols,sizeof(float)))==NULL)
			return -1;
	if ((cld_diags->std_t6_clear=(float **)malloc(cld_diags->nbrows*sizeof(float *)))==NULL)
		return -1;
	for (i=0;i<cld_diags->nbrows;i++)
		if ((cld_diags->std_t6_clear[i]=(float *)calloc(cld_diags->nbcols,sizeof(float)))==NULL)
			return -1;
	if ((cld_diags->avg_b7_clear=(float **)malloc(cld_diags->nbrows*sizeof(float *)))==NULL)
		return -1;
	for (i=0;i<cld_diags->nbrows;i++) 
		if ((cld_diags->avg_b7_clear[i]=(float *)calloc(cld_diags->nbcols,sizeof(float)))==NULL)
			return -1;
	if ((cld_diags->std_b7_clear=(float **)malloc(cld_diags->nbrows*sizeof(float *)))==NULL)
		return -1;
	for (i=0;i<cld_diags->nbrows;i++) 
		if ((cld_diags->std_b7_clear[i]=(float *)calloc(cld_diags->nbcols,sizeof(float)))==NULL)
			return -1;
	if ((cld_diags->airtemp_2m=(float **)malloc(cld_diags->nbrows*sizeof(float *)))==NULL)
		return -1;
	for (i=0;i<cld_diags->nbrows;i++) 
		if ((cld_diags->airtemp_2m[i]=(float *)calloc(cld_diags->nbcols,sizeof(float)))==NULL)
			return -1;
	if ((cld_diags->nb_t6_clear=(int **)malloc(cld_diags->nbrows*sizeof(int *)))==NULL)
		return -1;
	for (i=0;i<cld_diags->nbrows;i++) 
		if ((cld_diags->nb_t6_clear[i]=(int *)calloc(cld_diags->nbcols,sizeof(int)))==NULL)
			return -1;
	return 0;
}

int free_cld_diags(cld_diags_t *cld_diags) {
	
	int i;
	
	for (i=0;i<cld_diags->nbrows;i++) 
		free(cld_diags->avg_t6_clear[i]);
	free(cld_diags->avg_t6_clear);
	
	for (i=0;i<cld_diags->nbrows;i++) 
		free(cld_diags->std_t6_clear[i]);
	free(cld_diags->std_t6_clear);
	
	for (i=0;i<cld_diags->nbrows;i++) 
		free(cld_diags->avg_b7_clear[i]);
	free(cld_diags->avg_b7_clear);
	
	for (i=0;i<cld_diags->nbrows;i++) 
		free(cld_diags->std_b7_clear[i]);
	free(cld_diags->std_b7_clear);
	
	for (i=0;i<cld_diags->nbrows;i++) 
		free(cld_diags->airtemp_2m[i]);
	free(cld_diags->airtemp_2m);
	
	for (i=0;i<cld_diags->nbrows;i++) 
		free(cld_diags->nb_t6_clear[i]);
	free(cld_diags->nb_t6_clear);
	
	return 0;
}

void fill_cld_diags(cld_diags_t *cld_diags) {
/*
!Description: fill in missing values in the T6Clear grid based
on existing values (spatial interpolation). Missing values have been previously
set to -9999. A filling can be distincted from a original value by looking at the
standart deviation of the optical depth which is set to -9999 for a filling.

!Input parameters: 

	
	
!Output Parameters: 


!Revision History:
$LOG: process_TM.c,v1.0$

Original version: Thu Dec  7 15:46:05 EST 1995


!Developer Header:
		   	      
	      
!References and Credits: 


!Design Notes:   
!END****************************************************************************
*/
   int i,j,k,l,pass,count,band;
   float lastt6,dist,sumt6,sumdist,lastb7,sumb7;
   char missing_flag[300][300];
	int min_nb_values,n,max_distance;
   
	count=0;
	for (i=0;i<cld_diags->nbrows;i++) {
		for (j=0;j<cld_diags->nbcols;j++) {
			missing_flag[i][j]=1;
			if (cld_diags->avg_t6_clear[i][j]!=-9999.)  {
      		count++;
				lastt6=cld_diags->avg_t6_clear[i][j];
				lastb7=cld_diags->avg_b7_clear[i][j];
				missing_flag[i][j]=0;
			}
		}
	} 
	if (count==0)
      return;
	if (count==1) {
		for (i=0;i<cld_diags->nbrows;i++)
			for (j=0;j<cld_diags->nbcols;j++) {
				cld_diags->avg_t6_clear[i][j]=lastt6;
				cld_diags->avg_b7_clear[i][j]=lastb7;
			}
		return;
	}
	
	for (i=0;i<cld_diags->nbrows;i++) {
		for (j=0;j<cld_diags->nbcols;j++) {

/**
Look for at least 3 neighboring valid values within 4 GPs
**/
			min_nb_values=3;
			max_distance=4;	
			pass=0;
			while ((cld_diags->avg_t6_clear[i][j]==-9999.)&&(pass<max_distance)) {
				pass++;
				sumdist=0.;
				sumt6=0.;
				sumb7=0.;
				n=0;
				for (k=i-pass;k<=(i+pass);k++) {
					if ((k>=0)&&(k<cld_diags->nbrows)) {
						for (l=j-pass;l<=(j+pass);l++) {
							if ((l>=0)&&(l<cld_diags->nbcols)) {
								if (!missing_flag[k][l]) {
									dist=sqrt((k-i)*(k-i)+(l-j)*(l-j));
									sumdist += dist;
									sumt6 += (dist*cld_diags->avg_t6_clear[k][l]);
									sumb7 += (dist*cld_diags->avg_b7_clear[k][l]);
									n++;
								}
							}
						}
					}
				}
				if ((n>=min_nb_values)&&(sumdist!=0.)) {
					cld_diags->avg_t6_clear[i][j]=sumt6/sumdist;
					cld_diags->avg_b7_clear[i][j]=sumb7/sumdist;
				}
			}
/**
Look for at least 2 neighboring valid values within 6 GPs
**/
			min_nb_values=2;
			max_distance=6;	
			pass=0;
			while ((cld_diags->avg_t6_clear[i][j]==-9999.)&&(pass<max_distance)) {
				pass++;
				sumdist=0.;
				sumt6=0.;
				sumb7=0.;
				n=0;
				for (k=i-pass;k<=(i+pass);k++) {
					if ((k>=0)&&(k<cld_diags->nbrows)) {
						for (l=j-pass;l<=(j+pass);l++) {
							if ((l>=0)&&(l<cld_diags->nbcols)) {
								if (!missing_flag[k][l]) {
									dist=sqrt((k-i)*(k-i)+(l-j)*(l-j));
									sumdist += dist;
									sumt6 += (dist*cld_diags->avg_t6_clear[k][l]);
									sumb7 += (dist*cld_diags->avg_b7_clear[k][l]);
									n++;
								}
							}
						}
					}
				}
				if ((n>=min_nb_values)&&(sumdist!=0.)) {
					cld_diags->avg_t6_clear[i][j]=sumt6/sumdist;
					cld_diags->avg_b7_clear[i][j]=sumb7/sumdist;
				}
			}
/**
Look for at least 1 neighboring valid values within 10 GPs
**/
			min_nb_values=1;
			max_distance=10;	
			pass=0;
			while ((cld_diags->avg_t6_clear[i][j]==-9999.)&&(pass<max_distance)) {
				pass++;
				sumdist=0.;
				sumt6=0.;
				sumb7=0.;
				n=0;
				for (k=i-pass;k<=(i+pass);k++) {
					if ((k>=0)&&(k<cld_diags->nbrows)) {
						for (l=j-pass;l<=(j+pass);l++) {
							if ((l>=0)&&(l<cld_diags->nbcols)) {
								if (!missing_flag[k][l]) {
									dist=sqrt((k-i)*(k-i)+(l-j)*(l-j));
									sumdist += dist;
									sumt6 += (dist*cld_diags->avg_t6_clear[k][l]);
									sumb7 += (dist*cld_diags->avg_b7_clear[k][l]);
									n++;
								}
							}
						}
					}
				}
				if ((n>=min_nb_values)&&(sumdist!=0.)) {
					cld_diags->avg_t6_clear[i][j]=sumt6/sumdist;
					cld_diags->avg_b7_clear[i][j]=sumb7/sumdist;
				}
			}
			
			
		}
	}
}

int interpol_clddiags_1pixel(cld_diags_t *cld_diags, int img_line, int img_sample,float *inter_value) 
/* 
  Point order:

    0 ---- 1    +--> sample
    |      |    |
    |      |    v
    2 ---- 3   line


	inter_value[0] => t6_clear
	inter_value[1] => b7_clear
	inter_value[3] => airtemp_2m

 */

{
typedef struct {
  int l;                /* line number */
  int s;                /* sample number */
} Img_coord_int_t;
	
  Img_coord_int_t p[4];
  int i, n;
  float dl, ds, w;
  float sum[10], sum_w;

  int cell_half_height,cell_half_width;
	float aot_fill=-1.;

  	for (i=0;i<3;i++) 
    inter_value[i] = -9999.;

  cell_half_height = (cld_diags->cellheight + 1) / 2;
  cell_half_width = (cld_diags->cellwidth + 1) / 2;

  p[0].l = (img_line - cell_half_height) / cld_diags->cellheight;
  if (p[0].l<0)
  	p[0].l=0;
  p[2].l = p[0].l + 1;
  if (p[2].l >= cld_diags->nbrows) {
    p[2].l = cld_diags->nbrows - 1;
    if (p[0].l > 0) p[0].l--;
  }    
      
  p[1].l = p[0].l;
  p[3].l = p[2].l;

  p[0].s = (img_sample - cell_half_width) / cld_diags->cellwidth;
  if (p[0].s < 0)
  	p[0].s=0;
  p[1].s = p[0].s + 1;

  if (p[1].s >= cld_diags->nbcols) {
    p[1].s = cld_diags->nbcols - 1;
    if (p[0].s > 0) p[0].s--;
  }    

  p[2].s = p[0].s;
  p[3].s = p[1].s;

  n = 0;
  sum_w = 0.0;
  for (i=0;i<3;i++)
  	sum[i]=0.;
  for (i = 0; i < 4; i++) {
    if (p[i].l != -1  &&  p[i].s != -1) {
		if (cld_diags->avg_t6_clear[p[i].l][p[i].s] == -9999.)
			 continue;

      dl = fabs(img_line - cell_half_height) - (p[i].l * cld_diags->cellheight);
      dl = fabs(dl) / cld_diags->cellheight; 
      ds = fabs(img_sample - cell_half_width) - (p[i].s * cld_diags->cellwidth);
      ds = fabs(ds) / cld_diags->cellwidth; 
      w = (1.0 - dl) * (1.0 - ds);

      n++;
      sum_w += w;
      sum[0] += (cld_diags->avg_t6_clear[p[i].l][p[i].s] * w);
      sum[1] += (cld_diags->avg_b7_clear[p[i].l][p[i].s] * w);
      sum[2] += (cld_diags->airtemp_2m[p[i].l][p[i].s] * w);
    }
  }

  if ((n > 0)&&(sum_w>0)) {
  	for (i=0;i<3;i++) {
    inter_value[i] = sum[i] / sum_w;
   }
  }

  return 0;
}
