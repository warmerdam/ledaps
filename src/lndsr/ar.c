/***************************************************************
03-08-2005: add other gases to band 7 correction
***************************************************************/
#include "ar.h"
#include "const.h"
#include "error.h"
#include "sixs_runs.h"

#define AOT_MIN_NB_SAMPLES 100

#ifdef DEBUG_AR
extern FILE *fd_ar_diags;
extern int diags_il_ar;
#endif

#ifndef  HPUX
#define chand chand_
#define csalbr csalbr_
#endif
void chand(float *phi,float *muv,float *mus,float *tau_ray,float *actual_rho_ray);
void csalbr(float *tau_ray,float *actual_S_r);

int compute_aot(int band,float rho_toa,float rho_surf_est,float ts,float tv, float phi, float uoz, float uwv, float spres,sixs_tables_t *sixs_tables,float *aot);
int update_gridcell_atmos_coefs(int irow,int icol,atmos_t *atmos_coef,Ar_gridcell_t *ar_gridcell, sixs_tables_t *sixs_tables,int **line_ar,Lut_t *lut,int nband, int bkgd_aerosol);

bool Ar(int il_ar,Lut_t *lut, Img_coord_int_t *size_in, int ***line_in, bool mask_flag, 
        char **mask_line, char **ddv_line, int **line_ar, int **line_ar_stats, Ar_stats_t *ar_stats,
		  Ar_gridcell_t *ar_gridcell,sixs_tables_t *sixs_tables) 
{
/***
ddv_line contains results of cloud_screening when this routine is called
bit 2 = adjacent clouds 1=yes 0=no
bit 3 = fill value 1=fill 0=valid
bit 4 = land/water mask 1=land 0=water
bit 5 = cloud 0=clear 1=cloudy
bit 6 = cloud shadow 
bit 7 = snow

The DDV flag in ddv_line (bit 0) is updated in this routine

***/
  int is, il,i,j;
  int is_ar;
  int is_start, is_end;
  int ib;
  double sum_band[3],sum_band_sq[3];
  double sum_srefl,sum_srefl_sq;
  float rho_surf,a,b;
  short *collect_band[3],*collect_band7,tmp_short;
  int collect_nbsamps;
  
  int nb_all_pixs,nb_water_pixs,nb_fill_pixs,nb_cld_pixs,nb_cldshadow_pixs,nb_snow_pixs;
  float fraction_water,fraction_clouds,fraction_cldshadow,fraction_snow;
  bool is_fill;
  int n,water;
  float ndvi;

	float avg_band[3],std_band[3];
        float avg_srefl,std_srefl;
	float fts,ftv,phi;
	float uwv,uoz,spres;
	float avg_aot,std_aot;
 	int band_number,start_i;

	float T_h2o_b7,T_g_b7,rho7,rho4,MP;
	float rho6,rho1,rho3;
	float a_h2o_b7=-3.7338, b_h2o_b7=0.76348,c_h2o_b7=-.030233;
	float a_CO2_b7=0.0071958, b_CO2_b7=0.55665;
	float a_NO2_b7=0.0013383, b_NO2_b7=0.95109;
	float a_CH4_b7=0.030172, b_CH4_b7=0.79652;


	atmos_t atmos_coef_ar;
	float rho,tmpf;
	int nb_negative_red,nb_red_obs,ipt;

	for (ib=0;ib<3;ib++)
		if ((collect_band[ib]=(short *)malloc(lut->ar_region_size.s*lut->ar_region_size.l*sizeof(short)))==NULL)
			return false;
	if ((collect_band7=(short *)malloc(lut->ar_region_size.s*lut->ar_region_size.l*sizeof(short)))==NULL)
		return false;
/**
	Allocate memory for atmos_coef_ar struct used in filtering aot based on AC red band
**/
	if (allocate_mem_atmos_coeff(ar_gridcell->nbrows*ar_gridcell->nbcols,&atmos_coef_ar))
		return false;

  /* Do for each region along a line */

  for (is_start = 0, is_ar = 0; 
       is_start < size_in->s; 
       is_start += lut->ar_region_size.s, is_ar++) {

    is_end = is_start + lut->ar_region_size.s - 1;
    if (is_end >= size_in->s) is_end = size_in->s - 1;

    n = 0;
	 for (ib=0;ib<3;ib++) {
    	sum_band[ib] = 0.0;
    	sum_band_sq[ib] = 0.0;
	 }
    sum_srefl = 0.0;
    sum_srefl_sq = 0.0;

	 collect_nbsamps=0;
	 
    fts=ar_gridcell->line_sun_zen[is_ar];
    ftv=ar_gridcell->line_view_zen[is_ar];
    phi=ar_gridcell->line_rel_az[is_ar];
    uwv=ar_gridcell->line_wv[is_ar];
    uoz=ar_gridcell->line_ozone[is_ar];
    spres=ar_gridcell->line_spres[is_ar];

/**
compute wv transmittance for band 7
**/
	MP=(1./cos(fts/DEG)+1./cos(ftv/DEG));
	T_h2o_b7=log(MP*uwv);
	T_h2o_b7=a_h2o_b7+b_h2o_b7*T_h2o_b7+c_h2o_b7*T_h2o_b7*T_h2o_b7;
	T_h2o_b7=exp(-exp(T_h2o_b7));
	T_g_b7=T_h2o_b7;
	T_g_b7 /= (1.+a_CO2_b7*pow(MP,b_CO2_b7));
	T_g_b7 /= (1.+a_NO2_b7*pow(MP,b_NO2_b7));
	T_g_b7 /= (1.+a_CH4_b7*pow(MP,b_CH4_b7));


	nb_all_pixs=0;
    nb_water_pixs=0;
    nb_cld_pixs=0;
    nb_snow_pixs=0;
	nb_fill_pixs=0;
    for (il = 0; il < lut->ar_region_size.l; il++) {
      for (is = is_start; is < (is_end + 1); is++) {
		if (ddv_line[il][is]&0x08)
            is_fill = true;
		else
        	is_fill = false;
		nb_all_pixs++;
		if (is_fill)
			nb_fill_pixs++;
	if (!is_fill) {
		water = ((ddv_line[il][is] & 0x10)==0);
		if (water) {
			nb_water_pixs++;
        		is_fill= true; 
		}

/***
exclude clouds & snow pixels flagged by the external cloud mask (ACCA)
***/
    	if ( mask_flag ) {
           if ( mask_line[il][is]==lut->cloud_snow  || 
         	 mask_line[il][is]==lut->cloud_cloud   ) is_fill= true; 
        }
/***
exclude clouds, cloud shadow & snow pixels flagged by the internal cloud mask
***/
		if (((ddv_line[il][is] & 0x20)!=0)||((ddv_line[il][is] & 0x04)!=0)) { /* clouds or adjacent clouds */
			is_fill=true;
			nb_cld_pixs++;
		}
		if ((ddv_line[il][is] & 0x40)!=0) { /* cloud shadow */
			is_fill=true;
			nb_cldshadow_pixs++;
		}
		if ((ddv_line[il][is] & 0x80)!=0) { /* snow */
			is_fill=true;
			nb_snow_pixs++;
		}

	}  
    if (!is_fill) {

		
/* band 7 water vapor correction */	 
	 rho7=line_in[il][5][is]/10000.;
	 rho4=line_in[il][3][is]/10000.;
	 rho7 /= T_g_b7;  /* correct for water vapor and other gases*/
	 

   /* update sums if dark target, dark target if not water and 0.015 < rho7 < 0.05 */
	ddv_line[il][is] &= 0xfe;  /* set bit 0 to 0 */


   if ((rho7>0.015) && (rho4 >0.10) /* &&(rho7<0.05) */) {
	  n++;
	 for (ib=0;ib<3;ib++) {
	  sum_band[ib] += (line_in[il][ib][is]/10000.); 
	  sum_band_sq[ib] += (line_in[il][ib][is]/10000.)*(line_in[il][ib][is]/10000.);
	   collect_band[ib][collect_nbsamps]=line_in[il][ib][is];
	 }
	  sum_srefl += rho7; 
	  sum_srefl_sq += (rho7*rho7); 
	 collect_band7[collect_nbsamps]=(short)(rho7*10000.);
	 collect_nbsamps++;
	  if (rho7<0.05)
	  	ddv_line[il][is] |= 0x01; /* set bit 0 to 1 */
	}
      
	}
      }
    }
    
#ifdef DEBUG_AR
	if (fd_ar_diags!=NULL)
   		fprintf(fd_ar_diags,"%d %d %d",diags_il_ar,is_ar,collect_nbsamps);
#endif
    if (collect_nbsamps == 0) {
      line_ar[0][is_ar] = lut->aerosol_fill;
      line_ar[1][is_ar] = lut->aerosol_fill;
      line_ar[2][is_ar] = lut->aerosol_fill;
      line_ar_stats[0][is_ar] = 0;
      line_ar_stats[1][is_ar] = lut->in_fill;
      line_ar_stats[2][is_ar] = lut->in_fill;
      ar_stats->nfill++;
#ifdef DEBUG_AR
	  if (fd_ar_diags != NULL) {
      	for (ib=0;ib<3;ib++) 
         	fprintf(fd_ar_diags," %f %f",-1.,-1.);
      	fprintf(fd_ar_diags," %f %f",-1.,-1.);
      	fprintf(fd_ar_diags," %f %f %f %f %f %f %f\n",ar_gridcell->line_sun_zen[is_ar],
             ar_gridcell->line_view_zen[is_ar],ar_gridcell->line_rel_az[is_ar],
             ar_gridcell->line_wv[is_ar],ar_gridcell->line_ozone[is_ar],
             ar_gridcell->line_spres[is_ar],-1.);
	  }
#endif
    } else {

/**
		Sort collected observations
**/
/*    printf("Sort %d obs\n",collect_nbsamps); fflush(stdout); */
		for (i=0;i<(collect_nbsamps-1);i++) {
			for (j=i+1;j<collect_nbsamps;j++) {
				if (collect_band[0][j] < collect_band[0][i]) {
					for (ib=0;ib<3;ib++) {
						tmp_short=collect_band[ib][i];
						collect_band[ib][i]=collect_band[ib][j];
						collect_band[ib][j]=tmp_short;
					}
					tmp_short=collect_band7[i];
					collect_band7[i]=collect_band7[j];
					collect_band7[j]=tmp_short;
				}
			}
		}


 /*		collect_nbsamps=(int)(collect_nbsamps *0.1);  Take 10% of the observations */
		if (collect_nbsamps >= 2*AOT_MIN_NB_SAMPLES) {

/*
		start_i=0;
		while (start_i < (collect_nbsamps-AOT_MIN_NB_SAMPLES)) {
			sum_srefl=0.;
			sum_srefl_sq=0.;
			for (i=start_i;i<(start_i+AOT_MIN_NB_SAMPLES);i++) {
				sum_srefl += (collect_band7[i]/10000.);
				sum_srefl_sq += ((collect_band7[i]/10000.)*(collect_band7[i]/10000.));
			}
      		avg_srefl = (sum_srefl) / AOT_MIN_NB_SAMPLES; 
			std_srefl=sqrt(fabs(sum_srefl_sq-(sum_srefl*sum_srefl/AOT_MIN_NB_SAMPLES))/(AOT_MIN_NB_SAMPLES-1));
			printf("IL_AR=%d IS_AR=%d START_I=%d AVG_REFL=%f STD_REFL=%f\n",il_ar,is_ar,start_i,avg_srefl,std_srefl);
			
			start_i+=AOT_MIN_NB_SAMPLES;
		}
*/
		start_i=AOT_MIN_NB_SAMPLES;
		collect_nbsamps=AOT_MIN_NB_SAMPLES; /* Take the first AOT_MIN_NB_SAMPLES samples only */
		for (ib=0;ib<3;ib++) {
			sum_band[ib]=0.;
			sum_band_sq[ib]=0.;
			for (i=0;i<collect_nbsamps;i++) {
				sum_band[ib] += (collect_band[ib][i+start_i]/10000.);
				sum_band_sq[ib] += ((collect_band[ib][i+start_i]/10000.)*(collect_band[ib][i+start_i]/10000.));
			}
		}
		sum_srefl=0.;
		sum_srefl_sq=0.;
		for (i=0;i<collect_nbsamps;i++) {
			sum_srefl += (collect_band7[i+start_i]/10000.);
			sum_srefl_sq += ((collect_band7[i+start_i]/10000.)*(collect_band7[i+start_i]/10000.));
		}
		/* update stats line */
      avg_srefl = (sum_srefl) / collect_nbsamps; 
		for (ib=0;ib<3;ib++)
			avg_band[ib]=sum_band[ib]/collect_nbsamps; 
			
      if (collect_nbsamps>1) {
        std_srefl=((sum_srefl_sq)+collect_nbsamps*avg_srefl*avg_srefl-2.*avg_srefl*(sum_srefl))/(collect_nbsamps-1);
        if (std_srefl>0)
         std_srefl=sqrt(std_srefl);
        else 
          std_srefl=0;
			for (ib=0;ib<3;ib++) {
        		std_band[ib]=((sum_band_sq[ib])+collect_nbsamps*avg_band[ib]*avg_band[ib]-2.*avg_band[ib]*(sum_band[ib]))/(collect_nbsamps-1);
        		if (std_band[ib]>0)
         		std_band[ib]=sqrt(std_band[ib]);
        		else 
          		std_band[ib]=0;
			}
      } else {
        std_srefl=0.;
		  for (ib=0;ib<3;ib++) {
          		std_band[ib]=0;
      	}
		}
      line_ar_stats[0][is_ar] = collect_nbsamps;
      line_ar_stats[1][is_ar] = (int)(avg_srefl*10000.);
      line_ar_stats[2][is_ar] = 0;
/*
      line_ar_stats[2][is_ar] = (int)(((float)nb_water_pixs/nb_all_pixs)*1000.);
*/
	fraction_water=(float)nb_water_pixs/(nb_all_pixs-nb_fill_pixs);
	fraction_clouds=(float)nb_cld_pixs/(nb_all_pixs-nb_fill_pixs);
	fraction_cldshadow=(float)nb_cldshadow_pixs/(nb_all_pixs-nb_fill_pixs);
	fraction_snow=(float)nb_snow_pixs/(nb_all_pixs-nb_fill_pixs);
	if (nb_snow_pixs >= 5 ) {
	line_ar_stats[0][is_ar]=-100;
	/* resetting the number of usable observation if snow pixel present EFV March 2007*/
	}
#ifdef DEBUG_AR
	  if (fd_ar_diags != NULL) {
      	for (ib=0;ib<3;ib++) 
         	fprintf(fd_ar_diags," %f %f",avg_band[ib],std_band[ib]);
      	fprintf(fd_ar_diags," %f %f",avg_srefl,std_srefl);
      	fprintf(fd_ar_diags," %f %f %f %f %f %f %f %f %f %f",ar_gridcell->line_sun_zen[is_ar],
             ar_gridcell->line_view_zen[is_ar],ar_gridcell->line_rel_az[is_ar],
             ar_gridcell->line_wv[is_ar],ar_gridcell->line_ozone[is_ar],
             ar_gridcell->line_spres[is_ar],fraction_water,fraction_clouds,
             fraction_cldshadow,fraction_snow);
	  }
#endif
		
/**
	Compute AOT blue band
***/

/*    printf("Compute AOT\n"); fflush(stdout); */

	 if ((std_srefl <= 1.015) && (avg_srefl <= 0.15) && (nb_snow_pixs < 5 )&& (fraction_water < 0.3) && (fraction_clouds < 1e-10)) {
		rho_surf=0.33*avg_srefl;
				
	   compute_aot(0,avg_band[0],avg_band[2],fts,ftv,phi,uoz,uwv,spres,sixs_tables,&avg_aot);
/*			printf("HOW MANY REAL ONE?\n");*/

		std_aot=0.;
/* eric debug 3/07/2011 begin*/
/* eric commented after sucessfull debug 4/29/2011 begin*/
/*      avg_aot=0.01;    */
/* eric commented after sucessfull debug 4/29/2011 end*/
/* eric debug 3/07/2011 end*/
      	
      line_ar[0][is_ar] = (int)(avg_aot*1000.);
      line_ar_stats[2][is_ar] =(int)(avg_aot*1000.);
#ifdef DEBUG_AR
	  if (fd_ar_diags!=NULL)
      	fprintf(fd_ar_diags," %f",avg_aot);
#endif

/***
	Filter aot : Correct red band using retreived aot. if over 30% of the corrected refelctances are
    negative, reject aot.
***/
		if (update_gridcell_atmos_coefs(il_ar,is_ar,&atmos_coef_ar,ar_gridcell,sixs_tables,line_ar,lut,6, 0))
			return false;
		printf("is this going here ? \n");	
		ib=2; /*  test with red band */
		nb_red_obs=0;
		nb_negative_red=0;
		ipt=il_ar*lut->ar_size.s+is_ar;
    	for (il = 0; il < lut->ar_region_size.l; il++) {
      		for (is = is_start; is < (is_end + 1); is++) {
			if (!(ddv_line[il][is]&0x08)) {
				rho=(float)line_in[il][ib][is]/10000.;
     	                        rho7=(float)line_in[il][5][is]/10000.;
     	                        rho4=(float)line_in[il][3][is]/10000.;
     	                        rho6=(float)line_in[il][4][is]/10000.;
     	                        rho1=(float)line_in[il][0][is]/10000.;
	                        rho7 /= T_g_b7;  /* correct for water vapor and other gases*/
     			        rho=(rho/atmos_coef_ar.tgOG[ib][ipt]-atmos_coef_ar.rho_ra[ib][ipt]);
				rho /= (atmos_coef_ar.tgH2O[ib][ipt]*atmos_coef_ar.td_ra[ib][ipt]*atmos_coef_ar.tu_ra[ib][ipt]);
				rho /= (1.+atmos_coef_ar.S_ra[ib][ipt]*rho);
     			        rho4=(rho/atmos_coef_ar.tgOG[3][ipt]-atmos_coef_ar.rho_ra[3][ipt]);
				rho4 /= (atmos_coef_ar.tgH2O[3][ipt]*atmos_coef_ar.td_ra[3][ipt]*atmos_coef_ar.tu_ra[3][ipt]);
				rho4 /= (1.+atmos_coef_ar.S_ra[3][ipt]*rho4);
     			        rho6=(rho/atmos_coef_ar.tgOG[4][ipt]-atmos_coef_ar.rho_ra[4][ipt]);
				rho6 /= (atmos_coef_ar.tgH2O[4][ipt]*atmos_coef_ar.td_ra[4][ipt]*atmos_coef_ar.tu_ra[4][ipt]);
				rho6 /= (1.+atmos_coef_ar.S_ra[4][ipt]*rho6);
     			        rho1=(rho/atmos_coef_ar.tgOG[0][ipt]-atmos_coef_ar.rho_ra[0][ipt]);
				rho1 /= (atmos_coef_ar.tgH2O[0][ipt]*atmos_coef_ar.td_ra[0][ipt]*atmos_coef_ar.tu_ra[0][ipt]);
				rho1 /= (1.+atmos_coef_ar.S_ra[0][ipt]*rho1);
				nb_red_obs++;
			
				if ((rho < 0.) || (rho > rho7 )) /*eric introduced that to get rid of the salt pan */
					nb_negative_red++;
/*				if ((rho7 > 0.2) && (((rho4-rho6)/(rho4+rho6)) < 0.1) && (avg_aot > 0.3) && ((rho/rho1) > 2.0))
				        {
					printf(" got some negative red\n");
			                nb_negative_red++;
					}*/
			}
			}
		}
		printf(" percent negative %f \n",(float)nb_negative_red/(float)nb_red_obs);
/*		if (((float)nb_negative_red/(float)nb_red_obs) > 0.3) { Eric Change for test*/
		if (((float)nb_negative_red/(float)nb_red_obs) > 0.01) {
			line_ar[0][is_ar]=lut->aerosol_fill;
			line_ar_stats[0][is_ar]=-100;
			printf("I FOUND A BAD ONE\n");
#ifdef DEBUG_AR
	  		if (fd_ar_diags!=NULL)
      			fprintf(fd_ar_diags," -1. %f\n",((float)nb_negative_red/(float)nb_red_obs));
#endif
		} else {
#ifdef DEBUG_AR
	  		if (fd_ar_diags!=NULL)
      			fprintf(fd_ar_diags," %f %f\n",avg_aot,((float)nb_negative_red/(float)nb_red_obs));
#endif
		}




/*     printf("DONE\n"); fflush(stdout); */
/**
	Compute AOT green band
		rho_surf=avg_srefl*1.4;  / * estimated reflectance in band2 = 7/5*band7 * /
				
	   compute_aot(1,avg_band[1],rho_surf,fts,ftv,phi,uoz,uwv,spres,sixs_tables,&avg_aot);
		std_aot=0.;
      line_ar[1][is_ar] = (int)(avg_aot*1000.);
***/
/**
	Compute AOT red band
		rho_surf=avg_srefl*0.5;  / * estimated reflectance in band3 = band7/2 * /
				
	   compute_aot(2,avg_band[2],rho_surf,fts,ftv,phi,uoz,uwv,spres,sixs_tables,&avg_aot);
		std_aot=0.;
      line_ar[2][is_ar] = (int)(avg_aot*1000.);

***/
      if (ar_stats->first) {

        ar_stats->ar_min = ar_stats->ar_max = line_ar[0][is_ar];
        ar_stats->first = false;

      } else {

        if (line_ar[0][is_ar] < ar_stats->ar_min)
          ar_stats->ar_min = line_ar[0][is_ar];

        if (line_ar[0][is_ar] > ar_stats->ar_max)
          ar_stats->ar_max = line_ar[0][is_ar];
      }
 	  } else {
      	line_ar[0][is_ar] = lut->aerosol_fill;
      	line_ar[1][is_ar] = lut->aerosol_fill;
      	line_ar[2][is_ar] = lut->aerosol_fill;
#ifdef DEBUG_AR
	  if (fd_ar_diags)
      	fprintf(fd_ar_diags," -1.\n");
#endif

	  }
	  } else {
      line_ar[0][is_ar] = lut->aerosol_fill;
      line_ar[1][is_ar] = lut->aerosol_fill;
      line_ar[2][is_ar] = lut->aerosol_fill;
/*      printf("DOES IT HAPPEN HERE\n");*/
      line_ar_stats[0][is_ar] = 0;
      line_ar_stats[1][is_ar] = lut->in_fill;
      line_ar_stats[2][is_ar] = lut->in_fill;
      ar_stats->nfill++;
#ifdef DEBUG_AR
	  if (fd_ar_diags!=NULL) {
      	for (ib=0;ib<3;ib++)
         	fprintf(fd_ar_diags," %f %f",-1.,-1.);
      	fprintf(fd_ar_diags," %f %f",-1.,-1.);
      	fprintf(fd_ar_diags," %f %f %f %f %f %f %f\n",ar_gridcell->line_sun_zen[is_ar],
             ar_gridcell->line_view_zen[is_ar],ar_gridcell->line_rel_az[is_ar],
             ar_gridcell->line_wv[is_ar],ar_gridcell->line_ozone[is_ar],
             ar_gridcell->line_spres[is_ar],-1.);
	  }
#endif

	  }
    }
  }
	for (ib=0;ib<3;ib++)
		free(collect_band[ib]);
	free(collect_band7);

	if(free_mem_atmos_coeff(&atmos_coef_ar))
		return false;

  return true;
}


int compute_aot(int band,float toarhoblue,float toarhored,float ts,float tv, float phi, float uoz, float uwv, float spres,sixs_tables_t *sixs_tables,float *aot){
	int i,iaot,ioat;
	float minimum,temp,eratio;
	float slope;
	float aot2;
	float surrhoblue[SIXS_NB_AOT],surrhored[SIXS_NB_AOT],aot_tab[SIXS_NB_AOT];
	double coef;
        float temp1,temp2;
	float actual_rho_ray,actual_T_ray,actual_S_r;
	float mus,muv,tau_ray,ratio;
	float tau_ray_sealevel[7]={0.16511,0.08614,0.04716,0.01835,0.00113,0.00037}; /* index=5 => band 7 */
	
/* correct the blue band */	
	band=0;
	mus=cos(ts*RAD);
	muv=cos(tv*RAD);
	ratio=spres/1013.;
	tau_ray=tau_ray_sealevel[band]*ratio;
	
	chand(&phi,&muv,&mus,&tau_ray,&actual_rho_ray);

	actual_T_ray=((2./3.+mus)+(2./3.-mus)*exp(-tau_ray/mus))/(4./3.+tau_ray); /* downward */
	actual_T_ray *= ((2./3.+muv)+(2./3.-muv)*exp(-tau_ray/muv))/(4./3.+tau_ray); /* total */

	csalbr(&tau_ray,&actual_S_r);
	

	for (i=0;i<SIXS_NB_AOT;i++) {
	        surrhoblue[i]=toarhoblue/sixs_tables->T_g_og[band];
		surrhoblue[i]=surrhoblue[i]-(actual_rho_ray+(sixs_tables->rho_ra[band][i]-sixs_tables->rho_r[band]));
		surrhoblue[i]=surrhoblue[i]/actual_T_ray*sixs_tables->T_a[band][i]*sixs_tables->T_g_wv[band];
		surrhoblue[i]=surrhoblue[i]/((1.+(actual_S_r+(sixs_tables->S_ra[band][i]-sixs_tables->S_r[band]))*surrhoblue[i]));
	}
	
/* correct the red band */	
	band=2;
	mus=cos(ts*RAD);
	muv=cos(tv*RAD);
	ratio=spres/1013.;
	tau_ray=tau_ray_sealevel[band]*ratio;
	
	chand(&phi,&muv,&mus,&tau_ray,&actual_rho_ray);

	actual_T_ray=((2./3.+mus)+(2./3.-mus)*exp(-tau_ray/mus))/(4./3.+tau_ray); /* downward */
	actual_T_ray *= ((2./3.+muv)+(2./3.-muv)*exp(-tau_ray/muv))/(4./3.+tau_ray); /* total */

	csalbr(&tau_ray,&actual_S_r);

	for (i=0;i<SIXS_NB_AOT;i++) {
	        surrhored[i]=toarhored/sixs_tables->T_g_og[band];
		surrhored[i]=surrhored[i]-(actual_rho_ray+(sixs_tables->rho_ra[band][i]-sixs_tables->rho_r[band]));
		surrhored[i]=surrhored[i]/actual_T_ray*sixs_tables->T_a[band][i]*sixs_tables->T_g_wv[band];
		surrhored[i]=surrhored[i]/((1.+(actual_S_r+(sixs_tables->S_ra[band][i]-sixs_tables->S_r[band]))*surrhored[i]));
	}
	
/* compute aot based on empirical ratio */	
	minimum=9999999;
	eratio=0.66;
	iaot=-1;
	for (i=0;i<SIXS_NB_AOT;i++) {
	    if ( surrhoblue[i] > 0. ) {
	        temp=surrhoblue[i]-eratio*surrhored[i];
	        if (fabs(temp) < minimum ) {
	            minimum=temp;
	            iaot=i;
		}
	     }
	
	}
	if (iaot==-1)  {
	    *aot=0.01;
	    }
	    else
	    {
        	if (iaot==0) iaot=1;
	        temp1=surrhoblue[iaot-1]-eratio*surrhored[iaot-1];
	        temp2=surrhoblue[iaot]-eratio*surrhored[iaot];
		slope=(sixs_tables->aot_wavelength[1][iaot]-sixs_tables->aot_wavelength[1][iaot-1])/(temp2-temp1);
/* let's hope all the aot table are independent of wavelength but that need to be checked */		
		*aot=sixs_tables->aot_wavelength[1][iaot-1]-temp1/slope;
		
	    }
/*	&aot2=*aot;*/
        printf(" variables in compaot toarhored %f toarhoblue %f aot %f\n",toarhored,toarhoblue,*aot);
        if (*aot < 0.01)
                *aot=0.01;

	return 0;		
}

int ArInterp(Lut_t *lut, Img_coord_int_t *input_loc, int ***line_ar,int *inter_aot) 
/* 
  Point order:

    0 ---- 1    +--> sample
    |      |    |
    |      |    v
    2 ---- 3   line

 */

{
  Img_coord_int_t p[4];
  int i, n;
  float dl, ds, w;
  float sum[3], sum_w;
  Img_coord_int_t ar_region_half;

  inter_aot[0] = lut->aerosol_fill;
  inter_aot[1] = lut->aerosol_fill;
  inter_aot[2] = lut->aerosol_fill;

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
  sum[0]=sum[1]=sum[2] = sum_w = 0.0;
  for (i = 0; i < 4; i++) {
    if (p[i].l != -1  &&  p[i].s != -1) {
      if (line_ar[p[i].l][0][p[i].s] == lut->aerosol_fill) continue;

      dl = (input_loc->l - ar_region_half.l) - (p[i].l * lut->ar_region_size.l);
      dl = fabs(dl) / lut->ar_region_size.l;
      ds = (input_loc->s - ar_region_half.s) - (p[i].s * lut->ar_region_size.s);
      ds = fabs(ds) / lut->ar_region_size.s;
      w = (1.0 - dl) * (1.0 - ds);

      n++;
      sum_w += w;
      sum[0] += (line_ar[p[i].l][0][p[i].s] * w);
      sum[1] += (line_ar[p[i].l][1][p[i].s] * w);
      sum[2] += (line_ar[p[i].l][2][p[i].s] * w);
    }
  }

  if ((n > 0)&&(sum_w>0)) {
    inter_aot[0] = floor((sum[0] / sum_w) + 0.5);
    inter_aot[1] = floor((sum[1] / sum_w) + 0.5);
    inter_aot[2] = floor((sum[2] / sum_w) + 0.5);
  }

  return 0;
}

int Old_Fill_Ar_Gaps(Lut_t *lut, int ***line_ar, int ib) 
/* 
  Point order:

    0 ---- 1    +--> sample
    |      |    |
    |      |    v
    2 ---- 3   line

 */

{
  Img_coord_int_t *gaps_loc;
  int *gaps_neighbors,nb_gaps;
  int i, j, i_aot,j_aot,n,more_gaps;
  float dl, ds, w;
  float sum, sum_w;



  /**
  allocate memory for gaps location and number of neighbors
  **/
  gaps_loc=(Img_coord_int_t *)malloc(lut->ar_size.l*lut->ar_size.s*sizeof(Img_coord_int_t));
  if (gaps_loc == NULL)
   ERROR("failed to allocate memory for gaps_loc","Fill_Ar_Gaps");
  gaps_neighbors=(int *)malloc(lut->ar_size.l*lut->ar_size.s*sizeof(int));
  if (gaps_neighbors == NULL)
   ERROR("failed to allocate memory for gaps_neighbors","Fill_Ar_Gaps");

  more_gaps=1;
  while (more_gaps) {
   more_gaps=0;
  /**
  get location of gaps and number of their valid neighbors
  **/
  nb_gaps=0;
  for (i=0;i<lut->ar_size.l;i++)
    for (j=0;j<lut->ar_size.s;j++)
      if (line_ar[i][ib][j] == lut->aerosol_fill) {  /* check aot in the red band since it is used to fill the others*/
        gaps_loc[nb_gaps].l=i;
        gaps_loc[nb_gaps].s=j;
		  gaps_neighbors[nb_gaps]=0;
		  for (i_aot=(i-1);i_aot<=(i+1);i_aot++) 
			 if ((i_aot>=0)&&(i_aot<lut->ar_size.l))
		  		for (j_aot=(j-1);j_aot<=(j+1);j_aot++) 
			 		if ((j_aot>=0)&&(j_aot<lut->ar_size.s))
						if (line_ar[i_aot][ib][j_aot] != lut->aerosol_fill)
							gaps_neighbors[nb_gaps]++;
			nb_gaps++;
      }

  /**
  sort list of gaps by decreasing number of neighbors
  **/
  for (i=0;i<(nb_gaps-1);i++) 
    for (j=(i+1);j<nb_gaps;j++)
      if (gaps_neighbors[i]<gaps_neighbors[j]) {
        n=gaps_neighbors[j];
		  gaps_neighbors[j]=gaps_neighbors[i];
		  gaps_neighbors[i]=n;
		  n=gaps_loc[j].l;
		  gaps_loc[j].l=gaps_loc[i].l;
		  gaps_loc[i].l=n;
		  n=gaps_loc[j].s;
		  gaps_loc[j].s=gaps_loc[i].s;
		  gaps_loc[i].s=n;
	   }

  for (i=0;i<nb_gaps;i++) { 
	 n=0;
    sum=0.;
    sum_w=0.;
	 for (i_aot=(gaps_loc[i].l-1);i_aot<=(gaps_loc[i].l+1);i_aot++)
		if ((i_aot>=0)&&(i_aot<lut->ar_size.l))
	 		for (j_aot=(gaps_loc[i].s-1);j_aot<=(gaps_loc[i].s+1);j_aot++)
				if ((j_aot>=0)&&(j_aot<lut->ar_size.s)) 
					if (line_ar[i_aot][ib][j_aot] != lut->aerosol_fill) {
      				dl = fabs((gaps_loc[i].l - i_aot)/3.);
      				ds = fabs((gaps_loc[i].s - j_aot)/3.);
      				w = (1.0 - dl) * (1.0 - ds);
      				n++;
      				sum_w += w;
      				sum += (line_ar[i_aot][ib][j_aot] * w);
						
					}
    if (n>0) {
      line_ar[gaps_loc[i].l][ib][gaps_loc[i].s] = floor((sum / sum_w) + 0.5);
		if (line_ar[gaps_loc[i].l][ib][gaps_loc[i].s] < 20)
			line_ar[gaps_loc[i].l][ib][gaps_loc[i].s]=20;
    } else
	 	more_gaps=1;
  }
  
  }
  
  free(gaps_loc);
  free(gaps_neighbors);
  return 0;
}


int Fill_Ar_Gaps(Lut_t *lut, int ***line_ar, int ib) {
/*
!Description: fill in missing values in the OZONE grid based
on existing values (spatial interpolation). 

!Input parameters: 

	
	
!Output Parameters: 


!Revision History:

Original version: Thu Dec  7 15:46:05 EST 1995


!Developer Header:
		   	      
	      
!References and Credits: 


!Design Notes:   
!END****************************************************************************
*/
   int i,j,k,l,pass,count,last_value,More_Gaps,nbfills;
   float dist,sum_value,sum_dist;
   char **missing_flag,*missing_flag_array;
   int min_nb_values,n,max_distance;

	if ((missing_flag_array=(char *)calloc(lut->ar_size.l*lut->ar_size.s,sizeof(char)))==NULL)
		return 0;
	if ((missing_flag=(char **)calloc(lut->ar_size.l,sizeof(char *)))==NULL)
		return 0;
	for (i=0;i<lut->ar_size.l;i++)
		missing_flag[i]=&missing_flag_array[i*lut->ar_size.s];
/**
Start by counting valid values
if nb gaps = 0 do nothing
if nb gaps = 1 duplicate value everywhere
**/
	count=0;
	for (i=0;i<lut->ar_size.l;i++) {
		for (j=0;j<lut->ar_size.s;j++) {
			if (line_ar[i][ib][j] != lut->aerosol_fill)  {
      			count++;
				last_value=line_ar[i][ib][j];
			}
		}
	} 
	if (count==0)
      return 0;
	if (count==1) {
		for (i=0;i<lut->ar_size.l;i++)
			for (j=0;j<lut->ar_size.s;j++) {
				line_ar[i][ib][j]=last_value;
			}
		return 0;
	}
 	More_Gaps=1;
	nbfills=1;
	while ((More_Gaps)&&(nbfills != 0)) {
 		More_Gaps=0;
                nbfills=0;
		for (i=0;i<lut->ar_size.l;i++) {
			for (j=0;j<lut->ar_size.s;j++) {
				missing_flag[i][j]=0;
				if (line_ar[i][ib][j]==lut->aerosol_fill)  {
					missing_flag[i][j]=1;
					More_Gaps=1;
				}
			}
		}
		if (More_Gaps) {
	
			for (i=0;i<lut->ar_size.l;i++) {
				for (j=0;j<lut->ar_size.s;j++) {

/**
Look for at least 3 neighboring valid values within 4 GPs
smooth
**/
					min_nb_values=3;
					max_distance=3;	
					if (missing_flag[i][j]) {
						sum_dist=0.;
						sum_value=0.;
						n=0;
						for (k=i-max_distance;k<=(i+max_distance);k++) {
							if ((k>=0)&&(k<lut->ar_size.l)) {
								for (l=j-max_distance;l<=(j+max_distance);l++) {
									if ((l>=0)&&(l<lut->ar_size.s)) {
										if (!missing_flag[k][l]) {
											dist=sqrt((k-i)*(k-i)+(l-j)*(l-j));
											sum_dist += dist;
											sum_value += (dist*line_ar[k][ib][l]);
											n++;
										}
									}
								}
							}
						}
						if ((n>=min_nb_values)&&(sum_dist!=0.)) {
/*						printf("Gaps in %d %d  band %d was filled \n",i,j,ib); */
							line_ar[i][ib][j]=sum_value/sum_dist;
							nbfills=nbfills+1;
						}
						else
						{
/*						printf("Gaps in %d %d  band %d was NOT filled \n",i,j,ib);*/
						}
					}
				}
			}
		}
 	}
	if ((More_Gaps)&&(nbfills == 0)) {
	for (i=0;i<lut->ar_size.l;i++) {
		for (j=0;j<lut->ar_size.s;j++) {
			if (line_ar[i][ib][j] == lut->aerosol_fill)  line_ar[i][ib][j]=60;
			}
			}
			}
	   
 	return 0;
}
