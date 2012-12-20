/**
 * lndsr_util.c - Support functions pulled out of lndsr.c so they can be used in
 *                ee_lndsr.c and related modules.
 *
 * Closely derived from original lndsr.c and related code.
 *
 * Author: Frank Warmerdam <warmerdam@google.com>
 */

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

/* Prototypes */
#ifndef  HPUX
#define chand chand_
#define csalbr csalbr_
#endif
void chand(float *phi,float *muv,float *mus,float *tau_ray,float *actual_rho_ray);
void csalbr(float *tau_ray,float *actual_S_r);

/* TODO(warmerdam): for debugging */
static int check_pixel_count = 0;
static int check_pixel_x[100];
static int check_pixel_y[100];

static double check_ar_lat = -360.0;
static double check_ar_lon = 0.0;

void add_check_pixel( int is, int il ) {
    check_pixel_x[check_pixel_count] = is;
    check_pixel_y[check_pixel_count] = il;
    check_pixel_count++;
}

void adjust_check_pixels( int offset_s, int offset_l ) {
    int i;
    for( i=0; i < check_pixel_count; i++ )
    {
        check_pixel_x[i] -= offset_s;
        check_pixel_y[i] -= offset_l;
    }
}

int is_check_pixel( int is, int il ) {
    int i;

    for (i=0; i < check_pixel_count; i++) {
        if (is == check_pixel_x[i] && il == check_pixel_y[i]) {
            return i+1;
        }
    }

    return 0;
}

int is_ar_check_pixel( int is_ar, int il_ar ) {
    int i;

    for (i=0; i < check_pixel_count; i++) {
        if (is_ar == check_pixel_x[i]/40 && il_ar == check_pixel_y[i]/40 ) {
            return i+1;
        }
    }

    return 0;
}


int allocate_mem_atmos_coeff(int nbpts,atmos_t *atmos_coef) {
	int ib;
	if ((atmos_coef->computed=(int *)malloc(nbpts*sizeof(int)))==(int *)NULL)
		return -1;
	for (ib=0;ib<7;ib++) {
		if ((atmos_coef->tgOG[ib]=(float *)malloc(nbpts*sizeof(float)))==(float *)NULL)
			return -1;
		if ((atmos_coef->tgH2O[ib]=(float *)malloc(nbpts*sizeof(float)))==(float *)NULL)
			return -1;
		if ((atmos_coef->td_ra[ib]=(float *)malloc(nbpts*sizeof(float)))==(float *)NULL)
			return -1;
		if ((atmos_coef->tu_ra[ib]=(float *)malloc(nbpts*sizeof(float)))==(float *)NULL)
			return -1;
		if ((atmos_coef->rho_mol[ib]=(float *)malloc(nbpts*sizeof(float)))==(float *)NULL)
			return -1;
		if ((atmos_coef->rho_ra[ib]=(float *)malloc(nbpts*sizeof(float)))==(float *)NULL)
			return -1;
		if ((atmos_coef->td_da[ib]=(float *)malloc(nbpts*sizeof(float)))==(float *)NULL)
			return -1;
		if ((atmos_coef->tu_da[ib]=(float *)malloc(nbpts*sizeof(float)))==(float *)NULL)
			return -1;
		if ((atmos_coef->S_ra[ib]=(float *)malloc(nbpts*sizeof(float)))==(float *)NULL)
			return -1;
		if ((atmos_coef->td_r[ib]=(float *)malloc(nbpts*sizeof(float)))==(float *)NULL)
			return -1;
		if ((atmos_coef->tu_r[ib]=(float *)malloc(nbpts*sizeof(float)))==(float *)NULL)
			return -1;
		if ((atmos_coef->S_r[ib]=(float *)malloc(nbpts*sizeof(float)))==(float *)NULL)
			return -1;
		if ((atmos_coef->rho_r[ib]=(float *)malloc(nbpts*sizeof(float)))==(float *)NULL)
			return -1;
	}
	return 0;
}

int free_mem_atmos_coeff(atmos_t *atmos_coef) {
	int ib;
	free(atmos_coef->computed);
	for(ib=0;ib<7;ib++) {
		free(atmos_coef->tgOG[ib]);
		free(atmos_coef->tgH2O[ib]);
		free(atmos_coef->td_ra[ib]);
		free(atmos_coef->tu_ra[ib]);
		free(atmos_coef->rho_mol[ib]);
		free(atmos_coef->rho_ra[ib]);
		free(atmos_coef->td_da[ib]);
		free(atmos_coef->tu_da[ib]);
		free(atmos_coef->S_ra[ib]);
		free(atmos_coef->td_r[ib]);
		free(atmos_coef->tu_r[ib]);
		free(atmos_coef->S_r[ib]);
		free(atmos_coef->rho_r[ib]);
	}
	return 0;
}

float calcuoz(short jday,float flat) {
/******************************************************************************
!C
!Routine: calcuoz

!Description:  Gets ozone concentration for a particular day and 
               latitude, interpolating if necessary.

!Revision History:
 Original version:    Nazmi Z El Saleous and Eric Vermote

!Input Parameters:
        jday       Julian day
        lat        latitude (in degrees)
        
!Output Parameters:
        uoz        Ozone concentration (in cm-atm) 

!Return value:
        returns 0 (success)
               -1 (latitude is beyond coverage of data, and so
                  the returned ozone value is an approximation).
                
                
!References and Credits:
      Data is from:
    LONDON J.,BOJKOV R.D.,OLTMANS S. AND KELLEY J.I.,1976,
    ATLAS OF THE GLOBAL DISTRIBUTION OF TOTAL OZONE .JULY 1957-JUNE 1967
    NCAR TECHNICAL NOTE, NCAR/TN/113+STR,PP276
      
!Developers:
      Nazmi Z El Saleous
      Eric Vermote
      University of Maryland / Dept. of Geography
      nazmi.elsaleous@gsfc.nasa.gov

!Design Notes:  
      
!END
*******************************************************************************/
    float t,u,tmpf;
    int i1,i2,j1,j2,Minf,Msup,Latinf,Latsup;
   

/*
  CREPARTITION ZONALE PAR BANDE DE 10 DEG DE LATITUDE A PARTIR DE 80 SUD
  C            ---   OZONE   ---
*/
    float oz[12][17] = {
        {.315,.320,.315,.305,.300,.280,.260,.240,.240,.240,.250,.280,.320,.350,.375,.380,.380},
        {.280,.300,.300,.300,.280,.270,.260,.240,.240,.240,.260,.300,.340,.380,.400,.420,.420},
        {.280,.280,.280,.280,.280,.260,.250,.240,.250,.250,.270,.300,.340,.400,.420,.440,.440},
        {.280,.280,.280,.280,.280,.260,.250,.250,.250,.260,.280,.300,.340,.380,.420,.430,.430},
        {.280,.290,.300,.300,.280,.270,.260,.250,.250,.260,.270,.300,.320,.360,.380,.400,.400},
        {.280,.300,.300,.305,.300,.280,.260,.250,.250,.260,.260,.280,.310,.330,.360,.370,.370},
        {.290,.300,.315,.320,.305,.280,.260,.250,.240,.240,.260,.270,.290,.310,.320,.320,.320},
        {.300,.310,.320,.325,.320,.300,.270,.260,.240,.240,.250,.260,.280,.290,.300,.300,.290},
        {.300,.320,.325,.335,.320,.300,.280,.260,.240,.240,.240,.260,.270,.280,.280,.280,.280},
        {.320,.340,.350,.345,.330,.300,.280,.260,.240,.240,.240,.260,.260,.280,.280,.280,.280},
        {.360,.360,.360,.340,.320,.300,.280,.260,.240,.240,.240,.260,.280,.300,.310,.310,.300},
        {.340,.350,.340,.320,.310,.280,.260,.250,.240,.240,.240,.260,.300,.320,.330,.340,.330}};
/*     
       C      /Begin of interpolation/
       C      /Find Neighbours/
       C      / First loop for time/
*/
    if (fabs(flat)>=80.) {
        tmpf=.270;
        return -1;
    }
    Minf= (int) ((jday-15.)/30.5);
    if (jday < 15) 
        Minf=Minf-1;
    Latinf=(int) (flat/10.);
    if (flat < 0.) 
        Latinf=Latinf-1;
    t=((jday-15.)-(30.5*Minf))/30.5;
    u=(flat-10.*Latinf)/10.;
   
    Minf=(Minf+12)%12;
    Msup=(Minf+13)%12;
    Latsup=Latinf+1;
    i1=Minf;
    j1=Latinf+8;
    i2=Msup;
    j2=Latsup+8;
   
/*     /Now Calculate Uo3 at the given point Xlat,Jjulien/ */
    tmpf=(1.-t)*(1.-u)*oz[i1][j1]+t*(1.-u)*oz[i2][j1]+t*u*oz[i2][j2]+(1.-t)*u*oz[i1][j2];
    return tmpf;
}

int update_atmos_coefs(atmos_t *atmos_coef,Ar_gridcell_t *ar_gridcell, sixs_tables_t *sixs_tables,int ***line_ar,Lut_t *lut,int nband, int bkgd_aerosol) {
    int irow,icol;
    for (irow=0;irow<ar_gridcell->nbrows;irow++)
        for (icol=0;icol<ar_gridcell->nbcols;icol++)
            update_gridcell_atmos_coefs(irow,icol,atmos_coef,ar_gridcell, sixs_tables,line_ar[irow],lut,nband,bkgd_aerosol);
    return 0;
}

int update_gridcell_atmos_coefs(int irow,int icol,atmos_t *atmos_coef,Ar_gridcell_t *ar_gridcell, sixs_tables_t *sixs_tables,int **line_ar,Lut_t *lut,int nband, int bkgd_aerosol) {
    int ib,ipt,k;
    float mus,muv,phi,ratio_spres,tau_ray,aot550;
    double coef;
    float actual_rho_ray,actual_T_ray_up,actual_T_ray_down,actual_S_r;
    float rho_ray_P0,T_ray_up_P0,T_ray_down_P0,S_r_P0;
    float lamda[7]={486.,570.,660.,835.,1669.,0.,2207.};
    float tau_ray_sealevel[7]={0.16511,0.08614,0.04716,0.01835,0.00113,0.00037}; /* index=5 => band 7 */

    ipt=irow*ar_gridcell->nbcols+icol;	
    mus=cos(ar_gridcell->sun_zen[ipt]*RAD);
    muv=cos(ar_gridcell->view_zen[ipt]*RAD);
    phi=ar_gridcell->rel_az[ipt];
    ratio_spres=ar_gridcell->spres[ipt]/1013.;
/* eric is trying something */
/*                if ((irow == 53 )&&(icol == 83)) 
                  {
                  printf("let's see what is going on\n"); 
                  }
*/
    if (bkgd_aerosol) {
        atmos_coef->computed[ipt]=1;
        aot550=0.01;
    } else {
        if (line_ar[0][icol] != lut->aerosol_fill) {
            atmos_coef->computed[ipt]=1;
            aot550=((float)line_ar[0][icol]/1000.)*pow((550./lamda[0]),-1.);
        } else {
/*
  atmos_coef->computed[ipt]=0;
*/
                                
            atmos_coef->computed[ipt]=1;
            aot550=0.01;
        }
    }

    for (k=1;k<SIXS_NB_AOT;k++) {
        if (aot550 < sixs_tables->aot[k])
            break;
    }
    k--;
    if (k>=(SIXS_NB_AOT-1))
        k=SIXS_NB_AOT-2;
    coef=(aot550-sixs_tables->aot[k])/(sixs_tables->aot[k+1]-sixs_tables->aot[k]);

/*                printf("pressure ratio %d %d %f \n",irow,icol,ratio_spres);*/
    for (ib=0;ib < nband; ib++) {
        atmos_coef->tgOG[ib][ipt]=sixs_tables->T_g_og[ib];				
        atmos_coef->tgH2O[ib][ipt]=sixs_tables->T_g_wv[ib];				
        atmos_coef->td_ra[ib][ipt]=(1.-coef)*sixs_tables->T_ra_down[ib][k]+coef*sixs_tables->T_ra_down[ib][k+1];
        atmos_coef->tu_ra[ib][ipt]=(1.-coef)*sixs_tables->T_ra_up[ib][k]+coef*sixs_tables->T_ra_up[ib][k+1];
        atmos_coef->rho_mol[ib][ipt]=sixs_tables->rho_r[ib];				
        atmos_coef->rho_ra[ib][ipt]=(1.-coef)*sixs_tables->rho_ra[ib][k]+coef*sixs_tables->rho_ra[ib][k+1];
        atmos_coef->td_da[ib][ipt]=(1.-coef)*sixs_tables->T_a_down[ib][k]+coef*sixs_tables->T_a_down[ib][k+1];
        atmos_coef->tu_da[ib][ipt]=(1.-coef)*sixs_tables->T_a_up[ib][k]+coef*sixs_tables->T_a_up[ib][k+1];
        atmos_coef->S_ra[ib][ipt]=(1.-coef)*sixs_tables->S_ra[ib][k]+coef*sixs_tables->S_ra[ib][k+1];
/**
   compute DEM-based pressure correction for each grid point
**/
        tau_ray=tau_ray_sealevel[ib]*ratio_spres;
	
        chand(&phi,&muv,&mus,&tau_ray,&actual_rho_ray);

        actual_T_ray_down=((2./3.+mus)+(2./3.-mus)*exp(-tau_ray/mus))/(4./3.+tau_ray); /* downward */
        actual_T_ray_up = ((2./3.+muv)+(2./3.-muv)*exp(-tau_ray/muv))/(4./3.+tau_ray); /* upward */

        csalbr(&tau_ray,&actual_S_r);
						
        rho_ray_P0=sixs_tables->rho_r[ib];
        T_ray_down_P0=sixs_tables->T_r_down[ib];
        T_ray_up_P0=sixs_tables->T_r_up[ib];
        S_r_P0=sixs_tables->S_r[ib];

        atmos_coef->rho_ra[ib][ipt]=actual_rho_ray+(atmos_coef->rho_ra[ib][ipt]-rho_ray_P0); /* will need to correct for uwv/2 */
        atmos_coef->td_ra[ib][ipt] *= (actual_T_ray_down/T_ray_down_P0);
        atmos_coef->tu_ra[ib][ipt] *= (actual_T_ray_up/T_ray_up_P0);
        atmos_coef->S_ra[ib][ipt] = atmos_coef->S_ra[ib][ipt]-S_r_P0+actual_S_r;
        atmos_coef->td_r[ib][ipt] = actual_T_ray_down;
        atmos_coef->tu_r[ib][ipt]=actual_T_ray_up;
        atmos_coef->S_r[ib][ipt]=actual_S_r;
        atmos_coef->rho_r[ib][ipt]=actual_rho_ray;
					
        if ( is_ar_check_pixel(icol, irow) ) {
            printf( "   line_ar[0][%d]=%d, aot550=%.15g, k=%d, coef=%.15g\n", 
                    icol, line_ar[0][icol], aot550, k, coef );
            report_atmos_coef( atmos_coef, ib, ipt );
        }        
    }

    return 0;
}


void sun_angles (short jday,float gmt,float flat,float flon,float *ts,float *fs) {
      
      double mst,tst,tet,et,ha,delta;
      double dlat,dlon,amuzero,elev,az,caz,azim;
      
      double A1=.000075,A2=.001868,A3=.032077,A4=.014615,A5=.040849;
      double B1=.006918,B2=.399912,B3=.070257,B4=.006758;
      double B5=.000907,B6=.002697,B7=.001480;
      
      dlat=(double)flat*M_PI/180.;
      dlon=(double)flon*M_PI/180.;
		                  
/*      
      SOLAR POSITION (ZENITHAL ANGLE ThetaS,AZIMUTHAL ANGLE PhiS IN DEGREES)
      J IS THE DAY NUMBER IN THE YEAR
 
     MEAN SOLAR TIME (HEURE DECIMALE)
*/
      mst=gmt+(flon)/15.;
      tet=2.*M_PI*(double)jday/365.;

/*    TIME EQUATION (IN MN.DEC) */
      et=A1+A2*cos(tet)-A3*sin(tet)-A4*cos(2.*tet)-A5*sin(2.*tet);
      et=et*12.*60./M_PI;

/*     TRUE SOLAR TIME */

      tst=mst+et/60.;
      tst=(tst-12.);

/*     HOUR ANGLE */

      ha=tst*15.*M_PI/180.;

/*     SOLAR DECLINATION   (IN RADIAN) */

      delta=B1-B2*cos(tet)+B3*sin(tet)-B4*cos(2.*tet)+B5*sin(2.*tet)-
            B6*cos(3.*tet)+B7*sin(3.*tet);

/*    ELEVATION,AZIMUTH */

      amuzero=sin(dlat)*sin(delta)+cos(dlat)*cos(delta)*cos(ha);
      elev=asin(amuzero);
      az=cos(delta)*sin(ha)/cos(elev);
      if (az<-1.)
       az=-1;
      else if (az>1.)
       az=1.;
      caz=(-cos(dlat)*sin(delta)+sin(dlat)*cos(delta)*cos(ha))/cos(elev);
      azim=asin(az);
      if (caz < 0.) 
       azim=M_PI-azim;
      if ((caz > 0.) && (az < 0.)) 
       azim=2*M_PI+azim;
      azim=azim+M_PI;
      if (azim > (2.*M_PI)) 
       azim=azim-2.*M_PI;
      elev=elev*180./M_PI;

/*     CONVERSION IN DEGREES */

      *ts=90.-elev;
      *fs=azim*180./M_PI;
      return;
}


/*
 * make sure the SIXS application exists in the path supplied.
 */
#define MAX_SIXS_PATH_LEN 2048
#ifndef SIXS_APP
#define SIXS_APP   "sixsV1.0B"
#endif
static char sixs_path[MAX_SIXS_PATH_LEN] = "";
void set_sixs_path_from(const char *path)
{
    struct stat stbuf;
    char *cptr = (char *)NULL;
  
    /*
     * if passed in arg is path to a file, get the directory path
     */
    strncpy(sixs_path, path, MAX_SIXS_PATH_LEN);
    /*
     * if path ends in /, remove
     */
    for (cptr = &(sixs_path[strlen(sixs_path)-1]); *cptr == '/' && cptr >= &sixs_path[0]; --cptr)
    {
        *cptr = (char) 0;
    }
  
    if (stat(sixs_path, &stbuf) != 0)    /* make sure the path exists */
    {
        /*    fprintf(stderr, "find_file: can't stat %s\n", sixs_path);
              exit(EXIT_FAILURE);*/
    }
  
    /* if it's a file, find the end of the directory */
    cptr = (char *)NULL;
    if (! S_ISDIR(stbuf.st_mode))   /* not a directory */
    {
        cptr = rindex(sixs_path, '/');
        if (cptr == NULL)            /* current directory */
            sixs_path[0] = (char) 0;
    }
  
   
    /*
     * create path to sixs application file
     */
    if (cptr)
    {
        strcpy((cptr+1), SIXS_APP);
    }
    else
    {
        strcat(sixs_path, SIXS_APP);
    }
   
    /*
     * and make sure it exists
     */
    if (stat(sixs_path, &stbuf) != 0)
    {
        /* fprintf(stderr, "find_file: can't stat %s\n", sixs_path);
           exit(EXIT_FAILURE);*/
    }
}

char *get_sixs_path()
{
    if (strlen(sixs_path) < 1)
    {
        printf("%s\n%s\n%s\n",
               "NO SIXS_PATH defined!",
               "  YOU FORGOT TO CALL set_sixs_path_from",
               "  Cannot continue"
               );
        exit(-1);
    }
  
    return &sixs_path[0];
}


void report_timer(const char *stage) {
    static time_t start_time = 0; 
    static time_t last_time = 0;

    time_t this_time = time(NULL);

    if (start_time == 0) {
        start_time = last_time = this_time;
    }

    printf( "\nTIMER / %s / %ds / %ds / %s\n", 
            stage, 
            (int) (this_time - last_time), 
            (int) (this_time - start_time),
            ctime(&this_time) );
    last_time = this_time;
}

void report_atmos_coef( atmos_t *atmos_coef, int ib, int ic ) {
    if( ib == 0 && ic == 76)
    {
        printf( "" );
    }
    printf( "   atmos_coef[%d,%d]: tgOG=%.15g, tgH2O=%.15g,\n",
            ib, ic,
            atmos_coef->tgOG[ib][ic],
            atmos_coef->tgH2O[ib][ic]);
    printf( "          td_ra=%.15g, tu_ra=%.15g, rho_mol=%.15g,\n",   
            atmos_coef->td_ra[ib][ic],
            atmos_coef->tu_ra[ib][ic],
            atmos_coef->rho_mol[ib][ic] );
    printf( "          rho_ra=%.15g, td_da=%.15g, tu_da=%.15g,\n",   
            atmos_coef->rho_ra[ib][ic],
            atmos_coef->td_da[ib][ic],
            atmos_coef->tu_da[ib][ic] );
    printf( "          S_ra=%.15g, td_r=%.15g, tu_r=%.15g,\n",   
            atmos_coef->S_ra[ib][ic],
            atmos_coef->td_r[ib][ic],
            atmos_coef->tu_r[ib][ic] );
    printf( "          S_r=%.15g, rho_r=%.15g.\n",   
            atmos_coef->S_r[ib][ic],
            atmos_coef->rho_r[ib][ic] );
}
