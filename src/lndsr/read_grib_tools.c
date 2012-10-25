#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "read_grib_tools.h"

#define   LEAPYR(y) (!((y)%400) || (!((y)%4) && ((y)%100)))
#define MARCH 3

short getdoy(short year,short month,short day);

int read_grib_anc(t_ncep_ancillary *anc,int datatype) {
	FILE *fd;
	char where[10],tag[10],date[30];
	int i,grib_ret,ny,nx;
	short year,doy,month,day,hour,minute;
	float time,sec;
	int read_grib_array(FILE *input, char *what, char *where, int *nx, int *ny, float **narray);
	int read_grib_date(FILE *input, char *what, char *where, char *date);

	switch (datatype) {
		case TYPE_OZONE_DATA:
			strcpy(tag,OZONE_GRIBTAG);
			strcpy(where,OZONE_GRIBLEVEL);
			break;
		case TYPE_WV_DATA:
			strcpy(tag,WV_GRIBTAG);
			strcpy(where,WV_GRIBLEVEL);
			break;
		case TYPE_SP_DATA:
			strcpy(tag,SP_GRIBTAG);
			strcpy(where,SP_GRIBLEVEL);
			break;
		case TYPE_ATEMP_DATA:
			strcpy(tag,TMP_GRIBTAG);
			strcpy(where,TMP_GRIBLEVEL);
			break;
		case TYPE_UWND_DATA:
			strcpy(tag,UWD_GRIBTAG);
			strcpy(where,UWD_GRIBLEVEL);
			break;
		case TYPE_VWND_DATA:
			strcpy(tag,VWD_GRIBTAG);
			strcpy(where,VWD_GRIBLEVEL);
			break;
		default:
			return -1;
	}		
	anc->latmin=-90;
	anc->latmax=90;
	anc->lonmin=-180;
	anc->lonmax=180;
	anc->deltalat=1;
	anc->deltalon=1;
	
	anc->nbrows=-1;
	anc->nbcols=-1;
	anc->year=-1;
	anc->doy=-1;
	for (i=0;i<anc->nblayers;i++) {

		printf("reading file %s\n",anc->filename[i]);
		if ((fd=fopen(anc->filename[i],"rb")) != NULL) {
			read_grib_date(fd, tag, where, date);
			printf("date=%s\n",date);
			sscanf(date,"%4hd-%2hd-%2hdT%2hd:%2hd:%f",&year,&month,&day,&hour,&minute,&sec);
			if (anc->year == -1)
				anc->year=year;
			else if (anc->year != year) {
				fprintf(stderr,"ERROR: inconsistent year in %s\n",anc->filename[i]);
				return (-1);
			}
			doy=getdoy(year,month,day);
			if (anc->doy==-1)
				anc->doy=doy;
			else if (anc->doy != doy) {
				fprintf(stderr,"ERROR: inconsistent day in %s\n",anc->filename[i]);
				return (-1);
			}
			anc->time[i]=sec/3600.+ (float)minute/60.+(float)hour;
			printf("date=%04d-%02d-%02dT%02d:%02d:%09.6f   %03d %09.6f\n",year,month,day,hour,minute,sec,anc->doy,anc->time[i]);
			
			grib_ret=read_grib_array(fd, tag, where, &ny, &nx, &(anc->data[i]));
			if (anc->nbrows == -1)
				anc->nbrows = ny;
			else if (anc->nbrows != ny) {
				fprintf(stderr,"ERROR: inconsistent nbrows in %s\n",anc->filename[i]);
				return (-1);
			}
			if (anc->nbcols == -1)
				anc->nbcols = nx;
			else if (anc->nbcols != nx) {
				fprintf(stderr,"ERROR: inconsistent ncols in %s\n",anc->filename[i]);
				return (-1);
			}
			fclose(fd);
		} else
			return -1;
	}
	
	return 0;
}

int interpol_spatial_anc(t_ncep_ancillary *anc,float lat, float lon,float *value) {
/* 
  Point order:

    0 ---- 1    +--> sample
    |      |    |
    |      |    v
    2 ---- 3   line

 */

	typedef struct {
  	int l;                /* line number */
  	int s;                /* sample number */
	} Img_coord_int_t;

  Img_coord_int_t p[4];
  int i, j, n;
  float dl, ds, w;
  float sum[10], sum_w;


	p[0].l = (int)((anc->latmax - lat + anc->deltalat/2.)/anc->deltalat);
	p[0].l = (int)((anc->latmax - lat)/anc->deltalat);
	p[2].l = p[0].l + 1;
	if (p[2].l >= anc->nbrows) {
		p[2].l = anc->nbrows - 1;
		if (p[0].l > 0) p[0].l--;
	}
	p[1].l = p[0].l;
	p[3].l = p[2].l;

	p[0].s = (int)((lon - anc->lonmin - anc->deltalon/2.)/anc->deltalon);
	p[0].s = (int)((lon - anc->lonmin )/anc->deltalon);
	p[1].s = p[0].s + 1;

	if (p[1].s >= anc->nbcols) {
		p[1].s = anc->nbcols - 1;
		if (p[0].s > 0) p[0].s--;
	}    

	p[2].s = p[0].s;
	p[3].s = p[1].s;
	    

	n = 0;
	for (j=0;j<anc->nblayers;j++) 
		sum[j] = 0.0;
	sum_w = 0.0;
  for (i = 0; i < 4; i++) {
    if (p[i].l != -1  &&  p[i].s != -1) {
/*	 
      if ((line_ar[p[i].l][0][p[i].s] == lut->aerosol_fill)||
			 (line_ar[p[i].l][1][p[i].s] == lut->aerosol_fill)||
			 (line_ar[p[i].l][2][p[i].s] == lut->aerosol_fill)) continue;
*/
      dl = (anc->latmax-p[i].l * anc->deltalat)-lat;
      dl = fabs(dl) / anc->deltalat;
      ds = lon - (p[i].s * anc->deltalon+anc->lonmin);
      ds = fabs(ds) / anc->deltalon;
      w = (1.0 - dl) * (1.0 - ds);

      n++;
      sum_w += w;
	for (j=0;j<anc->nblayers;j++) 
      		sum[j] += (anc->data[j][p[i].l*anc->nbcols+p[i].s] * w);
    }
  }

  if (n > 0) {
 		for (j=0;j<anc->nblayers;j++) 
		{
			value[j]=sum[j] / sum_w;
/*			printf("Value %d %f\n",j,value[j]);
			printf("lat %f row %d col %d\n",lat,p[0].l,p[0].s);*/
			
                }
  }

  return 0;
}

/***
int interpol_spatial_anc(t_ncep_ancillary *anc,float lat, float lon,float *value) {
	int i,iline,isamp;
	double rowfactor,colfactor;
	float v1,v2,v3,v4;
	
	if(lat>90.0F||lat<-90.0F||lon>180.0F||lon<-180.0F)
		return(-1);

	iline = (int)(((anc->latmax - lat)/anc->deltalat)+0.5F);
	rowfactor = (double)((anc->latmax - ((float)(iline)*anc->deltalat)) - lat)/anc->deltalat;
	if((rowfactor)<0.0) {
		iline--;
		rowfactor=1.0-fabs(rowfactor);
	}

	/ * special wrap around case * /
	if (lon>(180.0F-(anc->deltalon*0.5F)))	lon=-180.0F-(180.0F-lon);

	isamp = (int)(((lon - anc->lonmin)/anc->deltalon) + 0.5F);
	colfactor = (double)(lon -(anc->lonmin + (isamp*anc->deltalon)))/anc->deltalon;
	if(colfactor<0.0){
		isamp--;
		/ * special wrap around case * /
		if(isamp<0) isamp=0;
	}
	for (i=0;i<anc->nblayers;i++) {
		v1 = anc->data[i][iline * anc->nbcols + isamp];
		v2 = anc->data[i][(iline + 1) * anc->nbcols + isamp];	
		v3 = anc->data[i][(iline + 1) * anc->nbcols + isamp+1];
		v4 = anc->data[i][iline * anc->nbcols + isamp+1];

		value[i]=(float)((1.0-rowfactor)*(1.0-colfactor)*v1 +  
    		rowfactor*(1.0-colfactor)*v2 +
    		rowfactor*colfactor*v3 + 
    		colfactor*(1.0-rowfactor)*v4);
	}
	
	return 0;
}
***/
int free_anc_data(t_ncep_ancillary *anc) {
	int i;
	for (i=0;i<anc->nblayers;i++)
		if (anc->data[i] !=NULL)
			free(anc->data[i]);
	return 0;
}

short getdoy(short year,short month,short day)
{
   int dayatmonth[]=
       {0, 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366 };
   int jdate;

   jdate = day + dayatmonth[month];
   if (!LEAPYR(year) && month >= MARCH)
      jdate--;
   return(jdate);
}
	
int getdaymonth(short year,short doy,short *month,short *day)
{
	int i;
   int daysinmonth[]=
       {0,31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
   short totaldays=0;

	if (LEAPYR(year))
		daysinmonth[2]=29;
		
	for (i=1;(i<=12)&&(totaldays < doy);i++)
		totaldays += daysinmonth[i];

	if ((totaldays < doy)||(doy<=0))
		return -1;
		
	totaldays -= daysinmonth[i-1];
	*month=(short)i-1;
	*day=doy-totaldays;

   return(0);
}

void print_anc_data(t_ncep_ancillary *anc, char* ancftype) {
  int i;
  printf("\n*********************************************\n");
  printf("**** anc file %s ***\n",ancftype);
  printf("**** source  = %d ***\n",anc->source);
  printf("**** nblayers= %d ***\n",anc->nblayers);
  printf("**** year,doy=(%d,%d) ***\n",anc->year,anc->doy);
  printf("**** timeres = %f ***\n",anc->timeres);
  printf("**** time    = ");
  for (i=0; i<anc->nblayers; i++)printf("%8.1f ",anc->time[i]);
  printf("****\n");
  printf("**** latmin = %f ***\n",anc->latmin);
  printf("**** latmax = %f ***\n",anc->latmax);
  printf("**** deltalat = %f ***\n",anc->deltalat);
  printf("**** lonmin = %f ***\n",anc->lonmin);
  printf("**** lonmax = %f ***\n",anc->lonmax);
  printf("**** deltalon = %f ***\n",anc->deltalon);
  printf("**** nbrows = %d ***\n",anc->nbrows);
  printf("**** nbcols = %d ***\n",anc->nbcols);
  printf("*********************************************\n");

}
