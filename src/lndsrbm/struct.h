#ifndef uchar
#define uchar unsigned char
#endif

/*****************************************************************/
/* Extra types */
/*****************************************************************/
   
typedef struct {
   long proj;
   int xsize,ysize;
   int center_x,center_y;
   short step_x,step_y,step_gx,step_gy;
   int nbvalues_x,nbvalues_y,nbvalues_gx,nbvalues_gy;
   float latmin,latmax,deltalat,lat0;
   float lonmin,lonmax,deltalon,lon0;
   double OriginLat,CentMerid;
   double CentLat,CentLon;
   double StdPar1,StdPar2;
   double LongPol,TrueScale;
   double x_pixelsize,y_pixelsize;
   double Factor;
   double center_easting,center_northing;
   char file27[256],file83[256];
   double proj_parm[15];
   long zonec;
   long outdatum;
   } t_project;
   
