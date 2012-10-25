#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

/*
 * cc -O0 -o LS_geoloc.o -c LS_geoloc.c 
 */

#define PI      3.14159265358979323846
#define HALF_PI 1.57079632679489661923
#define TWO_PI  6.28318530717958647692
#define EPSLN   1.0e-10
#define R2D     57.2957795131
#define D2R     1.745329251994328e-2
#define S2R     4.848136811095359e-6
#define SQUARE(x)       x * x  

#define DATMCT 20
static double major[20] = {6378206.4, 6378249.145, 6377397.155, 6378157.5,
                    6378388.0, 6378135.0, 6377276.3452, 6378145.0,
                    6378137.0, 6377563.396, 6377304.063, 6377340.189,
                    6378137.0, 6378155.0, 6378160.0, 6378245.0,
                    6378270.0, 6378166.0, 6378150.0, 6370997.0};

static double minor[20] = {6356583.8, 6356514.86955, 6356078.96284, 6356772.2,
                    6356911.94613, 6356750.519915, 6356075.4133,
                    6356759.769356, 6356752.31414, 6356256.91,
                    6356103.039, 6356034.448, 6356752.314245,
                    6356773.3205, 6356774.719, 6356863.0188,
                    6356794.343479, 6356784.283666, 6356768.337303,
                    6370997.0};
static double scale_factor   = 0.9996;
static double false_easting  = 0.0;
static double false_northing = 0.0;

static double r_major;          /* major axis                           */
static double r_minor;          /* minor axis                           */
static double scale_factor;     /* scale factor                         */
static double lon_center;       /* Center longitude (projection center) */
static double lat_origin;       /* center latitude                      */
static double e0,e1,e2,e3;      /* eccentricity constants               */
static double e,es,esp;         /* eccentricity constants               */
static double ml0;              /* small value m                        */
static long ind;		/* sphere flag value			*/

double pixel_size, sin_orien, cos_orien, ul_corner[2];

int LSsphdz(float coordinates[8], double *parm, double *radius, double corner[2])

/*long isph;		/ * spheroid code number also known as datum	*/
/*long zone;		/ * zone code                            	*/
/*double *parm;		/ * projection parameters			*/
/*double *radius;	/ * radius					*/


/* this is meant to be visible from other routines */

{

double r_major;
double r_minor;
int utminvint(double r_maj, double r_min, double scale_fact, long zone);
int utmforint(double r_maj, double r_min, double scale_fact, long zone);
int ret;
long isph;		/* spheroid code number also known as datum	*/
long zone;		/* zone code                            	*/

double t_major;		/* temporary major axis				*/
double t_minor;		/* temporary minor axis				*/
long jsph;		/* spheroid code number				*/
double orient;		/* orientation   				*/

zone = (long)coordinates[4];
isph = (long)coordinates[5];
orient     = coordinates[6];
pixel_size = coordinates[7];
sin_orien = sin(orient);
cos_orien = cos(orient);

ul_corner[0] = corner[0];
ul_corner[1] = corner[1];


if (isph < 0)
   {
   t_major = fabs(parm[0]);
   t_minor = fabs(parm[1]);
   
   if (t_major  > 0.0) 
     {
     if (t_minor > 1.0)
        {
        r_major = t_major;
  	r_minor = t_minor;
	*radius = t_major;
        } 
     else
     if (t_minor > 0.0)
        {
        r_major = t_major;
	*radius = t_major;
        r_minor = (sqrt(1.0 - t_minor)) * t_major; 
        }
     else
        {
        r_major = t_major;
	*radius = t_major;
        r_minor = t_major;
        }
     }
   else
   if (t_minor > 0.0)	/* t_major = 0 */
     {
     r_major = major[0];
     *radius = major[0];
     r_minor = minor[0];
     }
   else
     {
     r_major = major[DATMCT - 1];
     *radius = major[DATMCT - 1];
     r_minor = 6370997.0;
     }
  }
else		/* isph >= 0 */
  {
  jsph = abs(isph);
  if (jsph > 19)
     {
     printf("Invalid spheroid selection\n");
     isph = 1;
     jsph = 0;
     }
  r_major = major[jsph];
  r_minor = minor[jsph];
  *radius = major[DATMCT - 1];

  if (zone == 0) zone = 31L;

  
#ifdef INV
  ret = utminvint(r_major, r_minor, scale_factor, zone);
#else  
  ret = utmforint(r_major, r_minor, scale_factor, zone);
#endif  
  }
return(ret);
}


int utmforint(double r_maj, double r_min, double scale_fact, long zone)

/*double r_maj;			/ * major axis				*/
/*double r_min;			/ * minor axis				*/
/*double scale_fact;		/ * scale factor			*/
/*long zone;			/ * zone number				*/

/* this is meant to be private to this file */

{
double temp;			/* temporary variable			*/
double e0fn(double);
double e1fn(double);
double e2fn(double);
double e3fn(double);
double mlfn(double, double, double, double, double);

if ((abs(zone) < 1) || (abs(zone) > 60))
   {
   printf("Illegal zone number\n");
   return(11);
   }
/*
if (zone == 0)
   zone = 1;
*/
r_major = r_maj;
r_minor = r_min;
scale_factor = scale_fact;
lat_origin = 0.0;
lon_center = ((6 * abs(zone)) - 183) * D2R;
false_easting = 500000.0;
false_northing = (zone < 0) ? 10000000.0 : 0.0;

temp = r_minor / r_major;
es = 1.0 - SQUARE(temp);
e = sqrt(es);
e0 = e0fn(es);
e1 = e1fn(es);
e2 = e2fn(es);
e3 = e3fn(es);
ml0 = r_major * mlfn(e0, e1, e2, e3, lat_origin);
esp = es / (1.0 - es);

if (es < .00001)
   ind = 1;
else 
   ind = 0;

return(0);
}


int utminvint(double r_maj, double r_min, double scale_fact, long zone)

/*double r_maj;			/ * major axis				*/
/*double r_min;			/ * minor axis				*/
/*double scale_fact;		/ * scale factor			*/
/*long zone;			/ * zone number				*/

/* this is meant to be private to this file */

{
double temp;			/* temporary variables			*/
double e0fn(double);
double e1fn(double);
double e2fn(double);
double e3fn(double);
double mlfn(double, double, double, double, double);

if ((abs(zone) < 1) || (abs(zone) > 60))
   {
   printf("Illegal zone number\n");
   return(11);
   }

r_major = r_maj;
r_minor = r_min;
scale_factor = scale_fact;
lat_origin = 0.0;
lon_center = ((6 * abs(zone)) - 183) * D2R;
false_easting = 500000.0;
false_northing = (zone < 0) ? 10000000.0 : 0.0;

temp = r_minor / r_major;
es = 1.0 - SQUARE(temp);
e = sqrt(es);
e0 = e0fn(es);
e1 = e1fn(es);
e2 = e2fn(es);
e3 = e3fn(es);
ml0 = r_major * mlfn(e0, e1, e2, e3, lat_origin);
esp = es / (1.0 - es);

if (es < .00001)
   ind = 1;
else 
   ind = 0;

return(0);
}

int sign(double x)    { if (x < 0.0) return(-1); else return(1);}
double e0fn(double x) {return(1.0-0.25*x*(1.0+x/16.0*(3.0+1.25*x)));}
double e1fn(double x) {return(0.375*x*(1.0+0.25*x*(1.0+0.46875*x)));}
double e2fn(double x) {return(0.05859375*x*x*(1.0+0.75*x));}
double e3fn(double x) {return(x*x*x*(35.0/3072.0));}
double mlfn(double e0, double e1, double e2, double e3, double phi) 
                {return(e0*phi-e1*sin(2.0*phi)+e2*sin(4.0*phi)-e3*sin(6.0*phi));}

double e4fn(double x)
{
 double con;
 double com;
 con = 1.0 + x;
 com = 1.0 - x;
 return (sqrt((pow(con,con))*(pow(com,com))));
 }
 
double asinz (double con)
{
 if (fabs(con) > 1.0)
   {
   if (con > 1.0)
     con = 1.0;
   else
     con = -1.0;
   }
 return(asin(con));
}


int LSutminv(double s, double l, double *lon, double *lat)
/*double x;		/ * (I) X projection coordinate 		*/
/*double y;		/ * (I) Y projection coordinate 		*/
/*double *lon;		/ * (O) Longitude 				*/
/*double *lat;		/ * (O) Latitude 				*/

/* this is meant to be visible from other routines */

{
double x, y;
double con,phi;		/* temporary angles				*/
double delta_phi;	/* difference between longitudes		*/
long i;			/* counter variable				*/
double sin_phi, cos_phi, tan_phi;	/* sin cos and tangent values	*/
double c, cs, t, ts, n, r, d, ds;	/* temporary variables		*/
double f, h, g, temp;			/* temporary variables		*/
long max_iter = 6;			/* maximun number of iterations	*/
double asinz();

  double dl, dp, dy, dx;

  dl = (l + 0.5) * pixel_size;
  dp = (s + 0.5) * pixel_size;

  dy = (dp * sin_orien) - (dl * cos_orien);
  dx = (dp * cos_orien) + (dl * sin_orien);

  y = ul_corner[1] + dy;
  x = ul_corner[0] + dx;

/* fortran code for spherical form 
--------------------------------*/
if (ind != 0)
   {
   f = exp(x/(r_major * scale_factor));
   g = .5 * (f - 1/f);
   temp = lat_origin + y/(r_major * scale_factor);
   h = cos(temp);
   con = sqrt((1.0 - h * h)/(1.0 + g * g));
   *lat = asinz(con);
   if (temp < 0)
     *lat = -*lat;
   if ((g == 0) && (h == 0))
     {
     *lon = lon_center;
     return(0);
     }
   else
     {
     *lon = (atan2(g,h) + lon_center);
     return(0);
     }
   }

/* Inverse equations
  -----------------*/
x = x - false_easting;
y = y - false_northing;

con = (ml0 + y / scale_factor) / r_major;
phi = con;
for (i=0;;i++)
   {
   delta_phi = ((con + e1 * sin(2.0*phi) - e2 * sin(4.0*phi) + e3 * sin(6.0*phi))
               / e0) - phi;
   phi += delta_phi;
   if (fabs(delta_phi) <= EPSLN) break;
   if (i >= max_iter) 
      { 
      printf("Latitude failed to converge\n"); 
      return(95);
      }
   }
if (fabs(phi) < HALF_PI)
   {
   cos_phi = cos(phi);
   sin_phi = sin(phi);
   tan_phi = tan(phi);
   c    = esp * SQUARE(cos_phi);
   cs   = SQUARE(c);
   t    = SQUARE(tan_phi);
   ts   = SQUARE(t);
   con  = 1.0 - es * SQUARE(sin_phi); 
   n    = r_major / sqrt(con);
   r    = n * (1.0 - es) / con;
   d    = x / (n * scale_factor);
   ds   = SQUARE(d);
   *lat = phi - (n * tan_phi * ds / r) * (0.5 - ds / 24.0 * (5.0 + 3.0 * t + 
          10.0 * c - 4.0 * cs - 9.0 * esp - ds / 30.0 * (61.0 + 90.0 * t +
          298.0 * c + 45.0 * ts - 252.0 * esp - 3.0 * cs)));
   *lon = lon_center + (d * (1.0 - ds / 6.0 * (1.0 + 2.0 * t +
          c - ds / 20.0 * (5.0 - 2.0 * c + 28.0 * t - 3.0 * cs + 8.0 * esp +
          24.0 * ts))) / cos_phi);
   }
else
   {
   *lat = HALF_PI * sign(y);
   *lon = lon_center;
   }

  *lat *= R2D;
  *lon *= R2D;

return(0);
}


int LSutmfor(double *x, double *y, double lon, double lat)
/*double *x;		/ * (O) X projection coordinate 		*/
/*double *y;		/ * (O) Y projection coordinate 		*/
/*double lon;		/ * (I) Longitude 				*/
/*double lat;		/ * (I) Latitude 				*/

/* this is meant to be visible from other routines */

{

double delta_lon;	/* Delta longitude (Given longitude - center 	*/
double theta;		/* angle					*/
double delta_theta;	/* adjusted longitude				*/
double mlfn();		/* function to compute small m			*/
double sin_phi, cos_phi;/* sin and cos value				*/
double al, als;		/* temporary values				*/
double b;		/* temporary values				*/
double c, t, tq;	/* temporary values				*/
double con, n, ml;	/* cone constant, small m			*/
  
  double dl, dp, dy, dx;

  lat *= D2R;
  lon *= D2R;

/* Forward equations
  -----------------*/

delta_lon = lon - lon_center;
sin_phi = sin(lat);
cos_phi = cos(lat);

/* This part was in the fortran code and is for the spherical form 
----------------------------------------------------------------*/
if (ind != 0)
  {
  b = cos_phi * sin(delta_lon);
  if ((fabs(fabs(b) - 1.0)) < .0000000001)
     {
     printf("Point projects into infinity, tm-for\n");
     return(93);
     }
  else
     {
     *x = .5 * r_major * scale_factor * log((1.0 + b)/(1.0 - b));
     con = acos(cos_phi * cos(delta_lon)/sqrt(1.0 - b*b));
     if (lat < 0)
        con = - con;
     *y = r_major * scale_factor * (con - lat_origin); 
     return(0);
     }
  }

al  = cos_phi * delta_lon;
als = SQUARE(al);
c   = esp * SQUARE(cos_phi);
tq  = tan(lat);
t   = SQUARE(tq);
con = 1.0 - es * SQUARE(sin_phi);
n   = r_major / sqrt(con);
ml  = r_major * mlfn(e0, e1, e2, e3, lat);


*x  = scale_factor * n * al * (1.0 + als / 6.0 * (1.0 - t + c + als / 20.0 *
      (5.0 - 18.0 * t + SQUARE(t) + 72.0 * c - 58.0 * esp))) + false_easting;

*y  = scale_factor * (ml - ml0 + n * tq * (als * (0.5 + als / 24.0 *
      (5.0 - t + 9.0 * c + 4.0 * SQUARE(c) + als / 30.0 * (61.0 - 58.0 * t
      + SQUARE(t) + 600.0 * c - 330.0 * esp))))) + false_northing;

  *x -= ul_corner[0];
  *y -= ul_corner[1];

  dx = ( (*x) * cos_orien) + ( (*y) * sin_orien);
  dy = ( (*x) * sin_orien) - ( (*y) * cos_orien);

  *x = dx/pixel_size;
  *y = dy/pixel_size;

  *x -= 0.5;
  *y -= 0.5;

return(0);

}

