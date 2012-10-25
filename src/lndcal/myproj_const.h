/*
!C****************************************************************************

!File: myproj_const.h

!Description: Header file for declaring constants describing the map 
 projections, spheres, and map projection parameters.

!Revision History:
 Revision 1.0 2001/05/08
 Robert Wolfe
 Original Version.

 Revision 1.1 2001/07/24
 Robert Wolfe
 1) Fixed projection parameters for ISINUS and HOM
 2) Used "radius" for 19 projections where an ellipsoidal version is not 
    available.
 3) Changed the number of parameters from 16 to 15 to match the GCTP package.


!Team Unique Header:
  This software was developed by the MODIS Land Science Team Support 
  Group for the Labatory for Terrestrial Physics (Code 922) at the 
  National Aeronautics and Space Administration, Goddard Space Flight 
  Center, under NASA Task 92-012-00.

 ! References and Credits:

  ! MODIS Science Team Member:
      Christopher O. Justice
      MODIS Land Science Team           University of Maryland
      justice@hermes.geog.umd.edu       Dept. of Geography
      phone: 301-405-1600               1113 LeFrak Hall
                                        College Park, MD, 20742

  ! Developers:
      Robert E. Wolfe (Code 922)
      MODIS Land Team Support Group     Raytheon ITSS
      robert.e.wolfe.1@gsfc.nasa.gov    4400 Forbes Blvd.
      phone: 301-614-5508               Lanham, MD 20706  

 ! Design Notes:
   1. Constants are declared describing the map projections, spheres, and 
      map projection parameters.

!END****************************************************************************
*/

#ifndef PROJ_CONST_H
#define PROJ_CONST_H

#include "myproj.h"

/* Constants describing the supported spheres */

const Proj_sphere_t Proj_sphere[PROJ_NSPHERE] = { 
  {"Clarke 1866",           6378206.4,    6356583.8},  
  {"Clarke 1880",           6378249.145,  6356514.86955},
  {"Bessel",                6377397.155,  6356078.96284}, 
  {"International 1967",    6378157.5,    6356772.2}, 
  {"International 1909",    6378388.0,    6356911.94613}, 
  {"WGS 72",                6378135.0,    6356750.519915}, 
  {"Everest",               6377276.3452, 6356075.4133}, 
  {"WGS 66",                6378145.0,    6356759.769356}, 
  {"GRS 1980/WGS 84",       6378137.0,    6356752.31414}, 
  {"Airy",                  6377563.396,  6356256.91}, 
  {"Modified Everest",      6377304.063,  6356103.039},  
  {"Modified Airy",         6377340.189,  6356034.448}, 
  {"Walbeck",               6378137.0,    6356752.314245}, 
  {"Southeast Asia",        6378155.0,    6356773.3205}, 
  {"Australian National",   6378160.0,    6356774.719}, 
  {"Krassovsky",            6378245.0,    6356863.0188}, 
  {"Hough",                 6378270.0,    6356794.343479}, 
  {"Mercury 1960",          6378166.0,    6356784.283666}, 
  {"Modified Mercury 1968", 6378150.0,    6356768.337303}, 
  {"Sphere of Radius 6370997 meters", 6370997.0, 6370997.0}};

/* Constants describing the supported map projections */

const Proj_type_t Proj_type[PROJ_NPROJ] = {
  {PROJ_GEO,    "GEO",    "Geographic"}, 
  {PROJ_UTM,    "UTM",    "Universal Transverse Mercator (UTM)"}, 
  {PROJ_SPCS,   "SPCS",   "State Plane Coordinates"}, 
  {PROJ_ALBERS, "ALBERS", "Albers Conical Equal Area"}, 
  {PROJ_LAMCC,  "LAMCC",  "Lambert Conformal Conic"}, 
  {PROJ_MERCAT, "MERCAT", "Mercator"}, 
  {PROJ_PS ,    "PS",     "Polar Stereographic"}, 
  {PROJ_POLYC,  "POLYC",  "Polyconic"}, 
  {PROJ_EQUIDC, "EQUIDC", "Equidistant Conic"}, 
  {PROJ_TM,     "TM",     "Transverse Mercator"}, 
  {PROJ_STEREO, "STEREO", "Stereographic"}, 
  {PROJ_LAMAZ,  "LAMAZ",  "Lambert Azimuthal Equal Area"}, 
  {PROJ_AZMEQD, "AZMEQD", "Azimuthal Equidistant"}, 
  {PROJ_GNOMON, "GNOMON", "Gnomonic"}, 
  {PROJ_ORTHO,  "ORTHO",  "Orthographic"}, 
  {PROJ_GVNSP,  "GVNSP",  "General Vertical Near-Side Perspective"}, 
  {PROJ_SNSOID, "SNSOID", "Sinusoidal"}, 
  {PROJ_EQRECT, "EQRECT", "Equirectangular"}, 
  {PROJ_MILLER, "MILLER", "Miller Cylindrical"}, 
  {PROJ_VGRINT, "VGRINT", "Van der Grinten"}, 
  {PROJ_HOM,    "HOM",    "(Hotine) Oblique Mercator"}, 
  {PROJ_ROBIN,  "ROBIN",  "Robinson"}, 
  {PROJ_SOM,    "SOM",    "Space Oblique Mercator (SOM)"}, 
  {PROJ_ALASKA, "ALASKA", "Alaska Conformal"}, 
  {PROJ_GOODE,  "GOODE",  "Interrupted Goode Homolosine"}, 
  {PROJ_MOLL,   "MOLL",   "Mollweide"}, 
  {PROJ_IMOLL,  "IMOLL",  "Interrupted Mollweide"}, 
  {PROJ_HAMMER, "HAMMER", "Hammer"}, 
  {PROJ_WAGIV,  "WAGIV",  "Wagner IV"}, 
  {PROJ_WAGVII, "WAGVII", "Wagner VII"}, 
  {PROJ_OBEQA,  "OBEQA",  "Oblated Equal Area"}, 
  {PROJ_ISINUS, "ISINUS", "Integerized Sinusoidial"}};
  /*   {PROJ_ISINUS_OLD, "ISINUS", "Integerized Sinusoidial"} */

/* Constants describing the map projection parameters */

const Proj_param_type_t Proj_param_type[PROJ_PARAM_NTYPE] = {
  {PROJ_PARAM_NULL,           "NULL",           "(null value)"}, 
  {PROJ_PARAM_MAJOR_AXIS,     "MAJOR_AXIS",     "Major axis (meters)"}, 
  {PROJ_PARAM_MINOR_AXIS,     "MINOR_AXIS",     "Minor axis (meters)"},  
  {PROJ_PARAM_RADIUS,         "RADIUS",         "Radius (meters)"},  
  {PROJ_PARAM_FALSE_EASTING,  "FALSE_EASTING",  "False easting (meters)"}, 
  {PROJ_PARAM_FALSE_NORTHING, "FALSE_NORTHING", "False northing (meters)"}, 
  {PROJ_PARAM_LAT1,           "LAT1",           "1st standard parallel"}, 
  {PROJ_PARAM_LAT2,           "LAT2",           "2nd standard parallel"}, 
  {PROJ_PARAM_CENTER_LONG,    "CENTER_LONG",    "Center longitude"}, 
  {PROJ_PARAM_LAT_ORIGIN,     "LAT_ORIGIN",     "Latitude at origin"}, 
  {PROJ_PARAM_MODE,           "MODE",           
                                               "Initialization method to use"}, 
  {PROJ_PARAM_SCALE_FACTOR,   "SCALE_FACTOR",   "Scale factor"}, 
  {PROJ_PARAM_CENTER_LAT,     "CENTER_LAT",     "Center latitude"}, 
  {PROJ_PARAM_AZIMUTH,        "AZIMUTH",        "Azimuth"}, 
  {PROJ_PARAM_LON1,           "LON1",           "1st longitude"}, 
  {PROJ_PARAM_LON2,           "LON2",           "2nd longitude"}, 
  {PROJ_PARAM_LON_ORIGIN,     "LON_ORIGIN",     "Longitude at origin"}, 
  {PROJ_PARAM_H,              "H",              "Height above sphere"}, 
  {PROJ_PARAM_PATH,           "PATH",           "Path number"},  
  {PROJ_PARAM_SATNUM,         "SATNUM",         "Satellite number"}, 
  {PROJ_PARAM_ALF,            "ALF",            "Angle"},
  {PROJ_PARAM_TIME,           "TIME",           "Time"},
  {PROJ_PARAM_SAT_RATIO,      "SAT_RATIO",      
                            "Satellite ratio which specifies the start point"}, 
  {PROJ_PARAM_START,          "START",          
                                              "Start flag - beginning or end"}, 
  {PROJ_PARAM_SHAPE_M,        "SHAPE_M",        "Constant 'm'"}, 
  {PROJ_PARAM_SHAPE_N,        "SHAPE_N",        "Constant 'n'"}, 
  {PROJ_PARAM_ANGLE,          "ANGLE",          "Rotation angle"}, 
  {PROJ_PARAM_NZONE,          "NZONE",          
                                               "Number of longitudinal zones"}, 
  {PROJ_PARAM_JUSTIFY,        "JUSTIFY",        "Justify flag"}}; 

/* Constants describing the map projection parameters for each map projection */

#define PROJ_NPARAM (15)

/* Need special cases for utm: scale_factor always 0.9996           */
/*                         if (zone == 0) zone calculated with      */
/*                            outparam[0]: long, outparam[1]: lat   */
/*                            if outdatum < 0: clarke_1866)         */
/* Need special case for SPCS: ...)                                 */
 
const int Proj_param_value_type[PROJ_NPROJ + 2][PROJ_NPARAM] = {
  /* GEO */
  {PROJ_PARAM_MAJOR_AXIS,    PROJ_PARAM_MINOR_AXIS,     PROJ_PARAM_NULL,
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL, 
   PROJ_PARAM_FALSE_EASTING, PROJ_PARAM_FALSE_NORTHING, PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL}, 
  /* UTM */
  {PROJ_PARAM_MAJOR_AXIS,    PROJ_PARAM_MINOR_AXIS,     PROJ_PARAM_NULL,
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL, 
   PROJ_PARAM_FALSE_EASTING, PROJ_PARAM_FALSE_NORTHING, PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL},
  /* SPCS */
  {PROJ_PARAM_MAJOR_AXIS,    PROJ_PARAM_MINOR_AXIS,     PROJ_PARAM_NULL,
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL, 
   PROJ_PARAM_FALSE_EASTING, PROJ_PARAM_FALSE_NORTHING, PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL},
  /* ALBERS */
  {PROJ_PARAM_MAJOR_AXIS,    PROJ_PARAM_MINOR_AXIS,     PROJ_PARAM_LAT1,
   PROJ_PARAM_LAT2,          PROJ_PARAM_CENTER_LONG,    PROJ_PARAM_LAT_ORIGIN, 
   PROJ_PARAM_FALSE_EASTING, PROJ_PARAM_FALSE_NORTHING, PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL},
  /* LAMCC */
  {PROJ_PARAM_MAJOR_AXIS,    PROJ_PARAM_MINOR_AXIS,     PROJ_PARAM_LAT1,
   PROJ_PARAM_LAT2,          PROJ_PARAM_CENTER_LONG,    PROJ_PARAM_LAT_ORIGIN, 
   PROJ_PARAM_FALSE_EASTING, PROJ_PARAM_FALSE_NORTHING, PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL},
  /* MERCAT */
  {PROJ_PARAM_MAJOR_AXIS,    PROJ_PARAM_MINOR_AXIS,     PROJ_PARAM_NULL,
   PROJ_PARAM_NULL,          PROJ_PARAM_CENTER_LONG,    PROJ_PARAM_LAT1,
   PROJ_PARAM_FALSE_EASTING, PROJ_PARAM_FALSE_NORTHING, PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL},
  /* PS */
  {PROJ_PARAM_MAJOR_AXIS,    PROJ_PARAM_MINOR_AXIS,     PROJ_PARAM_NULL,
   PROJ_PARAM_NULL,          PROJ_PARAM_CENTER_LONG,    PROJ_PARAM_LAT1,
   PROJ_PARAM_FALSE_EASTING, PROJ_PARAM_FALSE_NORTHING, PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL},
  /* POLYC */
  {PROJ_PARAM_MAJOR_AXIS,    PROJ_PARAM_MINOR_AXIS,     PROJ_PARAM_NULL,
   PROJ_PARAM_NULL,          PROJ_PARAM_CENTER_LONG,    PROJ_PARAM_LAT_ORIGIN,
   PROJ_PARAM_FALSE_EASTING, PROJ_PARAM_FALSE_NORTHING, PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL},
  /* EQUIDC */
  {PROJ_PARAM_MAJOR_AXIS,    PROJ_PARAM_MINOR_AXIS,     PROJ_PARAM_LAT1,
   PROJ_PARAM_LAT2,          PROJ_PARAM_CENTER_LONG,    PROJ_PARAM_LAT_ORIGIN, 
   PROJ_PARAM_FALSE_EASTING, PROJ_PARAM_FALSE_NORTHING, PROJ_PARAM_MODE, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL},
  /* TM */
  {PROJ_PARAM_MAJOR_AXIS,    PROJ_PARAM_MINOR_AXIS,     PROJ_PARAM_SCALE_FACTOR,
   PROJ_PARAM_NULL,          PROJ_PARAM_CENTER_LONG,    PROJ_PARAM_LAT_ORIGIN,
   PROJ_PARAM_FALSE_EASTING, PROJ_PARAM_FALSE_NORTHING, PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL},
  /* STEREO */
  {PROJ_PARAM_RADIUS,        PROJ_PARAM_NULL,           PROJ_PARAM_NULL,
   PROJ_PARAM_NULL,          PROJ_PARAM_CENTER_LONG,    PROJ_PARAM_CENTER_LAT, 
   PROJ_PARAM_FALSE_EASTING, PROJ_PARAM_FALSE_NORTHING, PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL},
  /* LAMAZ */
  {PROJ_PARAM_RADIUS,        PROJ_PARAM_NULL,           PROJ_PARAM_NULL,
   PROJ_PARAM_NULL,          PROJ_PARAM_CENTER_LONG,    PROJ_PARAM_CENTER_LAT, 
   PROJ_PARAM_FALSE_EASTING, PROJ_PARAM_FALSE_NORTHING, PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL},
  /* AZMEQD */
  {PROJ_PARAM_RADIUS,        PROJ_PARAM_NULL,           PROJ_PARAM_NULL,
   PROJ_PARAM_NULL,          PROJ_PARAM_CENTER_LONG,    PROJ_PARAM_CENTER_LAT, 
   PROJ_PARAM_FALSE_EASTING, PROJ_PARAM_FALSE_NORTHING, PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL},
  /* GNOMON */
  {PROJ_PARAM_RADIUS,        PROJ_PARAM_NULL,           PROJ_PARAM_NULL,
   PROJ_PARAM_NULL,          PROJ_PARAM_CENTER_LONG,    PROJ_PARAM_CENTER_LAT, 
   PROJ_PARAM_FALSE_EASTING, PROJ_PARAM_FALSE_NORTHING, PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL},
  /* ORTHO */
  {PROJ_PARAM_RADIUS,        PROJ_PARAM_NULL,           PROJ_PARAM_NULL,
   PROJ_PARAM_NULL,          PROJ_PARAM_CENTER_LONG,    PROJ_PARAM_CENTER_LAT, 
   PROJ_PARAM_FALSE_EASTING, PROJ_PARAM_FALSE_NORTHING, PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL},
  /* GVNSP */
  {PROJ_PARAM_RADIUS,        PROJ_PARAM_NULL,           PROJ_PARAM_H,
   PROJ_PARAM_NULL,          PROJ_PARAM_CENTER_LONG,    PROJ_PARAM_CENTER_LAT, 
   PROJ_PARAM_FALSE_EASTING, PROJ_PARAM_FALSE_NORTHING, PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL},
  /* SNSOID */
  {PROJ_PARAM_RADIUS,        PROJ_PARAM_NULL,           PROJ_PARAM_NULL,
   PROJ_PARAM_NULL,          PROJ_PARAM_CENTER_LONG,    PROJ_PARAM_NULL, 
   PROJ_PARAM_FALSE_EASTING, PROJ_PARAM_FALSE_NORTHING, PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL},
  /* EQRECT */
  {PROJ_PARAM_RADIUS,        PROJ_PARAM_NULL,           PROJ_PARAM_NULL,
   PROJ_PARAM_NULL,          PROJ_PARAM_CENTER_LONG,    PROJ_PARAM_LAT1, 
   PROJ_PARAM_FALSE_EASTING, PROJ_PARAM_FALSE_NORTHING, PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL},
  /* MILLER */
  {PROJ_PARAM_RADIUS,        PROJ_PARAM_NULL,           PROJ_PARAM_NULL,
   PROJ_PARAM_NULL,          PROJ_PARAM_CENTER_LONG,    PROJ_PARAM_NULL, 
   PROJ_PARAM_FALSE_EASTING, PROJ_PARAM_FALSE_NORTHING, PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL},
  /* VGRINT */
  {PROJ_PARAM_RADIUS,        PROJ_PARAM_NULL,           PROJ_PARAM_NULL,
   PROJ_PARAM_NULL,          PROJ_PARAM_CENTER_LONG,    PROJ_PARAM_NULL, 
   PROJ_PARAM_FALSE_EASTING, PROJ_PARAM_FALSE_NORTHING, PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL},
  /* HOM */
  {PROJ_PARAM_MAJOR_AXIS,    PROJ_PARAM_MINOR_AXIS,     PROJ_PARAM_SCALE_FACTOR,
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_LAT_ORIGIN, 
   PROJ_PARAM_FALSE_EASTING, PROJ_PARAM_FALSE_NORTHING, PROJ_PARAM_LON1, 
   PROJ_PARAM_LAT1,          PROJ_PARAM_LON2,           PROJ_PARAM_LAT2, 
   PROJ_PARAM_MODE,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL},
  /* ROBIN */
  {PROJ_PARAM_RADIUS,        PROJ_PARAM_NULL,           PROJ_PARAM_NULL,
   PROJ_PARAM_NULL,          PROJ_PARAM_CENTER_LONG,    PROJ_PARAM_NULL, 
   PROJ_PARAM_FALSE_EASTING, PROJ_PARAM_FALSE_NORTHING, PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL},
  /* SOM */
  {PROJ_PARAM_MAJOR_AXIS,    PROJ_PARAM_MINOR_AXIS,     PROJ_PARAM_SATNUM,
   PROJ_PARAM_PATH,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL, 
   PROJ_PARAM_FALSE_EASTING, PROJ_PARAM_FALSE_NORTHING, PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL, 
   PROJ_PARAM_MODE,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL},
  /* ALASKA */
  {PROJ_PARAM_MAJOR_AXIS,    PROJ_PARAM_MINOR_AXIS,     PROJ_PARAM_NULL,
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL, 
   PROJ_PARAM_FALSE_EASTING, PROJ_PARAM_FALSE_NORTHING, PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL},
  /* GOODE */
  {PROJ_PARAM_RADIUS,        PROJ_PARAM_NULL,           PROJ_PARAM_NULL,
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL, 
   PROJ_PARAM_FALSE_EASTING, PROJ_PARAM_FALSE_NORTHING, PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL},
  /* MOLL */
  {PROJ_PARAM_RADIUS,        PROJ_PARAM_NULL,           PROJ_PARAM_NULL,
   PROJ_PARAM_NULL,          PROJ_PARAM_CENTER_LONG,    PROJ_PARAM_NULL, 
   PROJ_PARAM_FALSE_EASTING, PROJ_PARAM_FALSE_NORTHING, PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL},
  /* IMOLL */
  {PROJ_PARAM_RADIUS,        PROJ_PARAM_NULL,           PROJ_PARAM_NULL,
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL, 
   PROJ_PARAM_FALSE_EASTING, PROJ_PARAM_FALSE_NORTHING, PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL},
  /* HAMMER */
  {PROJ_PARAM_RADIUS,        PROJ_PARAM_NULL,           PROJ_PARAM_NULL,
   PROJ_PARAM_NULL,          PROJ_PARAM_CENTER_LONG,    PROJ_PARAM_NULL, 
   PROJ_PARAM_FALSE_EASTING, PROJ_PARAM_FALSE_NORTHING, PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL},
  /* WAGIV */
  {PROJ_PARAM_RADIUS,        PROJ_PARAM_NULL,           PROJ_PARAM_NULL,
   PROJ_PARAM_NULL,          PROJ_PARAM_CENTER_LONG,    PROJ_PARAM_NULL, 
   PROJ_PARAM_FALSE_EASTING, PROJ_PARAM_FALSE_NORTHING, PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL},
  /* WAGVII */
  {PROJ_PARAM_RADIUS,        PROJ_PARAM_NULL,           PROJ_PARAM_NULL,
   PROJ_PARAM_NULL,          PROJ_PARAM_CENTER_LONG,    PROJ_PARAM_NULL, 
   PROJ_PARAM_FALSE_EASTING, PROJ_PARAM_FALSE_NORTHING, PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL},
  /* OBEQA */
  {PROJ_PARAM_RADIUS,        PROJ_PARAM_NULL,           PROJ_PARAM_SHAPE_M,
   PROJ_PARAM_SHAPE_N,       PROJ_PARAM_CENTER_LONG,    PROJ_PARAM_CENTER_LAT, 
   PROJ_PARAM_FALSE_EASTING, PROJ_PARAM_FALSE_NORTHING, PROJ_PARAM_ANGLE, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL},
  /* ISINUS */
  {PROJ_PARAM_RADIUS,        PROJ_PARAM_NULL,           PROJ_PARAM_NULL,
   PROJ_PARAM_NULL,          PROJ_PARAM_CENTER_LONG,    PROJ_PARAM_NULL, 
   PROJ_PARAM_FALSE_EASTING, PROJ_PARAM_FALSE_NORTHING, PROJ_PARAM_NZONE, 
   PROJ_PARAM_NULL,          PROJ_PARAM_JUSTIFY,        PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL},
  /* ISINUS_OLD */
/*
  {PROJ_PARAM_RADIUS,        PROJ_PARAM_NULL,           PROJ_PARAM_NULL,
   PROJ_PARAM_NULL,          PROJ_PARAM_CENTER_LONG,    PROJ_PARAM_NULL, 
   PROJ_PARAM_FALSE_EASTING, PROJ_PARAM_FALSE_NORTHING, PROJ_PARAM_NZONE, 
   PROJ_PARAM_NULL,          PROJ_PARAM_JUSTIFY,        PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL},
*/
  /* HOM_ALT */
  {PROJ_PARAM_MAJOR_AXIS,    PROJ_PARAM_MINOR_AXIS,     PROJ_PARAM_SCALE_FACTOR,
   PROJ_PARAM_AZIMUTH,       PROJ_PARAM_LON_ORIGIN,     PROJ_PARAM_LAT_ORIGIN, 
   PROJ_PARAM_FALSE_EASTING, PROJ_PARAM_FALSE_NORTHING, PROJ_PARAM_NULL, 
   PROJ_PARAM_NULL,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL, 
   PROJ_PARAM_MODE,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL},
  /* SOM_ALT */
  {PROJ_PARAM_MAJOR_AXIS,    PROJ_PARAM_MINOR_AXIS,     PROJ_PARAM_SATNUM,
   PROJ_PARAM_ALF,           PROJ_PARAM_LON1,           PROJ_PARAM_NULL, 
   PROJ_PARAM_FALSE_EASTING, PROJ_PARAM_FALSE_NORTHING, PROJ_PARAM_TIME, 
   PROJ_PARAM_SAT_RATIO,     PROJ_PARAM_START,          PROJ_PARAM_NULL, 
   PROJ_PARAM_MODE,          PROJ_PARAM_NULL,           PROJ_PARAM_NULL}};

#endif
