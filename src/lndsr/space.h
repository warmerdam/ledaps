/*
!C****************************************************************************

!File: space.h

!Description: Header file for space.c - see space.c for more information.

!Revision History:
 Revision 1.0 2001/05/08
 Robert Wolfe
 Original Version.

 Revision 1.1 2001/07/24
 Robert Wolfe
 Changed the number of parameters from 16 to 15 to match the GCTP package.

 Revision 1.2 2002/03/02
 Robert Wolfe
 Added new type (Space_isin_t) for special ISIN input case.

 Revision 1.5 2002/11/27
 Robert Wolfe
 Added support for data that is not oriented to map north.

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
   1. Structures are declared for the 'map_coord', 'geo_coord', 'space_def' and 
      'space' data types.
   2. The number of projection parameters is set at 'NPROJ_PARAM'.
   3. GCTP stands for the General Cartographic Transformation Package.
  
!END****************************************************************************
*/

#ifndef SPACE_H
#define SPACE_H

#include "bool.h"

/* Constants */

#define NPROJ_PARAM (15)

typedef enum {SPACE_NOT_ISIN, SPACE_ISIN_NEST_1, SPACE_ISIN_NEST_2, 
              SPACE_ISIN_NEST_4} Space_isin_t;
#define SPACE_MAX_NEST (4)

/* Integer image coordinates data structure */

#ifndef IMG_COORD_INT_TYPE_DEFINED

#define IMG_COORD_INT_TYPE_DEFINED

typedef struct {
  int l;                /* line number */        
  int s;                /* sample number */
} Img_coord_int_t;

#endif

/* Floating point image coordinates data structure */

typedef struct {
  float l;               /* line number */
  float s;               /* sample number */
  bool is_fill;          /* fill value flag; 
                            'true' = is fill; 'false' = is not fill */
} Img_coord_float_t;

/* Structure to store map projection coordinates */

typedef struct {
  double x;             /* Map projection X coordinate (meters) */
  double y;             /* Map projection Y coordinate (meters) */
  bool is_fill;         /* Flag to indicate whether the point is a fill value;
                           'true' = fill; 'false' = not fill */
} Map_coord_t;


/* Structure to store geodetic coordinates */

typedef struct {
  double lon;           /* Geodetic longitude coordinate (radians) */ 
  double lat;           /* Geodetic latitude coordinate (radians) */ 
  bool is_fill;         /* Flag to indicate whether the point is a fill value;
                           'true' = fill; 'false' = not fill */
} Geo_coord_t;

/* Structure to store the space definition */
typedef struct {
  int proj_num;         /* GCTP map projection number */
  double proj_param[NPROJ_PARAM]; /* GCTP map projection parameters */
  float pixel_size;     /* Pixel size (meters) */
  Map_coord_t ul_corner;  /* Map projection coordinates of the upper left 
                             corner of the pixel in the upper left corner 
			     of the image */
  bool ul_corner_set;   /* Flag to indicate whether the upper left corner
                           has been set; 'true' = set; 'false' = not set */
  Img_coord_int_t img_size;  /* Image size (lines, samples) */
  int zone;             /* GCTP zone number */
  int sphere;           /* GCTP sphere number */
  bool zone_set;        /* Flag to indicate whether the zone has been set;
                           'true' = set; 'false' = not set */
  Space_isin_t isin_type;  /* Flag to indicate whether the projection is ISIN,
                              and if it is, the ISIN nesting */
  double orientation_angle;  /* Orientation of the image with respect to 
                                map north (radians).  A positive angle causes 
				the image to be rotated counter-clockwise. */
} Space_def_t;

/* Structure to store the space information */

typedef struct {
  Space_def_t def;       /* Space definition structure */
  long (*for_trans)(double lat, double lon, double *x, double *y);
                         /* Forward transformation function call */
  long (*inv_trans)(double x, double y, double *lat, double *lon);
                         /* Inverse transformation function call */
  double cos_orien;      /* Cosine of the orientation angle */
  double sin_orien;      /* Sine of the orientation angle */
} Space_t;

/* Prototypes */

Space_t *SetupSpace(Space_def_t *space_def);
bool ToSpace(Space_t *this, Geo_coord_t *geo, Img_coord_float_t *map);
bool FromSpace(Space_t *this, Img_coord_float_t *map, Geo_coord_t *geo);
bool FreeSpace(Space_t *this);
bool GetSpaceDefFile(char *file_name, Space_def_t *this);
bool PutSpaceDefHDF(Space_def_t *this, char *file_name, int nsds, 
                    char **sds_names, int *sds_types, char *grid_name);
bool GetSpaceDefHDF(Space_def_t *this, char *file_name, char *grid_name);
 
#endif
