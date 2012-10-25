/*
!C****************************************************************************

!File: space.c
  
!Description: Functions for mapping from an image space to latitude and
 longitude and visa versa.

!Revision History:
 Revision 1.0 2001/05/08
 Robert Wolfe
 Original Version.

 Revision 1.1 2002/03/02
 Robert Wolfe
 Added special handling for input ISIN case.

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
      phone: 301-614-5508               Lanham, MD 20770  
  
 ! Design Notes:
   1. The following functions handle the image space mapping:

       SetupSpace - Setup the space mapping.
       ToSpace - Map from geodetic to image coordinates.
       FromSpace - Map from image to geodetic coordinates. 
       FreeSpace - Frees space data structure memory.

   2. Geodetic coordinates are the geodetic latitude and longitude of the 
      point to be mapped. Geodetic cooridnates are in radians.
   3. Only one space mapping definition (instance) should be used at any 
      one time because the GCTP package is not re-entrant.
   4. Image coordinates are in pixels with the origin (0.0, 0.0) at the 
      center of the upper left corner pixel.  Sample coordinates are positive 
      to the right and line coordinates are positive downward.
   5. 'SetupSpace' must be called before any of the other routines.  'FreeSpace'
      should be to free the 'space' data structure.
   5. GCTP stands for the General Cartographic Transformation Package.

!END****************************************************************************
*/

#include <stdlib.h>
#include <math.h>
#include "space.h"
#include "hdf.h"
#include "mfhdf.h"
#include "HdfEosDef.h"
#include "proj.h"
#include "mystring.h"
#include "myproj.h"
#include "myproj_const.h"
#include "const.h"
#include "bool.h"
#include "myhdf.h"
#include "error.h"

/* Constants */

#define MAX_PROJ (99)  /* Maximum map projection number */
#define GCTP_OK (0)    /* Okay status return from the GCTP package */

/* Prototypes for initializing the GCTP projections */

void for_init(long proj_num, long zone, double *proj_param, long sphere,
              char *file27, char *file83, long *iflag, 
	      long (*for_trans[MAX_PROJ + 1])());
void inv_init(long proj_num, long zone, double *proj_param, long sphere,
              char *file27, char *file83, long *iflag, 
	      long (*inv_trans[MAX_PROJ + 1])());

/* Functions */

Space_t *SetupSpace(Space_def_t *space_def)
/* 
!C******************************************************************************

!Description: 'SetupSpace' sets up the 'space' data structure and initializes 
 the forward and inverse mapping.
 
!Input Parameters:
 space_def      space definition; the following fields are input:
                  pixel_size, ul_corner.x, ul_corner.y, img_size.l,
		  img_size.s, proj_num, zone, sphere, proj_param[*], isin_type

!Output Parameters:
 (returns)      'space' data structure or NULL when error occurs

!Team Unique Header:

 ! Design Notes:
   1. An error status is returned when:
       a. either image size dimension is zero or negative
       b. the pixel size is zero or negative
       c. the project number is less than zero or greater than 'MAX_PROJ'
       d. there is an error allocating memory
       e. an error occurs initializing either the forward or inverse map 
          projection.
   2. Error messages are handled with the 'RETURN_ERROR' macro.
   3. 'FreeSpace' should be called to deallocate the space structure.

!END****************************************************************************
*/
{
  Space_t *this;
  char file27[28] = "FILE27";
  char file83[28] = "FILE83";
  long (*for_trans[MAX_PROJ + 1])();
  long (*inv_trans[MAX_PROJ + 1])();
  int ip;
  long iflag;

  /* Verify some of the space definition parameters */
  
  if (space_def->img_size.l < 1) 
    RETURN_ERROR("invalid number of lines", "SetupSpace",  
                 (Space_t *)NULL);
  if (space_def->img_size.s < 1) 
    RETURN_ERROR("invalid number of samples per line", "SetupSpace",  
                 (Space_t *)NULL);
  if (space_def->pixel_size <= 0.0)
    RETURN_ERROR("invalid pixel size", "SetupSpace",  
                 (Space_t *)NULL);
  if (space_def->proj_num < 0  ||  space_def->proj_num > MAX_PROJ)
    RETURN_ERROR("invalid projection number", "SetupSpace",  
                 (Space_t *)NULL);

  /* Create the space data structure */

  this = (Space_t *)malloc(sizeof(Space_t));
  if (this == (Space_t *)NULL) 
    RETURN_ERROR("allocating Space structure", "SetupSpace", 
                 (Space_t *)NULL);

  /* Copy the space definition */

  this->def.pixel_size = space_def->pixel_size;
  this->def.ul_corner.x = space_def->ul_corner.x;
  this->def.ul_corner.y = space_def->ul_corner.y;
  this->def.img_size.l = space_def->img_size.l;
  this->def.img_size.s = space_def->img_size.s;
  this->def.proj_num = space_def->proj_num;
  this->def.zone = space_def->zone;
  this->def.sphere = space_def->sphere;
  this->def.isin_type = space_def->isin_type;
  this->def.orientation_angle = space_def->orientation_angle;
  for (ip = 0; ip < NPROJ_PARAM; ip++) 
    this->def.proj_param[ip] = space_def->proj_param[ip];

  /* Calculate the orientation cosine/sine */

  this->sin_orien = sin(space_def->orientation_angle);
  this->cos_orien = cos(space_def->orientation_angle);

  /* Setup the forward transform */

  for_init(this->def.proj_num, this->def.zone, this->def.proj_param, 
           this->def.sphere, file27, file83, &iflag, for_trans);
  if (iflag) {
    free(this);
    RETURN_ERROR("bad return from for_init", "SetupSpace", (Space_t *)NULL);
  }
  this->for_trans = for_trans[this->def.proj_num];

  /* Setup the inverse transform */

  inv_init(this->def.proj_num, this->def.zone, this->def.proj_param, 
           this->def.sphere, file27, file83, &iflag, inv_trans);
  if (iflag) {
    free(this);
    RETURN_ERROR("bad return from inv_init", "SetupSpace", (Space_t *)NULL);
  }
  this->inv_trans = inv_trans[this->def.proj_num];

  return this;
}


bool ToSpace(Space_t *this, Geo_coord_t *geo, Img_coord_float_t *img)
/* 
!C******************************************************************************

!Description: 'ToSpace' maps a point from geodetic to image coordinates.
 
!Input Parameters:
 this           'space' data structure; the following fields are input:
                   for_trans
 geo            geodetic coordinates (radians)

!Output Parameters:
 img            image space coordinates
 (returns)      status:
                  'true' = okay
		  'false' = error return

!Team Unique Header:

 ! Design Notes:
   1. 'SetupSpace' must be called before this routine is called.
   2. An error status is returned when:
       a. a fill value ('geo->is_fill' = 'true') is given
       b. an error occurs in the forward transformation from 
          coordinates to map projection coordinates.
   3. Error messages are handled with the 'RETURN_ERROR' macro.
   4. The image space coordinate is set to fill ('img->is_fill' = 'true') if 
      an error occurs.

!END****************************************************************************
*/
{
  Map_coord_t map;
  double dx, dy, dl, ds;

  img->is_fill = true;
  if (geo->is_fill)
    RETURN_ERROR("called with fill value", "ToSpace", false);

  if (this->for_trans(geo->lon, geo->lat, &map.x, &map.y) != GCTP_OK) 
    RETURN_ERROR("forward transform", "ToSpace", false);

  dx = map.x - this->def.ul_corner.x;
  dy = map.y - this->def.ul_corner.y;

  dl = (dx * this->sin_orien) - (dy * this->cos_orien);
  ds = (dx * this->cos_orien) + (dy * this->sin_orien);

  img->l = (dl / this->def.pixel_size) - 0.5;
  img->s = (ds / this->def.pixel_size) - 0.5;
  img->is_fill = false;

  return true;
}


bool FromSpace(Space_t *this, Img_coord_float_t *img, Geo_coord_t *geo)
/* 
!C******************************************************************************

!Description: 'FromSpace' maps a point from image to geodetic coordinates.
 
!Input Parameters:
 this           'space' data structure; the following fields are input:
                   inv_trans
 img            image space coordinates

!Output Parameters:
 geo            geodetic coordinates (radians)
 (returns)      status:
                  'true' = okay
		  'false' = error return

!Team Unique Header:

 ! Design Notes:
   1. 'SetupSpace' must be called before this routine is called.
   2. An error status is returned when:
       a. a fill value ('img->is_fill' = 'true') is given
       b. an error occurs in the inverse transformation from map projection 
          to geodetic coordinates.
   3. Error messages are handled with the 'RETURN_ERROR' macro.
   4. The image space coordinate is set to fill ('geo->is_fill' = 'true') if 
      an error occurs.

!END****************************************************************************
*/
{
  Map_coord_t map;
  double dl, ds, dy, dx;

  geo->is_fill = true;
  if (img->is_fill) 
    RETURN_ERROR("called with fill value", "FromSpace", false);

  dl = (img->l + 0.5) * this->def.pixel_size;
  ds = (img->s + 0.5) * this->def.pixel_size;

  dy = (ds * this->sin_orien) - (dl * this->cos_orien);
  dx = (ds * this->cos_orien) + (dl * this->sin_orien);

  map.y = this->def.ul_corner.y + dy;
  map.x = this->def.ul_corner.x + dx;

  if (this->inv_trans(map.x, map.y, &geo->lon, &geo->lat) != GCTP_OK) 
    RETURN_ERROR("inverse transform", "FromSpace", false);
  geo->is_fill = false;

  return true;
}

bool FreeSpace(Space_t *this)
/* 
!C******************************************************************************

!Description: 'FreeSpace' frees 'space' data structure memory.
 
!Input Parameters:
 this           'space' data structure

!Output Parameters:
 (returns)      status:
                  'true' = okay (always returned)

!Team Unique Header:

 ! Design Notes:
   1. 'SetupSpace' must be called before this routine is called.
   2. An error status is never returned.

!END****************************************************************************
*/
{
  if (this != (Space_t *)NULL) {
    free(this);
  }

  return true;
}


#define NLINE_MAX (20000)
#define NSAMP_MAX (20000)

typedef enum {
  SPACE_NULL = -1,
  SPACE_START = 0,
  SPACE_PROJ_NUM,
  SPACE_PROJ_PARAM,
  SPACE_PIXEL_SIZE,
  SPACE_UL_CORNER,
  SPACE_NSAMPLE,
  SPACE_NLINE,
  SPACE_ZONE,
  SPACE_SPHERE,
  SPACE_ORIEN_ANGLE,
  SPACE_END,
  SPACE_MAX
} Space_key_t;

Key_string_t Space_string[SPACE_MAX] = {
  {(int)SPACE_START,       "HEADER_FILE"},
  {(int)SPACE_PROJ_NUM,    "PROJECTION_NUMBER"},
  {(int)SPACE_PROJ_PARAM,  "PROJECTION_PARAMETERS"},
  {(int)SPACE_PIXEL_SIZE,  "PIXEL_SIZE"},
  {(int)SPACE_UL_CORNER,   "UPPER_LEFT_CORNER"},
  {(int)SPACE_NSAMPLE,     "NSAMPLE"},
  {(int)SPACE_NLINE,       "NLINE"},
  {(int)SPACE_ZONE,        "PROJECTION_ZONE"},
  {(int)SPACE_SPHERE,      "PROJECTION_SPHERE"},
  {(int)SPACE_ORIEN_ANGLE, "ORIENTATION"},
  {(int)SPACE_END,         "END"}
};

bool GetSpaceDefFile(char *file_name, Space_def_t *this) {
  char *error_string = (char *)NULL;
  FILE *fp;
  Key_t key;
  char temp[MAX_STR_LEN + 1];
  bool got_start, got_end;
  int len;
  char line[MAX_STR_LEN + 1];
  Space_key_t space_key;
  int i;
  double dval;

  int ip;
  double r_sphere, pixel_size_nest, z_nest;
  int n_nest;

  /* Populate the data structure */

  this->proj_num = -1;
  for (ip = 0; ip < NPROJ_PARAM; ip++)
    this->proj_param[ip] = 0.0;
  this->pixel_size = -1.0;
  this->ul_corner.x = -1.0;
  this->ul_corner.y = -1.0;
  this->ul_corner_set = false;
  this->img_size.l = -1;
  this->img_size.s = -1;
  this->zone = 0;
  this->zone_set = false;
  this->sphere = -1;
  this->isin_type = SPACE_NOT_ISIN;
  this->orientation_angle = 0.0;

  /* Open the header file */
  
  if ((fp = fopen(file_name, "r")) == (FILE *)NULL) {
    free(this);
    RETURN_ERROR("unable to open header file", "GetSpaeceDefFile", false);
  }

  /* Parse the header file */

  got_start = got_end = false;

  while((len = GetLine(fp, line)) > 0) {
    if (!StringParse(line, &key)) {
      sprintf(temp, "parsing header file; line = %s", line);
      error_string = temp;
      break;
    }
    if (key.len_key <= 0) continue;
    if (key.key[0] == '#') continue;
    space_key = (Space_key_t)
       KeyString(key.key, key.len_key, Space_string, 
		 (int)SPACE_NULL, (int)SPACE_MAX);
    if (space_key == SPACE_NULL) {
      key.key[key.len_key] = '\0';
      /* printf("skipping unknown key; key = %s", key.key); */
      continue;
    }
    if (!got_start) {
      if (space_key == SPACE_START) {
        if (key.nval != 0) {
	  error_string = "no value expected";
	  break;
	}
        got_start = true;
	continue;
      } else {
        error_string  = "no start key in file";
        break;
      }
    }

    /* Get the value for each keyword */

    switch (space_key) {

      case SPACE_PROJ_NUM:
        if (key.nval <= 0) {
	  error_string = "no projection number value";
	  break; 
	} else if (key.nval > 1) {
	  error_string = "too many projection number values";
	  break; 
	}
        if (key.len_value[0] < 1) {
	  error_string = "invalid projection number string";
	  break;
	}
	strncpy(temp, key.value[0], key.len_value[0]);
        temp[key.len_value[0]] = '\0';

	if (sscanf(temp, "%d", &this->proj_num) != 1) {
	  error_string = "converting projection number";
	  break;
	}
	if (this->proj_num < 0  ||  
	    (this->proj_num >= PROJ_NPROJ  &&  
	     this->proj_num != PROJ_ISINUS)) {
	  error_string = "projection number out of range";
	  break;
	}
        break;

      case SPACE_PROJ_PARAM:
        if (key.nval <= 0) {
	  error_string = "no projection parameter values";
	  break; 
	} else if (key.nval > NPROJ_PARAM) {
	  error_string = "too many projection parameter values";
	  break; 
	}
      for (ip = 0; ip < key.nval; ip++) {
        if (key.len_value[ip] < 1) {
          sprintf(temp, "invalid projection parameter string, parameter = %d", ip);
 	    error_string = temp;
	    break;
        }
	  strncpy(temp, key.value[ip], key.len_value[ip]);
          temp[key.len_value[ip]] = '\0';
	  if (sscanf(temp, "%lf", &this->proj_param[ip]) != 1) {
          sprintf(temp, "convering projection parameter values, parameter = %d", ip);
 	    error_string = temp;
	    break;
        }
	}
        break;

      case SPACE_PIXEL_SIZE:
        if (key.nval <= 0) {
	  error_string = "no pixel size value";
	  break; 
	} else if (key.nval > 1) {
	  error_string = "too many pixel size values";
	  break; 
	}
        if (key.len_value[0] < 1) {
	  error_string = "invalid pixel size string";
	  break;
	}
	strncpy(temp, key.value[0], key.len_value[0]);
        temp[key.len_value[0]] = '\0';

	if (sscanf(temp, "%f", &this->pixel_size) != 1) {
	  error_string = "converting pixel size";
	  break;
	}
	if (this->pixel_size < 5.0  ||  
	    this->pixel_size > 100.0) {
	  error_string = "pixel size out of range";
	  break;
	}
        break;

      case SPACE_UL_CORNER:
        if (key.nval <= 0) {
	  error_string = "no upper left corner values";
	  break; 
	} else if (key.nval > 2) {
	  error_string = "too many upper left corner values";
	  break; 
	} else if (key.nval == 1) {
	  error_string = "only one left corner value";
	  break; 
	}
      for (i = 0; i < key.nval; i++) {
        if (key.len_value[i] < 1) {
          sprintf(temp, "invalid upper left corner string, value = %d", i);
 	    error_string = temp;
	    break;
        }
	  strncpy(temp, key.value[i], key.len_value[i]);
          temp[key.len_value[i]] = '\0';
	  if (sscanf(temp, "%lf", &dval) != 1) {
          sprintf(temp, "convering upper left corner values, value = %d", i);
 	    error_string = temp;
	    break;
        }
        if (i == 0) this->ul_corner.x = dval;
        else this->ul_corner.y = dval;
	}
        this->ul_corner_set = true;
        break;

      case SPACE_NSAMPLE:
        if (key.nval <= 0) {
	  error_string = "no number of sample value";
	  break; 
	} else if (key.nval > 1) {
	  error_string = "too many number of sample values";
	  break; 
	}
        if (key.len_value[0] < 1) {
	  error_string = "invalid number of samples string";
	  break;
	}
	strncpy(temp, key.value[0], key.len_value[0]);
        temp[key.len_value[0]] = '\0';

	if (sscanf(temp, "%d", &this->img_size.s) != 1) {
	  error_string = "converting number of samples";
	  break;
	}
	if (this->img_size.s < 1  ||
	    this->img_size.s > NSAMP_MAX)
	  error_string = "sample number out of range";
        break;

      case SPACE_NLINE:
        if (key.nval <= 0) {
	  error_string = "no number of line value";
	  break; 
	} else if (key.nval > 1) {
	  error_string = "too many number of line values";
	  break; 
	}
        if (key.len_value[0] < 1) {
	  error_string = "invalid number of lines string";
	  break;
	}
	strncpy(temp, key.value[0], key.len_value[0]);
        temp[key.len_value[0]] = '\0';

	if (sscanf(temp, "%d", &this->img_size.l) != 1) {
	  error_string = "converting number of lines";
	  break;
	}
	if (this->img_size.l < 1  ||
	    this->img_size.l > NLINE_MAX)
	  error_string = "line number out of range";
        break;

      case SPACE_ZONE:
        if (key.nval <= 0) {
	  error_string = "no zone number value";
	  break; 
	} else if (key.nval > 1) {
	  error_string = "too many zone number values";
	  break; 
	}
        if (key.len_value[0] < 1) {
	  error_string = "invalid zone number string";
	  break;
	}
	strncpy(temp, key.value[0], key.len_value[0]);
        temp[key.len_value[0]] = '\0';

	if (sscanf(temp, "%d", &this->zone) != 1) {
	  error_string = "converting zone number";
	  break;
	}
      this->zone_set = true;
        break;

      case SPACE_END:
        if (key.nval != 0) {
	  error_string = "no value expected";
	  break; 
	}
        got_end = true;
        break;

      case SPACE_SPHERE:
        if (key.nval <= 0) {
	  error_string = "no projection sphere number value";
	  break; 
	} else if (key.nval > 1) {
	  error_string = "too many projection sphere number values";
	  break; 
	}
        if (key.len_value[0] < 1) {
	  error_string = "invalid projection sphere number string";
	  break;
	}
	strncpy(temp, key.value[0], key.len_value[0]);
        temp[key.len_value[0]] = '\0';

	if (sscanf(temp, "%d", &this->sphere) != 1) {
	  error_string = "converting projection sphere number";
	  break;
	}
	if (this->sphere < 0  ||  this->sphere >= PROJ_NSPHERE) {
	  error_string = "projection sphere number out of range";
	  break;
	}
        break;

      case SPACE_ORIEN_ANGLE:
        if (key.nval <= 0) {
	  error_string = "no orientation angle value";
	  break; 
	} else if (key.nval > 1) {
	  error_string = "too many orientation angle values";
	  break; 
	}
        if (key.len_value[0] < 1) {
	  error_string = "invalid orientation angle string";
	  break;
	}
	strncpy(temp, key.value[0], key.len_value[0]);
        temp[key.len_value[0]] = '\0';

	if (sscanf(temp, "%lf", &this->orientation_angle) != 1) {
	  error_string = "converting orientation angle";
	  break;
	}
      if (this->orientation_angle < -360.0  ||
          this->orientation_angle >  360.0) {
        error_string = "orientation angle out of range";
        break;
      }
      this->orientation_angle *= RAD;
        break;

      default:
        error_string = "key not implmented";

    }
    if (error_string != (char *)NULL) break;
    if (got_end) break;

  }

  /* Handle errors */

  if (error_string != (char *)NULL) 
    RETURN_ERROR(error_string, "GetHeaderInput", false);

  /* Make sure both the start and end of file exist */
  
  if (!got_start) 
    error_string = "no start key in header";
  else if (!got_end)
    error_string = "no end key in header";

  if (error_string != (char *)NULL)
    RETURN_ERROR(error_string, "GetHeaderInput", false);

  /* Check output space definition */

  if (this->proj_num < 0) {
    sprintf(temp, "projection number not given");
    error_string = temp;
  }

  if (this->pixel_size < 0.0) {
    sprintf(temp, "pixel size not given");
    error_string = temp;
  }

  if (!this->ul_corner_set) {
    sprintf(temp, "upper left corner not given");
    error_string = temp;
  }

  if ((this->img_size.l <= 0)  ||
      (this->img_size.s <= 0)) {
    sprintf(temp, "image size not given");
    error_string = temp;
  }

  if ((this->proj_param[0] <= 0.0)   &&  (this->sphere < 0)) {
    sprintf(temp, "either sphere or first projection "
            "parameter value must be given");
    error_string = temp;
  }

  if ((this->proj_num == PROJ_UTM)  && 
       !this->zone_set) { 
    sprintf(temp, "zone number not given for UTM");
    error_string = temp;
  }

  if (this->proj_num == 99) /* ISINUS => proj_num = 99 */
    this->isin_type= SPACE_ISIN_NEST_1;

  /* Comput nesting level for ISINUS projection */

  if (this->proj_num == PROJ_ISINUS) {
    this->isin_type= SPACE_ISIN_NEST_1;
    if (this->sphere >= 0)
      r_sphere = Proj_sphere[this->sphere].major_axis;
    else if (this->proj_param[0] > 0.0) 
      r_sphere = this->proj_param[0];
    else
      r_sphere = 6371007.181;
    if (this->proj_param[8] <= 0.0) {
      sprintf(temp, "number of input zones invalid in ISINUS"
                      " projection");
      error_string = temp;
    } else {
      pixel_size_nest = PI * r_sphere / this->proj_param[8];
      z_nest = this->pixel_size / pixel_size_nest;
      n_nest = (int)(z_nest + 0.5);
      if (n_nest == 1) 
        this->isin_type = SPACE_ISIN_NEST_1;
      else if (n_nest == 2) 
        this->isin_type = SPACE_ISIN_NEST_2;
      else if (n_nest == 4) 
        this->isin_type = SPACE_ISIN_NEST_4;
      else {
        sprintf(temp, "invalid ISIN nesting\n");
        error_string = temp;
      }
    }
  }

  if (error_string != (char *)NULL)
    RETURN_ERROR(error_string, "GetHeaderInput", false);

  return true;
}


/*

   HDF-EOS Metadata Example (from ncdump)

                :HDFEOSVersion = "HDFEOS_V2.4" ;
                :StructMetadata.0 = "GROUP=SwathStructure\n",
    "END_GROUP=SwathStructure\n",
    "GROUP=GridStructure\n",
    "\tGROUP=GRID_1\n",
    "\t\tGridName=\"MOD_Grid_Pntr_2d\"\n",
    "\t\tXDim=1200\n",
    "\t\tYDim=1200\n",
    "\t\tUpperLeftPointMtrs=(-8895604.158132,5559752.598833)\n",
    "\t\tLowerRightMtrs=(-7783653.638366,4447802.079066)\n",
    "\t\tProjection=GCTP_ISINUS\n",
    "\t\tProjParams=(6371007.181000,0,0,0,0,0,0,0,86400,0,1,0,0)\n",
    "\t\tSphereCode=-1\n",
    "\t\tGridOrigin=HDFE_GD_UL\n",
    "\t\tGROUP=Dimension\n",
    "\t\tEND_GROUP=Dimension\n",
    "\t\tGROUP=DataField\n",
    "\t\t\tOBJECT=DataField_1\n",
    "\t\t\t\tDataFieldName=\"num_observations\"\n",
    "\t\t\t\tDataType=DFNT_INT8\n",
    "\t\t\t\tDimList=(\"YDim\",\"XDim\")\n",
    "\t\t\tEND_OBJECT=DataField_1\n",
    "\t\t\tOBJECT=DataField_2\n",
    "\t\t\t\tDataFieldName=\"obscov_1\"\n",
    "\t\t\t\tDataType=DFNT_INT8\n",
    "\t\t\t\tDimList=(\"YDim\",\"XDim\")\n",
    "\t\t\tEND_OBJECT=DataField_2\n",
    "\t\tEND_GROUP=DataField\n",
    "\t\tGROUP=MergedFields\n",
    "\t\tEND_GROUP=MergedFields\n",
    "\tEND_GROUP=GRID_1\n",
    "END_GROUP=GridStructure\n",
    "GROUP=PointStructure\n",
    "END_GROUP=PointStructure\n",
    "END\n",
    "" ;
    
***/

#define SPACE_NTYPE_HDF (10)

struct {
  int32 type;		/* type values */
  char *name;		/* type name */
} space_hdf_type[SPACE_NTYPE_HDF] = {
    {DFNT_CHAR8,   "DFNT_CHAR8"},
    {DFNT_UCHAR8,  "DFNT_UCHAR8"},
    {DFNT_INT8,    "DFNT_INT8"},
    {DFNT_UINT8,   "DFNT_UINT8"},
    {DFNT_INT16,   "DFNT_INT16"},
    {DFNT_UINT16,  "DFNT_UINT16"},
    {DFNT_INT32,   "DFNT_INT32"},
    {DFNT_UINT32,  "DFNT_UINT32"},
    {DFNT_FLOAT32, "DFNT_FLOAT32"},
    {DFNT_FLOAT64, "DFNT_FLOAT64"}
};

#define SPACE_HDFEOS_VERSION ("HDFEOSVersion");
#define SPACE_STRUCT_METADATA ("StructMetadata.0");
#define SPACE_ORIENTATION_ANGLE_HDF ("OrientationAngle")
#define SPACE_PIXEL_SIZE_HDF ("PixelSize")
#define NPROJ_PARAM_HDFEOS (13)

bool AppendMeta(char *cbuf, int *ic, char *s) {
  int nc, i;

  if (ic < 0) return false;
  nc = strlen(s);
  if (nc <= 0) return false;
  if (*ic + nc > MYHDF_MAX_NATTR_VAL) return false;

  for (i = 0; i < nc; i++) {
    cbuf[*ic] = s[i];
    (*ic)++;
  }

  cbuf[*ic] = '\0';

  return true;
}

bool PutSpaceDefHDF(Space_def_t *this, char *file_name, int nsds, 
                    char **sds_names, int *sds_types, char *grid_name) {
  int32 sds_file_id;
  char struct_meta[MYHDF_MAX_NATTR_VAL];
  char cbuf[MYHDF_MAX_NATTR_VAL];
  char *hdfeos_version = "HDFEOS_V2.4";
  char *dim_names[2] = {"YDim", "XDim"};
  int ic;
  Map_coord_t lr_corner;
  double ds, dl, dy, dx;
  double sin_orien, cos_orien;
  Myhdf_attr_t attr;
  int isds;
  int ip;
  double f_fractional, f_integral;
  char *cproj;
  double dval[1];
  char *ctype;
  int it;
  int32 hdf_id;
  int32 vgroup_id[3];
  int32 sds_index, sds_id;

  /* Check inputs */

  if (nsds <= 0) 
    RETURN_ERROR("invalid number of SDSs", "PutSpaceDefHDF", false);
  for (isds = 0; isds < nsds; isds++)
    if (strlen(sds_names[isds]) < 1)  
      RETURN_ERROR("invalid SDS name", "PutSpaceDefHDF", false);
  if (strlen(grid_name) < 1)  
    RETURN_ERROR("invalid grid name [0]", "PutSpaceDefHDF", false);

  /* Put header */

  ic = 0;
  sprintf(cbuf, 
    "GROUP=SwathStructure\n" 
    "END_GROUP=SwathStructure\n" 
    "GROUP=GridStructure\n" 
    "\tGROUP=GRID_1\n");

  if (!AppendMeta(struct_meta, &ic, cbuf))
    RETURN_ERROR("appending to metadata string (head)", "PutSpaceDefHDF", 
                  false);

  /* Compute lower right corner */

  dl = this->img_size.l * this->pixel_size;
  ds = this->img_size.s * this->pixel_size;

  sin_orien = sin(this->orientation_angle);
  cos_orien = cos(this->orientation_angle);

  dy = (ds * sin_orien) - (dl * cos_orien);
  dx = (ds * cos_orien) + (dl * sin_orien);

  lr_corner.y = this->ul_corner.y + dy;
  lr_corner.x = this->ul_corner.x + dx;


  /* Get the projection name string */

  cproj = (char *)NULL;

  for (ip = 0; ip < PROJ_NPROJ; ip++) {
    if (this->proj_num == ip) { 
      cproj = Proj_type[ip].short_name;
      break;
    }
  }
  if (cproj == (char *)NULL)
    RETURN_ERROR("getting projection name string", "PutSpaceDefHDF", false);

  /* Put Grid information */

  sprintf(cbuf, 
    "\t\tGridName=\"%s\"\n" 
    "\t\tXDim=%d\n" 
    "\t\tYDim=%d\n" 
    "\t\tUpperLeftPointMtrs=(%.6f,%.6f)\n" 
    "\t\tLowerRightMtrs=(%.6f,%.6f)\n" 
    "\t\tProjection=GCTP_%s\n", 
    grid_name, 
    this->img_size.s, this->img_size.l, 
    this->ul_corner.x, this->ul_corner.y,
    lr_corner.x, lr_corner.y, 
    cproj);

  if (!AppendMeta(struct_meta, &ic, cbuf))
    RETURN_ERROR("appending to metadata string (grid information start)", 
                  "PutSpaceDefHDF", false);

  if (this->proj_num == PROJ_UTM  ||  this->proj_num == PROJ_SPCS) {

    sprintf(cbuf, "\t\tZoneCode=%d\n", this->zone);

    if (!AppendMeta(struct_meta, &ic, cbuf))
      RETURN_ERROR("appending to metadata string "
                   "(zone number)", 
                    "PutSpaceDefHDF", false);

  } else {

    sprintf(cbuf, "\t\tProjParams=(");
  
    if (!AppendMeta(struct_meta, &ic, cbuf))
      RETURN_ERROR("appending to metadata string "
                   "(grid projection parameters start)", 
                   "PutSpaceDefHDF", false);

    for (ip = 0; ip < NPROJ_PARAM_HDFEOS; ip++) {
      f_fractional = modf(this->proj_param[ip], &f_integral);
      if (fabs(f_fractional) < 0.5e-6) {
        if (ip < (NPROJ_PARAM_HDFEOS - 1)) 
          sprintf(cbuf, "%g,", this->proj_param[ip]);
        else 
          sprintf(cbuf, "%g)", this->proj_param[ip]);
      } else {
        if (ip < (NPROJ_PARAM_HDFEOS + 1)) 
          sprintf(cbuf, "%.6f,", this->proj_param[ip]);
        else 
          sprintf(cbuf, "%.6f)", this->proj_param[ip]);
      }

      if (!AppendMeta(struct_meta, &ic, cbuf))
        RETURN_ERROR("appending to metadata string "
                     "(individual grid projection parameters)", 
                      "PutSpaceDefHDF", false);
    }

    sprintf(cbuf, "\n");

    if (!AppendMeta(struct_meta, &ic, cbuf))
      RETURN_ERROR("appending to metadata string "
                   "(grid projection parameters end)", 
                    "PutSpaceDefHDF", false);

  }


  sprintf(cbuf, 
    "\t\tSphereCode=%d\n" 
    "\t\tGridOrigin=HDFE_GD_UL\n",
    this->sphere);

  if (!AppendMeta(struct_meta, &ic, cbuf))
    RETURN_ERROR("appending to metadata string (grid information end)", 
                  "PutSpaceDefHDF", false);

  /* Put SDS group */

  sprintf(cbuf, 
    "\t\tGROUP=Dimension\n" 
    "\t\tEND_GROUP=Dimension\n"
    "\t\tGROUP=DataField\n");

  if (!AppendMeta(struct_meta, &ic, cbuf))
    RETURN_ERROR("appending to metadata string (SDS group start)", 
                 "PutSpaceDefHDF", false);

  for (isds = 0; isds < nsds; isds++) {

    /* Get the hdf type name string */

    ctype = (char *)NULL;

    for (it = 0; it < SPACE_NTYPE_HDF; it++) {
      if (sds_types[isds] == space_hdf_type[it].type) { 
        ctype = space_hdf_type[it].name;
        break;
      }
    }
    if (ctype == (char *)NULL)
      RETURN_ERROR("getting hdf type name string", "PutSpaceDefHDF", false);

    sprintf(cbuf, 
      "\t\t\tOBJECT=DataField_%d\n"
      "\t\t\t\tDataFieldName=\"%s\"\n"
      "\t\t\t\tDataType=%s\n"
      "\t\t\t\tDimList=(\"%s\",\"%s\")\n"
      "\t\t\tEND_OBJECT=DataField_%d\n",
      isds, sds_names[isds], ctype, dim_names[0], dim_names[1], isds);

    if (!AppendMeta(struct_meta, &ic, cbuf))
      RETURN_ERROR("appending to metadata string (SDS group)", 
                   "PutSpaceDefHDF", false);
  }

  sprintf(cbuf, 
    "\t\tEND_GROUP=DataField\n" 
    "\t\tGROUP=MergedFields\n" 
    "\t\tEND_GROUP=MergedFields\n");

  if (!AppendMeta(struct_meta, &ic, cbuf))
    RETURN_ERROR("appending to metadata string (SDS group end)", 
                 "PutSpaceDefHDF", false);

  /* Put trailer */

  sprintf(cbuf, 
    "\tEND_GROUP=GRID_1\n"
    "END_GROUP=GridStructure\n"
    "GROUP=PointStructure\n"
    "END_GROUP=PointStructure\n"
    "END\n");

  if (!AppendMeta(struct_meta, &ic, cbuf))
    RETURN_ERROR("appending to metadata string (tail)", "PutSpaceDefHDF", 
                  false);

  /* Write file attributes */

  sds_file_id = SDstart((char *)file_name, DFACC_RDWR);
  if (sds_file_id == HDF_ERROR) 
    RETURN_ERROR("opening output file for SD access (1)", "PutSpaceDefHDF", 
                 false); 
  
  attr.type = DFNT_FLOAT64;
  attr.nval = 1;
  attr.name = SPACE_ORIENTATION_ANGLE_HDF;
  dval[0] = (double)this->orientation_angle * DEG;
  if (!PutAttrDouble(sds_file_id, &attr, dval))
    RETURN_ERROR("writing attribute (orientation angle)", 
                 "PutMetadata", false);

  attr.type = DFNT_FLOAT64;
  attr.nval = 1;
  attr.name = SPACE_PIXEL_SIZE_HDF;
  dval[0] = (double)this->pixel_size;
  if (!PutAttrDouble(sds_file_id, &attr, dval))
    RETURN_ERROR("writing attribute (pixel_size)", 
                 "PutMetadata", false);

  attr.type = DFNT_CHAR8;
  attr.nval = strlen(hdfeos_version);
  attr.name = SPACE_HDFEOS_VERSION;
  if (!PutAttrString(sds_file_id, &attr, hdfeos_version))
    RETURN_ERROR("writing attribute (hdfeos_version)", "PutSpaceDefHDF", false);

  attr.type = DFNT_CHAR8;
  attr.nval = strlen(struct_meta);
  attr.name = SPACE_STRUCT_METADATA;
  if (!PutAttrString(sds_file_id, &attr, struct_meta))
    RETURN_ERROR("writing attribute (struct_meta)", "PutSpaceDefHDF", false);

  attr.type = DFNT_CHAR8;
  attr.nval = strlen(struct_meta);
  attr.name = SPACE_STRUCT_METADATA;
  if (!PutAttrString(sds_file_id, &attr, struct_meta))
    RETURN_ERROR("writing attribute (struct_meta)", "PutSpaceDefHDF", false);

  if (SDend(sds_file_id) == HDF_ERROR) 
    RETURN_ERROR("ending SD access (1)", "PutSpaceDefHDF", false);

  /* Setup the HDF Vgroup */

  hdf_id = Hopen((char *)file_name, DFACC_RDWR, 0);
  if (hdf_id == HDF_ERROR) 
    RETURN_ERROR("opening HDF file (Hopen)", "PutSpaceDefHDF", false);

  /* Start the Vgroup access */

  if (Vstart(hdf_id) == HDF_ERROR) 
    RETURN_ERROR("starting Vgroup access (Vstart)", "PutSpaceDefHDF", false);

  /* Create Root Vgroup for Grid */

  vgroup_id[0] = Vattach(hdf_id, -1, "w");
  if (vgroup_id[0] == HDF_ERROR) 
    RETURN_ERROR("getting Grid Vgroup ID (Vattach)", 
                 "PutSpaceDefHDF", false);
  if (Vsetname(vgroup_id[0], grid_name) == HDF_ERROR) 
    RETURN_ERROR("setting Grid Vgroup name (Vsetname)", 
                 "PutSpaceDefHDF", false);
  if (Vsetclass(vgroup_id[0], "GRID") == HDF_ERROR) 
    RETURN_ERROR("setting Grid Vgroup class (Vsetclass)", 
                 "PutSpaceDefHDF", false);

  /* Create Data Fields Vgroup */

  vgroup_id[1] = Vattach(hdf_id, -1, "w");
  if (vgroup_id[1] == HDF_ERROR) 
    RETURN_ERROR("getting Data Fields Vgroup ID (Vattach)", 
                 "PutSpaceDefHDF", false);
  if (Vsetname(vgroup_id[1], "Data Fields") == HDF_ERROR) 
    RETURN_ERROR("setting Data Fields Vgroup name (Vsetname)", 
                 "PutSpaceDefHDF", false);
  if (Vsetclass(vgroup_id[1], "GRID Vgroup") == HDF_ERROR) 
    RETURN_ERROR("setting Data Fields Vgroup class (Vsetclass)", 
                 "PutSpaceDefHDF", false);
  if (Vinsert(vgroup_id[0], vgroup_id[1]) == HDF_ERROR) 
    RETURN_ERROR("inserting Data Fields Vgroup (Vinsert)", 
                 "PutSpaceDefHDF", false);

  /* Create Attributes Vgroup */

  vgroup_id[2] = Vattach(hdf_id, -1, "w");
  if (vgroup_id[2] == HDF_ERROR) 
    RETURN_ERROR("getting Attributes Vgroup ID (Vattach)", 
                 "PutSpaceDefHDF", false);
  if (Vsetname(vgroup_id[2], "Grid Attributes") == HDF_ERROR) 
    RETURN_ERROR("setting Attributes Vgroup name (Vsetname)", 
                 "PutSpaceDefHDF", false);
  if (Vsetclass(vgroup_id[2], "GRID Vgroup") == HDF_ERROR) 
    RETURN_ERROR("setting Attributes Vgroup class (Vsetclass)", 
                 "PutSpaceDefHDF", false);
  if (Vinsert(vgroup_id[0], vgroup_id[2]) == HDF_ERROR) 
    RETURN_ERROR("inserting Attributes Vgroup (Vinsert)", 
                 "PutSpaceDefHDF", false);

  /* Attach SDSs to Data Fields Vgroup */

  sds_file_id = SDstart((char *)file_name, DFACC_RDWR);
  if (sds_file_id == HDF_ERROR) 
    RETURN_ERROR("opening output file for SD access (2)", 
                 "PutSpaceDefHDF", false);

  for (isds = 0; isds < nsds; isds++) {
    sds_index = SDnametoindex(sds_file_id, sds_names[isds]);
    if (sds_index == HDF_ERROR) 
      RETURN_ERROR("getting SDS index (SDnametoindex)", 
                   "PutSpaceDefHDF", false);
    sds_id = SDselect(sds_file_id, sds_index);
    if (sds_id == HDF_ERROR) 
      RETURN_ERROR("getting SDS ID (SDselect)", "PutSpaceDefHDF", false);
    if (Vaddtagref(vgroup_id[1], DFTAG_NDG, SDidtoref(sds_id)) == 
        HDF_ERROR) 
      RETURN_ERROR("adding reference tag to SDS (Vaddtagref)", 
                   "PutSpaceDefHDF", false);
    if (SDendaccess(sds_id) == HDF_ERROR) 
      RETURN_ERROR("ending access to SDS", "PutSpaceDefHDF", false);
  }
  
  if (SDend(sds_file_id) == HDF_ERROR) 
    RETURN_ERROR("ending SD access (2)", "PutSpaceDefHDF", false);

  /* Detach Vgroups */

  if (Vdetach(vgroup_id[0]) == HDF_ERROR) 
    RETURN_ERROR("detaching from Grid Vgroup (Vdetach)", 
                 "PutSpaceDefHDF", false);
  if (Vdetach(vgroup_id[1]) == HDF_ERROR) 
    RETURN_ERROR("detaching from Data Fields Vgroup (Vdetach)", 
                 "PutSpaceDefHDF", false);
  if (Vdetach(vgroup_id[2]) == HDF_ERROR) 
    RETURN_ERROR("detaching from Attributes Vgroup (Vdetach)", 
                 "PutSpaceDefHDF", false);

  /* Close access */

  if (Vend(hdf_id) == HDF_ERROR) 
    RETURN_ERROR("end Vgroup access (Vend)", "PutSpaceDefHDF", false);
  if (Hclose(hdf_id) == HDF_ERROR) 
    RETURN_ERROR("end HDF access (Hclose)", "PutSpaceDefHDF", false);

  return true;
}

bool GetSpaceDefHDF(Space_def_t *this, char *file_name, char *grid_name) {
  int32 gd_id;
  int32 gd_file_id;
  int32 xdim_size, ydim_size;
  float64 ul_corner[2], lr_corner[2];
  int32 proj_num, zone, sphere;
  float64 proj_param[NPROJ_PARAM_HDFEOS];
  int ip;
  double r_sphere, pixel_size_nest, z_nest;
  int n_nest;
  int32 sds_file_id;
  double dval[MYHDF_MAX_NATTR_VAL];
  Myhdf_attr_t attr;

  /* Populate the data structure */

  this->proj_num = -1;
  for (ip = 0; ip < NPROJ_PARAM; ip++)
    this->proj_param[ip] = 0.0;
  this->pixel_size = -1.0;
  this->ul_corner.x = -1.0;
  this->ul_corner.y = -1.0;
  this->ul_corner_set = false;
  this->img_size.l = -1;
  this->img_size.s = -1;
  this->zone = 0;
  this->zone_set = false;
  this->sphere = -1;
  this->isin_type = SPACE_NOT_ISIN;
  this->orientation_angle = 0.0;

  /* Open the HDF-EOS file. */

  printf("HDF Filename: %s\n",file_name);
  gd_file_id = GDopen((char *)file_name, DFACC_READ);
  if (gd_file_id == HDF_ERROR) 
    RETURN_ERROR("opening HDF-EOS file (GDopen)", "GetSpaceDefHDF", false);

  /* Attach to the grid */
  printf("the grid name is %s\n",grid_name);
  gd_id = GDattach(gd_file_id, grid_name);
  if (gd_id == HDF_ERROR) 
    RETURN_ERROR("attaching to HDF-EOS grid (GDattach)", 
                 "GetSpaceDefHDF", false);

  /* Get Grid information */

  if (GDgridinfo(gd_id, &xdim_size, &ydim_size,
                 ul_corner, lr_corner) == HDF_ERROR)
    RETURN_ERROR("getting HDF-EOS grid information (GDgridinfo)", 
                 "GetSpaceDefHDF", false);

  this->img_size.l = ydim_size;
  this->img_size.s = xdim_size;
  this->ul_corner.x = ul_corner[0];
  this->ul_corner.y = ul_corner[1];
  this->ul_corner_set = true;

  if(GDprojinfo(gd_id, &proj_num, &zone, &sphere, proj_param) == HDF_ERROR)
    RETURN_ERROR("getting HDF-EOS map projection information (GDprojinfo)", 
                 "GetSpaceDefHDF", false);

  this->proj_num = proj_num;
  if (this->proj_num == PROJ_UTM  ||  this->proj_num == PROJ_SPCS) {
    this->zone = zone;
    this->zone_set = true;
  }
  this->sphere = sphere;
  for (ip = 0; ip < NPROJ_PARAM_HDFEOS; ip++)
    this->proj_param[ip] = proj_param[ip];

  /* Comput nesting level for ISINUS projection */

  if (this->proj_num == PROJ_ISINUS) {
    this->isin_type= SPACE_ISIN_NEST_1;
    if (this->sphere >= 0)
      r_sphere = Proj_sphere[this->sphere].major_axis;
    else if (this->proj_param[0] > 0.0) 
      r_sphere = this->proj_param[0];
    else
      r_sphere = 6371007.181;
    if (this->proj_param[8] <= 0.0) {
      RETURN_ERROR("number of input zones invalid in ISINUS projection",
                   "GetSpaceDefHDF", false);
    } else {
      pixel_size_nest = PI * r_sphere / this->proj_param[8];
      z_nest = this->pixel_size / pixel_size_nest;
      n_nest = (int)(z_nest + 0.5);
      if (n_nest == 1) 
        this->isin_type = SPACE_ISIN_NEST_1;
      else if (n_nest == 2) 
        this->isin_type = SPACE_ISIN_NEST_2;
      else if (n_nest == 4) 
        this->isin_type = SPACE_ISIN_NEST_4;
      else 
        RETURN_ERROR("invalid ISIN nesting", "GetSpaceDefHDF", false);
    }
  }

  /* Close the HDF-EOS file. */

  if (GDclose(gd_file_id) == HDF_ERROR) 
    RETURN_ERROR("closing HDF-EOS file (GDclose)", "GetSpaceDefHDF", false);

  /* Read file attributes */

  sds_file_id = SDstart((char *)file_name, DFACC_READ);
  if (sds_file_id == HDF_ERROR) 
    RETURN_ERROR("opening output file for SD access", 
                 "GetSpaceDefHDF", false); 
  
  attr.type = DFNT_FLOAT64;
  attr.nval = 1;
  attr.name = SPACE_ORIENTATION_ANGLE_HDF;
  if (!GetAttrDouble(sds_file_id, &attr, dval))
    RETURN_ERROR("reading attribute (orientation angle)", 
                 "PutMetadata", false);
  if (attr.nval != 1)
    RETURN_ERROR("invalid number of values (orientation angle)",
                 "PutMetadata", false);
  this->orientation_angle = dval[0] * RAD;

  attr.type = DFNT_FLOAT64;
  attr.nval = 1;
  attr.name = SPACE_PIXEL_SIZE_HDF;
  if (!GetAttrDouble(sds_file_id, &attr, dval))
    RETURN_ERROR("reading attribute (pixel_size)", 
                 "PutMetadata", false);
  if (attr.nval != 1)
    RETURN_ERROR("invalid number of values (pixel size)",
                 "PutMetadata", false);
  this->pixel_size = (float)dval[0];

  if (SDend(sds_file_id) == HDF_ERROR) 
    RETURN_ERROR("ending SD access", "GetSpaceDefHDF", false);

  return true;
}
