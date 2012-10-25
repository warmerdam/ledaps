/*
!C****************************************************************************

!File: lut.h

!Description: Header file for 'input.c' - see 'input.c' for more information.

!Revision History:
 Revision 1.0 2001/05/08
 Robert Wolfe
 Original Version.

!Team Unique Header:
  This software was developed by the MODIS Land Science Team Support 
  Group for the Laboratory for Terrestrial Physics (Code 922) at the 
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
   1. Structure is declared for the 'input' data type.
  
!END****************************************************************************
*/

#ifndef LUT_H
#define LUT_H

#include "lndsr.h"
#include "input.h"
#include "bool.h"

#define NBAND_SR_LUT (3)

/* Structure for the 'lut' data type */

typedef struct {
  char *file_name;         /* Lookup table file name */
  int nband;               /* Number of bands */
  int in_fill;             /* Input fill value */
  int output_fill;         /* Output fill value */
  int aerosol_fill;        /* Aerosol fill value */
  int out_satu;            /* Output saturation value (Feng, 3/23/09) */
  Img_coord_int_t ar_region_size;  
                           /* Size of the aerosol retreval region */
  Img_coord_int_t ar_size;  
                           /* Size of the aerosol retreval image */
  int min_valid_sr;        /* Minimum valid surface reflectance */
  int max_valid_sr;        /* Maximum valid surface reflectance */
  Input_meta_t meta;       /* Input metadata */
  int cloud_fill;          /* cloud_mask_fill_value */
  int cloud_land;          /* cloud_mask_water_value */
  int cloud_water;         /* cloud_mask_water_value */
  int cloud_shadow;        /* cloud_mask_shadow_value */
  int cloud_snow;          /* cloud_mask_snow_value */
  int cloud_cloud;         /* cloud_mask_cloud_value */
  char* long_name_prefix;  /* long name prefix (append band num) */
  char* units;             /* units */
  double scale_factor;     /* scale factor */
  double scale_factor_err; /* scale factor error */
  double add_offset;       /* add offset */
  double add_offset_err;   /* add offset error */
  double calibrated_nt;    /* calibrated nt */
} Lut_t;

/* Prototypes */

Lut_t *GetLut(char *file_name, int nband, Input_meta_t *input_meta,
              Img_coord_int_t *input_size);
bool FreeLut(Lut_t *this);

#endif
