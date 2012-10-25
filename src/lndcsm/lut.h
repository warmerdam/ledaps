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

#include "lndcsm.h"
#include "input.h"
#include "bool.h"

#define NBAND_CSM_LUT (1)

/* Structure for the 'lut' data type */

typedef struct {
  char *file_name;          /* Lookup table file name */
  int nband;                /* Number of bands */
  int in_fill;              /* Input fill value */
  int output_fill;          /* Output fill value */
  int csm_fill;             /* CSM output fill value */
  int csm_land;             /* CSM output land (non-cloud/snow) value */
  int csm_cloud;            /* CSM output cloud value */
  int csm_snow;             /* CSM output snow/icd value */
  int csm_water;            /* CSM output water value */
  int csm_shadow;           /* CSM output shadow value */
  float refl_scale_factor;  /* reflective scale factor */
  float therm_scale_factor; /* thermal scale factor */
  float celc2kelvin;        /* celcius to kelvin */
  float top_distance;       /* Distance of the top line from the bottom line
                               in unscaled brightness units */
  float keep_ratio;         /* Ratio of values to keep */
  Input_meta_t meta;        /* Input metadata */

  char* long_name;          /* long name */
  char* units;              /* units */
  int valid_range[2];       /* valid range */
  double scale_factor;      /* scale factor */
  double scale_factor_err;  /* scale factor error */
  double add_offset;        /* add offset */
  double add_offset_err;    /* add offset error */
  double calibrated_nt;     /* calibrated nt */
} Lut_t;

/* Prototypes */

Lut_t *GetLut(char *file_name, int nband, Input_meta_t *input_meta);
bool FreeLut(Lut_t *this);

#endif
