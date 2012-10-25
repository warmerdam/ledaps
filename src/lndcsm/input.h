/*
!C****************************************************************************

!File: input.h

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

#ifndef INPUT_H
#define INPUT_H

#include <stdlib.h>
#include <stdio.h>
#include "lndcsm.h"
#include "myhdf.h"
#include "const.h"
#include "date.h"

typedef struct {
  Provider_t provider;
                           /* Data provider type */
  Sat_t sat;               /* Satellite */
  Inst_t inst;             /* Instrument */
  Date_t acq_date;         /* Acqsition date/time (scene center) */
  Date_t prod_date;        /* Production date (must be available for ETM) */
  float sun_zen;           /* Solar zenith angle (radians; scene center) */
  float sun_az;            /* Solar azimuth angle (radians; scene center) */
  Wrs_t wrs_sys;	   /* WRS system */
  int ipath;               /* WRS path number */
  int irow;                /* WRS row number */
  int fill;                /* Fill value */
  int iband[NBAND_REFL_MAX];
                           /* Band numbers */
} Input_meta_t;

/* Structure for the 'input' data type */

typedef struct {
  char *file_name;         /* Input image file name */
  bool open;               /* Open file flag; open = true */
  Input_meta_t meta;       /* Input metadata */
  int nband;               /* Number of input image files (bands) */
  Img_coord_int_t size;    /* Input file size */
  int32 sds_file_id;       /* SDS file id */
  Myhdf_sds_t sds[NBAND_REFL_MAX];
                           /* SDS data structures */
  int16 *buf[NBAND_REFL_MAX];
                           /* Input data buffer (one line of data) */
} Input_t;

/* Prototypes */

Input_t *OpenInput(char *file_name);
bool GetInputLine(Input_t *this, int iband, int iline, int *line);
bool CloseInput(Input_t *this);
bool FreeInput(Input_t *this);
bool InputMetaCopy(Input_meta_t *this, int nband, Input_meta_t *copy);
bool GetInputMeta(Input_t *this);

#endif
