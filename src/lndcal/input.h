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
#define MYHDF_H
#include "lndcal.h"
#include "bool.h"
#include "const.h"
#include "date.h"
#include "param.h"

#define ANGLE_FILL (-999.0)
#define WRS_FILL (-1)
#define GAIN_BIAS_FILL (-999.0)

/* Input file type definition */

typedef enum {
  INPUT_TYPE_NULL = -1,
  INPUT_TYPE_BINARY = 0, 
  INPUT_TYPE_BINARY_2BYTE_BIG,
  INPUT_TYPE_BINARY_2BYTE_LITTLE,
  INPUT_TYPE_GEOTIFF,
  INPUT_TYPE_MAX
} Input_type_t;

/* Structure for the 'input metadata' data type */

typedef struct {
  Provider_t provider;
                           /* Data provider type */
  Sat_t sat;               /* Satellite */
  Inst_t inst;             /* Instrument */
  Date_t acq_date;         /* Acqsition date/time (scene center) */
  bool time_fill;          /* Acqsition time fill; true = fill value (0h) */
  Date_t prod_date;        /* Production date (must be available for ETM) */
  float sun_zen;           /* Solar zenith angle (radians; scene center) */
  float sun_az;            /* Solar azimuth angle (radians; scene center) */
  Wrs_t wrs_sys;	   /* WRS system */
  int ipath;               /* WRS path number */
  int irow;                /* WRS row number */
  unsigned char fill;      /* Fill value */
  int iband[NBAND_REFL_MAX];
                           /* Band numbers */
  int iband_th;            /* Thermal Band number= (6) */
  Gain_t gain_set[NBAND_REFL_MAX]; 
                           /* Band gain settings (ETM+ only) */
  Gain_t gain_setting_th; 
                           /* Band gain settings Thermal (ETM+ only) */
  float gain[NBAND_REFL_MAX];
                           /* Band gain (MSS and TM only) */
  float bias[NBAND_REFL_MAX]; 
                           /* Band bias (MSS and TM only) */
  float gain_th;           /* Thermal band gain (MSS and TM only) */
  float bias_th;           /* Thermal band bias (MSS and TM only) */
} Input_meta_t;

/* Structure for the 'input' data type */

typedef struct {
  char *file_header_name;  /* Input image header file name */
  Input_type_t file_type;  /* Type of the input image files */
  Input_meta_t meta;       /* Input metadata */
  int nband;               /* Number of input image files (bands) */
  int nband_th;            /* Number of thermal input image files (0 or 1) */
  Img_coord_int_t size;    /* Input file size */
  Img_coord_int_t size_th; /* Input (thermal) file size */
  char *file_name[NBAND_REFL_MAX];  
                           /* Name of the input image files */
  char *file_name_th;      /* Name of the thermal input image files */
  bool open[NBAND_REFL_MAX]; 
                           /* Flag to indicate whether the specific input file 
			      is open for access; 'true' = open, 
			     'false' = not open */
  bool open_th;            /* thermal open flag */
  FILE *fp_bin[NBAND_REFL_MAX];  
                           /* File pointer for binary files */
  FILE *fp_bin_th;         /* File pointer for thremal binary file */
  void *fp_tiff[NBAND_REFL_MAX];
                           /* File pointer for TIFF file */
  void *fp_tiff_th;        /* File pointer for thermal TIFF file */
  void *buf[NBAND_REFL_MAX];
                           /* Input data buffer (one line of data) */
  void *buf_th;            /* Input (thermal) data buffer (one line of data) */
  int short_flag;          /* =0 8-bit input,  =1 16-bit input */
  int swap_flag;           /* =0 no byte swap, =1 perform byte swap */
  bool dnout;              /* dn output flag */
  float dn_map[4];         /* map from(0,1) -> (2,3) */
} Input_t;

/* Prototypes */

Input_t *OpenInput(char *file_header_name, Param_t *param);
bool GetInputLine(Input_t *this, int iband, int iline, unsigned char *line);
bool GetInputLineTh(Input_t *this, int iline, unsigned char *line);
bool CloseInput(Input_t *this);
bool FreeInput(Input_t *this);
bool InputMetaCopy(Input_meta_t *this, int nband, Input_meta_t *copy);
bool GetHeaderInput(Input_t *this, char *file_header_name, Param_t *param);

#endif
