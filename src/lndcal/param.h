/*
!C****************************************************************************

!File: param.h

!Description: Header file for 'param.c' - see 'param.c' for more information.

!Revision History:
 Revision 1.0 2000/11/07
 Robert Wolfe
 Original Version.

 Revision 1.1 2000/12/13
 Sadashiva Devadiga
 Modified to accept parameters from command line or file.

 Revision 1.2 2001/05/08
 Sadashiva Devadiga
 Cleanup.

 Revision 1.3 2002/05/10
 Robert Wolfe
 Added separate output SDS name.

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
      phone: 301-614-5508               Lanham, MD 20706  

      Sadashiva Devadiga (Code 922)
      MODIS Land Team Support Group     SSAI
      devadiga@ltpmail.gsfc.nasa.gov    5900 Princess Garden Pkway, #300
      phone: 301-614-5549               Lanham, MD 20706

 ! Design Notes:
   1. Structures are declared for the 'input_space_type' and 'param' data types.
   2. The acronym SDS stands for Science Data Set.
  
!END****************************************************************************
*/

#ifndef PARAM_H
#define PARAM_H

#include "bool.h"

/* Parameter data structure type definition */

typedef struct {
  char *param_file_name;         /* Parameter file name                */
  char *input_header_file_name;  /* Input image header file name       */
  char *lut_file_name;           /* Lookup table file name             */
  char *work_order_file_name;    /* Work Order (wo) file name          */
  char *gold_file_name;          /* G-old file name                    */
  char *gold_2003_name;          /* G-old file name                    */
  char *gnew_file_name;          /* G-new file name                    */
  char *output_file_name;        /* Output image HDF file name         */
  char *output_therm_file_name;  /* Output thermal image HDF file name */
  char *PGEVersion;              /* PGE Version                        */
  char *ProcessVersion;          /* Process Version                    */
  float dn_map[4];               /* map from(0,1) -> (2,3)             */
  bool dnout;                    /* dn output flag                     */
  bool work_order_flag;          /* work order input flag              */
  bool gold_flag;                /* G-old input flag                   */
  bool gold_2003_flag;           /* G-old input flag                   */
  bool gnew_flag;                /* G-new input flag                   */
  bool ETM_GB;                   /* ETM_GB flag                        */
  bool RE_CAL;                   /* re-calibration flag  (fgao added)  */             
  int  est_gainbias;             /* estimate Gain/bias?  1=yes         */
} Param_t;

/* Prototypes */

Param_t *GetParam(int argc, const char **argv);
bool FreeParam(Param_t *this);
bool setETM_GB(Param_t *this);

#endif
