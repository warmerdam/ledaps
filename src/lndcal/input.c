/*
!C****************************************************************************

!File: input.c
  
!Description: Functions reading data from the input data file.

!Revision History:
 Revision 1.0 2001/05/08
 Robert Wolfe
 Original Version.

 Revision 1.1 2002/05/02
 Robert Wolfe
 Added handling for SDS's with ranks greater than 2.

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
   1. The following public functions handle the input data:

	OpenInput - Setup 'input' data structure and open file for access.
	CloseInput - Close the input file.
	FreeInput - Free the 'input' data structure memory.

   2. 'OpenInput' must be called before any of the other routines.  
   3. 'FreeInput' should be used to free the 'input' data structure.
   4. The only input file type supported is HDF.

!END****************************************************************************
*/

#include <stdlib.h>
#include "input.h"
#include "util.h"
#include "error.h"
#include "mystring.h"
#include "const.h"
#include "date.h"
#include "geotiffio.h"
#include "xtiffio.h"
#define INPUT_FILL (0)
typedef enum {
  HEADER_NULL = -1,
  HEADER_START = 0,
  HEADER_FILE_TYPE,
  HEADER_PROIVDER,
  HEADER_SAT,
  HEADER_INST,
  HEADER_ACQ_DATE,
  HEADER_ACQ_TIME,
  HEADER_PROD_DATE,
  HEADER_SUN_ZEN,
  HEADER_SUN_AZ,
  HEADER_WRS_SYS,
  HEADER_WRS_PATH,
  HEADER_WRS_ROW,
  HEADER_NBAND,
  HEADER_BANDS,
  HEADER_GAIN_SET,
  HEADER_GAIN_SETTINGS_TH,
  HEADER_GAIN,
  HEADER_BIAS,
  HEADER_NSAMPLE,
  HEADER_NLINE,
  HEADER_FILES,
  HEADER_NBAND_TH,
  HEADER_BANDS_TH,
  HEADER_FILES_TH,
  HEADER_NSAMPLE_TH,
  HEADER_NLINE_TH,
  HEADER_GAIN_TH,
  HEADER_BIAS_TH,
  HEADER_END,
  HEADER_MAX
} Header_key_t;

Key_string_t Header_string[HEADER_MAX] = {
  {(int)HEADER_START,     "HEADER_FILE"},
  {(int)HEADER_FILE_TYPE, "FILE_TYPE"},
  {(int)HEADER_PROIVDER,  "DATA_PROVIDER"},
  {(int)HEADER_SAT,       "SATTELITE"},
  {(int)HEADER_INST,      "INSTRUMENT"},
  {(int)HEADER_ACQ_DATE,  "ACQUISITION_DATE"},
  {(int)HEADER_ACQ_TIME,  "ACQUISITION_TIME"},
  {(int)HEADER_PROD_DATE, "PRODUCTION_DATE"},
  {(int)HEADER_SUN_ZEN,   "SOLAR_ZENITH"},
  {(int)HEADER_SUN_AZ,    "SOLAR_AZIMUTH"},
  {(int)HEADER_WRS_SYS,   "WRS_SYSTEM"},
  {(int)HEADER_WRS_PATH,  "WRS_PATH"},
  {(int)HEADER_WRS_ROW,   "WRS_ROW"},
  {(int)HEADER_NBAND,     "NBAND"},
  {(int)HEADER_BANDS,     "BANDS"},
  {(int)HEADER_GAIN_SET,  "GAIN_SETTINGS"},
  {(int)HEADER_GAIN_SETTINGS_TH, "GAIN_SETTINGS_TH"},
  {(int)HEADER_GAIN,      "GAIN"},
  {(int)HEADER_BIAS,      "BIAS"},
  {(int)HEADER_NSAMPLE,   "NSAMPLE"},
  {(int)HEADER_NLINE,     "NLINE"},
  {(int)HEADER_FILES,     "FILE_NAMES"},
  {(int)HEADER_NBAND_TH,  "NBAND_TH"},
  {(int)HEADER_BANDS_TH,  "BANDS_TH"},
  {(int)HEADER_FILES_TH,  "FILE_NAMES_TH"},
  {(int)HEADER_NSAMPLE_TH,"NSAMPLE_TH"},
  {(int)HEADER_NLINE_TH,  "NLINE_TH"},
  {(int)HEADER_GAIN_TH,   "GAIN_TH"},
  {(int)HEADER_BIAS_TH,   "BIAS_TH"},
  {(int)HEADER_END,       "END"}
};

Key_string_t Input_type_string[INPUT_TYPE_MAX] = {
  {(int)INPUT_TYPE_BINARY,              "BINARY"},
  {(int)INPUT_TYPE_BINARY_2BYTE_BIG,    "BINARY_2BYTE_BIG"},
  {(int)INPUT_TYPE_BINARY_2BYTE_LITTLE, "BINARY_2BYTE_LITTLE"},
  {(int)INPUT_TYPE_GEOTIFF,             "GEOTIFF"}
};

/* Functions */

Input_t *OpenInput(char *file_header_name, Param_t *param)
/* 
!C******************************************************************************

!Description: 'OpenInput' sets up the 'input' data structure, opens the
 input file for read access.
 
!Input Parameters:
 file_name      input file name
 sds_name       name of sds to be read
 iband          band number for application of band offset
 rank           rank of the input SDS
 dim            dimension flags; the line and sample dimensions are 
                indicated by a -1 and -2 in this array, respectively;
		the index in the other dimensions are indicated by a value 
		of zero or greater
 input_space_type  input space type; either 'SWATH_SPACE' for L1 and L2 input 
                data or 'GRID_SPACE' for L2G, L3 and L4 input data

!Output Parameters:
 dim            dimension flags
 (returns)      'input' data structure or NULL when an error occurs

!Team Unique Header:

 ! Design Notes:
   1. When 'OpenInput' returns, the file is open for HDF access and the 
      SDS is open for access.
   2. For an input space type of 'GRID_SPACE', a band number of -1 should be 
      given.
   3. The only input HDF data types supported are CHAR8, UNIT8, INT16 and
      UNIT16.
   4. An error status is returned when:
       a. the SDS rank is less than 2 or greater than 'MYHDF_MAX_RANK'
       b. the band number is less than -1 or greater than or equal to
          'NBAND_OFFSET'
       c. either none or more than one dimension is given for the line 
          or sample dimension
       d. an invalid dimension field is given
       e. duplicating strings is not successful
       f. errors occur when opening the input HDF file
       g. errors occur when reading SDS dimensions or attributes
       h. errors occur when opening the SDS for read access
       i. the given SDS rank or dimensions do not match the input file
       j. for a input space type of SWATH_SPACE, the dimensions of a swath 
          are not 1, 2 or 4 times the nominal size of a MODIS swath
       k. for a input space type of SWATH_SPACE, the number of lines is not 
          an integral multiple of the size of the scan at the given resolution
       l. memory allocation is not successful
       m. an invalid input data type is not supported.
   5. Error messages are handled with the 'RETURN_ERROR' macro.
   6. 'FreeInput' should be called to deallocate memory used by the 
      'input' data structures.
   7. 'CloseInput' should be called after all of the data is written and 
      before the 'input' data structure memory is released.

!END****************************************************************************
*/
{
  Input_t *this;
  char *error_string = (char *)NULL;
  int ib;
  TIFF *fp_tiff;
  TIFF *fp_tiff_th;
  uint8 *buf_uint8;
  uint8 *buf_uint8_th;
  uint16 *buf_uint16;
  uint16 *buf_uint16_th;

  /* Create the Input data structure */
  
  this = (Input_t *)malloc(sizeof(Input_t));
  if (this == (Input_t *)NULL) 
    RETURN_ERROR("allocating Input data structure", "OpenInput", 
                 (Input_t *)NULL);

  /* Initialize and get input from header file */

  if (!GetHeaderInput(this, file_header_name, param)) { 
    free(this);
    RETURN_ERROR("getting input from header file", 
                 "OpenInput", (Input_t *)NULL);
  }

  /* Open files for access */

  this->short_flag= this->file_type==INPUT_TYPE_BINARY_2BYTE_BIG  ||
                    this->file_type==INPUT_TYPE_BINARY_2BYTE_LITTLE ? 1 : 0;

  this->swap_flag= 0;
  if ( big_endian() ) {
    if ( this->file_type == INPUT_TYPE_BINARY_2BYTE_LITTLE )this->swap_flag= 1;
  } else {
    if ( this->file_type == INPUT_TYPE_BINARY_2BYTE_BIG )   this->swap_flag= 1;
  }

  if (
      this->file_type == INPUT_TYPE_BINARY              ||
      this->file_type == INPUT_TYPE_BINARY_2BYTE_BIG    ||
      this->file_type == INPUT_TYPE_BINARY_2BYTE_LITTLE
     ) {

    for (ib = 0; ib < this->nband; ib++) {
      this->fp_bin[ib] = fopen(this->file_name[ib], "r");
      if (this->fp_bin[ib] == (FILE *)NULL) 
        {
        error_string = "opening binary file";
        break;
        }
      this->open[ib] = true;
    }
    if ( this->nband_th == 1 )
      {
      this->fp_bin_th = fopen(this->file_name_th, "r");
      if (this->fp_bin_th == (FILE *)NULL) 
        error_string = "opening thermal binary file";
      else
        this->open_th = true;
      }

  } else if (this->file_type == INPUT_TYPE_GEOTIFF) {

    for (ib = 0; ib < this->nband; ib++) {
      fp_tiff = XTIFFOpen(this->file_name[ib], "r");
      if (fp_tiff == (TIFF *)NULL) {
        error_string = "opening TIFF file";
        break;
      }
      this->fp_tiff[ib] = (void *)fp_tiff;
      this->open[ib] = true;
    }
    if ( this->nband_th == 1 )
      {
      fp_tiff_th = XTIFFOpen(this->file_name[ib], "r");
      if (fp_tiff_th == (TIFF *)NULL) 
        error_string = "opening thermal TIFF file";
      this->fp_tiff_th = (void *)fp_tiff_th;
      this->open_th = true;

      }
  } else 
    error_string = "invalid file type";

  /* Allocate input buffer (if needed) */

  if (error_string == (char *)NULL) {

    if ( this->file_type == INPUT_TYPE_BINARY_2BYTE_BIG    ||
         this->file_type == INPUT_TYPE_BINARY_2BYTE_LITTLE    ){
      if (sizeof(uint16) != sizeof(unsigned short int)) {
        for (ib = 0; ib < this->nband; ib++) {
          buf_uint16 = (uint16 *)calloc(this->size.s, sizeof(uint16));
          if (buf_uint16 == (uint16 *)NULL) {
            error_string = "allocating input buffer";
            break;
  	  }
  	  this->buf[ib] = (void *)buf_uint16;
        }
        if ( this->nband_th == 1 )
          {
          buf_uint16_th = (uint16 *)calloc(this->size.s, sizeof(uint16));
  	  this->buf_th = (void *)buf_uint16_th;
          }
        }
      }
    else
      {
      if (sizeof(uint8) != sizeof(unsigned char)) {
        for (ib = 0; ib < this->nband; ib++) {
          buf_uint8 = (uint8 *)calloc(this->size.s, sizeof(uint8));
          if (buf_uint8 == (uint8 *)NULL) {
            error_string = "allocating input buffer";
            break;
  	  }
  	  this->buf[ib] = (void *)buf_uint8;
        }
        if ( this->nband_th == 1 )
          {
          buf_uint8_th = (uint8 *)calloc(this->size.s, sizeof(uint8));
  	  this->buf_th = (void *)buf_uint8_th;
          }
        }
      }
    }

  if (error_string != (char *)NULL) {
    for (ib = 0; ib < this->nband; ib++) {
      if (this->file_name[ib] != (char *)NULL) free(this->file_name[ib]);
      if (this->buf[ib] != (void *)NULL) free(this->buf[ib]);
      if (this->open[ib]) {
        if ( this->file_type == INPUT_TYPE_BINARY              ||
             this->file_type == INPUT_TYPE_BINARY_2BYTE_BIG    ||
             this->file_type == INPUT_TYPE_BINARY_2BYTE_LITTLE    )
	  fclose(this->fp_bin[ib]);
        else if (this->file_type == INPUT_TYPE_GEOTIFF) {
	  fp_tiff = (TIFF *)this->fp_tiff[ib];
          XTIFFClose(fp_tiff);
        }
        this->open[ib] = false;
      }
    }
    if (this->file_name_th != (char *)NULL) free(this->file_name_th);
    if (this->buf_th != (void *)NULL) free(this->buf_th);
    if ( this->file_type == INPUT_TYPE_BINARY              ||
         this->file_type == INPUT_TYPE_BINARY_2BYTE_BIG    ||
         this->file_type == INPUT_TYPE_BINARY_2BYTE_LITTLE    )
      fclose(this->fp_bin_th);  
    else if (this->file_type == INPUT_TYPE_GEOTIFF) 
      {
      fp_tiff_th = (TIFF *)this->fp_tiff_th;
      XTIFFClose(fp_tiff_th);
      }
    this->open_th = false;
    free(this->file_header_name);
    free(this);
    RETURN_ERROR(error_string, "OpenInput", (Input_t *)NULL);
  }

  return this;
}

bool GetInputLine(Input_t *this, int iband, int iline, unsigned char *line) 
{
  long loc;
  void *buf_void;
  int is;
  TIFF *fp_tiff;
  uint8 *buf_uint8;
  uint16 *buf_uint16;

  if (this == (Input_t *)NULL) 
    RETURN_ERROR("invalid input structure", "GetInputLine", false);
  if (iband < 0  ||  iband >= this->nband) 
    RETURN_ERROR("band index out of range", "GetInputLine", false);
  if (iline < 0  ||  iline >= this->size.l) 
    RETURN_ERROR("line index out of range", "GetInputLine", false);
  if (!this->open[iband])
    RETURN_ERROR("band not open", "GetInputLine", false);

  if ( this->file_type == INPUT_TYPE_BINARY_2BYTE_BIG    ||
       this->file_type == INPUT_TYPE_BINARY_2BYTE_LITTLE    ){
    if (sizeof(uint16) == sizeof(unsigned short int))
      buf_void = (void *)line;
    else {
      buf_void = this->buf[iband];
      buf_uint16 = (uint16 *)this->buf[iband];
    }
  }
 else
  {
    if (sizeof(uint8) == sizeof(unsigned char))
      buf_void = (void *)line;
    else {
      buf_void = this->buf[iband];
      buf_uint8 = (uint8 *)this->buf[iband];
    }
  }

  if (this->file_type == INPUT_TYPE_BINARY) {
    loc = (long)iline * this->size.s * sizeof(uint8);
    if (fseek(this->fp_bin[iband], loc, SEEK_SET))
      RETURN_ERROR("error seeking line (binary)", "GetInputLine", false);
    if (fread(buf_void, sizeof(uint8), (size_t)this->size.s, 
              this->fp_bin[iband]) != (size_t)this->size.s)
      RETURN_ERROR("error reading line (binary)", "GetInputLine", false);
   } else if (this->file_type == INPUT_TYPE_BINARY_2BYTE_BIG    ||
             this->file_type == INPUT_TYPE_BINARY_2BYTE_LITTLE    ){
   loc = (long)iline * this->size.s * sizeof(uint16);
    if (fseek(this->fp_bin[iband], loc, SEEK_SET))
      RETURN_ERROR("error seeking line (binary)", "GetInputLine", false);
    if (fread(buf_void, sizeof(uint16), (size_t)this->size.s, 
              this->fp_bin[iband]) != (size_t)this->size.s)
      RETURN_ERROR("error reading line (binary)", "GetInputLine", false);
  } else if (this->file_type == INPUT_TYPE_GEOTIFF) {
    fp_tiff = (TIFF *)this->fp_tiff[iband];
    if (!TIFFReadScanline(fp_tiff, buf_void, iline, 0)) 
      RETURN_ERROR("error reading line (TIFF)", "GetInputLine", false);
  }

  if ( this->file_type == INPUT_TYPE_BINARY_2BYTE_BIG    ||
       this->file_type == INPUT_TYPE_BINARY_2BYTE_LITTLE    ){
    if (sizeof(uint16) != sizeof(unsigned short int))
      for (is = 0; is < this->size.s; is++) 
        line[is] = (unsigned char)buf_uint16[is];
    }
  else
    {
    if (sizeof(uint8) != sizeof(unsigned char))
      for (is = 0; is < this->size.s; is++) 
        line[is] = (unsigned char)buf_uint8[is];
    }

  return true;
}
bool GetInputLineTh(Input_t *this, int iline, unsigned char *line) 
{
  long loc;
  void *buf_void;
  int is;
  TIFF *fp_tiff_th;
  uint8 *buf_uint8;
  uint16 *buf_uint16;

  if (this == (Input_t *)NULL) 
    RETURN_ERROR("invalid input structure", "GetInputLine", false);
  if ( this->nband_th < 1 ) 
    RETURN_ERROR("no thermal input band", "GetInputLine", false);
  if (iline < 0  ||  iline >= this->size_th.l) 
    RETURN_ERROR("line index out of range", "GetInputLine", false);
  if (!this->open_th)
    RETURN_ERROR("band not open", "GetInputLine", false);

  if ( this->file_type == INPUT_TYPE_BINARY_2BYTE_BIG    ||
       this->file_type == INPUT_TYPE_BINARY_2BYTE_LITTLE    ){
    if (sizeof(uint16) == sizeof(unsigned short int))
      buf_void = (void *)line;
    else {
      buf_void = this->buf_th;
      buf_uint16 = (uint16 *)this->buf_th;
      }
    }
   else
    {
    if (sizeof(uint8) == sizeof(unsigned char))
      buf_void = (void *)line;
    else {
      buf_void = this->buf_th;
      buf_uint8 = (uint8 *)this->buf_th;
      }
    }

  if (this->file_type == INPUT_TYPE_BINARY) {
    loc = (long)iline * this->size_th.s * sizeof(uint8);
    if (fseek(this->fp_bin_th, loc, SEEK_SET))
      RETURN_ERROR("error seeking line (binary)", "GetInputLine", false);
    if (fread(buf_void, sizeof(uint8), (size_t)this->size_th.s, 
              this->fp_bin_th) != (size_t)this->size_th.s)
      RETURN_ERROR("error reading line (binary)", "GetInputLine", false);
  } else if (this->file_type == INPUT_TYPE_BINARY_2BYTE_BIG    ||
             this->file_type == INPUT_TYPE_BINARY_2BYTE_LITTLE    ){
    loc = (long)iline * this->size_th.s * sizeof(uint16);
    if (fseek(this->fp_bin_th, loc, SEEK_SET))
      RETURN_ERROR("error seeking line (binary)", "GetInputLine", false);
    if (fread(buf_void, sizeof(uint16), (size_t)this->size_th.s, 
              this->fp_bin_th) != (size_t)this->size_th.s)
      RETURN_ERROR("error reading line (binary)", "GetInputLine", false);
  } else if (this->file_type == INPUT_TYPE_GEOTIFF) {
    fp_tiff_th = (TIFF *)this->fp_tiff_th;
    if (!TIFFReadScanline(fp_tiff_th, buf_void, iline, 0)) 
      RETURN_ERROR("error reading line (TIFF)", "GetInputLine", false);
  }

  if ( this->file_type == INPUT_TYPE_BINARY_2BYTE_BIG    ||
       this->file_type == INPUT_TYPE_BINARY_2BYTE_LITTLE    ){
    if (sizeof(uint16) != sizeof(unsigned short int))
      for (is = 0; is < this->size_th.s; is++) 
        line[is] = (unsigned char)buf_uint16[is];
    }
  else
    {
    if (sizeof(uint8) != sizeof(unsigned char))
      for (is = 0; is < this->size_th.s; is++) 
        line[is] = (unsigned char)buf_uint8[is];
    }

  return true;
}


bool CloseInput(Input_t *this)
/* 
!C******************************************************************************

!Description: 'CloseInput' ends SDS access and closes the input file.
 
!Input Parameters:
 this           'input' data structure; the following fields are input:
                   open, sds.id, sds_file_id

!Output Parameters:
 this           'input' data structure; the following fields are modified:
                   open
 (returns)      status:
                  'true' = okay
		  'false' = error return

!Team Unique Header:

 ! Design Notes:
   1. An error status is returned when:
       a. the file is not open for access
       b. an error occurs when closing access to the SDS.
   2. Error messages are handled with the 'RETURN_ERROR' macro.
   3. 'OpenInput' must be called before this routine is called.
   4. 'FreeInput' should be called to deallocate memory used by the 
      'input' data structure.

!END****************************************************************************
*/
{
  int ib;
  bool none_open;
  TIFF *fp_tiff;
  TIFF *fp_tiff_th;

  if (this == (Input_t *)NULL) 
    RETURN_ERROR("invalid input structure", "CloseInput", false);

  none_open = true;
  for (ib = 0; ib < this->nband; ib++) {
    if (this->open[ib]) {
      none_open = false;
      if (this->file_type == INPUT_TYPE_BINARY)
        fclose(this->fp_bin[ib]);
      else if (this->file_type == INPUT_TYPE_GEOTIFF) {
        fp_tiff = this->fp_tiff[ib];
        XTIFFClose(fp_tiff);
      }
      this->open[ib] = false;
    }
  }
  /*** now close the thermal file ***/
  if (this->open_th) 
    {
    if (this->file_type == INPUT_TYPE_BINARY)
      fclose(this->fp_bin_th);
    else if (this->file_type == INPUT_TYPE_GEOTIFF) 
      {
      fp_tiff_th = this->fp_tiff_th;
      XTIFFClose(fp_tiff_th);
      }
    this->open_th = false;
    }

  if (none_open)
    RETURN_ERROR("no files open", "CloseInput", false);

  return true;
}


bool FreeInput(Input_t *this)
/* 
!C******************************************************************************

!Description: 'FreeInput' frees the 'input' data structure memory.
 
!Input Parameters:
 this           'input' data structure; the following fields are input:
                   sds.rank, sds.dim[*].name, sds.name, file_name

!Output Parameters:
 (returns)      status:
                  'true' = okay (always returned)

!Team Unique Header:

 ! Design Notes:
   1. 'OpenInput' and 'CloseInput' must be called before this routine is called.
   2. An error status is never returned.

!END****************************************************************************
*/
{
  int ib;

  if (this != (Input_t *)NULL) {
    for (ib = 0; ib < this->nband; ib++) {
      if (this->file_name[ib] != (char *)NULL) free(this->file_name[ib]);
      if (this->buf[ib] != (void *)NULL) free(this->buf[ib]);
    }
    if (this->file_header_name != (char *)NULL) free(this->file_header_name);
    free(this);
  }

  return true;
}

bool InputMetaCopy(Input_meta_t *this, int nband, Input_meta_t *copy) 
{
  int ib;

  if (this == (Input_meta_t *)NULL) 
    RETURN_ERROR("invalid input structure", "InputMetaCopy", false);

  copy->provider = this->provider;
  copy->sat = this->sat;
  copy->inst = this->inst;
  if (!DateCopy(&this->acq_date, &copy->acq_date)) 
    RETURN_ERROR("copying acquisition date/time", "InputMetaCopy", false);
  if (!DateCopy(&this->prod_date, &copy->prod_date)) 
    RETURN_ERROR("copying production date/time", "InputMetaCopy", false);
  copy->time_fill = this->time_fill;
  copy->sun_zen = this->sun_zen;
  copy->sun_az = this->sun_az;
  copy->wrs_sys = this->wrs_sys;
  copy->ipath = this->ipath;
  copy->irow = this->irow;
  copy->fill = this->fill;

  for (ib = 0; ib < nband; ib++) {
    copy->iband[ib] = this->iband[ib];
    copy->gain_set[ib] = this->gain_set[ib];
    copy->gain[ib] = this->gain[ib];
    copy->bias[ib] = this->bias[ib];
  }
  copy->gain_th = this->gain_th;
  copy->bias_th = this->bias_th;
  copy->iband_th =this->iband_th;
  return true;
}

#define DATE_STRING_LEN (10)
#define TIME_STRING_LEN (15)
#define NLINE_MAX (20000)
#define NSAMP_MAX (20000)

bool GetHeaderInput(Input_t *this, char *file_header_name, Param_t *param) {
  char *error_string = (char *)NULL;
  int ib;
  FILE *fp;
  Key_t key;
  int nband_iband, nband_gain_set, nband_gain, nband_bias, nband_file_name;
  char acq_date[DATE_STRING_LEN + 1];
  char prod_date[DATE_STRING_LEN + 1];
  char acq_time[TIME_STRING_LEN + 1];
  char temp[MAX_STR_LEN + 1];
  bool got_start, got_end;
  int len;
  char line[MAX_STR_LEN + 1];
  Header_key_t header_key;
  int i;
  bool flag_gain_set, flag_gain, flag_bias;

  /* Populate the data structure */

  this->file_header_name = DupString(file_header_name);
  if (this->file_header_name == (char *)NULL) {
    free(this);
    RETURN_ERROR("duplicating file name", "GetHeaderInput", false);
  }

  this->file_type = INPUT_TYPE_NULL;
  this->meta.provider = PROVIDER_NULL;
  this->meta.sat = SAT_NULL;
  this->meta.inst = INST_NULL;
  this->meta.acq_date.fill = true;
  this->meta.time_fill = true;
  this->meta.prod_date.fill = true;
  this->meta.sun_zen = ANGLE_FILL;
  this->meta.sun_az = ANGLE_FILL;
  this->meta.wrs_sys = (Wrs_t)WRS_FILL;
  this->meta.ipath = -1;
  this->meta.irow = -1;
  this->meta.fill = INPUT_FILL;
  this->nband = 0;
  this->nband_th = 0;
  this->size.s = this->size.l = -1;
  for (ib = 0; ib < NBAND_REFL_MAX; ib++) {
    this->meta.iband[ib] = -1;
    this->meta.gain_set[ib] = GAIN_NULL;
    this->meta.gain[ib] = GAIN_BIAS_FILL;
    this->meta.bias[ib] = GAIN_BIAS_FILL;
    this->file_name[ib] = (char *)NULL;
    this->open[ib] = false;
    this->fp_bin[ib] = (FILE *)NULL;
    this->fp_tiff[ib] = (void *)NULL;
    this->buf[ib] = (void *)NULL;
  }
  this->nband_th = 0;
  this->open_th = false;
  this->meta.gain_th = GAIN_BIAS_FILL;
  this->meta.bias_th = GAIN_BIAS_FILL;
  this->dnout= false;
  this->dn_map[0]= 0.0;
  this->dn_map[1]= 0.0;
  this->dn_map[2]= 0.0;
  this->dn_map[3]= 0.0;
  this->file_name_th  = (char *)NULL;
  this->fp_bin_th  = (FILE *)NULL;
  this->fp_tiff_th  = (void *)NULL;
  this->buf_th  = (void *)NULL;

  /* Open the header file */
  
  if ((fp = fopen(this->file_header_name, "r")) == (FILE *)NULL) {
    free(this->file_header_name);
    free(this);
    RETURN_ERROR("unable to open header file", "GetHeaderInput", false);
  }

  /* Parse the header file */

  got_start = got_end = false;
  nband_iband = nband_gain_set = -1;
  nband_gain = nband_bias = nband_file_name = -1;
  acq_date[0] = acq_time[0] = '\0';
  prod_date[0] = '\0';
  flag_gain_set = flag_gain = flag_bias = false;

  while((len = GetLine(fp, line)) > 0) {
    if (!StringParse(line, &key)) {
      sprintf(temp, "parsing header file; line = %s", line);
      error_string = temp;
      break;
    }
    if (key.len_key <= 0) continue;
    if (key.key[0] == '#') continue;
    header_key = (Header_key_t)
       KeyString(key.key, key.len_key, Header_string, 
		 (int)HEADER_NULL, (int)HEADER_MAX);

    if (header_key == HEADER_NULL) {
      key.key[key.len_key] = '\0';
    }
    if (!got_start) {
      if (header_key == HEADER_START) {
        if (key.nval != 0) {
	  error_string = "no value expected";
	  break;
	}
        got_start = true;
	continue;
      } else {
        error_string  = "no start key in header";
        break;
      }
    }

    /* Get the value for each keyword */

    switch (header_key) {

      case HEADER_FILE_TYPE:
        if (key.nval <= 0) {
	  error_string = "no file type value";
	  break; 
	} else if (key.nval > 1) {
	  error_string = "too many file type values";
	  break; 
	}
        this->file_type = (Input_type_t)
	  KeyString(key.value[0], key.len_value[0],
	            Input_type_string, 
		    (int)INPUT_TYPE_NULL, (int)INPUT_TYPE_MAX);
        if (this->file_type == INPUT_TYPE_NULL) {
	  key.value[0][key.len_value[0]] = '\0';
          sprintf(temp, "invalid file type; value = %s", key.value[0]);
          error_string = temp;
	}
        break;

      case HEADER_PROIVDER:
        if (key.nval <= 0) {
	  error_string = "no provider value";
	  break; 
	} else if (key.nval > 1) {
	  error_string = "too many provider values";
	  break; 
	}
        this->meta.provider = (Provider_t)
	  KeyString(key.value[0], key.len_value[0],
	            Provider_string, (int)PROVIDER_NULL, (int)PROVIDER_MAX);
        if (this->meta.provider == PROVIDER_NULL) {
	  key.value[0][key.len_value[0]] = '\0';
          sprintf(temp, "invalid data provider; value = %s", key.value[0]);
          error_string = temp;
	}
        break;

      case HEADER_SAT:
        if (key.nval <= 0) {
	  error_string = "no satellite value";
	  break; 
	} else if (key.nval > 1) {
	  error_string = "too many satellite values";
	  break; 
	}
        this->meta.sat = (Sat_t)
	  KeyString(key.value[0], key.len_value[0],
	            Sat_string, (int)SAT_NULL, (int)SAT_MAX);
        if (this->meta.sat == SAT_NULL) {
	  key.value[0][key.len_value[0]] = '\0';
          sprintf(temp, "invalid satellite; value = %s", key.value[0]);
          error_string = temp;
	}
        break;

      case HEADER_INST:
        if (key.nval <= 0) {
	  error_string = "no instrument value";
	  break; 
	} else if (key.nval > 1) {
	  error_string = "too many instrument values";
	  break; 
	}
        this->meta.inst = (Inst_t)
	  KeyString(key.value[0], key.len_value[0],
	            Inst_string, (int)INST_NULL, (int)INST_MAX);
        if (this->meta.inst == INST_NULL) {
	  key.value[0][key.len_value[0]] = '\0';
          sprintf(temp, "invalid satellite; value = %s", key.value[0]);
          error_string = temp;
	}
        break;

      case HEADER_ACQ_DATE:
        if (key.nval <= 0) {
	  error_string = "no acquisition date value";
	  break; 
	} else if (key.nval > 1) {
	  error_string = "too many acquisition date values";
	  break; 
	}
        if (key.len_value[0] != DATE_STRING_LEN) {
	  error_string = "invalid acquisition date string";
	  break;
	}
	strncpy(acq_date, key.value[0], key.len_value[0]);
	acq_date[key.len_value[0]] = '\0';
        break;

      case HEADER_ACQ_TIME:
        if (key.nval <= 0) {
	  error_string = "no acquisition time value";
	  break; 
	} else if (key.nval > 1) {
	  error_string = "too many acquisition time values";
	  break; 
	}
        if (key.len_value[0] < 8  ||
	    key.len_value[0] > TIME_STRING_LEN) {
	  error_string = "invalid time string";
	  break;
	}
	strncpy(acq_time, key.value[0], key.len_value[0]);
	acq_time[key.len_value[0]] = '\0';
        break;

      case HEADER_PROD_DATE:
        if (key.nval <= 0) {
	  error_string = "no production date value";
	  break; 
	} else if (key.nval > 1) {
	  error_string = "too many production date values";
	  break; 
	}
        if (key.len_value[0] >0 && 
            key.len_value[0] != DATE_STRING_LEN) {
	  error_string = "invalid production date string";
	  break; 
	}
        if (key.len_value[0] >0 ){
  	  strncpy(prod_date, key.value[0], key.len_value[0]);
	  prod_date[key.len_value[0]] = '\0';
	}
        break;

      case HEADER_SUN_ZEN:
        if (key.nval <= 0) {
	  error_string = "no solar zenith angle value";
	  break; 
	} else if (key.nval > 1) {
	  error_string = "too many solar zenith angle values";
	  break; 
	}
        if (key.len_value[0] < 1) {
	  error_string = "invalid solar zenith angle string";
	  break;
	}
	strncpy(temp, key.value[0], key.len_value[0]);
        temp[key.len_value[0]] = '\0';

	if (sscanf(temp, "%f", &this->meta.sun_zen) != 1) {
	  error_string = "converting solar zenith angle";
	  break;
	}
	if (this->meta.sun_zen < -90.0  ||  
	    this->meta.sun_zen >  90.0) {
	  error_string = "solar zenith angle out of range";
	  break;
	}
	this->meta.sun_zen *= RAD;
        break;

      case HEADER_SUN_AZ:
        if (key.nval <= 0) {
	  error_string = "no solar azimuth angle value";
	  break; 
	} else if (key.nval > 1) {
	  error_string = "too many solar azimuth angle values";
	  break; 
	}
        if (key.len_value[0] < 1) {
	  error_string = "invalid solar azimuth angle string";
	  break;
	}
	strncpy(temp, key.value[0], key.len_value[0]);
        temp[key.len_value[0]] = '\0';

	if (sscanf(temp, "%f", &this->meta.sun_az) != 1) {
	  error_string = "converting solar azimuth angle";
	  break;
	}
	if (this->meta.sun_az < -360.0  ||  
	    this->meta.sun_az >  360.0) {
	  error_string = "solar azimuth angle out of range";
	  break;
	}
	this->meta.sun_az *= RAD;
        break;

      case HEADER_WRS_SYS:
        if (key.nval <= 0) {
	  error_string = "no WRS system value";
	  break; 
	} else if (key.nval > 1) {
	  error_string = "too many WRS system values";
	  break; 
	}
        this->meta.wrs_sys = (Wrs_t)
	  KeyString(key.value[0], key.len_value[0],
	            Wrs_string, (int)WRS_NULL, (int)WRS_MAX);
        if (this->meta.wrs_sys == WRS_NULL) {
	  key.value[0][key.len_value[0]] = '\0';
          sprintf(temp, "invalid WRS system; value = %s", key.value[0]);
          error_string = temp;
	}
        break;

      case HEADER_WRS_PATH:
        if (key.nval <= 0) {
	  error_string = "no WRS path value";
	  break; 
	} else if (key.nval > 1) {
	  error_string = "too many WRS path values";
	  break; 
	}
        if (key.len_value[0] < 1) {
	  error_string = "invalid WRS path string";
	  break;
	}
	strncpy(temp, key.value[0], key.len_value[0]);
        temp[key.len_value[0]] = '\0';

	if (sscanf(temp, "%d", &this->meta.ipath) != 1) {
	  error_string = "converting WRS path";
	  break;
	}
	if (this->meta.ipath < 1)
	  error_string = "WRS path number out of range";
        break;

      case HEADER_WRS_ROW:
        if (key.nval <= 0) {
	  error_string = "no WRS row value";
	  break; 
	} else if (key.nval > 1) {
	  error_string = "too many WRS row values";
	  break; 
	}
        if (key.len_value[0] < 1) {
	  error_string = "invalid WRS row string";
	  break;
	}
	strncpy(temp, key.value[0], key.len_value[0]);
        temp[key.len_value[0]] = '\0';

	if (sscanf(temp, "%d", &this->meta.irow) != 1) {
	  error_string = "converting WRS row";
	  break;
	}
	if (this->meta.irow < 1)
	  error_string = "WRS row number out of range";
        break;

      case HEADER_NBAND:
        if (key.nval <= 0) {
	  error_string = "no number of bands value";
	  break; 
	} else if (key.nval > 1) {
	  error_string = "too many number of bands values";
	  break; 
	}
        if (key.len_value[0] < 1) {
	  error_string = "invalid number of bands string";
	  break;
	}
	strncpy(temp, key.value[0], key.len_value[0]);
        temp[key.len_value[0]] = '\0';

	if (sscanf(temp, "%d", &this->nband) != 1) {
	  error_string = "converting number of bands";
	  break;
	}
	if (this->nband < 1  ||
	    this->nband > NBAND_REFL_MAX)
	  error_string = "number of bands out of range";
        break;

      case HEADER_BANDS:
        if (key.nval <= 0) {
	  error_string = "no band number value";
	  break; 
	} else if (key.nval > NBAND_REFL_MAX) {
	  error_string = "too many band number values";
	  break; 
	}
	for (i = 0; i < key.nval; i++) {
          if (key.len_value[i] < 1) {
	    error_string = "invalid band number string";
	    break;
	  }
	  strncpy(temp, key.value[i], key.len_value[i]);
	  temp[key.len_value[i]] = '\0';

	  if (sscanf(temp, "%d", &this->meta.iband[i]) != 1) {
	    error_string = "converting band number";
	    break;
	  }
	  if (this->meta.iband[i] < 1) {
	    error_string = "band number out of range";
	    break;
	  }
	}
	if (error_string != (char *)NULL) break;
	nband_iband = key.nval;
        break;

      case HEADER_GAIN_SET:
        if (key.nval <= 0) {
	  error_string = "no gain setting value";
	  break; 
	} else if (key.nval > NBAND_REFL_MAX) {
	  error_string = "too many gain setting values";
	  break; 
	}
	for (i = 0; i < key.nval; i++) {
          this->meta.gain_set[i] = (Gain_t)
	    KeyString(key.value[i], key.len_value[i],
	              Gain_string, (int)GAIN_NULL, (int)GAIN_MAX);
          if (this->meta.inst == GAIN_NULL) {
            error_string = "invalid gain setting string";
	    break;
	  }
	}
	if (error_string != (char *)NULL) break;
	nband_gain_set = key.nval;
	flag_gain_set = true;
        break;

      case HEADER_GAIN:
        if (key.nval <= 0) {
	  error_string = "no gain values";
	  break; 
	} else if (key.nval > NBAND_REFL_MAX) {
	  error_string = "too many gain values";
	  break; 
	}
	for (i = 0; i < key.nval; i++) {
          if (key.len_value[i] < 1) {
	    error_string = "invalid gain string";
	    break;
	  }
	  strncpy(temp, key.value[i], key.len_value[i]);
          temp[key.len_value[i]] = '\0';

	  if (sscanf(temp, "%g", &this->meta.gain[i]) != 1) {
	    error_string = "converting gain";
	    break;
	  }
        }
	if (error_string != (char *)NULL) break;
	nband_gain = key.nval;
	flag_gain = true;
        break;

      case HEADER_BIAS:
        if (key.nval <= 0) {
	  error_string = "no bais values";
	  break; 
	} else if (key.nval > NBAND_REFL_MAX) {
	  error_string = "too many gain values";
	  break; 
	}
	for (i = 0; i < key.nval; i++) {
          if (key.len_value[i] < 1) {
	    error_string = "invalid bias string";
	    break;
	  }
	  strncpy(temp, key.value[i], key.len_value[i]);
          temp[key.len_value[i]] = '\0';

	  if (sscanf(temp, "%g", &this->meta.bias[i]) != 1) {
	    error_string = "converting bias";
	    break;
	  }
        }
	if (error_string != (char *)NULL) break;
	nband_bias = key.nval;
	flag_bias = true;
        break;

      case HEADER_NSAMPLE:
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

	if (sscanf(temp, "%d", &this->size.s) != 1) {
	  error_string = "converting number of samples";
	  break;
	}
	if (this->size.s < 1  ||
	    this->size.s > NSAMP_MAX)
	  error_string = "sample number out of range";
        break;

      case HEADER_NLINE:
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

	if (sscanf(temp, "%d", &this->size.l) != 1) {
	  error_string = "converting number of lines";
	  break;
	}
	if (this->size.l < 1  ||
	    this->size.l > NLINE_MAX)
	  error_string = "line number out of range";
        break;

      case HEADER_FILES:
        if (key.nval <= 0) {
	  error_string = "no file name value";
	  break; 
	} else if (key.nval > NBAND_REFL_MAX) {
	  error_string = "too many file name values";
	  break; 
	}
	for (i = 0; i < key.nval; i++) {
          if (key.len_value[i] < 1) {
	    error_string = "invalid file name string";
	    break;
	  }
	  strncpy(temp, key.value[i], key.len_value[i]);
          temp[key.len_value[i]] = '\0';

	  this->file_name[i] = DupString(temp);
          if (this->file_name[i] == (char *)NULL) {
            error_string = "converting file name string";
	    break;
	  }
	}
	if (error_string != (char *)NULL) break;
        nband_file_name = key.nval;
        break;

/******************************************************************/
/*********************** get thermal values ***********************/
/******************************************************************/

      case HEADER_NSAMPLE_TH:
        if (key.nval <= 0) {
	  error_string = "no number of thermal sample value";
	  break; 
	} else if (key.nval > 1) {
	  error_string = "too many number of thermal sample values";
	  break; 
	}
        if (key.len_value[0] < 1) {
	  error_string = "invalid number of thermal samples string";
	  break;
	}
	strncpy(temp, key.value[0], key.len_value[0]);
        temp[key.len_value[0]] = '\0';

	if (sscanf(temp, "%d", &this->size_th.s) != 1) {
	  error_string = "converting number of thermal samples";
	  break;
	}
	if (this->size_th.s < 1  ||
	    this->size_th.s > NSAMP_MAX)
	  error_string = "thermal sample number out of range";
        break;

      case HEADER_NLINE_TH:
        if (key.nval <= 0) {
	  error_string = "no number of thermal line value";
	  break; 
	} else if (key.nval > 1) {
	  error_string = "too many number of thremal line values";
	  break; 
	}
        if (key.len_value[0] < 1) {
	  error_string = "invalid number of thermal lines string";
	  break;
	}
	strncpy(temp, key.value[0], key.len_value[0]);
        temp[key.len_value[0]] = '\0';

	if (sscanf(temp, "%d", &this->size_th.l) != 1) {
	  error_string = "converting number of thremal lines";
	  break;
	}
	if (this->size_th.l < 1  ||
	    this->size_th.l > NLINE_MAX)
	  error_string = "thremal line number out of range";
        break;

      case HEADER_FILES_TH:
        if (key.nval <= 0) {
	  error_string = "no thermal file name value";
	  break; 
	} else if (key.nval > 1) {
	  error_string = "too many thermal file name values";
	  break; 
	}
          if (key.len_value[0] < 1) {
	    error_string = "invalid thermal file name string";
	    break;
	  }
	  strncpy(temp, key.value[0], key.len_value[0]);
          temp[key.len_value[0]] = '\0';

	  this->file_name_th = DupString(temp);
          if (this->file_name_th == (char *)NULL) {
            error_string = "converting thermal file name string";
	    break;
	  }

	if (error_string != (char *)NULL) break;
        break;

      case HEADER_NBAND_TH:
        if (key.nval <= 0) {
	  error_string = "no number of thermal bands value";
	  break; 
	} else if (key.nval > 1) {
	  error_string = "too many number of thermal bands values";
	  break; 
	}
        if (key.len_value[0] < 1) {
	  error_string = "invalid number of thermal bands string";
	  break;
	}
	strncpy(temp, key.value[0], key.len_value[0]);
        temp[key.len_value[0]] = '\0';

	if (sscanf(temp, "%d", &this->nband_th) != 1) {
	  error_string = "converting number of thermal bands";
	  break;
	}
	if (this->nband < 1  ||
	    this->nband > NBAND_REFL_MAX)
	  error_string = "number of thermal bands out of range";
        break;

      case HEADER_BANDS_TH:
        if (key.nval <= 0) {
	  error_string = "no thermal band number value";
	  break; 
	} else if (key.nval > 1) {
	  error_string = "too many thermal band number values";
	  break; 
	}
          if (key.len_value[0] < 1) {
	    error_string = "invalid thermal band number string";
	    break;
	  }
	  strncpy(temp, key.value[0], key.len_value[0]);
	  temp[key.len_value[0]] = '\0';

	  if (sscanf(temp, "%d", &this->meta.iband_th) != 1) {
	    error_string = "converting thermal band number";
	    break;
	  }
	  if (this->meta.iband_th < 1) {
	    error_string = "thermal band number out of range";
	    break;
	  }

	if (error_string != (char *)NULL) break;
        break;

      case HEADER_GAIN_SETTINGS_TH:
        if (key.nval <= 0) {
	  error_string = "no thermal gain values";
	  break; 
	} else if (key.nval > 1) {
	  error_string = "too many GAIN SETTING TH values";
	  break; 
	}

          if (key.len_value[0] < 1) {
	    error_string = "invalid thermal gain setting string";
	    break;
	  }
	  strncpy(temp, key.value[0], key.len_value[0]);
          temp[key.len_value[0]] = '\0';

          this->meta.gain_setting_th = (Gain_t)
	    KeyString(key.value[0], key.len_value[0],
	              Gain_string, (int)GAIN_NULL, (int)GAIN_MAX);
		      		      
	if (error_string != (char *)NULL) break;
        break;
	
      case HEADER_GAIN_TH:
        if (key.nval <= 0) {
	  error_string = "no thermal gain values";
	  break; 
	} else if (key.nval > 1) {
	  error_string = "too many thermal gain values";
	  break; 
	}

          if (key.len_value[0] < 1) {
	    error_string = "invalid thermal gain string";
	    break;
	  }
	  strncpy(temp, key.value[0], key.len_value[0]);
          temp[key.len_value[0]] = '\0';

/*        printf("gain_th string=(%s)\n",temp); */

	  if (sscanf(temp, "%g", &this->meta.gain_th) != 1) {
	    error_string = "converting thermal gain";
	    break;
	  }

	if (error_string != (char *)NULL) break;
        break;

      case HEADER_BIAS_TH:
        if (key.nval <= 0) {
	  error_string = "no bais values";
	  break; 
	} else if (key.nval > NBAND_REFL_MAX) {
	  error_string = "too many thermal bias values";
	  break; 
	}

          if (key.len_value[0] < 1) {
	    error_string = "invalid thermal bias string";
	    break;
	  }
	  strncpy(temp, key.value[0], key.len_value[0]);
          temp[key.len_value[0]] = '\0';

/*        printf("bias_th string=(%s)\n",temp); */
	  
	  if (sscanf(temp, "%g", &this->meta.bias_th) != 1) {
	    error_string = "converting thermal bias";
	    break;
	  }

	if (error_string != (char *)NULL) break;
        break;

/******************************************************************/
/******************** end of get thermal values *******************/
/******************************************************************/
      case HEADER_END:
        if (key.nval != 0) {
	  error_string = "no value expected";
	  break; 
	}
        got_end = true;
        break;

      default:
	/*
        error_string = "key not implmented";
	*/
	break;
    }
    if (error_string != (char *)NULL) break;
    if (got_end) break;
  }

  /* Handle errors */

  if (error_string != (char *)NULL) {
    for (ib = 0; ib < NBAND_REFL_MAX; ib++)
      if (this->file_name[ib] != (char *)NULL) free(this->file_name[ib]);
    free(this->file_header_name);
    RETURN_ERROR(error_string, "GetHeaderInput", false);
  }

  /* Make sure both the start and end of file exist */
  
  if (!got_start) 
    error_string = "no start key in header";
  else if (!got_end)
    error_string = "no end key in header";

  if (error_string != (char *)NULL) {
    for (ib = 0; ib < NBAND_REFL_MAX; ib++)
      if (this->file_name[ib] != (char *)NULL) free(this->file_name[ib]);
    free(this->file_header_name);
    RETURN_ERROR(error_string, "GetHeaderInput", false);
  }

  /* Handle null values */
  
  if (this->file_type == INPUT_TYPE_NULL) 
    error_string = "no file type given";
  else if (this->meta.provider == PROVIDER_NULL) 
    error_string = "no data provider given";
  else if (this->meta.sat == SAT_NULL) 
    error_string = "no satellite given";
  else if (this->meta.inst == INST_NULL) 
    error_string = "no instrument type given";
  else if (acq_date[0] == '\0') 
    error_string = "no acquisition date given";
  else if (this->meta.inst == INST_ETM  &&  prod_date[0] == '\0') 
    error_string = "no production date given";
  else if (this->meta.sun_zen == ANGLE_FILL) 
    error_string = "no solar zenith angle given";
  else if (this->meta.sun_az == ANGLE_FILL) 
    error_string = "no solar azimuth angle given";
  else if (this->meta.wrs_sys == WRS_FILL) 
    error_string = "no WRS system given";
  else if (this->meta.ipath == WRS_FILL) 
    error_string = "no number of paths given";
  else if (this->meta.irow == WRS_FILL) 
    error_string = "no number of rows given";
  else if (this->nband < 1) 
    error_string = "no number of bands given";
  else if (this->size.s < 1) 
    error_string = "no number of samples given";
  else if (this->size.l < 1) 
    error_string = "no number of lines given";
  else if (this->nband != nband_iband)
    error_string = "inconsistant number of values (band numbers)";
  else if (flag_gain_set  &&  this->nband != nband_gain_set)
    error_string = "inconsistant number of values (gain settings)";
  else if (flag_gain  &&  this->nband != nband_gain)
    error_string = "inconsistant number of values (gains)";
  else if (flag_bias  &&  this->nband != nband_bias)
    error_string = "inconsistant number of values (biases)";
  else if (this->nband != nband_file_name)
    error_string = "inconsistant number of values (file names)";

  if (error_string != (char *)NULL) {
    for (ib = 0; ib < NBAND_REFL_MAX; ib++)
      if (this->file_name[ib] != (char *)NULL) free(this->file_name[ib]);
    free(this->file_header_name);
    RETURN_ERROR(error_string, "GetHeaderInput", false);
  }

  for (ib = 0; ib < this->nband;    ib++) {
    if (this->meta.iband[ib] < 1) {
      error_string = "no band number given";
      break;
    }
    if (this->meta.inst == INST_ETM && !param->ETM_GB) {
      if (this->meta.gain_set[ib] == GAIN_NULL) {
        error_string = "no gain setting given";
        break;
      }
      if (this->meta.gain[ib] != GAIN_BIAS_FILL) {
        error_string = "no gain allowed";
        break;
      }
      if (this->meta.bias[ib] != GAIN_BIAS_FILL) {
        error_string = "no bias allowed";
        break;
      }
    } else if (!param->ETM_GB) {
      if (this->meta.gain_set[ib] != GAIN_NULL) {
        error_string = "no gain setting allowed";
        break;
      }
      if (this->meta.gain[ib] == GAIN_BIAS_FILL) {
        error_string = "no gain given";
        break;
      }
      if (this->meta.bias[ib] == GAIN_BIAS_FILL) {
        error_string = "no bias given";
        break;
      }
    }
    if (this->meta.inst == INST_ETM && param->ETM_GB) {
    this->meta.gain_th= this->meta.gain_th;
    }
    
    if (this->file_name[ib] == (char *)NULL) {
      error_string = "no file name given";
      break;
    }
  }
  
  if (error_string != (char *)NULL) {
    for (ib = 0; ib < this->nband;    ib++)
      if (this->file_name[ib] != (char *)NULL) free(this->file_name[ib]);
    free(this->file_header_name);
    RETURN_ERROR(error_string, "GetHeaderInput", false);
  }

  /* Check WRS path/rows */
  
  if (this->meta.wrs_sys == WRS_1) {
    if (this->meta.ipath > 251)
      error_string = "WRS path number out of range";
    else if (this->meta.irow > 248)
      error_string = "WRS row number out of range";
  } else if (this->meta.wrs_sys == WRS_2) {
    if (this->meta.ipath > 233)
      error_string = "WRS path number out of range";
    else if (this->meta.irow > 248)
      error_string = "WRS row number out of range";
  } else
    error_string = "invalid WRS system";

  if (error_string != (char *)NULL) {
    for (ib = 0; ib < this->nband;    ib++)
      if (this->file_name[ib] != (char *)NULL) free(this->file_name[ib]);
    free(this->file_header_name);
    RETURN_ERROR(error_string, "GetHeaderInput", false);
  }

  /* Check satellite/instrument combination */
  
  if (this->meta.inst == INST_MSS) {
    if (this->meta.sat != SAT_LANDSAT_1  &&  
        this->meta.sat != SAT_LANDSAT_2  &&
        this->meta.sat != SAT_LANDSAT_3  &&  
	this->meta.sat != SAT_LANDSAT_4  &&
        this->meta.sat != SAT_LANDSAT_5)
      error_string = "invalid insturment/satellite combination";
  } else if (this->meta.inst == INST_TM) {
    if (this->meta.sat != SAT_LANDSAT_4  &&  
        this->meta.sat != SAT_LANDSAT_5)
      error_string = "invalid insturment/satellite combination";
  } else if (this->meta.inst == INST_ETM) {
    if (this->meta.sat != SAT_LANDSAT_7)
      error_string = "invalid insturment/satellite combination";
  } else
    error_string = "invalid instrument type";

  if (error_string != (char *)NULL) {
    for (ib = 0; ib < this->nband;    ib++) 
      if (this->file_name[ib] != (char *)NULL) free(this->file_name[ib]);
    free(this->file_header_name);
    RETURN_ERROR(error_string, "GetHeaderInput", false);
  }

  /* Check band numbers */

  if (this->meta.inst == INST_MSS) {
    for (ib = 0; ib < this->nband; ib++)
      if (this->meta.iband[ib] > 4) {
        error_string = "band number out of range";
	break;
      }
  } else if (this->meta.inst == INST_TM  ||
             this->meta.inst == INST_ETM) {
    for (ib = 0; ib < this->nband; ib++)
      if (this->meta.iband[ib] == 6  ||  this->meta.iband[ib] > 7) {
        error_string = "band number out of range";
	break;
      }
  }
  
  if (error_string != (char *)NULL) {
    for (ib = 0; ib < this->nband;    ib++)
      if (this->file_name[ib] != (char *)NULL) free(this->file_name[ib]);
    free(this->file_header_name);
    RETURN_ERROR(error_string, "GetHeaderInput", false);
  }

  /* Convert the date/time values */

  if (acq_time[0] != '\0') 
    this->meta.time_fill = false;
  else
    strcpy(acq_time, "00:00:00");
    
  if (sprintf(temp, "%sT%sZ", acq_date, acq_time) < 0) 
    error_string = "converting acquisition date/time";
  else if (!DateInit(&this->meta.acq_date, temp, DATE_FORMAT_DATEA_TIME))
    error_string = "converting acquisition date/time";

  if (error_string == (char *)NULL  &&  
      prod_date[0] != '\0') {
    if (!DateInit(&this->meta.prod_date, prod_date, DATE_FORMAT_DATEA))
      error_string = "converting production date/time";
  }

  if (error_string != (char *)NULL) {
    for (ib = 0; ib < this->nband;    ib++)
      if (this->file_name[ib] != (char *)NULL) free(this->file_name[ib]);
    free(this->file_header_name);
    RETURN_ERROR(error_string, "GetHeaderInput", false);
  }

  return true;
}
