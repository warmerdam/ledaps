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
	FreeOutput - Free the 'input' data structure memory.

   2. 'OpenInput' must be called before any of the other routines.  
   3. 'FreeInput' should be used to free the 'input' data structure.
   4. The only input file type supported is HDF.

!END****************************************************************************
*/

#include <stdlib.h>
#include <math.h>
#include "input.h"
#include "error.h"
#include "mystring.h"
#include "myhdf.h"
#include "const.h"

#define SDS_PREFIX ("band")
#define SDS_CSM ("cloud_snow_mask")
#define INPUT_FILL (-9999)

#define INPUT_PROVIDER ("DataProvider")
#define INPUT_SAT ("Satellite")
#define INPUT_INST ("Instrument")
#define INPUT_ACQ_DATE ("AcquisitionDate")
#define INPUT_PROD_DATE ("Level1ProductionDate")
#define INPUT_SUN_ZEN ("SolarZenith")
#define INPUT_SUN_AZ ("SolarAzimuth")
#define INPUT_WRS_SYS ("WRS_System")
#define INPUT_WRS_PATH ("WRS_Path")
#define INPUT_WRS_ROW ("WRS_Row")
#define INPUT_NBAND ("NumberOfBands")
#define INPUT_BANDS ("BandNumbers")

#define N_LSAT_WRS1_ROWS  (251)
#define N_LSAT_WRS1_PATHS (233)
#define N_LSAT_WRS2_ROWS  (248)
#define N_LSAT_WRS2_PATHS (233)

/* Functions */

Input_t *OpenInput(char *file_name)
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
  char sds_name[40];
  int ir;
  bool sds_open[NBAND_REFL_MAX],qa_sds_open;
  Myhdf_dim_t *dim[2];
  int ib,read_qa;
  int16 *buf;

  /* Create the Input data structure */

  this = (Input_t *)malloc(sizeof(Input_t));
  if (this == (Input_t *)NULL) 
    RETURN_ERROR("allocating Input data structure", "OpenInput", 
                 (Input_t *)NULL);

  /* Populate the data structure */

  this->file_name = DupString(file_name);
  if (this->file_name == (char *)NULL) {
    free(this);
    RETURN_ERROR("duplicating file name", "OpenInput", (Input_t *)NULL);
  }

  /* Open file for SD access */

  this->sds_file_id = SDstart((char *)file_name, DFACC_RDONLY);
  if (this->sds_file_id == HDF_ERROR) {
    free(this->file_name);
    free(this);  
    RETURN_ERROR("opening input file", "OpenInput", 
                 (Input_t *)NULL); 
  }
  this->open = true;

  /* Get the input metadata */

  if (!GetInputMeta(this)) {
    free(this->file_name);
    free(this);  
    RETURN_ERROR("getting input metadata", "OpenInput", 
                 (Input_t *)NULL); 
  }

  if (this->nband == 1)
    read_qa=0;  /* do not read qa if thermal file */
  else
    read_qa=1;
  /* Get SDS information and start SDS access */

  for (ib = 0; ib < this->nband; ib++) {
    this->sds[ib].name = (char *)NULL;
    this->sds[ib].dim[0].name = (char *)NULL;
    this->sds[ib].dim[1].name = (char *)NULL;
    sds_open[ib] = false;
    this->buf[ib] = (int16 *)NULL;
  }
  this->qa_sds.name = (char *)NULL;
  this->qa_sds.dim[0].name = (char *)NULL;
  this->qa_sds.dim[1].name = (char *)NULL;
  qa_sds_open = false;
  this->qabuf = (int8 *)NULL;

  for (ib = 0; ib < this->nband; ib++) {

    if (sprintf(sds_name, "%s%d", SDS_PREFIX, this->meta.iband[ib]) < 0) {
      error_string = "creating SDS name";
      break;
    }

    this->sds[ib].name = DupString(sds_name);
    if (this->sds[ib].name == (char *)NULL) {
      error_string = "setting SDS name";
      break;
    }

    if (!GetSDSInfo(this->sds_file_id, &this->sds[ib])) {
      error_string = "getting sds info";
      break;
    }
    sds_open[ib] = true;

    /* Check rank */

    if (this->sds[ib].rank != 2) {
      error_string = "invalid rank";
      break;
    }

    /* Check SDS type */

    if (this->sds[ib].type != DFNT_INT16) {
      error_string = "invalid number type";
      break;
    }

    /* Get dimensions */

    for (ir = 0; ir < this->sds[ib].rank; ir++) {
      dim[ir] = &this->sds[ib].dim[ir];
      if (!GetSDSDimInfo(this->sds[ib].id, dim[ir], ir)) {
        error_string = "getting dimensions";
        break;
      }
    }

    if (error_string != (char *)NULL) break;

    /* Save and check line and sample dimensions */

    if (ib == 0) {
      this->size.l = dim[0]->nval;
      this->size.s = dim[1]->nval;
    } else {
      if (this->size.l != dim[0]->nval) {
        error_string = "all line dimensions do not match";
        break;
      }
      if (this->size.s != dim[1]->nval) {
        error_string = "all sample dimensions do not match";
        break;
      }
    }

  }


  if (read_qa) {
    strcpy(sds_name, "lndcal_QA");
    this->qa_sds.name = DupString(sds_name);
    if (this->qa_sds.name == (char *)NULL) {
      error_string = "setting qa SDS name";
    }

    if (!GetSDSInfo(this->sds_file_id, &this->qa_sds)) {
      error_string = "getting qa sds info";
    }
    qa_sds_open = true;

    /* Check rank */

    if (this->qa_sds.rank != 2) {
      error_string = "invalid qa rank";
    }

    /* Check SDS type */
    if (this->qa_sds.type != DFNT_UINT8) {
      error_string = "invalid qa number type";
    }

    /* Get dimensions */
    for (ir = 0; ir < this->qa_sds.rank; ir++) {
      dim[ir] = &this->qa_sds.dim[ir];
      if (!GetSDSDimInfo(this->qa_sds.id, dim[ir], ir)) {
        error_string = "getting qa dimensions";
        break;
      }
    }

    /* Save and check line and sample dimensions */

     if (this->size.l != dim[0]->nval) {
       error_string = "all line dimensions do not match";
     }
     if (this->size.s != dim[1]->nval) {
       error_string = "all sample dimensions do not match";
     }
    /* Allocate QA buffer */
   this->qabuf=(int8 *)calloc((size_t) this->size.s, sizeof(int8));
    if (this->qabuf == (int8 *)NULL)
      error_string = "allocating qa input buffer";
  }

  /* Allocate input buffer */

  if (sizeof(int16) != sizeof(int)) {
    buf = (int16 *)calloc((size_t)(this->size.s * this->nband), sizeof(int16));
    if (buf == (int16 *)NULL)
      error_string = "allocating input buffer";
    else {
      this->buf[0] = buf;
      for (ib = 1; ib < this->nband; ib++)
        this->buf[ib] = this->buf[ib - 1] + this->size.s;
    }
  }

  if (error_string != (char *)NULL) {
    for (ib = 0; ib < this->nband; ib++) {
      for (ir = 0; ir < this->sds[ib].rank; ir++) {
        if (this->sds[ib].dim[ir].name != (char *)NULL)
          free(this->sds[ib].dim[ir].name);
      }
      if (sds_open[ib])
        SDendaccess(this->sds[ib].id);
      if (this->sds[ib].name != (char *)NULL)
        free(this->sds[ib].name);
      if (this->buf[ib] != (int16 *)NULL)
        free(this->buf[ib]);

    }
	if (qa_sds_open)
	    SDendaccess(this->qa_sds.id);
    for (ir = 0; ir < this->qa_sds.rank; ir++) {
	    if (this->qa_sds.dim[ir].name != (char *)NULL)
          free(this->qa_sds.dim[ir].name);
    }
    if (this->qa_sds.name != (char *)NULL)
        free(this->qa_sds.name);
    if (this->qabuf != (int8 *)NULL)
        free(this->qabuf);

    SDend(this->sds_file_id);
    free(this->file_name);
    free(this);
    RETURN_ERROR(error_string, "OpenInput", (Input_t *)NULL);
  }

  return this;
}

InputMask_t *OpenInputMask(char *file_name)
{
  InputMask_t *this;
  char *error_string = (char *)NULL;
  int ir;
  bool sds_open[NBAND_REFL_MAX];
  Myhdf_dim_t *dim[2];
  uint8 *buf;

  /* Create the Input data structure */

  this = (InputMask_t *)malloc(sizeof(InputMask_t));
  if (this == (InputMask_t *)NULL) 
    RETURN_ERROR("allocating Input data structure", "OpenInput", 
                 (InputMask_t *)NULL);

  /* Populate the data structure */

  this->file_name = DupString(file_name);
  if (this->file_name == (char *)NULL) {
    free(this);
    RETURN_ERROR("duplicating file name", "OpenInput", (InputMask_t *)NULL);
  }
  /* Open file for SD access */

  this->sds_file_id = SDstart((char *)file_name, DFACC_RDONLY);
  if (this->sds_file_id == HDF_ERROR) {
    free(this->file_name);
    free(this);  
    RETURN_ERROR("opening input file", "OpenInput", 
                 (InputMask_t *)NULL); 
  }
  this->open = true;

  /* Get SDS information and start SDS access */

  this->sds[0].name = (char *)NULL;
  this->sds[0].dim[0].name = (char *)NULL;
  this->sds[0].dim[1].name = (char *)NULL;
  sds_open[0] = false;
  this->buf[0] = (uint8 *)NULL;

  this->sds[0].name = DupString(SDS_CSM);
  if (this->sds[0].name == (char *)NULL) {
    error_string = "setting SDS name";
  }

  if (!GetSDSInfo(this->sds_file_id, &this->sds[0])) {
    error_string = "getting sds info";
  }
  sds_open[0] = true;

    /* Check rank */

  if (this->sds[0].rank != 2) {
    error_string = "invalid rank";
  }

    /* Check SDS type */

  if (this->sds[0].type != DFNT_UINT8) {
    error_string = "invalid number type";
  }

    /* Get dimensions */

  for (ir = 0; ir < this->sds[0].rank; ir++) {
    dim[ir] = &this->sds[0].dim[ir];
    if (!GetSDSDimInfo(this->sds[0].id, dim[ir], ir)) {
      error_string = "getting dimensions";
      break;
    }
  }

    /* Save and check line and sample dimensions */

    this->size.l = dim[0]->nval;
    this->size.s = dim[1]->nval;

  /* Allocate input buffer */

  if (sizeof(uint8) != sizeof(unsigned char)) {
    buf = (uint8 *)calloc((size_t)(this->size.s), sizeof(uint8));
    if (buf == (uint8 *)NULL)
      error_string = "allocating input buffer";
    else {
      this->buf[0] = buf;
    }
  }

  if (error_string != (char *)NULL) {
    for (ir = 0; ir < this->sds[0].rank; ir++) {
      if (this->sds[0].dim[ir].name != (char *)NULL)
        free(this->sds[0].dim[ir].name);
      }
    if (sds_open[0])
      SDendaccess(this->sds[0].id);
    if (this->sds[0].name != (char *)NULL)
      free(this->sds[0].name);
    if (this->buf[0] != (uint8 *)NULL)
      free(this->buf[0]);

    SDend(this->sds_file_id);
    free(this->file_name);
    free(this);
    RETURN_ERROR(error_string, "OpenInputMask", (InputMask_t *)NULL);
  }

  return this;
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

  if (!this->open)
    RETURN_ERROR("file not open", "CloseInput", false);

  for (ib = 0; ib < this->nband; ib++) {
    if (SDendaccess(this->sds[ib].id) == HDF_ERROR) 
      RETURN_ERROR("ending sds access", "CloseInput", false);
  }

  if (SDendaccess(this->qa_sds.id) == HDF_ERROR)
      RETURN_ERROR("ending qa sds access", "CloseInput", false);

  SDend(this->sds_file_id);
  this->open = false;

  return true;
}

bool CloseInputMask(InputMask_t *this)
{
  if (!this->open)
    RETURN_ERROR("file not open", "CloseInputMask", false);

  if (SDendaccess(this->sds[0].id) == HDF_ERROR) 
    RETURN_ERROR("ending sds access", "CloseInputMask", false);

  SDend(this->sds_file_id);
  this->open = false;

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
  int ib, ir;

  if (this != (Input_t *)NULL) {
    for (ib = 0; ib < this->nband; ib++) {
      for (ir = 0; ir < this->sds[ib].rank; ir++) {
        if (this->sds[ib].dim[ir].name != (char *)NULL) 
          free(this->sds[ib].dim[ir].name);
      }
      if (this->sds[ib].name != (char *)NULL) 
        free(this->sds[ib].name);
    }
    for (ir = 0; ir < this->qa_sds.rank; ir++) {
      if (this->qa_sds.dim[ir].name != (char *)NULL) 
        free(this->qa_sds.dim[ir].name);
    }
    if (this->qa_sds.name != (char *)NULL) 
      free(this->qa_sds.name);

    if (this->buf[0] != (int16 *)NULL)
      free(this->buf[0]);
    if (this->qabuf != (int8 *)NULL)
      free(this->qabuf);

    if (this->file_name != (char *)NULL) free(this->file_name);
    free(this);
  }

  return true;
}

bool FreeInputMask(InputMask_t *this)
{
  int ir;

  if (this != (InputMask_t *)NULL) {

      for (ir = 0; ir < this->sds[0].rank; ir++) {
        if (this->sds[0].dim[ir].name != (char *)NULL) 
          free(this->sds[0].dim[ir].name);
      }
      if (this->sds[0].name != (char *)NULL) 
        free(this->sds[0].name);

    if (this->buf[0] != (uint8 *)NULL)
      free(this->buf[0]);
    if (this->file_name != (char *)NULL) free(this->file_name);
    free(this);
  }

  return true;
}

bool GetInputLine(Input_t *this, int iband, int iline, int *line)
/* 
!C******************************************************************************

!Description: 'WriteOutput' writes a line of data to the output file.
 
!Input Parameters:
 this           'output' data structure; the following fields are input:
                   open, size, sds.id
 iline          output line number
 buf            buffer of data to be written

!Output Parameters:
 this           'output' data structure; the following fields are modified:
 (returns)      status:
                  'true' = okay
		  'false' = error return

!Team Unique Header:

 ! Design Notes:
   1. An error status is returned when:
       a. the file is not open for access
       b. the line number is invalid (< 0; >= 'this->size.l')
       b. an error occurs when writting to the SDS.
   2. Error messages are handled with the 'RETURN_ERROR' macro.
   3. 'OutputFile' must be called before this routine is called.

!END****************************************************************************
*/
{
  int32 start[MYHDF_MAX_RANK], nval[MYHDF_MAX_RANK];
  void *buf;
  int is;

  /* Check the parameters */

  if (this == (Input_t *)NULL) 
    RETURN_ERROR("invalid input structure", "GetIntputLine", false);
  if (!this->open)
    RETURN_ERROR("file not open", "GetInputLine", false);
  if (iband < 0  ||  iband >= this->nband)
    RETURN_ERROR("invalid band number", "GetInputLine", false);
  if (iline < 0  ||  iline >= this->size.l)
    RETURN_ERROR("invalid line number", "GetInputLine", false);

  if (sizeof(int16) == sizeof(int))
    buf = (void *)line;
  else
    buf = (void *)this->buf[iband];

  /* Read the data */

  start[0] = iline;
  start[1] = 0;
  nval[0] = 1;
  nval[1] = this->size.s;

  if (SDreaddata(this->sds[iband].id, start, NULL, nval, buf) == HDF_ERROR)
    RETURN_ERROR("reading input", "GetInputLine", false);

  if (sizeof(int16) != sizeof(int)) {
    for (is = 0; is < this->size.s; is++) 
      line[is] = (int)this->buf[iband][is];
  }

  return true;
}

bool GetInputQALine(Input_t *this, int iline, int8 *line)
{
  int32 start[MYHDF_MAX_RANK], nval[MYHDF_MAX_RANK];
  void *buf;
  int is;

  /* Check the parameters */

  if (this == (Input_t *)NULL) 
    RETURN_ERROR("invalid input structure", "GetIntputLine", false);
  if (!this->open)
    RETURN_ERROR("file not open", "GetInputLine", false);
  if (iline < 0  ||  iline >= this->size.l)
    RETURN_ERROR("invalid line number", "GetInputLine", false);

  buf = (void *)line;

  /* Read the data */

  start[0] = iline;
  start[1] = 0;
  nval[0] = 1;
  nval[1] = this->size.s;

  if (SDreaddata(this->qa_sds.id, start, NULL, nval, buf) == HDF_ERROR)
    RETURN_ERROR("reading input", "GetInputLine", false);

  return true;
}

bool GetInputMaskLine(InputMask_t *this, int iline, char *line)
{
  int32 start[MYHDF_MAX_RANK], nval[MYHDF_MAX_RANK];
  void *buf;
  int is;

  /* Check the parameters */

  if (this == (InputMask_t *)NULL) 
    RETURN_ERROR("invalid input structure", "GetIntputLine", false);
  if (!this->open)
    RETURN_ERROR("file not open", "GetInputLine", false);
  if (iline < 0  ||  iline >= this->size.l)
    RETURN_ERROR("invalid line number", "GetInputLine", false);

  if (sizeof(uint8) == sizeof(unsigned char)) 
    buf = (void *)line;
  else
    buf = (void *)this->buf[0];

  /* Read the data */

  start[0] = iline;
  start[1] = 0;
  nval[0] = 1;
  nval[1] = this->size.s;


  if (SDreaddata(this->sds[0].id, start, NULL, nval, buf) == HDF_ERROR)
    RETURN_ERROR("reading input", "GetInputMaskLine", false);

  if (sizeof(uint8) != sizeof(unsigned char)) {
    for (is = 0; is < this->size.s; is++) 
      line[is] = (int)this->buf[0][is];
  }

  return true;
}

bool GetInputMeta(Input_t *this) 
{
  Myhdf_attr_t attr;
  double dval[NBAND_REFL_MAX];
  char string[MAX_STR_LEN + 1];
  char date[MAX_DATE_LEN + 1];
  int ib;
  Input_meta_t *meta;
  char *error_string;

  /* Check the parameters */

  if (!this->open)
    RETURN_ERROR("file not open", "GetInputMeta", false);

  meta = &this->meta;
  meta->fill = INPUT_FILL;

  /* Read the metadata */

  attr.type = DFNT_CHAR8;
  attr.nval = MAX_STR_LEN;
  attr.name = INPUT_PROVIDER;
  if (!GetAttrString(this->sds_file_id, &attr, string))
    RETURN_ERROR("reading attribute (data provider)", "GetInputMeta", false);
  meta->provider = (Provider_t)KeyString(string, attr.nval, Provider_string, 
		                         (int)PROVIDER_NULL, (int)PROVIDER_MAX);
  if (meta->provider == PROVIDER_NULL)
    RETURN_ERROR("invalid input metadata (data provider)", 
                 "GetInputMeta", false);

  attr.type = DFNT_CHAR8;
  attr.nval = MAX_STR_LEN;
  attr.name = INPUT_SAT;
  if (!GetAttrString(this->sds_file_id, &attr, string))
    RETURN_ERROR("reading attribute (instrument)", "GetInputMeta", false);
  meta->sat = (Sat_t)KeyString(string, attr.nval, Sat_string, 
		               (int)SAT_NULL, (int)SAT_MAX);
  if (meta->sat == SAT_NULL)
    RETURN_ERROR("invalid input metadata (satellite)", "GetInputMeta", false);

  attr.type = DFNT_CHAR8;
  attr.nval = MAX_STR_LEN;
  attr.name = INPUT_INST;
  if (!GetAttrString(this->sds_file_id, &attr, string))
    RETURN_ERROR("reading attribute (instrument)", "GetInputMeta", false);
  meta->inst = (Inst_t)KeyString(string, attr.nval, Inst_string, 
		                 (int)INST_NULL, (int)INST_MAX);
  if (meta->inst == INST_NULL)
    RETURN_ERROR("invalid input metadata (instrument)", "GetInputMeta", false);

  attr.type = DFNT_CHAR8;
  attr.nval = MAX_DATE_LEN;
  attr.name = INPUT_ACQ_DATE;
  if (!GetAttrString(this->sds_file_id, &attr, date))
    RETURN_ERROR("reading attribute (acquisition date)", "GetInputMeta", false);
  if (!DateInit(&meta->acq_date, date, DATE_FORMAT_DATEA_TIME))
    RETURN_ERROR("converting acquisition date", "GetInputMeta", false);

  attr.type = DFNT_CHAR8;
  attr.nval = MAX_DATE_LEN;
  attr.name = INPUT_PROD_DATE;
  if (!GetAttrString(this->sds_file_id, &attr, date))
    RETURN_ERROR("reading attribute (production date)", "GetInputMeta", false);
  if (!DateInit(&meta->prod_date, date, DATE_FORMAT_DATEA_TIME))
    RETURN_ERROR("converting production date", "GetInputMeta", false);

  attr.type = DFNT_FLOAT32;
  attr.nval = 1;
  attr.name = INPUT_SUN_ZEN;
  if (!GetAttrDouble(this->sds_file_id, &attr, dval))
    RETURN_ERROR("reading attribute (solar zenith)", "GetInputMeta", false);
  if (attr.nval != 1) 
    RETURN_ERROR("invalid number of values (solar zenith)", 
                  "GetInputMeta", false);
  if (dval[0] < -90.0  ||  dval[0] > 90.0)
    RETURN_ERROR("solar zenith angle out of range", "GetInputMeta", false);
  meta->sun_zen = (float)(dval[0] * RAD);

  attr.type = DFNT_FLOAT32;
  attr.nval = 1;
  attr.name = INPUT_SUN_AZ;
  if (!GetAttrDouble(this->sds_file_id, &attr, dval))
    RETURN_ERROR("reading attribute (solar azimuth)", "GetInputMeta", false);
  if (attr.nval != 1) 
    RETURN_ERROR("invalid number of values (solar azimuth)", 
                 "GetInputMeta", false);
  if (dval[0] < -360.0  ||  dval[0] > 360.0)
    RETURN_ERROR("solar azimuth angle out of range", "GetInputMeta", false);
  meta->sun_az = (float)(dval[0] * RAD);

  attr.type = DFNT_CHAR8;
  attr.nval = MAX_STR_LEN;
  attr.name = INPUT_WRS_SYS;
  if (!GetAttrString(this->sds_file_id, &attr, string))
    RETURN_ERROR("reading attribute (WRS system)", "GetInputMeta", false);
  meta->wrs_sys = (Wrs_t)KeyString(string, attr.nval, Wrs_string, 
		                       (int)WRS_NULL, (int)WRS_MAX);
  if (meta->wrs_sys == WRS_NULL)
    RETURN_ERROR("invalid input metadata (WRS sytem)", "GetInputMeta", false);

  attr.type = DFNT_INT16;
  attr.nval = 1;
  attr.name = INPUT_WRS_PATH;
  if (!GetAttrDouble(this->sds_file_id, &attr, dval))
    RETURN_ERROR("reading attribute (WRS path)", "GetInputMeta", false);
  if (attr.nval != 1) 
    RETURN_ERROR("invalid number of values (WRS path)", "GetInputMeta", false);
  meta->ipath = (int)floor(dval[0] + 0.5);
  if (meta->ipath < 1) 
    RETURN_ERROR("WRS path out of range", "GetInputMeta", false);

  attr.type = DFNT_INT16;
  attr.nval = 1;
  attr.name = INPUT_WRS_ROW;
  if (!GetAttrDouble(this->sds_file_id, &attr, dval))
    RETURN_ERROR("reading attribute (WRS row)", "GetInputMeta", false);
  if (attr.nval != 1) 
    RETURN_ERROR("invalid number of values (WRS row)", "GetInputMeta", false);
  meta->irow = (int)floor(dval[0] + 0.5);
  if (meta->irow < 1) 
    RETURN_ERROR("WRS path out of range", "GetInputMeta", false);

  attr.type = DFNT_INT8;
  attr.nval = 1;
  attr.name = INPUT_NBAND;
  if (!GetAttrDouble(this->sds_file_id, &attr, dval))
    RETURN_ERROR("reading attribute (number of bands)", "GetInputMeta", false);
  if (attr.nval != 1) 
    RETURN_ERROR("invalid number of values (number of bands)", 
                 "GetInputMeta", false);
  this->nband = (int)floor(dval[0] + 0.5);
  if (this->nband < 1  ||  this->nband > NBAND_REFL_MAX) 
    RETURN_ERROR("number of bands out of range", "GetInputMeta", false);

  attr.type = DFNT_INT8;
  attr.nval = this->nband;
  attr.name = INPUT_BANDS;
  if (!GetAttrDouble(this->sds_file_id, &attr, dval))
    RETURN_ERROR("reading attribute (band numbers)", "GetInputMeta", false);
  if (attr.nval != this->nband) 
    RETURN_ERROR("invalid number of values (band numbers)", 
                 "GetInputMeta", false);
  for (ib = 0; ib < this->nband; ib++) {
    meta->iband[ib] = (int)floor(dval[ib] + 0.5);
    if (meta->iband[ib] < 1)
      RETURN_ERROR("band number out of range", "GetInputMeta", false);
  }

  /* Check WRS path/rows */

  error_string = (char *)NULL;

  if (meta->wrs_sys == WRS_1) {
    if (meta->ipath > N_LSAT_WRS1_PATHS)
      error_string = "WRS path number out of range";
    else if (meta->irow > N_LSAT_WRS1_ROWS)
      error_string = "WRS row number out of range";
  } else if (meta->wrs_sys == WRS_2) {
    if (meta->ipath > N_LSAT_WRS2_PATHS)
      error_string = "WRS path number out of range";
    else if (meta->irow > N_LSAT_WRS2_ROWS)
      error_string = "WRS row number out of range";
  } else
    error_string = "invalid WRS system";

  if (error_string != (char *)NULL)
    RETURN_ERROR(error_string, "GetInputMeta", false);

  /* Check satellite/instrument combination */
  
  if (meta->inst == INST_MSS) {
    if (meta->sat != SAT_LANDSAT_1  &&  
        meta->sat != SAT_LANDSAT_2  &&
        meta->sat != SAT_LANDSAT_3  &&  
	  meta->sat != SAT_LANDSAT_4  &&
        meta->sat != SAT_LANDSAT_5)
      error_string = "invalid insturment/satellite combination";
  } else if (meta->inst == INST_TM) {
    if (meta->sat != SAT_LANDSAT_4  &&  
        meta->sat != SAT_LANDSAT_5)
      error_string = "invalid insturment/satellite combination";
  } else if (meta->inst == INST_ETM) {
    if (meta->sat != SAT_LANDSAT_7)
      error_string = "invalid insturment/satellite combination";
  } else
    error_string = "invalid instrument type";

  if (error_string != (char *)NULL)
    RETURN_ERROR(error_string, "GetInputMeta", false);

  /* Check band numbers */

  if (meta->inst == INST_MSS) {
    for (ib = 0; ib < this->nband; ib++)
      if (meta->iband[ib] > 4) {
        error_string = "band number out of range";
	break;
      }
  } else if (meta->inst == INST_TM  ||
             meta->inst == INST_ETM) {
    for (ib = 0; ib < this->nband; ib++)
      if ( meta->iband[ib] > 7) {
        error_string = "band number out of range";
	break;
      }
  } else
    error_string = "invalid instrument type";

  if (error_string != (char *)NULL)
    RETURN_ERROR(error_string, "GetInputMeta", false);

  return true;
}

bool InputMetaCopy(Input_meta_t *this, int nband, Input_meta_t *copy) 
{
  int ib;

  if (this == (Input_meta_t *)NULL) 
    RETURN_ERROR("invalid input structure (original)", "InputMetaCopy", false);
  if (copy == (Input_meta_t *)NULL) 
    RETURN_ERROR("invalid input structure (copy)", "InputMetaCopy", false);
  if (nband < 1  ||  nband > NBAND_REFL_MAX)
    RETURN_ERROR("invalid number of bands", "InputMetaCopy", false);

  copy->provider = this->provider;
  copy->sat = this->sat;
  copy->inst = this->inst;
  if (!DateCopy(&this->acq_date, &copy->acq_date)) 
    RETURN_ERROR("copying acquisition date/time", "InputMetaCopy", false);
  if (!DateCopy(&this->prod_date, &copy->prod_date)) 
    RETURN_ERROR("copying production date/time", "InputMetaCopy", false);
  copy->sun_zen = this->sun_zen;
  copy->sun_az = this->sun_az;
  copy->wrs_sys = this->wrs_sys;
  copy->ipath = this->ipath;
  copy->irow = this->irow;
  copy->fill = this->fill;

  for (ib = 0; ib < nband; ib++)
    copy->iband[ib] = this->iband[ib];

  return true;
}
