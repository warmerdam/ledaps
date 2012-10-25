/*
!C****************************************************************************

!File: output.c
  
!Description: Functions creating and writting data to the product output file.

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
   1. The following public functions handle the output files:

        CreateOutput - Create new output file.
	OutputFile - Setup 'output' data structure.
	CloseOutput - Close the output file.
	FreeOutput - Free the 'output' data structure memory.
	PutMetadata - Write the output product metadata.
	WriteOutput - Write a line of data to the output product file.

   2. 'OutputFile' must be called before any of the other routines (except for 
      'CreateOutput').  
   3. 'FreeOutput' should be used to free the 'output' data structure.
   4. The only output file type supported is HDF.

!END****************************************************************************
*/

#include <stdlib.h>
#include "output.h"
#include "input.h"
#include "names.h"
#include "error.h"

#define OUTPUT_PROVIDER ("DataProvider")
#define OUTPUT_SAT ("Satellite")
#define OUTPUT_INST ("Instrument")
#define OUTPUT_ACQ_DATE ("AcquisitionDate")
#define OUTPUT_L1_PROD_DATE ("Level1ProductionDate")
#define OUTPUT_SUN_ZEN ("SolarZenith")
#define OUTPUT_SUN_AZ ("SolarAzimuth")
#define OUTPUT_WRS_SYS ("WRS_System")
#define OUTPUT_WRS_PATH ("WRS_Path")
#define OUTPUT_WRS_ROW ("WRS_Row")
#define OUTPUT_NBAND ("NumberOfBands")
#define OUTPUT_BANDS ("BandNumbers")
#define OUTPUT_HIGH_COUNT ("high_count")
#define OUTPUT_LOW_COUNT ("low_count")
#define OUTPUT_SHORT_NAME ("ShortName")
#define OUTPUT_LOCAL_GRAN_ID ("LocalGranuleID")
#define OUTPUT_PROD_DATE ("ProductionDate")
#define OUTPUT_PGEVERSION ("PGEVersion")
#define OUTPUT_PROCESSVERSION ("ProcessVersion")
#define DIM_CSM_LINE ("grid_line")
#define DIM_CSM_SAMP ("grid_sample")

#define OUTPUT_LONG_NAME        ("long_name")
#define OUTPUT_UNITS            ("units")
#define OUTPUT_VALID_RANGE      ("valid_range")
#define OUTPUT_FILL_VALUE       ("_FillValue")
#define OUTPUT_SCALE_FACTOR     ("scale_factor")
#define OUTPUT_ADD_OFFSET       ("add_offset")
#define OUTPUT_SCALE_FACTOR_ERR ("scale_factor_err")
#define OUTPUT_ADD_OFFSET_ERR   ("add_offset_err")
#define OUTPUT_CALIBRATED_NT    ("calibrated_nt")

#define OUTPUT_FILL_VALUE       ("_FillValue")
#define OUTPUT_LAND_VALUE       ("_LandValue")
#define OUTPUT_CLOUD_VALUE      ("_CloudValue")
#define OUTPUT_SNOW_VALUE       ("_SnowValue")
#define OUTPUT_SHADOW_VALUE     ("_ShadowValue")
#define OUTPUT_WATER_VALUE      ("_WaterValue")

bool CreateOutput(char *file_name)
/* 
!C******************************************************************************

!Description: 'CreateOuptut' creates a new output file.
 
!Input Parameters:
 file_name      output file name

!Output Parameters:
 (returns)      status:
                  'true' = okay
		  'false' = error return

!Team Unique Header:

 ! Design Notes:
   1. An error status is returned when:
       a. the creation of the output file failes.
   2. Error messages are handled with the 'RETURN_ERROR' macro.
   3. The output file is in HDF format and closed after it is created.

!END****************************************************************************
*/
{
  int32 hdf_file_id;

  /* Create the file with HDF open */

  hdf_file_id = Hopen(file_name, DFACC_CREATE, DEF_NDDS); 
  if(hdf_file_id == HDF_ERROR) {
    RETURN_ERROR("creating output file", "CreateOutput", false); 
  }

  /* Close the file */

  Hclose(hdf_file_id);

  return true;
}

/* 
!C******************************************************************************
Output_t *OpenOutput(char *file_name, Img_coord_int_t *size)
!C******************************************************************************

!Description: 'OutputFile' sets up the 'output' data structure, opens the
 output file for write access, and creates the output Science Data Set (SDS).
 
!Input Parameters:thi
 file_name      output file name
 nband          number of bands (SDSs) to be created
 sds_name       names of SDSs to be created
 size           SDS size

!Output Parameters:
 (returns)      'output' data structure or NULL when an error occurs

!Team Unique Header:

 ! Design Notes:
   1. When 'OutputFile' returns sucessfully, the file is open for HDF access 
      and the SDS is open for access.
   2. An error status is returned when:
       a. the output image dimensions are zero or negative
       b. an invalid output data type is passed
       c. memory allocation is not successful
       d. duplicating strings is not successful
       e. errors occur when opening the output HDF file
       f. errors occur when creating and initializing the SDS.
   3. Error messages are handled with the 'RETURN_ERROR' macro.
   4. 'FreeOutput' should be called to deallocate memory used by the 
      'output' data structures.
   5. 'CloseFile' should be called after all of the data is written and 
      before the 'output' data structure memory is released.

!END****************************************************************************
*/
Output_t *OpenOutput(char *file_name, Img_coord_int_t *size)
{
  Output_t *this;
  char *error_string = (char *)NULL;
  Myhdf_dim_t *dim[MYHDF_MAX_RANK];
  Myhdf_sds_t *sds;
  int ir, ib;

  int nband_csm_setup;
  char *sds_name_csm[NBAND_CSM] = {"cloud_snow_mask"};

  /* Check parameters */
  
  if (size->l < 1)
    RETURN_ERROR("invalid number of output lines", 
                 "OpenOutput", (Output_t *)NULL);

  if (size->s < 1)
    RETURN_ERROR("invalid number of samples per output line", 
                 "OpenOutput", (Output_t *)NULL);

  /* Create the Output data structure */

  this = (Output_t *)malloc(sizeof(Output_t));
  if (this == (Output_t *)NULL) 
    RETURN_ERROR("allocating Output data structure", "OpenOutput", 
                 (Output_t *)NULL);

  /* Populate the data structure */

  this->file_name = DupString(file_name);
  if (this->file_name == (char *)NULL) {
    free(this);
    RETURN_ERROR("duplicating file name", "OpenOutput", (Output_t *)NULL);
  }

  this->open = false;
  this->size.l = size->l;
  this->size.s = size->s;

  for (ib = 0; ib < NBAND_CSM; ib++) {
    this->sds_csm[ib].name = (char *)NULL;
    this->sds_csm[ib].dim[0].name = (char *)NULL;
    this->sds_csm[ib].dim[1].name = (char *)NULL;
  }

  /* Open file for SD access */

  this->sds_file_id = SDstart((char *)file_name, DFACC_RDWR);
  if (this->sds_file_id == HDF_ERROR) {
    free(this->file_name);
    free(this);  
    RETURN_ERROR("opening output file for SD access", "OpenOutput", 
                 (Output_t *)NULL); 
  }
  this->open = true;

  /* Set up the KT SDSs */

/*****************************************************************/
/****            new            new            new            ****/
/*****************************************************************/

  nband_csm_setup = 0;
  for (ib = 0; ib < NBAND_CSM; ib++) {
    sds = &this->sds_csm[ib];
    sds->type = DFNT_UINT8;
    sds->rank = 2;
    sds->name = DupString(sds_name_csm[ib]);
    if (sds->name == (char *)NULL) {
      error_string = "duplicating sds name";
      break;
    }

    dim[0] = &sds->dim[0];
    dim[1] = &sds->dim[1];

    dim[0]->nval = this->size.l;
    dim[1]->nval = this->size.s;

    dim[0]->type = dim[1]->type = sds->type;

    dim[0]->name = DupString(DIM_CSM_LINE);
    if (dim[0]->name == (char *)NULL) {
      error_string = "duplicating dim name (l)";
      break;
    }

    dim[1]->name = DupString(DIM_CSM_SAMP);
    if (dim[1]->name == (char *)NULL) {
      error_string = "duplicating dim name (s)";
      break;
    }

    if (!PutSDSInfo(this->sds_file_id, sds)) {
      error_string = "setting up the SDS";
      break;
    }

    for (ir = 0; ir < sds->rank; ir++) {
      if (!PutSDSDimInfo(sds->id, dim[ir], ir)) {
        error_string = "setting up the dimension";
        break;
      }
    }
    if (error_string != (char *)NULL) break;

    nband_csm_setup++;
  }

  if (error_string == (char *)NULL) {
    if (sizeof(uint8) != sizeof(unsigned char)) {
      for (ib = 0; ib < NBAND_CSM; ib++) {
        this->buf_csm[ib] = (uint8 *)calloc(this->size.s, sizeof(uint8));
        if (this->buf_csm[ib] == (uint8 *)NULL) {
          error_string = "allocating input buffer";
          break;
	}
      }
    }
  }



  if (error_string != (char *)NULL) {

    for (ib = 0; ib < NBAND_CSM; ib++) {
      sds = &this->sds_csm[ib];
      if (sds->name != (char *)NULL) free(sds->name);
      if (ib < nband_csm_setup) SDendaccess(sds->id);
      if ( sizeof(uint8) != sizeof(unsigned char) && 
           this->buf_csm[ib] != (uint8 *)NULL        ) free(this->buf_csm[ib]);
      if (sds->dim[0].name != (char *)NULL) free(sds->dim[0].name);
      if (sds->dim[1].name != (char *)NULL) free(sds->dim[1].name);
    }    
    SDend(this->sds_file_id);
    free(this->file_name);
    free(this);  
    RETURN_ERROR(error_string, "OpenOutput", (Output_t *)NULL); 

  }

  return this;
}


/* 
!C******************************************************************************
bool CloseOutput(Output_t *this)
!C******************************************************************************

!Description: 'CloseOutput' ends SDS access and closes the output file.
 
!Input Parameters:
 this           'output' data structure; the following fields are input:
                   open, sds.id, sds_file_id

!Output Parameters:
 this           'output' data structure; the following fields are modified:
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
   3. 'OutputFile' must be called before this routine is called.
   4. 'FreeOutput' should be called to deallocate memory used by the 
      'output' data structures.

!END****************************************************************************
*/
bool CloseOutput(Output_t *this)
{
  int ib;

  if (!this->open)
    RETURN_ERROR("file not open", "CloseOutput", false);

  for (ib = 0; ib < NBAND_CSM; ib++) 
    if (SDendaccess(this->sds_csm[ib].id) == HDF_ERROR) 
      RETURN_ERROR("ending sds access", "CloseOutput", false);
  SDend(this->sds_file_id);
  this->open = false;

  return true;
}


/* 
!C******************************************************************************
bool FreeOutput(Output_t *this)
!C******************************************************************************

!Description: 'FreeOutput' frees the 'output' data structure memory.
 
!Input Parameters:
 this           'output' data structure; the following fields are input:
                   sds.rank, dim[*].name, sds.name, file_name

!Output Parameters:
 this           'output' data structure; the following fields are modified:
                   dim[*].name, sds.name, file_name
 (returns)      status:
                  'true' = okay (always returned)

!Team Unique Header:

 ! Design Notes:
   1. 'OutputFile' must be called before this routine is called.
   2. An error status is never returned.

!END****************************************************************************
*/
bool FreeOutput(Output_t *this)
{
  int ir, ib;
  Myhdf_sds_t *sds;

  if (this->open) 
    RETURN_ERROR("file still open", "FreeOutput", false);

  if (this != (Output_t *)NULL) {

    for (ib = 0; ib < NBAND_CSM; ib++) {
      sds = &this->sds_csm[ib];

      if (sds->name != (char *)NULL) 
        free(sds->name);
      for (ir = 0; ir < sds->rank; ir++) {
        if (sds->dim[ir].name != (char *)NULL) 
          free(sds->dim[ir].name);
      }

      if ( sizeof(uint8) != sizeof(unsigned char) && 
          this->buf_csm[ib] != (uint8 *)NULL          )free(this->buf_csm[ib]);
    }

    if (this->file_name != (char *)NULL) free(this->file_name);
    free(this);
  }

  return true;
}

/* 
!C******************************************************************************
bool PutOutputLine(Output_t *this, int iband, int iline, int *line)
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
bool PutOutputLine(Output_t *this, int iband, int iline, unsigned char *line) 
{
  int32 start[MYHDF_MAX_RANK], nval[MYHDF_MAX_RANK];
  void *buf;
  int is;

  /* Check the parameters */

  if (this == (Output_t *)NULL) 
    RETURN_ERROR("invalid input structure", "PutOutputLine", false);
  if (!this->open)
    RETURN_ERROR("file not open", "PutOutputLine", false);
  if (iband < 0  ||  iband >= NBAND_CSM)
    RETURN_ERROR("invalid band number", "PutOutputLine", false);
  if (iline < 0  ||  iline >= this->size.l)
    RETURN_ERROR("invalid line number", "PutOutputLine", false);

  if (sizeof(uint8) == sizeof(unsigned char))
    buf = (void *)line;
  else {
    for (is = 0; is < this->size.s; is++) 
      this->buf_csm[iband][is] = (uint8)line[is];
    buf = (void *)this->buf_csm[iband];
  }

  /* Write the data */

  start[0] = iline;
  start[1] = 0;
  nval[0] = 1;
  nval[1] = this->size.s;

  if (SDwritedata(this->sds_csm[iband].id, start, NULL, nval, 
                  buf) == HDF_ERROR)
    RETURN_ERROR("writing output", "PutOutputLine", false);

  return true;
}


/* 
!C****************************************************************************** 
bool GetOutputLine(Output_t *this, int iband, int iline, int *line)
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
bool GetOutputLine(Output_t *this, int iband, int iline, int *line)
{
  int32 start[MYHDF_MAX_RANK], nval[MYHDF_MAX_RANK];
  void *buf;
  int is;

  /* Check the parameters */

  if (this == (Output_t *)NULL) 
    RETURN_ERROR("invalid input structure", "GetOutputLine", false);
  if (!this->open)
    RETURN_ERROR("file not open", "GetOutputLine", false);
  if (iband < 0  ||  iband >= NBAND_CSM)
    RETURN_ERROR("invalid band number", "GetOutputLine", false);
  if (iline < 0  ||  iline >= this->size.l)
    RETURN_ERROR("invalid line number", "GetOutputLine", false);


  if (sizeof(uint8) == sizeof(unsigned char))
    buf = (void *)line;
  else
    buf = (void *)this->buf_csm[iband];

  /* Read the data */

  start[0] = iline;
  start[1] = 0;
  nval[0] = 1;
  nval[1] = this->size.s;

  if (SDreaddata(this->sds_csm[iband].id, start, NULL, nval, buf) == HDF_ERROR)
    RETURN_ERROR("reading input", "GetInputLine", false);

  if (sizeof(int16) != sizeof(int)) {
    for (is = 0; is < this->size.s; is++)
      line[is] = (int)this->buf_csm[iband][is];
  }

  return true;
} 


/* 
!C******************************************************************************
bool PutMetadata(Output_t *this, Input_meta_t *meta, Param_t *param)
!C******************************************************************************

!Description: 'PutMetadata' writes metadata to the output file.
 
!Input Parameters:
 this           'output' data structure; the following fields are input:
                   open
 geoloc         'geoloc' data structure; the following fields are input:
                   (none)
 input          'input' data structure;  the following fields are input:
                   (none)

!Output Parameters:
 this           'output' data structure; the following fields are modified:
                   (none)
 (returns)      status:
                  'true' = okay
		  'false' = error return

!Team Unique Header:

 ! Design Notes:
   1. This routine is currently a 'stub' and will be implemented later.
   2. An error status is returned when:
       a. the file is not open for access.
   3. Error messages are handled with the 'RETURN_ERROR' macro.
   4. 'OutputFile' must be called before this routine is called.

!END****************************************************************************
*/
bool PutMetadata(Output_t *this, Input_meta_t *meta, Lut_t *lut,  Param_t *param)
{
  Myhdf_attr_t attr;
  char *string;
  char date[MAX_DATE_LEN + 1];
  double dval[2];
  char short_name[10], local_granule_id[100], production_date[MAX_DATE_LEN + 1];
  char pge_ver[10], process_ver[10];
  int ib;
  const int nband=1;

  /* Check the parameters */

  if (!this->open)
    RETURN_ERROR("file not open", "PutMetadata", false);

  /* Write the metadata */

  attr.id = -1;

  if (meta->provider != PROVIDER_NULL) {
    string = Provider_string[meta->provider].string;
    attr.type = DFNT_CHAR8;
    attr.nval = strlen(string);
    attr.name = OUTPUT_PROVIDER;
    if (!PutAttrString(this->sds_file_id, &attr, string))
      RETURN_ERROR("writing attribute (data provider)", "PutMetadata", false);
  }

  if (meta->sat != SAT_NULL) {
    string = Sat_string[meta->sat].string;
    attr.type = DFNT_CHAR8;
    attr.nval = strlen(string);
    attr.name = OUTPUT_SAT;
    if (!PutAttrString(this->sds_file_id, &attr, string))
      RETURN_ERROR("writing attribute (satellite)", "PutMetadata", false);
  }

  if (meta->inst != INST_NULL) {
    string = Inst_string[meta->inst].string;
    attr.type = DFNT_CHAR8;
    attr.nval = strlen(string);
    attr.name = OUTPUT_INST;
    if (!PutAttrString(this->sds_file_id, &attr, string))
      RETURN_ERROR("writing attribute (instrument)", "PutMetadata", false);
  }

  if (!FormatDate(&meta->acq_date, DATE_FORMAT_DATEA_TIME, date))
    RETURN_ERROR("formating acquisition date", "PutMetadata", false);
  attr.type = DFNT_CHAR8;
  attr.nval = strlen(date);
  attr.name = OUTPUT_ACQ_DATE;
  if (!PutAttrString(this->sds_file_id, &attr, date))
    RETURN_ERROR("writing attribute (acquisition date)", "PutMetadata", false);

  if (!FormatDate(&meta->prod_date, DATE_FORMAT_DATEA_TIME, date))
    RETURN_ERROR("formating production date", "PutMetadata", false);
  attr.type = DFNT_CHAR8;
  attr.nval = strlen(date);
  attr.name = OUTPUT_L1_PROD_DATE;
  if (!PutAttrString(this->sds_file_id, &attr, date))
    RETURN_ERROR("writing attribute (production date)", "PutMetadata", false);

  attr.type = DFNT_FLOAT32;
  attr.nval = 1;
  attr.name = OUTPUT_SUN_ZEN;
  dval[0] = (double)meta->sun_zen * DEG;
  if (!PutAttrDouble(this->sds_file_id, &attr, dval))
    RETURN_ERROR("writing attribute (solar zenith)", "PutMetadata", false);

  attr.type = DFNT_FLOAT32;
  attr.nval = 1;
  attr.name = OUTPUT_SUN_AZ;
  dval[0] = (double)meta->sun_az * DEG;
  if (!PutAttrDouble(this->sds_file_id, &attr, dval))
    RETURN_ERROR("writing attribute (solar azimuth)", "PutMetadata", false);

  if (meta->wrs_sys != WRS_NULL) {

    string = Wrs_string[meta->wrs_sys].string;
    attr.type = DFNT_CHAR8;
    attr.nval = strlen(string);
    attr.name = OUTPUT_WRS_SYS;
    if (!PutAttrString(this->sds_file_id, &attr, string))
      RETURN_ERROR("writing attribute (WRS system)", "PutMetadata", false);

    attr.type = DFNT_INT16;
    attr.nval = 1;
    attr.name = OUTPUT_WRS_PATH;
    dval[0] = (double)meta->ipath;
    if (!PutAttrDouble(this->sds_file_id, &attr, dval))
      RETURN_ERROR("writing attribute (WRS path)", "PutMetadata", false);

    attr.type = DFNT_INT16;
    attr.nval = 1;
    attr.name = OUTPUT_WRS_ROW;
    dval[0] = (double)meta->irow;
    if (!PutAttrDouble(this->sds_file_id, &attr, dval))
      RETURN_ERROR("writing attribute (WRS row)", "PutMetadata", false);

  }

  attr.type = DFNT_INT8;
  attr.nval = nband;
  attr.name = OUTPUT_NBAND;
  dval[0] = (double)1;
  if (!PutAttrDouble(this->sds_file_id, &attr, dval))
    RETURN_ERROR("writing attribute (number of bands)", "PutMetadata", false);

  attr.type = DFNT_INT8;
  attr.nval = nband;
  attr.name = OUTPUT_BANDS;
  for (ib = 0; ib < nband; ib++)
    dval[ib] = (double)1;
  if (!PutAttrDouble(this->sds_file_id, &attr, dval))
    RETURN_ERROR("writing attribute (band numbers)", "PutMetadata", false);

  /* Get the short name, local granule id and production date/time */
  
  if (!Names(meta->sat, meta->inst, "CSM", &meta->acq_date, 
             meta->wrs_sys, meta->ipath, meta->irow, 
	     short_name, local_granule_id, production_date))
    RETURN_ERROR("creating the short name and local granule id", 
                 "PutMetadata", false);

  attr.type = DFNT_CHAR8;
  attr.nval = strlen(short_name);
  attr.name = OUTPUT_SHORT_NAME;
  if (!PutAttrString(this->sds_file_id, &attr, short_name))
    RETURN_ERROR("writing attribute (short name)", "PutMetadata", false);

  attr.type = DFNT_CHAR8;
  attr.nval = strlen(local_granule_id);
  attr.name = OUTPUT_LOCAL_GRAN_ID;
  if (!PutAttrString(this->sds_file_id, &attr, local_granule_id))
    RETURN_ERROR("writing attribute (local granule id)", "PutMetadata", false);

  attr.type = DFNT_CHAR8;
  attr.nval = strlen(production_date);
  attr.name = OUTPUT_PROD_DATE;
  if (!PutAttrString(this->sds_file_id, &attr, production_date))
    RETURN_ERROR("writing attribute (production date)", "PutMetadata", false);

  if (sprintf(pge_ver, "%s", param->PGEVersion) < 0) RETURN_ERROR("creating PGEVersion","PutMetadata", false);
  attr.type = DFNT_CHAR8;
  attr.nval = strlen(pge_ver);
  attr.name = OUTPUT_PGEVERSION;
  if (!PutAttrString(this->sds_file_id, &attr, pge_ver))
    RETURN_ERROR("writing attribute (PGE Version)", "PutMetadata", false);

  if (sprintf(process_ver, "%s", param->ProcessVersion) < 0) RETURN_ERROR("creating ProcessVersion","PutMetadata", false);
  attr.type = DFNT_CHAR8;
  attr.nval = strlen(process_ver);
  attr.name = OUTPUT_PROCESSVERSION;
  if (!PutAttrString(this->sds_file_id, &attr, process_ver))
    RETURN_ERROR("writing attribute (Process Version)", "PutMetadata", false);

  /* now write out the per sds attributes */

  attr.type = DFNT_CHAR8;
  attr.nval = strlen(lut->long_name);
  attr.name = OUTPUT_LONG_NAME;
  if (!PutAttrString(this->sds_csm[0].id, &attr, lut->long_name))
    RETURN_ERROR("writing attribute (long name)", "PutMetadata", false);

  attr.type = DFNT_CHAR8;
  attr.nval = strlen(lut->units);
  attr.name = OUTPUT_UNITS;
  if (!PutAttrString(this->sds_csm[0].id, &attr, lut->units))
    RETURN_ERROR("writing attribute (units ref)", "PutMetadata", false);

  attr.type = DFNT_INT16;
  attr.nval = 2;
  attr.name = OUTPUT_VALID_RANGE;
  dval[0] = (double)lut->valid_range[0];
  dval[1] = (double)lut->valid_range[1];
  if (!PutAttrDouble(this->sds_csm[0].id, &attr, dval))
    RETURN_ERROR("writing attribute (valid range ref)","PutMetadata",false);

  attr.type = DFNT_UINT8;
  attr.nval = 1;
  attr.name = OUTPUT_FILL_VALUE;
  dval[0] = (double)lut->output_fill;
  if (!PutAttrDouble(this->sds_csm[0].id, &attr, dval))
    RETURN_ERROR("writing attribute (valid range ref)","PutMetadata",false);

  attr.type = DFNT_FLOAT64;
  attr.nval = 1;
  attr.name = OUTPUT_SCALE_FACTOR;
  dval[0] = (double)lut->scale_factor;
  if (!PutAttrDouble(this->sds_csm[0].id, &attr, dval))
    RETURN_ERROR("writing attribute (scale factor ref)", "PutMetadata", false);

  attr.type = DFNT_FLOAT64;
  attr.nval = 1;
  attr.name = OUTPUT_ADD_OFFSET;
  dval[0] = (double)lut->add_offset;
  if (!PutAttrDouble(this->sds_csm[0].id, &attr, dval))
    RETURN_ERROR("writing attribute (add offset ref)", "PutMetadata", false);

  attr.type = DFNT_FLOAT64;
  attr.nval = 1;
  attr.name = OUTPUT_SCALE_FACTOR_ERR;
  dval[0] = (double)lut->scale_factor_err;
  if (!PutAttrDouble(this->sds_csm[0].id, &attr, dval))
    RETURN_ERROR("writing attribute (scale factor err ref)", "PutMetadata", false);

  attr.type = DFNT_FLOAT64;
  attr.nval = 1;
  attr.name = OUTPUT_ADD_OFFSET_ERR;
  dval[0] = (double)lut->add_offset_err;
  if (!PutAttrDouble(this->sds_csm[0].id, &attr, dval))
    RETURN_ERROR("writing attribute (add offset err ref)", "PutMetadata", false);

  attr.type = DFNT_FLOAT32;
  attr.nval = 1;
  attr.name = OUTPUT_CALIBRATED_NT;
  dval[0] = (double)lut->calibrated_nt;
  if (!PutAttrDouble(this->sds_csm[0].id, &attr, dval))
    RETURN_ERROR("writing attribute (calibrated nt ref)","PutMetadata",false);

  attr.type = DFNT_UINT8;
  attr.nval = 1;
  attr.name = OUTPUT_LAND_VALUE;
  dval[0] = (double)lut->csm_land;
  if (!PutAttrDouble(this->sds_csm[0].id, &attr, dval))
    RETURN_ERROR("writing attribute (land value)","PutMetadata",false);

  attr.type = DFNT_UINT8;
  attr.nval = 1;
  attr.name = OUTPUT_CLOUD_VALUE;
  dval[0] = (double)lut->csm_cloud;
  if (!PutAttrDouble(this->sds_csm[0].id, &attr, dval))
    RETURN_ERROR("writing attribute (cloud value)","PutMetadata",false);

  attr.type = DFNT_UINT8;
  attr.nval = 1;
  attr.name = OUTPUT_SNOW_VALUE;
  dval[0] = (double)lut->csm_snow;
  if (!PutAttrDouble(this->sds_csm[0].id, &attr, dval))
    RETURN_ERROR("writing attribute (snow value)","PutMetadata",false);

  attr.type = DFNT_UINT8;
  attr.nval = 1;
  attr.name = OUTPUT_SHADOW_VALUE;
  dval[0] = (double)lut->csm_shadow;
  if (!PutAttrDouble(this->sds_csm[0].id, &attr, dval))
    RETURN_ERROR("writing attribute (shadow value)","PutMetadata",false);

  attr.type = DFNT_UINT8;
  attr.nval = 1;
  attr.name = OUTPUT_WATER_VALUE;
  dval[0] = (double)lut->csm_water;
  if (!PutAttrDouble(this->sds_csm[0].id, &attr, dval))
    RETURN_ERROR("writing attribute (water value)","PutMetadata",false);


  return true;
}
