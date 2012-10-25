/*
!C****************************************************************************

!File: param.c
  
!Description: Functions for accepting parameters from the command line or 
 a file.

!Revision History:
 Revision 1.0 2000/11/07
 Robert Wolfe
 Original Version.

 Revision 1.1 2000/12/13
 Sadashiva Devadiga
 Modified to accept parameters from command line or file.

 Revision 1.2 2001/05/08
 Sadashiva Devadiga
 Added checks for required parameters.

 Revision 1.3 2002/03/02
 Robert Wolfe
 Added special handling for input ISINUS case.

 Revision 1.4 2002/05/10
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
      phone: 301-614-5508               Lanham, MD 20770  

      Sadashiva Devadiga (Code 922)
      MODIS Land Team Support Group     SSAI
      devadiga@ltpmail.gsfc.nasa.gov    5900 Princess Garden Pkway, #300
      phone: 301-614-5549               Lanham, MD 20706
  
 ! Design Notes:
   1. The following public functions handle the input data:

	GetParam - Setup 'param' data structure and populate with user
	           parameters.
	FreeParam - Free the 'param' data structure memory.

   2. 'GetParam' must be called before 'FreeParam'.  
   3. 'FreeParam' should be used to free the 'param' data structure.

!END****************************************************************************
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "param.h"
#include "mystring.h"
#include "error.h"

#define LUT_FILE_NAME ("lndcsm_lut_NULL.txt")


typedef enum {
  PARAM_NULL = -1,
  PARAM_START = 0,
  PARAM_REF_FILE,
  PARAM_THERM_FILE,
  PARAM_CSM_FILE,
  PARAM_CSM_FILTER,
  PARAM_PGEVERSION,
  PARAM_PROCESSVERSION,
  PARAM_END,
  PARAM_MAX
} Param_key_t;

Key_string_t Param_string[PARAM_MAX] = {
  {(int)PARAM_START,     "PARAMETER_FILE"},
  {(int)PARAM_REF_FILE,  "REF_FILE"},
  {(int)PARAM_THERM_FILE,"THERM_FILE"},
  {(int)PARAM_CSM_FILE,  "CSM_FILE"},
  {(int)PARAM_CSM_FILTER,"CSM_FILTER"},
  {(int)PARAM_PGEVERSION,  "PGEVersion"},
  {(int)PARAM_PROCESSVERSION,  "ProcessVersion"},
  {(int)PARAM_END,       "END"}
};

/* Functions */

Param_t *GetParam(int argc, const char **argv)
/* 
!C******************************************************************************

!Description: 'GetParam' sets up the 'param' data structure and populate with user
 parameters, either from the command line or from a parameter file.
 
!Input Parameters:
 argc           number of command line arguments
 argv           command line argument list

!Output Parameters:
 (returns)      'param' data structure or NULL when an error occurs

!Team Unique Header:

 ! Design Notes:
   1. An error status is returned when:
       a. memory allocation is not successful
       b. an error is returned from the ReadCmdLine function
       c. certain required parameters are invalid or not entered:
            input file name, output file name, geolocation file name,
	    SDS name, output space projection number, 
	    output space pixel size, output space upper left corner, 
	    output space image size, either output space sphere or projection
	    parameter number 0
	    output space zone number not given for UTM
       d. the input is 'GRID_SPACE' then and certain required parameters 
          are invalid or not entered:
	    input space projection number, input space pixel size, 
	    input space upper left corner, either input space sphere or 
	    projection parameter number 0
       e. the input is 'GRID_SPACE' and the geolocation file name is given
       f. the resampling kernel is 'USER_KERNEL' and the kernel file name 
          is invalid or not entered.
   2. Error of type 'a' are handled with the 'RETURN_ERROR' macro and 
      the others are handled by writting the error messages to 'stderr' and 
      then printing the usage information.
   3. 'FreeParam' should be called to deallocate memory used by the 
      'param' data structures.

!END****************************************************************************
*/
{
  Param_t *this;
  char *error_string = (char *)NULL;
  FILE *fp;
  Key_t key;
  int len;
  char line[MAX_STR_LEN + 1];
  char temp[MAX_STR_LEN + 1];
  Param_key_t param_key;
  char *param_file_name;
  int itemp;
  bool got_start, got_end;

  if (argc < 2) 
    RETURN_ERROR("no command line parameter", "GetParam", 
                 (Param_t *)NULL);
  if (argc > 3) 
    RETURN_ERROR("too many command line parameters", "GetParam", 
                 (Param_t *)NULL);
  if (strlen(argv[1]) < 1)
    RETURN_ERROR("no paramter file name", "GetParam", 
                 (Param_t *)NULL);
  param_file_name = (char *)argv[1];

  /* Open the parameter file */
  
  if ((fp = fopen(param_file_name, "r")) == (FILE *)NULL)
    RETURN_ERROR("unable to open parameter file", "GetParam", 
                 (Param_t *)NULL);

  /* Create the Param data structure */

  this = (Param_t *)malloc(sizeof(Param_t));
  if (this == (Param_t *)NULL) {
    fclose(fp);
    RETURN_ERROR("allocating Input structure", "GetParam", 
                 (Param_t *)NULL);
  }

  /* set default parameters */

  this->param_file_name  = (char *)NULL;
  this->input_file_name  = (char *)NULL;
  this->input_therm_file_name  = (char *)NULL;
  this->lut_file_name    = (char *)NULL;
  this->output_file_name = (char *)NULL;
  this->PGEVersion              = (char *)NULL;
  this->ProcessVersion          = (char *)NULL;
  this->therm_flag  = false;
  this->apply_kernel =  0;
  this->ksize =         0;
  this->sieve_thresh =  0;

  /* Populate the data structure */

  this->param_file_name = DupString(param_file_name);
  if (this->param_file_name == (char *)NULL)
    error_string = "duplicating parameter file name";

  if (error_string == (char *)NULL) {
    this->lut_file_name = DupString(LUT_FILE_NAME);
    if (this->lut_file_name == (char *)NULL)
      error_string = "duplicating lookup table file name";
  }

  if (error_string != (char *)NULL) {
    if (this->param_file_name != (char *)NULL) 
      free(this->param_file_name);
    if (this->lut_file_name != (char *)NULL) 
      free(this->lut_file_name);
    FreeParam(this);
    RETURN_ERROR(error_string, "GetParam", (Param_t *)NULL);
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

    param_key = (Param_key_t)
       KeyString(key.key, key.len_key, Param_string, 
		 (int)PARAM_NULL, (int)PARAM_MAX);
    if (param_key == PARAM_NULL) {
      key.key[key.len_key] = '\0';
      sprintf(temp, "invalid key; key = %s", key.key);
      error_string = temp;
      break;
    }
    if (!got_start) {
      if (param_key == PARAM_START) {
        if (key.nval != 0) {
	  error_string = "no value expected (start key)";
	  break;
	}
        got_start = true;
	continue;
      } else {
        error_string  = "no start key in parameter file";
        break;
      }
    }

    /* Get the value for each keyword */

    switch (param_key) {

      case PARAM_REF_FILE:
        if (key.nval <= 0) {
	  error_string = "no reflectance file name";
	  break; 
	} else if (key.nval > 1) {
	  error_string = "too many reflectance file names";
	  break; 
	}
	if (key.len_value[0] < 1) {
	  error_string = "no reflectance file name";
	  break;
	}
	key.value[0][key.len_value[0]] = '\0';
        this->input_file_name = DupString(key.value[0]);
	if (this->input_file_name == (char *)NULL) {
	  error_string = "duplicating reflectance file name";
	  break;
        }
        break;

      case PARAM_THERM_FILE:
        if (key.nval <= 0) {
          printf("warning %s \n","no band6 file name"); 
	  break; 
	} else if (key.nval > 1) {
	  error_string = "too many band6 file names";
          printf("error_string=(%s)\n",error_string); 
	  break; 
	}
	if (key.len_value[0] < 1) {
          printf("warning %s \n","no band6 file name"); 
	  break;
	}
	key.value[0][key.len_value[0]] = '\0';
        this->input_therm_file_name = DupString(key.value[0]);
        this->therm_flag  = true;
	if (this->input_file_name == (char *)NULL) {
	  error_string = "duplicating band6 file name";
	  break;
        }
        break;

      case PARAM_CSM_FILE:
        if (key.nval <= 0) {
	  error_string = "no CLoud Mask file name";
	  break; 
	} else if (key.nval > 1) {
	  error_string = "too many Cloud Mask file names";
	  break; 
	}
	if (key.len_value[0] < 1) {
	  error_string = "no Cloud Mask file name";
	  break;
	}
	key.value[0][key.len_value[0]] = '\0';
        this->output_file_name = DupString(key.value[0]);
	if (this->output_file_name == (char *)NULL) {
	  error_string = "duplicating Cloud Mask file name";
	  break;
        }
        break;

      case PARAM_CSM_FILTER:
        if (key.nval <= 0) {
	  break; 
	} else if (key.nval > 3) {
	  error_string = "too many filter values";
	  break; 
	}
	if (key.len_value[0] < 1) {
	  break;
	}
	if (key.nval > 0)
	  {
          key.value[0][key.len_value[0]] = '\0';
          this->ksize = atoi(key.value[0]);
          this->apply_kernel = 1;
	  }
	if (key.nval > 1)
	  {
          key.value[1][key.len_value[1]] = '\0';
          itemp = atoi(key.value[1]);
          if ( itemp > this->apply_kernel )this->apply_kernel = itemp;
	  }
	if (key.nval > 2)
	  {
          key.value[2][key.len_value[2]] = '\0';
          this->sieve_thresh= atoi(key.value[2]);
	  }
        break;

      case PARAM_PGEVERSION:
        if (key.nval <= 0) {
          error_string = "no PGEVersion number";
          break;
        } else if (key.nval > 1) {
          error_string = "too many PGEVersion numbers";
          break;
        }
        if (key.len_value[0] < 1) {
          error_string = "no PGEVersion number";
          break;
        }
        key.value[0][key.len_value[0]] = '\0';
        this->PGEVersion = DupString(key.value[0]);
        if (this->PGEVersion == (char *)NULL) {
          error_string = "duplicating PGEVersion number";
          break;
        }
        break;

      case PARAM_PROCESSVERSION:
        if (key.nval <= 0) {
          error_string = "no ProcessVersion number";
          break;
        } else if (key.nval > 1) {
          error_string = "too many ProcessVersion numbers";
          break;
        }
        if (key.len_value[0] < 1) {
          error_string = "no ProcessVersion number";
          break;
        }
        key.value[0][key.len_value[0]] = '\0';
        this->ProcessVersion = DupString(key.value[0]);
        if (this->ProcessVersion == (char *)NULL) {
          error_string = "duplicating ProcessVersion number";
          break;
        }
        break;

      case PARAM_END:
        if (key.nval != 0) {
	  error_string = "no value expected (end key)";
	  break; 
	}
        got_end = true;
        break;

      default:
        error_string = "key not implmented";

    }
    if (error_string != (char *)NULL) break;
    if (got_end) break;

  }

  /* Close the parameter file */

  fclose(fp);

  if (error_string == (char *)NULL) {
    if (!got_start) 
      error_string = "no start key in header";
    else if (!got_end)
      error_string = "no end key in header";
  }

  /* Handle null values */
  
  if (error_string == (char *)NULL) {
    if (this->input_file_name == (char *)NULL) 
      error_string = "no input file name given";
    if (this->lut_file_name == (char *)NULL) 
      error_string = "no lookup table file name given";
    if (this->output_file_name == (char *)NULL) 
      error_string = "no output file name given";
    if (this->PGEVersion == (char *)NULL)
      error_string = "no PGE Version given";
    if (this->ProcessVersion == (char *)NULL)
      error_string = "no Process Version given";
  }

  /* Handle errors */

  if (error_string != (char *)NULL) {
    if (this->param_file_name != (char *)NULL) 
      free(this->param_file_name);
    if (this->input_file_name != (char *)NULL) 
      free(this->input_file_name);
    if (this->lut_file_name != (char *)NULL) 
      free(this->lut_file_name);
    if (this->output_file_name != (char *)NULL) 
      free(this->output_file_name);
    if (this->PGEVersion != (char *)NULL)
      free(this->PGEVersion);
    if (this->ProcessVersion != (char *)NULL)
      free(this->ProcessVersion);
    if (this != (Param_t *)NULL) free(this);
    RETURN_ERROR(error_string, "GetParam", (Param_t *)NULL);
  }
  
  return this;
}


bool FreeParam(Param_t *this)
/* 
!C******************************************************************************

!Description: 'FreeParam' frees the 'param' data structure memory.
 
!Input Parameters:
 this           'param' data structure; the following fields are input:
                   input_file_name, output_file_name, geoloc_file_name, 
		   input_sds_name, output_sds_name, kernel_file_name

!Output Parameters:
 (returns)      status:
                  'true' = okay (always returned)

!Team Unique Header:

 ! Design Notes:
   1. 'GetParam' must be called before this routine is called.
   2. An error status is never returned.

!END****************************************************************************
*/
{
  if (this != (Param_t *)NULL) {
    if (this->param_file_name        != (char *)NULL) free(this->param_file_name);
    if (this->input_file_name        != (char *)NULL) free(this->input_file_name);
    if (this->input_therm_file_name  != (char *)NULL) free(this->input_therm_file_name);
    if (this->lut_file_name          != (char *)NULL) free(this->lut_file_name);
    if (this->output_file_name       != (char *)NULL) free(this->output_file_name);
    free(this);
  }
  return true;
}                            
