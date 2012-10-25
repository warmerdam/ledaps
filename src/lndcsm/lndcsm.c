#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "keyvalue.h"
#include "lndcsm.h"
#include "keyvalue.h"
#include "const.h"
#include "param.h"
#include "input.h"
#include "lut.h"
#include "output.h"
#include "csm.h"
#include "bool.h"
#include "space.h"
#include "error.h"
#define NSDS 1
/* Type definitions */

typedef enum {FAILURE = 0, SUCCESS = 1} Status_t;

/* Functions */

int main (int argc, const char **argv) {
  Param_t *param;
  Input_t *input;
  Input_t *input_th;
  Lut_t *lut;
  Output_t *output;
  int i,odometer_flag=0;

  /* Space info */
  Space_def_t space_def;
  Space_t *space;
  char *grid_name = "Grid";
  char *sds_names[NSDS];
  int sds_types[NSDS];

  for (i=1; i<argc; i++)if ( !strcmp(argv[i],"-o") )odometer_flag=1;
  param = GetParam(argc, argv);
  if (param == (Param_t *)NULL) ERROR("getting runtime parameters", "main");

  /* Open input file */

  input = OpenInput(param->input_file_name);
  if (input == (Input_t *)NULL) ERROR("bad input file", "main");

  /* Get space definition from HDF */

  if (!GetSpaceDefHDF(&space_def, param->input_file_name, grid_name))
    ERROR("getting definition from file", "main");

  /* Setup Space */

  space = SetupSpace(&space_def);
  if (space == (Space_t *)NULL) ERROR("setting up space", "main");

  /* Open input thermal file */
  if ( param->therm_flag )
    {
    input_th = OpenInput(param->input_therm_file_name);
    if (input_th == (Input_t *)NULL) ERROR("bad thermal input file", "main");

  if ( input_th->size.s != input->size.s || input_th->size.l != input->size.l )
     ERROR("mismatch size on inputs (reflective vs thermal)", "main");
    }

  /* Get Lookup table */
  lut = GetLut(param->lut_file_name, input->nband, &input->meta);
  if (lut == (Lut_t *)NULL) ERROR("bad lut file", "main");
  
  /* Create and open output file */

  if (!CreateOutput(param->output_file_name))
    ERROR("creating output file", "main");

  output = OpenOutput( param->output_file_name, &input->size );
  if (output == (Output_t *)NULL) ERROR("opening output file", "main");

  /******************************************************************/
  /*                         call snow mask                         */
  /******************************************************************/

  if ( !CloudMask( 
                  input 
                , input_th
                , output
                , lut
                , param
                , odometer_flag
		  ) )
     ERROR("In Cloud/Snow Mask Computation","main");

  if (!PutMetadata(output, &input->meta,lut, param) )
    ERROR("writing the metadata", "main");


  if (!CloseInput(input)) ERROR("closing input file", "main");

  if (!CloseOutput(output)) ERROR("closing output file", "main");

  sds_names[0] = output->sds_csm[0].name;
  sds_types[0] = output->sds_csm[0].type;

  if (!PutSpaceDefHDF(&space_def, param->output_file_name, 1, 
                      sds_names, sds_types, grid_name))
    ERROR("putting space metadata in HDF file", "main");

  if (!FreeInput(input)) 
    ERROR("freeing input file stucture", "main");

  if (!FreeLut(lut)) 
    ERROR("freeing lut file stucture", "main");

  if (!FreeOutput(output)) 
    ERROR("freeing output file stucture", "main");

  /* All done */

  exit(EXIT_SUCCESS);

  return (EXIT_SUCCESS);
}
