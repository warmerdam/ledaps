
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "names.h"
#include "error.h"

#define VERSION_ID (20)

const char *sat_names[SAT_MAX] = {"1", "2", "3", "4", "5", "7"};
const char *inst_names[INST_MAX] = {"M", "T", "E"};
const char *wrs_names[WRS_MAX] = {"1", "2"};

bool Names(Sat_t sat, Inst_t inst, char *product_id, 
           Date_t *acq_date, Wrs_t wrs_sys, 
           int ipath, int irow, char *short_name, 
	   char *local_granule_id, char *production_date) {

  char acq_date_string[MAX_DATE_LEN + 1];
  char prod_date_string[14];
  time_t tp;
  struct tm *tm;
  int version_id;

  if (sat == SAT_NULL  ||  sat < SAT_LANDSAT_1 || sat >= SAT_MAX)
    RETURN_ERROR("invalid parameter (satellite)", "Names", false);
  if (inst == INST_NULL  ||  inst < INST_MSS  ||  inst >= INST_MAX)
    RETURN_ERROR("invalid parameter (instrument)", "Names", false);
  if (wrs_sys == WRS_NULL  ||  wrs_sys < WRS_1  || wrs_sys >= WRS_MAX) 
    RETURN_ERROR("invalid parameter (wrs system)", "Names", false);
  if (ipath < 0  ||  ipath > 999) 
    RETURN_ERROR("invalid parameter (path number)", "Names", false);
  if (irow < 0  ||  irow > 999) 
    RETURN_ERROR("invalid parameter (row number)", "Names", false);

  if (sprintf(short_name, "L%s%s%s", sat_names[sat], 
              inst_names[inst], product_id) < 0)
    RETURN_ERROR("creating short name", "Names", false);

  if (time(&tp) == -1)
    RETURN_ERROR("getting time", "Names", false);
    
  tm = gmtime(&tp);
  if (tm == (struct tm *)NULL)
    RETURN_ERROR("converting time to UTC", "Names", false);

  if (strftime(prod_date_string, 14, "%Y%j%H%M%S", tm) == 0)
    RETURN_ERROR("formating production date/time (1)", "Names", false);

  if (strftime(production_date, (MAX_DATE_LEN + 1), 
               "%Y-%m-%dT%H:%M:%SZ", tm) == 0)
    RETURN_ERROR("formating production date/time (2)", "Names", false);

  if (!FormatDate(acq_date, DATE_FORMAT_DATEB, acq_date_string))
    RETURN_ERROR("formatting acquisition date", "Names", false);
  acq_date_string[4] = '\0';

  version_id = VERSION_ID;  /* place holder */

  if (sprintf(local_granule_id, 
              "%s.a%4s%3s.w%1sp%03dr%03d.%03d.%13s.hdf",
              short_name, acq_date_string, &acq_date_string[5], 
	      wrs_names[wrs_sys], ipath, irow, 
	      version_id, prod_date_string) < 0)
    RETURN_ERROR("creating local granule id", "Names", false);

  return true;
}
