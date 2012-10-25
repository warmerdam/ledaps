#ifndef NAMES_H
#define NAMES_H

#include "lndcal.h"
#include "date.h"
#include "bool.h"

bool Names(Sat_t sat, Inst_t inst, char *product_id, 
           Date_t *acq_date, Wrs_t wrs_sys, 
           int ipath, int irow, char *short_name, 
	   char *local_granule_id, char *production_date);
	   
#endif
