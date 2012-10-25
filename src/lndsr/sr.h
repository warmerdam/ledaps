#ifndef SR_H
#define SR_H

#include "lndsr.h"
#include "bool.h"
#include "lut.h"
#include "error.h"

typedef struct {
  bool first[NBAND_SR_MAX];
  int sr_min[NBAND_SR_MAX];
  int sr_max[NBAND_SR_MAX];
  long nfill[NBAND_SR_MAX];
  long nout_range[NBAND_SR_MAX];
} Sr_stats_t;


bool Sr(Lut_t *lut, int nsamp, int il, int **line_in, bool mask_flag, 
        char* mask_line, int **line_out, Sr_stats_t *sr_stats);

bool computeBounds(Geo_bounds_t *bounds, Space_t *space, int nps, int nls);
#endif
