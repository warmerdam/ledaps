#ifndef CLOUDS_H
#define CLOUDS_H

#include "lndsr.h"
#include "bool.h"
#include "lut.h"
#include "sixs_runs.h"

#define CLDDIAGS_CELLHEIGHT_1KM 40
#define CLDDIAGS_CELLWIDTH_1KM 40
#define CLDDIAGS_CELLHEIGHT_5KM 160
#define CLDDIAGS_CELLWIDTH_5KM 160
#define CLDDIAGS_CELLHEIGHT_10KM 330
#define CLDDIAGS_CELLWIDTH_10KM 330

typedef struct cld_diags_t {
	int nbrows,nbcols,cellheight,cellwidth;
	float **avg_t6_clear,**std_t6_clear,**avg_b7_clear,**std_b7_clear,**airtemp_2m;
	int **nb_t6_clear;
}cld_diags_t;



int allocate_cld_diags(struct cld_diags_t *cld_diags,int cell_height, int cell_width, int scene_height, int scene_width);
int free_cld_diags(struct cld_diags_t *cld_diags);
void fill_cld_diags(cld_diags_t *cld_diags);
int interpol_clddiags_1pixel(cld_diags_t *cld_diags, int img_line, int img_sample,float *inter_value);

bool cloud_detection_pass1(Lut_t *lut, int nsamp, int il, int **line_in, int8 *qa_line, int *b6_line,float *atemp_line,
		  cld_diags_t *cld_diags);
bool cloud_detection_pass2(Lut_t *lut, int nsamp, int il, int **line_in, int8 *qa_line, int *b6_line,
		  cld_diags_t *cld_diags,char *ddv_line);
bool cast_cloud_shadow(Lut_t *lut, int nsamp, int il_start, int ***line_in, int **b6_line,
		  cld_diags_t *cld_diags,char ***cloud_buf, Ar_gridcell_t *ar_gridcell,float pixel_size,float adjust_north);
bool dilate_cloud_mask(Lut_t *lut, int nsamp, char ***cloud_buf, int dilate_dist);
bool dilate_shadow_mask(Lut_t *lut, int nsamp, char ***cloud_buf, int dilate_dist);

#endif
