#ifndef EE_LNDSR_MAIN_H
#define EE_LNDSR_MAIN_H
int ee_lndsr_main(
    const char *metadata,
    const char *sixs_filename,

    int full_input_size_s, int full_input_size_l,
    int input_offset_s, int input_offset_l,
    int input_size_s, int input_size_l, int nband,
    short *band1, short *band2, short *band3,
    short *band4, short *band5, short *band7,
    const short *thermal, const int8 *lndcal_QA,
    const char *csm_mask,
    short *atmos_opacity, short *lndsr_QA,

    int utm_zone,
    double ul_x, double ul_y, double pixel_size,

    int ar_size_s, int ar_size_l,
    const short *aerosol,
    const float *anc_SP,
    const float *anc_WV,
    const float *anc_ATEMP,
    const float *anc_O3,
    const float *anc_dem);

#endif /* EE_LNDSR_MAIN_H */
