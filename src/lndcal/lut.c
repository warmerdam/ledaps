#include "stdlib.h"
#include "lut.h"
#include "input.h"
#include "mystring.h"
#include "error.h"
#include "const.h"

#define OUTPUT_FILL (-9999)
#define OUTPUT_SATU (20000)
#define QA_FILL     (1)
#define QA_SAT      (2)
#define LONG_NAME_PREFIX_REF ("band %d reflectance")
#define UNITS_REF            ("reflectance")
#define LONG_NAME_PREFIX_TH  ("band %d temperature")
#define UNITS_TH             ("temperature (celcius)")
#define VALID_MIN_REF        (-100)
#define VALID_MAX_REF        (16000)
#define VALID_MIN_TH         (-7000)
#define VALID_MAX_TH         (7000)
#define SCALE_FACTOR_REF     (0.0001)
#define SCALE_FACTOR_TH      (0.01)
#define ADD_OFFSET_REF       (0.0)
#define ADD_OFFSET_TH        (0.0)
#define SCALE_FACTOR_ERR_REF (0.0)
#define SCALE_FACTOR_ERR_TH  (0.0)
#define ADD_OFFSET_ERR_REF   (0.0)
#define ADD_OFFSET_ERR_TH    (0.0)
#define CALIBRATED_NT_REF    (DFNT_FLOAT32)
#define CALIBRATED_NT_TH     (DFNT_FLOAT32)
#define RECAL_YEAR  (2003)
#define RECAL_MONTH (05)
#define RECAL_DAY   (05)
#define RECAL_DOY   (125)
#define NORECAL_YEAR  (2007)
#define NORECAL_DOY   (92)

/* Landsat 7 ETM+ Minimum and Maximum Spectral Radiance (LMIN and LMAX)
   Units = watts / (meter squared * steradian * micrometer)
   Table A is for data produced before July 1, 2000 and 
   Table B is for data produced on or after
   Reference:
 http://lpwww.gsfc.nasa.gov/IAS/handbook/handbook_htmls/chaper11/chapter11.html
   (web site visited on August 22, 2002)
*/


#define TABLE_B_DATE ("2000-07-01")
#define QCALMIN_LPGS (1)
#define QCALMIN_NLAPS (0)
#define QCALMAX (255)

const float etm_lmin_low_a[8] = {
  -6.2, -6.0, -4.5, -4.5, -1.0, 0.0, -0.35, -5.0
};

const float etm_lmax_low_a[8] = {
  297.5, 303.4, 235.5, 235.0, 47.70, 17.04, 16.60, 244.00
};

const float etm_lmin_high_a[8] = {
  -6.2, -6.0, -4.5, -4.5, -1.0, 3.2, -0.35, -5.0
};

const float etm_lmax_high_a[8] = {
  194.3, 202.4, 158.6, 157.5, 31.76, 12.65, 10.932, 158.40
};

const float etm_lmin_low_b[8] = {
  -6.2, -6.4, -5.0, -5.1, -1.0, 0.0, -0.35, -4.7
};

const float etm_lmax_low_b[8] = {
  293.7, 300.9, 234.4, 241.1, 47.57, 17.04, 16.54, 243.1
};

const float etm_lmin_high_b[8] = {
   -6.2, -6.4, -5.0, -5.1, -1.0, 3.2, -0.35, -4.7
};

const float etm_lmax_high_b[8] = {
   191.6, 196.5, 152.9, 157.4, 31.06, 12.65, 10.80, 158.3
};

const float mss_lmin_l1[4] = {  0.0,  0.0,  0.0,  0.0 };
const float mss_lmin_l2[4] = {  0.8,  0.6,  0.6,  0.4 };
const float mss_lmin_l3[4] = {  0.4,  0.3,  0.3,  0.1 };
const float mss_lmin_l4[4] = {  0.4,  0.4,  0.5,  0.4 };
const float mss_lmin_l5[4] = {  0.3,  0.3,  0.5,  0.3 };

const float mss_lmax_l1[4] = { 24.8, 20.0, 17.6, 15.3 };
const float mss_lmax_l2[4] = { 26.3, 17.6, 15.2, 13.0 };
const float mss_lmax_l3[4] = { 25.9, 17.9, 14.9, 12.8 };
const float mss_lmax_l4[4] = { 23.8, 16.4, 14.2, 11.6 };
const float mss_lmax_l5[4] = { 26.8, 17.9, 14.8, 12.3 };

const float delta= 0.00001;

/* Relative sun-earth distance variation (Ref. L7 handbook) */

#define NDSUN (25)
typedef struct {
  int doy;         /* Day of year */
  float dsun;      /* Relative sun-earth distance */
} Dsun_table_t;

const Dsun_table_t dsun_table[NDSUN] = {
  {  1, 0.9832}, { 15, 0.9836}, { 32, 0.9853}, { 46, 0.9878}, { 60, 0.9909},
  { 74, 0.9945}, { 91, 0.9993}, {106, 1.0033}, {121, 1.0076}, {135, 1.0109}, 
  {152, 1.0140}, {166, 1.0158}, {182, 1.0167}, {196, 1.0165}, {213, 1.0149},
  {227, 1.0128}, {242, 1.0092}, {258, 1.0057}, {274, 1.0011}, {288, 0.9972},
  {305, 0.9925}, {319, 0.9892}, {335, 0.9860}, {349, 0.9843}, {366, 0.9833}
};


/* Solar Spectral Irradiances 
   Units = watts / (meter squared * micrometers)
   Reference:
     ETM+ values are from the Landsat 7 User Handbook (see above) 
     TM and MSS values from EOSAT Landsat Technical Notes, August 1986 
     
     New version of L4 and L5 given Solar Irradinaces in  (mW/m^2-nm) to Jeff Masek from Brian Markham on 11/3/06 : 
>L7: 1969      1840      1551      1044    225.7      82.07     1368
>L5: 1957      1826      1554      1036    215.0      80.67  -
>L4: 1957      1825      1557      1033    214.9      80.72  -
*/

const float esun_etm[8] = {
  1969.000, 1840.000, 1551.000, 1044.000, 225.700, -1.0, 82.07, 1368.000
};

/*
const float esun_tm_4[7] = {
  1958.0, 1828.0, 1559.0, 1045.0, 219.1, -1.0, 74.57
};
*/
const float esun_tm_4[7] = {
  1957.0, 1825.0, 1557.0, 1033.0, 214.9, -1.0, 80.72
};

/*
const float esun_tm_5[7] = {
  1957.0, 1829.0, 1557.0, 1047.0, 219.3, -1.0, 74.57
};
*/

const float esun_tm_5[7] = {
  1957.0, 1826.0, 1554.0, 1036.0, 215.0, -1.0, 80.67
};

const float esun_mss_1[4] = {
  1852.0, 1584.0, 1276.0, 904.0
};

const float esun_mss_2[4] = {
  1856.0, 1559.0, 1269.0, 906.0
};

const float esun_mss_3[4] = {
  1860.0, 1571.0, 1289.0, 910.0
};

const float esun_mss_4[4] = {
  1851.0, 1593.0, 1260.0, 878.0
};

const float esun_mss_5[4] = {
  1849.0, 1595.0, 1253.0, 870.0
};

/* Landsat 4/5 TM Bandwidths
   Units = micrometers
   Reference: EOSAT Fast Format Document Version B, Effective Dec. 1, 1993
*/

const float band_width_tm_4[7] = {
  0.066, 0.081, 0.069, 0.129, 0.216, 1.00, 0.250
};

const float band_width_tm_5[7] = {
  0.066, 0.082, 0.067, 0.128, 0.217, 1.00, 0.252
};

const float K1_tm_4=  671.62;
const float K2_tm_4= 1284.30;
const float K1_tm_5=  607.76;
const float K2_tm_5= 1260.56;
const float K1_etm =  666.09;
const float K2_etm = 1282.71;

/* This program uses gains & biases in Landsat 7 radiance units of [W/(m^2 sr micron)].*/

Lut_t *GetLut(Param_t *param, int nband, Input_t *input) {
  Lut_t *this;
  Date_t date_table_b;
  float lmin[8], lmax[8];
  float lmin_th=0.0, lmax_th=0.0;
  float qcalmin, qcalmax;
  int ib, iband;
  int iband6= 6-1;
  double diff;
  int jdoy, i;
  float dsun;
  char msgbuf[1024];
  char *file_name= param->lut_file_name;
  bool after_2003;
  Input_meta_t *input_meta= &(input->meta);
  

  /* Create the lookup table data structure */

  this = (Lut_t *)malloc(sizeof(Lut_t));
  if (this == (Lut_t *)NULL) 
    RETURN_ERROR("allocating Input data structure", "OpenInput", 
                 (Lut_t *)NULL);

  /* Populate the data structure */

  this->in_fill = 0;
  this->qa_fill = QA_FILL;
  this->qa_sat =  QA_SAT;
  this->file_name = DupString(file_name);
  if (this->file_name == (char *)NULL) {
    free(this);
    RETURN_ERROR("duplicating file name", "GetLut", (Lut_t *)NULL);
  }

  if (!InputMetaCopy(input_meta, nband, &this->meta)) {
    free(this->file_name);
    free(this);
    RETURN_ERROR("copying input metadata", "GetLut", (Lut_t *)NULL);
  }
  this->out_fill = OUTPUT_FILL;
  this->out_satu = OUTPUT_SATU;

  /****************************************/
  /*        Check the dates for MSS       */
  /****************************************/

  if (input_meta->inst == INST_MSS) {

/*
Notes:  values for Landsats 4,5 only valid for data acquired after 4/1/83 and 11/9/84, respectively
Values for data *processed* during 1972-1978 may be different for Landsat-1-3
*/							

/*
    printf("***  acq date=(%2.2d/%2.2d/%4.4d) doy=(%3.3d) ***\n"
   ,  input_meta->acq_date.month
   ,  input_meta->acq_date.day
   ,  input_meta->acq_date.year
   ,  input_meta->acq_date.doy);  fflush(stdout);
    printf("*** prod date=(%2.2d/%2.2d/%4.4d) doy=(%3.3d) ***\n"
   ,  input_meta->prod_date.month
   ,  input_meta->prod_date.day
   ,  input_meta->prod_date.year
   ,  input_meta->prod_date.doy);  fflush(stdout);
 */

     if ( input_meta->sat == SAT_LANDSAT_1 ||
          input_meta->sat == SAT_LANDSAT_2 ||
          input_meta->sat == SAT_LANDSAT_3     ){
        if ( input_meta->prod_date.year > 1971 && input_meta->prod_date.year < 1979 ){
          sprintf(msgbuf,"Landsat_%1.1d production year=%4.4d valid only between 1972 and 1978"  
	  , ( (int)input_meta->sat +1 ),input_meta->prod_date.year );
          RETURN_ERROR(msgbuf, "GetLut", false);
          }
        }
     if ( input_meta->sat == SAT_LANDSAT_4 ){
       if ( input_meta->acq_date.year < 1983 || 
          ( input_meta->acq_date.year == 1983 &&
            input_meta->acq_date.doy <= 91 ) ) {
          sprintf(msgbuf,"Landsat_%1.1d acq date (%2.2d/%2.2d/%4.4d) invalid before 4/1/83"  
	  , ( (int)input_meta->sat +1 )
          ,  input_meta->acq_date.month
          ,  input_meta->acq_date.day
          ,  input_meta->acq_date.year );
          RETURN_ERROR(msgbuf, "GetLut", false);
          }
        }
     if ( input_meta->sat == SAT_LANDSAT_5 ){
       if ( input_meta->acq_date.year < 1984 || 
          ( input_meta->acq_date.year == 1984 &&
            input_meta->acq_date.doy <= 314 ) ) {
          sprintf(msgbuf,"Landsat_%1.1d acq date (%2.2d/%2.2d/%4.4d) invalid before 11/19/84"  
	  , ( (int)input_meta->sat +1 )
          ,  input_meta->acq_date.month
          ,  input_meta->acq_date.day
          ,  input_meta->acq_date.year );
          RETURN_ERROR(msgbuf, "GetLut", false);
          }
         }
    }
  
  /* Compute the gains and biases for ETM */
  
  if (input_meta->inst == INST_ETM) {

    /* Get the proper LMIN/LMAX based on the production data and 
       individual band gain settting */

    if (!DateInit(&date_table_b, TABLE_B_DATE, DATE_FORMAT_DATEA)) {
      free(this->file_name);
      free(this);
      RETURN_ERROR("converting table date", "GetLut", false);
    }
    if (!DateDiff(&date_table_b, &input_meta->prod_date, &diff)) {
      free(this->file_name);
      free(this);
      RETURN_ERROR("comparing dates", "GetLut", false);
    }
    if (diff > 0.0) {  /* Dates before table change (Table a) */
      for (ib = 0; ib < nband; ib++) {
        iband = input_meta->iband[ib] - 1;
	if (input_meta->gain_set[ib] == GAIN_LOW) {
   	  lmin[ib] = etm_lmin_low_a[iband];
   	  lmax[ib] = etm_lmax_low_a[iband];
	  if (ib == 0) {
            if (input_meta->gain_setting_th == GAIN_LOW) {
 	      lmin_th = etm_lmin_low_a[iband6];
	      lmax_th = etm_lmax_low_a[iband6];
	      }
	    else
	      {
 	      lmin_th = etm_lmin_high_a[iband6];
	      lmax_th = etm_lmax_high_a[iband6];
	      }
	  }
	} else {
   	  lmin[ib] = etm_lmin_high_a[iband];
   	  lmax[ib] = etm_lmax_high_a[iband];
	  if (ib == 0) {

            if (input_meta->gain_setting_th == GAIN_LOW) {
  	      lmin_th = etm_lmin_low_a[iband6];
	      lmax_th = etm_lmax_low_a[iband6];
	      }
	    else
	      {
  	      lmin_th = etm_lmin_high_a[iband6];
	      lmax_th = etm_lmax_high_a[iband6];
	      }
	  }
	}
      }
    } else {             /* Dates after table change (Table b) */
      for (ib = 0; ib < nband; ib++) {
        iband = input_meta->iband[ib] - 1;
	if (input_meta->gain_set[ib] == GAIN_LOW) {
   	  lmin[ib] = etm_lmin_low_b[iband];
   	  lmax[ib] = etm_lmax_low_b[iband];
	  if (ib == 0) {
            if (input_meta->gain_setting_th == GAIN_LOW) {
  	      lmin_th = etm_lmin_low_b[iband6];
	      lmax_th = etm_lmax_low_b[iband6];
	      }
	    else
	      {
  	      lmin_th = etm_lmin_high_b[iband6];
	      lmax_th = etm_lmax_high_b[iband6];
	      }
	  }
	} else {
   	  lmin[ib] = etm_lmin_high_b[iband];
   	  lmax[ib] = etm_lmax_high_b[iband];
	  if (ib == 0) {

            if (input_meta->gain_setting_th == GAIN_LOW) {
  	      lmin_th = etm_lmin_low_b[iband6];
	      lmax_th = etm_lmax_low_b[iband6];
	      }
	    else
	      {
  	      lmin_th = etm_lmin_high_b[iband6];
	      lmax_th = etm_lmax_high_b[iband6];
	      }
	  }
	}
      }
    }

    /* Compute the gain and bais (see L7 data handbook) */

    qcalmin = (float)QCALMIN_LPGS;  /* Assume LPGS (for now!) */
    qcalmax = (float)QCALMAX;
        
    if ( !param->ETM_GB )
      {
      for (ib = 0; ib < nband; ib++) {
        this->meta.gain[ib] = (lmax[ib] - lmin[ib]) / (qcalmax - qcalmin);
        this->meta.bias[ib] = lmin[ib] - (this->meta.gain[ib] * qcalmin);
        }
      this->meta.gain_th = (lmax_th - lmin_th) / (qcalmax - qcalmin);
      this->meta.bias_th = lmin_th - (this->meta.gain_th * qcalmin);
      }
    else
      {
      this->meta.gain_th = this->meta.gain_th;
      this->meta.bias_th = this->meta.bias_th;
      }

  }
  
  /* Compute the coefficients for the reflectance */
  
  for (ib = 0; ib < nband; ib++) {
  
    
    iband = input_meta->iband[ib] - 1;
    switch (input_meta->sat) {
      case SAT_LANDSAT_1:
        this->esun[ib] = esun_mss_1[iband];
        if (input_meta->inst == INST_MSS && param->RE_CAL ) {  
	/*
	printf("sat=%1d band=%1d Lmax=%4.1f Lmin=%4.1f Lmax-Lmin=%4.1f (Lmax-Lmin)/254.=%f gb=(%f,%4.1f)\n"
	, input_meta->sat+1,ib+1
	,   mss_lmax_l1[ib], mss_lmin_l1[ib]
	,  (mss_lmax_l1[ib]- mss_lmin_l1[ib])
	,  (mss_lmax_l1[ib]- mss_lmin_l1[ib])/254.
	,  (mss_lmax_l1[ib]- mss_lmin_l1[ib])/25.4
	,  10.0*mss_lmin_l2[ib]
	);
	*/
	
          if ( fabs(this->meta.gain[ib]-1.0)<delta && fabs(this->meta.bias[ib])<delta ) {
	    this->meta.gain[ib]= 10.0*(mss_lmax_l1[ib]-mss_lmin_l1[ib])/254.0;
	    this->meta.bias[ib]= 10.0*mss_lmin_l1[ib];
	    }
	  else
	    {
	    if ( ib==3 ) { /* replace gain for band 4 */ 
	      this->meta.gain[ib]= 0.6024;
	      this->meta.gain[ib]= 10.0*(mss_lmax_l1[ib]-mss_lmin_l1[ib])/254.0;
	      }
	    }
	  }
        break;
      case SAT_LANDSAT_2:
        this->esun[ib] = esun_mss_2[iband];
        if (input_meta->inst == INST_MSS && param->RE_CAL) {  
	if (0)
	printf("sat=%1d band=%1d Lmax=%4.1f Lmin=%4.1f Lmax-Lmin=%4.1f (Lmax-Lmin)/254.=%f gb=(%f,%4.1f)\n"
	, input_meta->sat+1,ib+1
	,   mss_lmax_l2[ib], mss_lmin_l2[ib]
	,  (mss_lmax_l2[ib]- mss_lmin_l2[ib])
	,  (mss_lmax_l2[ib]- mss_lmin_l2[ib])/254.
	,  (mss_lmax_l2[ib]- mss_lmin_l2[ib])/25.4
	,  10.0*mss_lmin_l2[ib]
	);
	  if ( fabs(this->meta.gain[ib]-1.0)<delta && fabs(this->meta.bias[ib])<delta ) {
	    this->meta.gain[ib]= 10.0*(mss_lmax_l2[ib]-mss_lmin_l2[ib])/254.0;
	    this->meta.bias[ib]= 10.0*mss_lmin_l2[ib];
	    }
	  }
        break;
      case SAT_LANDSAT_3:
        this->esun[ib] = esun_mss_3[iband];
        if (input_meta->inst == INST_MSS && param->RE_CAL) {  
	  if ( fabs(this->meta.gain[ib]-1.0)<delta && fabs(this->meta.bias[ib])<delta ) {
	    this->meta.gain[ib]= 10.0*(mss_lmax_l3[ib]-mss_lmin_l3[ib])/254.0;
	    this->meta.bias[ib]= 10.0*mss_lmin_l3[ib];
	    }
	  }
        break;
      case SAT_LANDSAT_4:
        if (input_meta->inst == INST_TM) {
          this->esun[ib] = esun_tm_4[iband];
          if ( ib==0 ){
            this->K1 = K1_tm_4;
            this->K2 = K2_tm_4;
  	  }
	} else
          this->esun[ib] = esun_mss_4[iband];
        if (input_meta->inst == INST_MSS && param->RE_CAL) {  
	  if ( fabs(this->meta.gain[ib]-1.0)<delta && fabs(this->meta.bias[ib])<delta ) {
	    this->meta.gain[ib]= 10.0*(mss_lmax_l4[ib]-mss_lmin_l4[ib])/254.0;
	    this->meta.bias[ib]= 10.0*mss_lmin_l4[ib];
	    }
	  }
        break;
      case SAT_LANDSAT_5:
        if (input_meta->inst == INST_TM) {
          this->esun[ib] = esun_tm_5[iband];
          if ( ib==0 ){
            this->K1 = K1_tm_5;
            this->K2 = K2_tm_5;
          }
	} else
          this->esun[ib] = esun_mss_5[iband];
        if (input_meta->inst == INST_MSS && param->RE_CAL) {  
	  if ( fabs(this->meta.gain[ib]-1.0)<delta && fabs(this->meta.bias[ib])<delta ) {
	    this->meta.gain[ib]= 10.0*(mss_lmax_l5[ib]-mss_lmin_l5[ib])/254.0;
	    this->meta.bias[ib]= 10.0*mss_lmin_l5[ib];
	    }
	  }
        break;
      case SAT_LANDSAT_7:
        this->esun[ib] = esun_etm[iband];
        this->K1 = K1_etm;
        this->K2 = K2_etm;
        break;
      case SAT_NULL: break;
      case SAT_MAX: break;
    }

    if ( input->file_type == INPUT_TYPE_BINARY_2BYTE_BIG    ||
         input->file_type == INPUT_TYPE_BINARY_2BYTE_LITTLE  ) {
      this->meta.gain[ib] = 0.01;
      this->meta.bias[ib] = 0.0;
      if ( ib==0 ) {
        this->meta.gain_th = 0.01;
        this->meta.bias_th = 0.0;
      }
    }
    if (input_meta->inst == INST_MSS) {
    if(0)printf("*** MSS LUT, band=%d gain/bias=(%f,%f) ***\n",ib
     ,this->meta.gain[ib] 
     ,this->meta.bias[ib] ); fflush(stdout);
     /*
      this->meta.gain[ib] = 1.0;
      this->meta.bias[ib] = 0.0;
      */
    }
  }

  this->cos_sun_zen = cos(input_meta->sun_zen);
  jdoy = input_meta->acq_date.doy;

  for (i = 0; i < (NDSUN - 1); i++) {
    if (jdoy >= dsun_table[i].doy  &&  jdoy <= dsun_table[i+1].doy)
      break;
  }

  if (i >= (NDSUN - 1)) {
    free(this->file_name);
    free(this);
    RETURN_ERROR("finding sun-earth distance", "GetLut", false);
  }

  dsun = dsun_table[i].dsun + 
         ((jdoy - dsun_table[i].doy) * 
          ((dsun_table[i + 1].dsun - dsun_table[i].dsun) /
	   (dsun_table[i + 1].doy - dsun_table[i].doy)));
  this->dsun2 = dsun * dsun;

  if (input_meta->inst == INST_MSS) {
    for (ib = 0; ib < nband; ib++) {
    this->refl_conv[ib] = (PI * this->dsun2) / (this->esun[ib] * this->cos_sun_zen);  
    }
  }



  this->long_name_prefix_ref = DupString(LONG_NAME_PREFIX_REF);
  if (this->long_name_prefix_ref == (char *)NULL) {
    free(this);
    RETURN_ERROR("duplicating long name prefix ref", "GetLut", (Lut_t *)NULL);
  }

  this->units_ref = DupString(UNITS_REF);
  if (this->units_ref == (char *)NULL) {
    free(this);
    RETURN_ERROR("duplicating ref units", "GetLut", (Lut_t *)NULL);
  }

  this->valid_range_ref[0]=     VALID_MIN_REF;
  this->valid_range_ref[1]=     VALID_MAX_REF;
  this->scale_factor_ref=       SCALE_FACTOR_REF;
  this->scale_factor_err_ref=   SCALE_FACTOR_ERR_REF;
  this->add_offset_ref=         ADD_OFFSET_REF;
  this->add_offset_err_ref=     ADD_OFFSET_ERR_REF;
  this->calibrated_nt_ref=      CALIBRATED_NT_REF;

  this->long_name_prefix_th = DupString(LONG_NAME_PREFIX_TH);
  if (this->long_name_prefix_th == (char *)NULL) {
    free(this);
    RETURN_ERROR("duplicating long name prefix th", "GetLut", (Lut_t *)NULL);
  }

  this->units_th = DupString(UNITS_TH);
  if (this->units_th == (char *)NULL) {
    free(this);
    RETURN_ERROR("duplicating th units", "GetLut", (Lut_t *)NULL);
  }

  this->valid_range_th[0]=     VALID_MIN_TH;
  this->valid_range_th[1]=     VALID_MAX_TH;
  this->scale_factor_th=       SCALE_FACTOR_TH;
  this->scale_factor_err_th=   SCALE_FACTOR_ERR_TH;
  this->add_offset_th=         ADD_OFFSET_TH;
  this->add_offset_err_th=     ADD_OFFSET_ERR_TH;
  this->calibrated_nt_th=      CALIBRATED_NT_TH;
  this->recal_flag=            false;
  this->work_order_flag=       param->work_order_flag;
  this->work_order=            (Param_wo_t*)NULL;
  this->gnew_flag= false;
  this->gold_flag= false;

  if ( input_meta->sat == SAT_LANDSAT_5 )this->recal_flag= true;
    
  if ( input_meta->sat == SAT_LANDSAT_5 )
    {
    
      if ( input_meta->prod_date.year < RECAL_YEAR || 
	   (input_meta->prod_date.year == RECAL_YEAR &&
	    input_meta->prod_date.doy < RECAL_DOY ) )
	after_2003= false;
      else
	after_2003 = true;
 
      /* read old and new calibration coefficients ONLY if RE_CAL=yes, Feng (9/24/09) */
      if( param->RE_CAL) { 
	if ( param->gnew_flag )
	  this->gnew= get_gnew(param->gnew_file_name,&input_meta->acq_date);
	else
	  {
	    ERROR("*** NO GNEW File supplied ***","GetLut");
	  }
	if ( this->gnew != (Gains_t*)NULL )
	  this->gnew_flag= true;
	else
	  {
	    sprintf(msgbuf,"attempt to get gnew from file \"%s\"",
		    param->gnew_file_name);
	    ERROR(msgbuf,"GetLut");
	  }
	if ( after_2003  )
	  {       
	    if ( param->gold_2003_flag )
	      this->gold= get_gnew(param->gold_2003_name,&input_meta->acq_date);
	    else
	      {
		ERROR("*** NO GOLD_2003 File supplied ***","GetLut");
	      }
	    if ( this->gold != (Gains_t*)NULL )
	      this->gold_flag= true;
	    else
	      {
		sprintf(msgbuf,"attempt to get gold from file \"%s\"",
			param->gold_2003_name);
		ERROR(msgbuf,"GetLut");
	      }
	  }
	else
	  {
	    if ( param->gold_flag )
	      this->gold= get_gold(param->gold_file_name,&input_meta->acq_date);
	    else
	      {
		ERROR("*** NO GOLD File supplied ***","GetLut");
	      }
	    if ( this->gold != (Gains_t*)NULL )
	      this->gold_flag= true;
	    else
	      {
		sprintf(msgbuf,"attempt to get gold from file \"%s\"",
			param->gold_file_name);
		ERROR(msgbuf,"GetLut");
	      }
	  }

	/* remove date rule, but use input flag from "RE_CAL" */  
	/*if ( ( input_meta->provider == PROVIDER_EDC  ||
	  input_meta->provider == PROVIDER_EROS ||
	  input_meta->provider == PROVIDER_UMD ) &&
	  ( input_meta->prod_date.year > NORECAL_YEAR || 
	  ( input_meta->prod_date.year == NORECAL_YEAR &&
	  input_meta->prod_date.doy >= NORECAL_DOY ) ) )
	  this->recal_flag= false;
      
	  if ( input_meta->provider != PROVIDER_EDC   &&
	  input_meta->provider != PROVIDER_EROS  &&
	  input_meta->provider != PROVIDER_UMD )
	  this->recal_flag= false;*/
      }
    }
  
  if ( input_meta->sat != SAT_LANDSAT_5 )this->recal_flag= false;
  /* do not apply re-calibration if lndcal*.txt says no - fgao in 3/23/2009 */ 
  if ( !param->RE_CAL )this->recal_flag= false; 

  return this;
}

bool FreeLut(Lut_t *this) {

  if (this != (Lut_t *)NULL) {
    if (this->file_name != (char *)NULL) free(this->file_name);
    free(this);
  }

  return true;
}

