#include "cal.h"
#include "const.h"
#define nint(A)(A<0?(int)(A-0.5):(int)(A+0.5))
/*
#define PRINT_X (1543)
#define PRINT_Y (  58)*/
#define PRINT_X (   0)
#define PRINT_Y (   0)

bool Cal(Lut_t *lut, int iband, Input_t *input, unsigned char *line_in, 
         int *line_out, int *line_out_qa, Cal_stats_t *cal_stats, int iy) {
  int is,val,jband;
  float gain, bias, rad, ref_conv, ref, fval;
  int nsamp= input->size.s;
  int ifill= input->short_flag ? FILL_VAL[iband]: (int)lut->in_fill;
  float scale;
  float Gold, Gnew, Alpha, Grescale;
  jband= iband==5 ? 6 : iband;

  if ( input->dnout )
    scale=(input->dn_map[3]-input->dn_map[2]+1)/
          (input->dn_map[1]-input->dn_map[0]+1);

  if ( lut->recal_flag && lut->gnew_flag && 
     ( lut->gold_flag || lut->work_order_flag ) )
    {
    jband= iband==5 ? 6 : iband;
    if (  lut->work_order_flag )
      {
      Alpha= lut->work_order->final_gain[jband];
      Grescale= lut->work_order->DN_to_Radiance_gain[jband];
      Gold=  Alpha / Grescale;
      }
    else
      Gold= lut->gold->gains[jband];
    Gnew= lut->gnew->gains[jband];
    }
  else 
    {
    Gold= Gnew= 1.0;
    }

  gain = lut->meta.gain[iband];
  bias = ( lut->meta.bias[iband] * ( Gold / Gnew ));

  ref_conv = (PI * lut->dsun2) / (lut->esun[iband] * lut->cos_sun_zen);
  
  if ( iy==0 )
    {
    printf("*** band=%1d gain=%f bias=%f ref_conv=%f=(PI*%f)/(%f*%f) ***\n",iband+1,gain,bias,ref_conv
    ,lut->dsun2 ,lut->esun[iband],   lut->cos_sun_zen
    );
    fflush(stdout);
    }

  for (is = 0; is < nsamp; is++) {

    val= getValue((unsigned char *)line_in, is, input->short_flag, input->swap_flag );
    /*if ( PRINT_X!=0 && PRINT_Y!=0 && is==PRINT_X && iy==PRINT_Y  ) 
      printf("val = %d ifill=%d qa_val=%d qa_fill=%d\n", val, ifill, line_out_qa[is], lut->qa_fill);*/

    if (val == ifill || line_out_qa[is]==lut->qa_fill ) {
      line_out[is] = lut->out_fill;
      cal_stats->nfill[iband]++;
      continue;
    }

    /* for saturated pixels, added by Feng (3/23/09) */
    if (val == SATU_VAL[iband]) {
      line_out[is] = lut->out_satu;
      continue;
    }

    cal_stats->nvalid[iband]++;
    fval= (float)val;

/*  if (  cal_stats->nvalid[iband]== 1  )  */
    if ( PRINT_X!=0 && PRINT_Y!=0 && is==PRINT_X && iy==PRINT_Y  ) 
      {
      jband= iband==5 ? 6 : iband;
      printf("%1d xy=(%4d,%4d) input_dn=%f ",jband,is,iy,fval);
      }
      
    if ( lut->gnew_flag && ( lut->gold_flag || lut->work_order_flag ) )
      fval*= ( Gold / Gnew );

    if ( lut->gnew_flag && ( lut->gold_flag || lut->work_order_flag ) )
    if ( PRINT_X!=0 && PRINT_Y!=0 && is==PRINT_X && iy==PRINT_Y  ) 
      {
      printf(" Gold=%f Gnew=%f Gold/Gnew=%f\n",Gold,Gnew,Gold/Gnew); fflush(stdout);
      }

    rad = (gain * fval) + bias;
    ref = rad * ref_conv;
    
    if ( PRINT_X!=0 && PRINT_Y!=0 && is==PRINT_X && iy==PRINT_Y  ) 
      {
      printf("  fval=%f gain=%f bias*(Gold/Gnew)=%f rad=%f ref_conv=%f ref=%f\n",fval,gain,bias,rad,ref_conv,ref);
      fflush(stdout);
      }

    if ( input->dnout )
      line_out[is]=
        (int)(input->dn_map[2]+(((float)val-input->dn_map[0])*scale)+0.5);
    else
      line_out[is] = (int)(ref * 10000.0) + 0.5;

    if (cal_stats->first[iband]) {
      
      cal_stats->idn_min[iband] = val;
      cal_stats->idn_max[iband] = val;

      cal_stats->rad_min[iband] = rad;
      cal_stats->rad_max[iband] = rad;

      cal_stats->ref_min[iband] = ref;
      cal_stats->ref_max[iband] = ref;

      cal_stats->iref_min[iband] = line_out[is];
      cal_stats->iref_max[iband] = line_out[is];

      cal_stats->first[iband] = false;

    } else {

      if (val < cal_stats->idn_min[iband]) 
        cal_stats->idn_min[iband] = val;
      if (val > cal_stats->idn_max[iband]) 
        cal_stats->idn_max[iband] = val;

      if (rad < cal_stats->rad_min[iband]) cal_stats->rad_min[iband] = rad;
      if (rad > cal_stats->rad_max[iband]) cal_stats->rad_max[iband] = rad;

      if (ref < cal_stats->ref_min[iband]) cal_stats->ref_min[iband] = ref;
      if (ref > cal_stats->ref_max[iband]) cal_stats->ref_max[iband] = ref;

      if (line_out[is] < cal_stats->iref_min[iband]) 
        cal_stats->iref_min[iband] = line_out[is];
      if (line_out[is] > cal_stats->iref_max[iband]) 
        cal_stats->iref_max[iband] = line_out[is];
    }
  }

  return true;
}

bool Cal6(Lut_t *lut, Input_t *input, unsigned char *line_in, int *line_out, 
         int *line_out_qa, Cal_stats6_t *cal_stats, int iy) {
  int is, val;
  float gain, bias, rad, temp;
  int nsamp= input->size_th.s;
  int ifill= input->short_flag ? FILL_VAL6 : (int)lut->in_fill;
  float scale;
  

  if ( input->dnout )
    scale=(input->dn_map[3]-input->dn_map[2]+1)/(input->dn_map[1]-input->dn_map[0]+1);

  gain = lut->meta.gain_th;
  bias = lut->meta.bias_th;
  
  if ( iy==0 )
    {
    printf("*** band=%1d gain=%f bias=%f ***\n",6,gain,bias);
    fflush(stdout);
    }

  for (is = 0; is < nsamp; is++) {
    val= getValue((unsigned char *)line_in, is, input->short_flag, input->swap_flag );
    if (val == ifill || line_out_qa[is]==lut->qa_fill ) {
      line_out[is] = lut->out_fill;
      cal_stats->nfill++;
      continue;
    }
    cal_stats->nvalid++;
 
/* Put 100.0 (scale_factor_th) and 273.15 (th_zero_celcius_in_degrees_kelvin)
   into lut. */

    rad = (gain * (float)val) + bias;
    temp = ( lut->K2 / log(1.0 + (lut->K1/rad) ) );
    temp -= 273.15;

    if ( input->dnout )
      {
      line_out[is]=(int)(input->dn_map[2]+(((float)val-input->dn_map[0])*scale)+0.5);
      }
    else
      line_out[is] = (int)(temp * 100.0 + 0.5);

    if (cal_stats->first) {
      
      cal_stats->idn_min = val;
      cal_stats->idn_max = val;

      cal_stats->rad_min = rad;
      cal_stats->rad_max = rad;

      cal_stats->temp_min = temp;
      cal_stats->temp_max = temp;

      cal_stats->itemp_min = line_out[is];
      cal_stats->itemp_max = line_out[is];

      cal_stats->first = false;

    } else {

      if (val < (int)cal_stats->idn_min) 
        cal_stats->idn_min = val;
      if (val > cal_stats->idn_max) 
        cal_stats->idn_max = val;

      if (rad < cal_stats->rad_min) cal_stats->rad_min = rad;
      if (rad > cal_stats->rad_max) cal_stats->rad_max = rad;

      if (temp < cal_stats->temp_min) cal_stats->temp_min = temp;
      if (temp > cal_stats->temp_max) cal_stats->temp_max = temp;

      if (line_out[is] < cal_stats->itemp_min) 
        cal_stats->itemp_min = line_out[is];
      if (line_out[is] > cal_stats->itemp_max) 
        cal_stats->itemp_max = line_out[is];
    }
  }

  return true;
}
/*************************************************************************
 *** this program returns the correct value (as an int) after swapping ***
 *** if if indicated                                                   ***
 *************************************************************************/
 int getValue(unsigned char* line_in, int ind, int short_flag, int swap_flag)
 {
 unsigned short int* line_in2_byte= (unsigned short int*)line_in;
 unsigned short int i2_temp;
 unsigned short int out_value;
 char* inc= (char*)&i2_temp;
 char* outc= (char*)&out_value;
 if ( short_flag )
   {
   if ( swap_flag ) 
     {
     i2_temp= (int)line_in2_byte[ ind ];
     outc[0]= inc[1];
     outc[1]= inc[0];
     }
   else
     out_value= (int)line_in2_byte[ ind ];
   }
 else
   out_value= (int)line_in[ ind ];
 return (int)out_value;
}
