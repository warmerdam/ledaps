#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h> 
#include <time.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "csm.h"
#include "const.h"
#include "error.h"
#include "util.h"
#define WRITE_SIEVE 0
#define LOG_FLAG 1
#define HIST_LOG_FLAG 1
#define CLMASK_LOG ("cloud_mask.log")
FILE* mout;      /* most of the log messages go here                         */
Pbuf_t *pst;     /* Print message log */
char pb[256];
bool CloudMask(  
 	         Input_t *input         /* input file for reflective bands   */
               , Input_t *input_th      /* input file for thermal band       */
               , Output_t *output       /* output file                       */
               , Lut_t *lut             /* lut file                          */
               , Param_t *param         /* paramater file                    */
               , int odometer_flag      /* odometer flag                     */
              )
{
/*--------------------------------------------------------------------------!*/
/*-                               constants                                -!*/
/*--------------------------------------------------------------------------!*/
 const double deg2rad=  0.01745329251994330 ; 
/*--------------------------------------------------------------------------!*/
/*-                            one byte consts                             -!*/
/*--------------------------------------------------------------------------!*/
 const unsigned char b255= 255;
 const unsigned char b254= 254;
 const unsigned char b200= 200;
 const unsigned char b125= 125;
 const unsigned char b55 =  55;
 const unsigned char b1  =   1;
 const unsigned char b0  =   0;
/*--------------------------------------------------------------------------!*/
/*-                            histogram consts                            -!*/
/*--------------------------------------------------------------------------!*/
 const float cloud_mean_temp = 295.;
 const int HMIN= 1;
 const int HMAX= 1024;
 const int HISTSIZ= HMAX-HMIN+1;
 const int bfact= 1024;
/*--------------------------------------------------------------------------!*/
/*-                              cloud masks                               -!*/
/*--------------------------------------------------------------------------!*/
 short int* clmaski2;
 unsigned char* clmask;
 unsigned char* clmaskb;
 unsigned char* clmaskb2;
 unsigned char* bmask_img;
 unsigned char* imask_img;
 unsigned char* omask_img;
/*--------------------------------------------------------------------------!*/
/*-                          cloud mask pointers                           -!*/
/*--------------------------------------------------------------------------!*/
 unsigned char* clmaskb2p;
 unsigned char* clmaskb2c;
 unsigned char* clmaskb2m;
 unsigned char maskvalue;
 unsigned char* sive_cflag;
 unsigned int* sive_flag;
/*--------------------------------------------------------------------------!*/
/*-                                  ints                                  -!*/
/*--------------------------------------------------------------------------!*/
   int i						
     , ix,iy
     , totpix
     , xhalf
     , yhalf
     , cloud[2*2]
     , b3tab 
     , ndsitab 
     , hotpixels 
     , count1 
     , count2 
     , count3 
     , count4 
     , count4x
     , count9 
     , count10
     , count42					           /* int */
     , b6delete
     , count6
     , snowtab 
     , landtab 
     , index 
     , index2 
     , bins
     , ibin
     , bins1
     , bins2
     , b6histogram[HISTSIZ]
     , b6histogram2[HISTSIZ]
     , histmin, histmax, histnum
     , cloudvalues 
     , cloudvalues2 
     , cloudval254
     , cloudval255
     , iclpercentorig                                      /* int */
     , cloudcut0125
     , cloudcut05
     , cloudcut025
     , cloudcut175
     , cloudsum 
     , hindex
     , filter_used
     , isum_cloud
     , isum_snow
     , ixp
     , ixm
     , iyp
     , iym
     , filt_new
     , new_snow
     , nps
     , nls
     , p
     , ib
     , ihi
     , *therm_line
     , *line_in_buf
     , *line_in[NBAND_REFL_MAX]
     , ik
     , ival
     , oval
     , yy,xx,yyy,xxx
     , ksize
     , bufptr=0
       ;                                                   /* int */
/*--------------------------------------------------------------------------!*/
/*-                                  doubles                               -!*/
/*--------------------------------------------------------------------------!*/
 double *ib2buf
      , *ib3buf
      , *ib4buf
      , *ib5buf
      , *ib6buf
    ;
/*--------------------------------------------------------------------------!*/
/*-                                  floats                                -!*/
/*--------------------------------------------------------------------------!*/
 float  *b6tempeture
      , sun_el
      , b3_thresh
      , b56_thresh
      , c_cover_tot
      , c_cover[2*2]
      , b2rflect
      , b3rflect
      , b4rflect
      , b5rflect 
      , b42ratio
      , b43ratio
      , b45ratio                                            /* float */
      , ndsi
      , clpercent
      , clpercentorig
      , clpercent254
      , clpercentorig1 
      , clpercentpass1
      , snpercent
      , thermeffecta
      , thermeffectb
      , clpercenta 
      , clpercentb 
      , result[6]
      , result1[6]
      , result2[6]
      , result_temp[6]
      , b6min
      , b6max
      , b6max_thresh1_diff 
      , b6min1
      , b6max1
      , b6min2
      , b6max2
      , b6sd
      , shiftfact
      , shift
      , b6mean
      , b6mean1
      , b6mean2
      , b6thresh1
      , b6thresh2
      , b6thresh3
      , b6thresh_max
      , b6thresh1_save 
      , b6thresh2_save 
      , desert_index                                        /* float */
      , mean_std[2]
      ;
/*--------------------------------------------------------------------------!*/
/*-                             threshold values                           -!*/
/*--------------------------------------------------------------------------!*/
  const float di_threshold = 0.5;
  const float b3thresh =0.08;
  const float b6_thresh =300.0;
  const float b45_thresh = 1.0; /*  1.075  */
  const float b42_thresh = 2.0;
  const float b43_thresh = 2.0;
  const float thresh_b56_hi = 225.0;
  const float thresh_b56_lo = 210.0;
  const float thresh_NDSI_max= 0.7;
  const float thresh_NDSI_min = -0.25;
  const float thresh_b3_lower = .07; /* 0.06; */
  const float thermal_effect_high= 35.0;
  const float thermal_effect_low = 25.;
  const float thresh_NDSI_snow = 0.8;
  const float snow_thresh = 1.0;
  const float cloud_254_thresh= 0.2;
/*--------------------------------------------------------------------------!*/
/*-                         image and virbuf structs                       -!*/
/*--------------------------------------------------------------------------!*/
 vbuf_t b6clouds;
 vbuf_t b6clouds200;
 vbuf_t b6clouds254;
 imgfile_t b6tempeture_File;
 imgfile_t clmask_File;
 imgfile_t sive_File;
 imgfile_t clmaskb_File;
 imgfile_t clmaskb2_File;
/*--------------------------------------------------------------------------!*/
/*-                            output flag values                          -!*/
/*--------------------------------------------------------------------------!*/
 unsigned char CLSTAT_S= lut->csm_snow;  /* snow                           -!*/
 unsigned char CLSTAT_L= lut->csm_land;  /* land  ( non-snow )             -!*/
 unsigned char CLSTAT_F= lut->csm_fill;  /* fill                           -!*/
 unsigned char CLSTAT_C= lut->csm_cloud; /* cloud                          -!*/
 int FILL_VALUE=         lut->in_fill;   /* -9999.000;                     -!*/
 float refl_scale_factor=lut->refl_scale_factor;   /* 10000.0;             -!*/
 float th_scale_factor=  lut->therm_scale_factor;  /*   100.0;             -!*/
 double th_zero_celcius_in_degrees_kelvin = (double)lut->celc2kelvin; /* 273 */
/*--------------------------------------------------------------------------!*/
/*-             do the initalization for the l7 cloud algorithm            -!*/
/*--------------------------------------------------------------------------!*/
 time_t clock;
 char date_time[21];
 bool therm_flag;
/*--------------------------------------------------------------------------!*/
 therm_flag=( input_th!=(Input_t *)NULL );

 if ( LOG_FLAG ) mout= fopen(CLMASK_LOG,"w"); else  mout= stdout;
 pst= InitPr(mout, pb, 79, '*');  
 cenb(pst,"*");  cenb(pst,"Cloud Mask Generation Log");
 clock = time( NULL ) ;
 strftime(date_time, 20,"%Y:%m:%d %H:%M:%S", localtime( &clock ) ) ;
 date_time[ 19 ]=(char)0;
 sprintf(pb,"created : %s",date_time); cenb(pst,pb); cenb(pst," ");
/*--------------------------------------------------------------------------!*/
/*-                       write sun and path/row to log                    -!*/
/*--------------------------------------------------------------------------!*/
 sprintf(pb,"wrs=(%3.3d,%3.3d)",input->meta.ipath,input->meta.irow); cen(pst);
 sun_el=      ( 90.0- ( input->meta.sun_zen/deg2rad ) );
 sprintf(pb,"doy=%d sun_el=%f (degrees) zenith (radians)=%f",
   input->meta.acq_date.doy,sun_el,input->meta.sun_zen); cen(pst); 
/*--------------------------------------------------------------------------!*/
/*-                               b3 threshold                             -!*/
/*--------------------------------------------------------------------------!*/
 b3_thresh = .08;
/*--------------------------------------------------------------------------!*/
/*-                           initalize image size                         -!*/
/*--------------------------------------------------------------------------!*/
 nps= input->size.s;
 nls= input->size.l;
 xhalf= nps / 2;
 yhalf= nls / 2;
/*--------------------------------------------------------------------------!*/
/*-            open for write (read/write) temporary image files           -!*/
/*--------------------------------------------------------------------------!*/
 if(!img_openw(&b6tempeture_File,"b6tempeture_temp.img",nps,nls,1,sizeof(float)
   ) )  ERROR("opening b6tempeture","csm");
 if(!img_openw( &clmaskb_File,  "clmaskb_temp.img",  nps, nls,1,sizeof(char) 
    ) )  ERROR("opening clmaskb","csm");
 if(!img_openw( &clmaskb2_File, "clmaskb2_temp.img", nps, nls,1,sizeof(char) 
    ) )  ERROR("opening clmaskb2","csm");
 if(!img_openw( &clmask_File, "clmask.img", nps, nls,1,sizeof(char)  ) ) 
   ERROR("opening clmask","csm");

/*--------------------------------------------------------------------------!*/
/*-                     initalize (open) virtual buffers                   -!*/
/*--------------------------------------------------------------------------!*/
 if ( !virinit( &b6clouds, "b6clouds.temp.file",bfact ) )
   ERROR("Error opening virtual buffer b6clouds","csm" );
 if ( ! virinit( &b6clouds254, "b6clouds254.temp.file",bfact ) )
   ERROR("Error opening virtual buffer b6clouds254","csm" );
/*--------------------------------------------------------------------------!*/
/*-                             initialize cloud                           -!*/
/*--------------------------------------------------------------------------!*/
  memset(cloud,0,sizeof(int)*2*2);
/*--------------------------------------------------------------------------!*/
/*-                         Allocate memory for lines                      -!*/
/*--------------------------------------------------------------------------!*/
  line_in_buf = (int *)calloc((size_t)(input->size.s * input->nband), 
                              sizeof(int));

  if ( line_in_buf==(int*)NULL )ERROR("allocating input line buffer","csm");
  line_in[0] = line_in_buf;
  for (ib = 1; ib < input->nband; ib++)
    line_in[ib] = line_in[ib - 1] + input->size.s;

 b6tempeture  = (float*)malloc ( nps* sizeof(float) );
 if ( !b6tempeture )ERROR(" allocate b6tempeture failed ","csm"); 
 clmaski2 =  (short int*)malloc ( nps* sizeof(short int) );
 if ( !clmaski2 )ERROR(" allocate clmaski2 failed ","csm"); 
 clmask=  (unsigned char*)malloc ( nps* sizeof(unsigned char) ); 
 if ( !clmask )ERROR(" allocate clmask failed ","csm"); 
 clmaskb=  (unsigned char*)malloc ( nps* sizeof(unsigned char) ); 
 if ( !clmaskb )ERROR(" allocate clmaskb failed ","csm"); 
 clmaskb2=  (unsigned char*)malloc ( nps*4* sizeof(unsigned char) ); 
 if ( !clmaskb2 )ERROR(" allocate clmaskb2 failed ","csm"); 
 ib2buf= (double*)malloc ( nps* sizeof(double) ); 
 if ( !ib2buf )ERROR(" allocate ib2buf failed ","csm"); 
 ib3buf= (double*)malloc ( nps* sizeof(double) ); 
 if ( !ib3buf )ERROR(" allocate ib3buf failed ","csm"); 
 ib4buf= (double*)malloc ( nps* sizeof(double) ); 
 if ( !ib4buf )ERROR(" allocate ib4buf failed ","csm"); 
 ib5buf= (double*)malloc ( nps* sizeof(double) ); 
 if ( !ib5buf )ERROR(" allocate ib5buf failed ","csm"); 
 ib6buf= (double*)malloc ( nps* sizeof(double) ); 
 if ( !ib6buf )ERROR(" allocate ib6buf failed ","csm"); 

 therm_line= (int*)malloc ( nps* sizeof(int) ); 
 if ( !therm_line )ERROR(" allocate therm_line failed ","csm"); 
 memset(therm_line,0,sizeof(int)*nps);


/*--------------------------------------------------------------------------!*/
      
/*--------------------------------------------------------------------------!*/
/*--------------------------------------------------------------------------!*/
/*-                        allocate the kernel buffers                     -!*/
/*--------------------------------------------------------------------------!*/
 if ( param->apply_kernel>0 )
   {
   bmask_img=  (unsigned char*)malloc ( 2*nps*nls* sizeof(unsigned char) ); 
   if ( !bmask_img )ERROR("alloacate both lmask img (kernel) failed ","csm"); 
   memset(bmask_img,0,2*nps*nls*sizeof(unsigned char));
   imask_img= & bmask_img[ ((bufptr++)%2)*nps*nls ];
   omask_img= & bmask_img[ ((bufptr  )%2)*nps*nls ]; 
   }

/*--------------------------------------------------------------------------!*/
/*-                         allocate the sieve buffers                     -!*/
/*--------------------------------------------------------------------------!*/
 if ( param->sieve_thresh>0 )
   {
   if ( WRITE_SIEVE )
     if(!img_openw( &sive_File, "sive.tif", nps, nls,1,sizeof(char)) ) 
       ERROR("opening sive_File","csm");

   sive_flag=(unsigned int*)malloc( nps*nls*sizeof(unsigned int) ); 
   if( sive_flag==(unsigned int*)NULL)ERROR("alloc sive flag","csm");
   memset((char*)sive_flag,0,(nps*nls)*sizeof(unsigned int));
 
   sive_cflag=(unsigned char*)malloc( nps*nls*sizeof(unsigned char) ); 
   if( sive_cflag==(unsigned char*)NULL)ERROR("alloc sive flag","csm");
   memset((char*)sive_cflag,0,(nps*nls)*sizeof(unsigned char));
   }

/*--------------------------------------------------------------------------!*/
/*-                         initalize some variables                       -!*/
/*--------------------------------------------------------------------------!*/
 b3tab = 0;
 ndsitab = 0;
 hotpixels = 0;
 count2= 0;
 count3= 0;
 count4= 0;
 count4x= 0;
 count42= 0;
 count1= 0;
 b6delete = 0;
 count6 = 0;
 ndsitab = 0;
 snowtab = 0;
/*--------------------------------------------------------------------------!*/
/*-                       main loop for each line (iy)                     -!*/
/*--------------------------------------------------------------------------!*/
 for (iy=0; iy<nls; iy++ )
   {
/*--------------------------------------------------------------------------!*/
/*-                  read in next image line for all 5 bands               -!*/
/*--------------------------------------------------------------------------!*/
   for (ib = 0; ib < input->nband; ib++) 
     {
     if (!GetInputLine(input, ib, iy, line_in[ib]))
        ERROR("reading input data for a line", "csm");
     }

   if ( therm_flag ){
     if (!GetInputLine(input_th, 0, iy, therm_line))
      ERROR("reading input data for a line", "csm");
   }

   if ( ( iy==0 || iy ==(nls-1) || iy%100==0 ) && odometer_flag )
     { 
     printf("--- main loop reading in line %d of %d --- \r",iy,nls-1); 
     fflush(stdout); 
     }

/*--------------------------------------------------------------------------!*/
/*-                      main loop for each sample (ix)                    -!*/
/*--------------------------------------------------------------------------!*/
   for ( ix=0; ix<nps; ix++)
     {
     clmaski2[ix]=           0;
     clmaskb[ix]=           b0;
/*--------------------------------------------------------------------------!*/
/*- Output from lndcal is in scaled reflectance * 10000, so convert to     -!*/
/*- unscaled reflectance                                                   -!*/
/*--------------------------------------------------------------------------!*/
     ib2buf[ix]= ( (float)line_in[1][ix] ) / refl_scale_factor;  /* band 2 */
     ib3buf[ix]= ( (float)line_in[2][ix] ) / refl_scale_factor;  /* band 3 */
     ib4buf[ix]= ( (float)line_in[3][ix] ) / refl_scale_factor;  /* band 4 */
     ib5buf[ix]= ( (float)line_in[4][ix] ) / refl_scale_factor;  /* band 5 */
     ib6buf[ix]= ( (float)therm_line[ix] ) / th_scale_factor;    /* band 6 */
     b6tempeture[ix]= (float)(th_zero_celcius_in_degrees_kelvin + ib6buf[ix]); 
/*--------------------------------------------------------------------------!*/
/*-                          initalize fill value                          -!*/
/*--------------------------------------------------------------------------!*/
     clmask[ix]=  (  line_in[1][ix] == FILL_VALUE ||
                     line_in[2][ix] == FILL_VALUE ||
                     line_in[3][ix] == FILL_VALUE ||
                     line_in[4][ix] == FILL_VALUE ||
                     therm_line[ix] == FILL_VALUE    )?  CLSTAT_F : CLSTAT_L;

     if (  line_in[1][ix] == FILL_VALUE )
       {
       ib2buf[ix]=      0.0;
       ib3buf[ix]=      0.0;
       ib4buf[ix]=      0.0;
       ib5buf[ix]=      0.0;
       ib6buf[ix]=      0.0;
       b6tempeture[ix]= 0.0;
       }
    
     }/* for ( ix=0; ix<nps; ix++) */

   if ( therm_flag ) {

   for ( ix=0; ix<nps; ix++)
     {
     b3rflect = ib3buf[ix]; 

/*--------------------------------------------------------------------------!*/
/*- band 2 _threshance threshold - screens out bright targets.             -!*/
/*- set at .12, changed to .15 for the mexico image, to .16 for 34-38      -!*/
/*--------------------------------------------------------------------------!*/

     if ( b3rflect > b3_thresh ) 
       {
       b2rflect = ib2buf[ix];
       b5rflect = ib5buf[ix];

       ndsi = (b2rflect - b5rflect) / (b2rflect + b5rflect );
       count1 = count1 + 1;
/*--------------------------------------------------------------------------!*/
/*- normalized difference snow index. less than .4? snow eliminated,       -!*/
/*- but not ice.                                                           -!*/
/*--------------------------------------------------------------------------!*/
       if (ndsi < thresh_NDSI_max && ndsi > thresh_NDSI_min)
         {
         count2 = count2 + 1;

         if(b6tempeture[ix] < b6_thresh ) 
           {
           b6delete++;
           b56_thresh = ( 1.0 - b5rflect ) * b6tempeture[ix];

/*--------------------------------------------------------------------------!*/
/*band 6/5 composite. eliminates ice - originally set at 225 (thresh_b56_hi)!*/
/*--------------------------------------------------------------------------!*/
         if( b56_thresh < thresh_b56_hi ) 
	   {
           count6++;
/*--------------------------------------------------------------------------!*/
/*- band 4/3 ratio. (equals about 1.0 for clouds), eliminates bright veg   -!*/
/*- and soil.                                                              -!*/
/*- 2.2 worked well except for cirrus, changed to 2.6 36/34 for cirrus     -!*/
/*- type pixel location                                                    -!*/
/*--------------------------------------------------------------------------!*/
             b4rflect = ib4buf[ix];
             b43ratio = b4rflect / b3rflect;

             if ( b43ratio < b43_thresh )
               {
               count3 = count3 + 1;
/*--------------------------------------------------------------------------!*/
/*- 4/2 ratio for senescing vegetation, chlorophyll absorbtion lacking     -!*/
/*--------------------------------------------------------------------------!*/
               b42ratio = b4rflect / b2rflect;

               if ( b42ratio < b42_thresh ) 
                 {
                 count42 = count42 + 1;
                 b45ratio = b4rflect / b5rflect;
/*--------------------------------------------------------------------------!*/
/*-  band 4/5 ratio. eliminates rocks and desert                           -!*/
/*--------------------------------------------------------------------------!*/
                 if ( b45ratio > b45_thresh ) 
                   {
                   count4 = count4 + 1;
                   clmaski2[ix] = 255;
                   clmaskb[ix] = b255;

                   if ( b56_thresh < thresh_b56_lo )
                     {
                     count4x = count4x + 1;
                     clmaskb[ix] = b254;


                     }                 /* b56_thresh*/
                   }                  /* b45 ratio endif*/
                 }                   /* b42 endif*/
               }                    /* b43 endif*/
             }                     /* b56 endif*/
           else
	     {
             if ( b5rflect < 0.08 )clmaskb[ix] = b55;
	     }
           }                      /* b6 endif*/
         else
           {
           clmaski2[ix]= 55;
           clmaskb[ix] = b55;
           hotpixels = hotpixels + 1;
           }
         }
       else 
         {
         if ( ndsi > thresh_NDSI_snow )
           {
           snowtab = snowtab + 1;
           clmask[ix]= CLSTAT_S;
           }
         clmaski2[ix]= 55;
         clmaskb[ix] = b55;
         ndsitab = ndsitab + 1;
         }
       }
     else 
       {
       if ( b3rflect < thresh_b3_lower )
         {
         clmaski2[ix]= 55;
         clmaskb[ix] = b55;
         b3tab = b3tab + 1;
         }
       }                           /* b3_thresh endif*/
     if( clmaskb[ix]>=b254 )cloud[(iy/yhalf)*2+ix/xhalf]++; 
     }                             /* for ( ix=0; ix<nps; ix++)*/
/*--------------------------------------------------------------------------!*/
/*-                   write b6tempeture and clmask line                    -!*/
/*--------------------------------------------------------------------------!*/
   if ( !put_line(&b6tempeture_File, (char*)b6tempeture, iy ) )
       ERROR("putline b6tempeture file","csm" );
   if ( !put_line(&clmaskb_File,   (char*)clmaskb,   iy ) )
       ERROR("putline clmaskb file","csm" );
   if ( !put_line(&clmask_File, (char*)clmask, iy ) )
     ERROR("putline clmask file","csm" );

   } else {
    if (!PutOutputLine(output, CLMASK, iy, clmask))
      ERROR("writing output data for a line (CLMASK)", "csm");
   }

   }                  /* enddo    // iy // endl clmaskb*/
 if ( odometer_flag )printf("\n");

 if ( !therm_flag ){
   if ( odometer_flag )printf("*** freeing it ***\n");
   goto FREE_IT;
 }

/*--------------------------------------------------------------------------!*/
/*-                      close all input image files                       -!*/
/*-        compute quadrant cloud cover percent for log file output        -!*/
/*--------------------------------------------------------------------------!*/
 c_cover_tot= 0.0;
 for ( iy=0; iy<2 ; iy++ )
   for ( ix=0; ix<2 ; ix++ )
     {
     c_cover[iy*2+ix]= ((float)cloud[iy*2+ix])*100.0/((float)(xhalf*yhalf));
     c_cover_tot+= c_cover[iy*2+ix]/4.0;
     }
 totpix= nps * nls;
/*--------------------------------------------------------------------------!*/
/*-                          first pass complete                           -!*/
/*--------------------------------------------------------------------------!*/
 landtab = b3tab + ndsitab + hotpixels ;
 cloudvalues = 0;
/*--------------------------------------------------------------------------!*/
/*-              index = where(clmaski2 eq 255, cloudvalues)               -!*/
/*--------------------------------------------------------------------------!*/
 cloudval254= 0;
 cloudvalues= 0;

 for (iy=0; iy<nls; iy++ )
   {
   if ( !get_line(&b6tempeture_File,(unsigned char*)b6tempeture, iy ) )
     ERROR("getline from b6tempeture","csm" );
   if ( !get_line(&clmaskb_File,  (unsigned char*)clmaskb,   iy ) )
         ERROR("getline from clmaskb","csm" );

   for ( ix=0; ix<nps; ix++)
     {
     if ( clmaskb[ix] == b254 )
       {
       virput( &b6clouds254, b6tempeture[ix] );
       cloudval254= cloudval254 + 1;
       }
     if ( clmaskb[ix] == b255 || clmaskb[ix] == b254   )
       {
       virput( &b6clouds, b6tempeture[ix] );
       cloudvalues++;
       }
     }

   if ( ( iy==0 || iy ==(nls-1) || iy%100==0 ) && odometer_flag )
     {
     printf("--- counting cloudvalues line %d of %d --- \r",iy,nls-1); 
     fflush(stdout); 
     }
   }
 if ( odometer_flag )printf("\n");
 virflush( &b6clouds );
 virflush( &b6clouds254 );

 clpercentorig1=  (float)cloudvalues / (float)totpix * 100.0;
 clpercent254=    (float)cloudval254 / (float)totpix * 100.0;
 clpercentorig =  c_cover_tot;
 iclpercentorig = nint( clpercentorig );

 cenb(pst," "); cenb(pst,"-"); cenb(pst," ");
 sprintf(pb,"clpercentorig1&2=%f,%f",clpercentorig,clpercentorig1);   pr(pst);
 sprintf(pb,"acca scores= (%5.3f,%5.3f,%5.3f,%5.3f) Total=%d",
        c_cover[0],c_cover[1],c_cover[2],c_cover[3],iclpercentorig);  pr(pst);

 sprintf(pb,"-> cloud count before thermal filter: %d,%f",
         cloudvalues, clpercentorig );                 pr(pst);
 sprintf(pb,"-> clpercent254= %f",clpercent254);       pr(pst);  cenb(pst," ");

/*--------------------------------------------------------------------------!*/
 sprintf(pb,"desert threshold:    %f",di_threshold);       pr(pst);
 sprintf(pb,"band 3 threshold:    %f",b3thresh);           pr(pst);
 sprintf(pb,"NDSI minimum:        %f",thresh_NDSI_min);    pr(pst);
 sprintf(pb,"NDSI maximum:        %f",thresh_NDSI_max);    pr(pst);
 sprintf(pb,"b6 threshold:        %f",b6_thresh);          pr(pst);
 sprintf(pb,"b43 threshold:       %f",b43_thresh);         pr(pst);
 sprintf(pb,"b42 threshold:       %f",b42_thresh);         pr(pst);
 sprintf(pb,"b45 threshold:       %f",b45_thresh);         pr(pst);
 sprintf(pb,"b56 hi threshold:    %f",thresh_b56_hi);      pr(pst);
 sprintf(pb,"b56 lo threshold:    %f",thresh_b56_lo);      pr(pst);
 sprintf(pb,"b3 lower threshold:  %f",thresh_b3_lower);    pr(pst);
 sprintf(pb,"thermal effect high: %f",thermal_effect_high);pr(pst);
 sprintf(pb,"thermal effect low:  %f",thermal_effect_low); pr(pst);
 sprintf(pb,"NDSI snow thresh:    %f",thresh_NDSI_snow);   pr(pst);
 cenb(pst," ");
 cenb(pst,"-"); cenb(pst,"PASS 2 Thresholds"); 
 cenb(pst,"-"); 
 cenb(pst," ");
 sprintf(pb,"snow thresh :        %f", snow_thresh);       pr(pst);
 sprintf(pb,"cloud_mean_temp:     %f", cloud_mean_temp);   pr(pst);
 sprintf(pb,"cloud 254 threshold: %f", cloud_254_thresh);  pr(pst);
/*--------------------------------------------------------------------------!*/
 cenb(pst,"-"); cenb(pst," ");
 cenb(pst,"  REAL LANDSAT DATA  ");
 cenb(pst," -------------------- ");
 cenb(pst," PIXEL FILTER TALLIES "); cenb(pst," ");
 sprintf(pb,"-> total image pixels:       %8d",totpix);   pr(pst);
 sprintf(pb,"-> band 3 reflectance count: %8d",count1);   pr(pst);
 sprintf(pb,"-> band 5/2 NDSI count:      %8d",count2);   pr(pst);
 sprintf(pb,"-> band 6 < 300k count:      %8d",b6delete); pr(pst);
 sprintf(pb,"-> band 5/6 threshold:       %8d",count6);   pr(pst);
 sprintf(pb,"-> band 43 ratio count:      %8d",count3);   pr(pst);
 sprintf(pb,"-> band 42 ratio count:      %8d",count42);  pr(pst);
 sprintf(pb,"-> band 4/5 threshold count: %8d",count4);   pr(pst);
 sprintf(pb,"-> NDSI pixels:              %8d",ndsitab);  pr(pst);
 sprintf(pb,"-> band 3 pixels eliminated: %8d",b3tab);    pr(pst);
 sprintf(pb,"-> hot pixels eliminated:    %8d",hotpixels);pr(pst);
 sprintf(pb,"-> snow pixels eliminated:   %8d",snowtab);  pr(pst);
 sprintf(pb,"-> unambiguous land total:   %8d",landtab);  pr(pst);
/*--------------------------------------------------------------------------!*/
/*- compute the cloud mean values and standard deviation but only if       -!*/
/*- clouds were at half of percent level.			           -!*/
/*---------------------------------------------------------------------------*/
/*- compute cloud mean values and standard deviation but only if clouds    -!*/
/*- were at the half of percent level.                                     -!*/
/*---------------------------------------------------------------------------*/
   desert_index= count42 > 0 ? (float)count4 / (float)count42 : 0.0 ;
   snpercent =  ( (float)snowtab / (float)totpix ) * 100.0;
/*--------------------------------------------------------------------------!*/
/*- If snow is present in the image then ice is too. Use the more          -!*/
/*- conservative cloud mask (254) for developing the training statistics   -!*/
/*- or for determining the cloud count if we're over desert.               -!*/
/*--------------------------------------------------------------------------!*/
 if( snpercent > 1.0 || desert_index <= di_threshold )
   {
   vir_reinit( &b6clouds );
   cloudval255= 0;
   cloudvalues= 0;
   for (iy=0; iy<nls; iy++ )
     {
     if ( !get_line(&clmaskb_File,  (unsigned char*)clmaskb,   iy ) )
       ERROR("getline from clmaskb","csm" );
     if ( !get_line(&b6tempeture_File,(unsigned char*)b6tempeture, iy ) )
       ERROR("getline from b6tempeture","csm" );

     for ( ix=0; ix<nps; ix++)
       {
       if ( clmaskb[ix] == b255 )
	 {
         cloudval255++;
         clmaskb[ix]= b0;
         } 
       if ( clmaskb[ix] == b254 )
	 {
         cloudvalues++;
         clmaskb[ix]=b255;
         virput( &b6clouds, b6tempeture[ix] );
	 }
       }
     if ( !put_line(&clmaskb_File, (char*)clmaskb, iy) )
       ERROR("putline clmaskb file","csm" );

     if ( ( iy==0 || iy ==(nls-1) || iy%100==0 ) && odometer_flag )
       {
       printf("--- desert cloud count pass line %d of %d --- \r",iy,nls-1); 
       fflush(stdout); 
       }

    }
    if ( odometer_flag )printf("\n");
    sprintf(pb,"cloudcount255=%d   cloudcount254=%d",cloudval255,cloudvalues);
    pr(pst);  cenb(pst," "); cenb(pst,"Conservative statistics employed!");
    cenb(pst," "); cenb(pst,"-"); cenb(pst," "); 
    sprintf(pb,"Snow (Ice) present %6.3f Desert Region - index:  %6.3f",
                snpercent, desert_index);
    pr(pst); cenb(pst," "); cenb(pst,"-"); cenb(pst," "); 
    }
  else
    {
    cenb(pst,"-"); cenb(pst," ");
    cenb(pst,"if( snpercent>1.0 || desert_index<=di_threshold )ELSE");
    cenb(pst," "); cenb(pst,"-"); 
    vir_reinit( &b6clouds );
    cloudvalues= 0;
    for (iy=0; iy<nls; iy++ )
      {
      if ( !get_line(&clmaskb_File,  (unsigned char*)clmaskb,   iy ) )
        ERROR("getline from clmaskb","csm" );
       if ( !get_line(&b6tempeture_File,(unsigned char*)b6tempeture, iy ) )
         ERROR("getline from b6tempeture","csm" );
      for ( ix=0; ix<nps; ix++)
        if( clmaskb[ix]==b254 || clmaskb[ix]==b255 )
          {
          clmaskb[ix]=b255;
          cloudvalues++;
          virput( &b6clouds, b6tempeture[ix] );
          }

     if ( !put_line(&clmaskb_File, (char*)clmaskb, iy) )
       ERROR("putline clmaskb file","csm" );

      if ( ( iy==0 || iy ==(nls-1) || iy%100==0 ) && odometer_flag )
        {
        printf("--- pass 2 cloudcount line %d of %d --- \r",iy,nls-1);
        fflush(stdout); 
        }
     }
    if ( odometer_flag )printf("\n");
    cenb(pst,"-"); cenb(pst," ");
    sprintf(pb,"Snow (Ice) absent %f%% moderate statistics employed!",
            snpercent);  pr(pst); cenb(pst," "); cenb(pst,"-");
    }

 virflush( &b6clouds );
 
 if ( cloudvalues >  0)
   {
   moment( &b6clouds, cloudvalues, result );
   b6mean = result[0];
   sprintf(pb,"cloud thermal mean     %f", result[0] ); pr(pst);
   sprintf(pb,"cloud thermal variance %f", result[1] ); pr(pst);
   sprintf(pb,"cloud thermal skew -->>%f", result[2] ); pr(pst);
   sprintf(pb,"cloud kurtosis         %f", result[3] ); pr(pst);
   sprintf(pb,"cloud meanabsdev       %f", result[4] ); pr(pst);
   sprintf(pb,"cloud stdv             %f", result[5] ); pr(pst);
   compute_std( &b6clouds, cloudvalues, mean_std );
   }
  else
    {
    sprintf(pb,
    "b6 mean to be zeroed because there is no band 6 mean, count=%d ",
    cloudvalues);    pr(pst);  cenb(pst,"-");
    b6mean =0.0;
    }
 sprintf(pb,"-> totpix= %d",totpix);  pr(pst);
 sprintf(pb,"-> desert_index= %f clpercent254= %f",desert_index,clpercent254);
 pr(pst);
 sprintf(pb,"  thresholds are ( desert_index > %f )",di_threshold);   pr(pst);
 sprintf(pb,"              && ( clpercent254 > %f )",cloud_254_thresh);
 pr(pst);
 sprintf(pb,"              && ( b6mean < %f       )",cloud_mean_temp);
 pr(pst);
 sprintf(pb,"-> snowtab= %d snpercent %f",snowtab,snpercent); pr(pst);

 if ( clpercent254 > cloud_254_thresh && 
      desert_index > di_threshold     && 
      b6mean       < cloud_mean_temp  )
   {
   cenb(pst,"-"); cenb(pst," ");
   sprintf(pb,"Entered main pass2 section");             pr(pst);
   sprintf(pb,"clpercent254=%f > cloud_254_thresh=%f",
   clpercent254,cloud_254_thresh);                       pr(pst);
   sprintf(pb,"desert_index=%f > di_threshhold=%f",desert_index,di_threshold);
   pr(pst);
   sprintf(pb,"b6mean=%f       < cloud_mean_temp=%f",b6mean,cloud_mean_temp);
   pr(pst);   cenb(pst," "); cenb(pst,"-"); 
/*--------------------------------------------------------------------------!*/
/*- this .95 thermal rule was born form 190 44 (JK)                        -!*/
/*--------------------------------------------------------------------------!*/
/*- compute cloud signature standard deviation and histogram the cloud     -!*/
/*- temperatures.                                                          -!*/
/*--------------------------------------------------------------------------!*/
   index= cloudvalues;
/*--------------------------------------------------------------------------!*/
/*- compute histogram                                                      -!*/
/*--------------------------------------------------------------------------!*/
   if ( index > 0 )
     {
     compute_std( &b6clouds, cloudvalues, mean_std );
     b6sd= mean_std[1];

     histo( &b6clouds, index, b6histogram, HMIN,  HMAX, &b6min, &b6max ) ;
     bins= nint( b6max ) - nint( b6min ) + 1; /* number of bins used*/

     if ( HIST_LOG_FLAG ) 
       {
       cenb(pst,"first histogram");
       sprintf(pb,"first histogram (2) nbins=%d ",bins);     
       for(ihi=nint(b6min)-1; ihi<nint(b6max)+1; ihi++)
         {
         sprintf(pb,"hist[%3.3i]=%10d",ihi,b6histogram[ihi]); 
         pr(pst);
         }
       }
     sprintf(pb,"cloud values:            %d    ",cloudvalues);   pr(pst);
     sprintf(pb,"cloud maximum:           %6.2f ",b6max);         pr(pst);
     sprintf(pb,"cloud minimum:           %6.2f ",b6min);         pr(pst);
     sprintf(pb,"standard deviation:      %6.2f ",mean_std[1]);   pr(pst);
     sprintf(pb,"number of histogram bins: %d   ",bins); pr(pst);cenb(pst," ");
     }
  else
    bins=0;
/*--------------------------------------------------------------------------!*/
/*- examine the thermal skew and set the skew shift factor accordingly.    -!*/
/*- The skew shift factor is limited to between -1 and 1.                  -!*/
/*--------------------------------------------------------------------------!*/
   shiftfact = result[2];
   sprintf(pb," thermal threshold shift factor:%6.2f ",shiftfact); pr(pst);

   if ( result[2] > 1.0  ) shiftfact=  1.0;
   if ( result[2] < -1.0 ) shiftfact= -1.0;

   sprintf(pb," thermal threshold new shift factor: %6.2f",shiftfact); pr(pst);
   shift = b6sd * shiftfact;
   sprintf(pb," threshold skew shift (shiftfact * b6sd): %6.2f",shift);pr(pst);

/*--------------------------------------------------------------------------!*/
/*- calculate 3 cloud temperature thresholds, 1.25%, 2.5%, and 17.5% below -!*/
/*- the max.                                                               -!*/
/*- first determine the total pixel count at each percentage level.        -!*/
/*--------------------------------------------------------------------------!*/
   cloudcut0125 = (int)( (float)cloudvalues * 0.0125 );
   cloudcut025  = (int)( (float)cloudvalues * 0.025 );
   cloudcut175  = (int)( (float)cloudvalues * 0.175 );

 cenb(pst," "); cenb(pst,"-"); cenb(pst,"Threshold Calculations");
 cenb(pst,"-");  cenb(pst," "); 
 sprintf(pb,"cloudvalues= %d totpix=%d",cloudvalues,totpix); pr(pst);
 sprintf(pb,"1.25%% of cloud pixels: %d",cloudcut0125);      pr(pst);
 sprintf(pb,"2.50%% of cloud pixels: %d",cloudcut025);       pr(pst);
 sprintf(pb,"17.50%% of cloud pixels: %d",cloudcut175);      pr(pst);

/*--------------------------------------------------------------------------!*/
/*- now determine the actual histogram bin number represented by each of   -!*/
/*- the 3 %'s.                                                             -!*/
/*--------------------------------------------------------------------------!*/
 cloudsum = 0;
 if ( bins>0 ) 
   {
   for (i=1; i<bins+1; i++)
     { 
     hindex= nint(b6max)+ 1 - i;
     cloudsum = cloudsum + b6histogram[ hindex-1 ] ;
     if ( cloudsum >= cloudcut025 )
       {
       b6thresh1 = b6max - i;
       ibin=i;
       break; /* goto gotit;*/
       }
     } 
/*--------------------------------------------------------------------------!*/
/*-                                 gotit:                                 -!*/
/*--------------------------------------------------------------------------!*/

   sprintf(pb,"2.5%% bin number:%d (%d) cloud pixels excluded:%d",
           hindex,ibin,cloudsum); 
   pr(pst);
   cloudsum = 0;
   for (i=1; i<bins+1; i++)
     {
     hindex= nint(b6max)+ 1 - i;
     cloudsum = cloudsum + b6histogram[ hindex-1 ] ;
     if (cloudsum >= cloudcut175 )
       {
       b6thresh2 = b6max - i;
       ibin=i;
       break; /* goto gotit2 */
       }
     }
/*--------------------------------------------------------------------------!*/
/*-                                gotit2:                                 -!*/
/*--------------------------------------------------------------------------!*/
   sprintf(pb,"17.5%% bin number: %d (%d)cloud pixels excluded:%d",
           hindex,ibin,cloudsum);                                    pr(pst);
   cloudsum = 0;
   for (i=1; i<bins+1; i++)
     {
     hindex= nint(b6max)+ 1 - i;
     cloudsum = cloudsum + b6histogram[ hindex-1 ];
     if (cloudsum >= cloudcut0125)
       {
       b6thresh3 = b6max - i;
       ibin=i;
       break; /* goto gotit3 */
       }
     }
   b6thresh_max= b6thresh3;
/*--------------------------------------------------------------------------!*/
/*-                                gotit3:                                 -!*/
/*--------------------------------------------------------------------------!*/
/*--------------------------------------------------------------------------!*/
/*-                       end b6thresh computations                        -!*/
/*--------------------------------------------------------------------------!*/
   sprintf(pb,"1.25%% bin number: %d (%d) cloud pixels excluded: %d",
           hindex,ibin,cloudsum);                pr(pst);  cenb(pst," "); 
   sprintf(pb,"band 6 threshold 97.5  %%: %f",b6thresh1);  pr(pst);
   sprintf(pb,"band 6 threshold 82.5  %%: %f",b6thresh2);  pr(pst);
   sprintf(pb,"band 6 threshold 98.75 %%: %f",b6thresh3);  pr(pst);
   pr(pst);   cenb(pst," "); cenb(pst,"-"); cenb(pst,"second pass results");
   cenb(pst,"-");  cenb(pst," "); 
   }
/*--------------------------------------------------------------------------!*/
/*- the b6skew threshold was borne from the 169/28 image. The old          -!*/
/*- algorithmhad the mean + 2SD as the threshold. This yielded an enormous -!*/
/*- threshold which  classified the entire image as cloud (large neg skew) -!*/
/*--------------------------------------------------------------------------!*/
/*--------------------------------------------------------------------------!*/
/*- exploit the cloud distribution skew characteristics but only if the    -!*/
/*- kew is positive.  change the two working thresholds (2.5%, 17.5%) by   -!*/
/*- the computed skew shift.                                               -!*/
/*--------------------------------------------------------------------------!*/

 if( shiftfact > 0.0 ) 
   {
   printf("*** Threshold Skew Enabled  ***\n");
   b6thresh1_save = b6thresh1;
   b6thresh2_save = b6thresh2;
   b6thresh1 = b6thresh1 + shift;
   b6thresh2 = b6thresh2 + shift;

/*--------------------------------------------------------------------------!*/
/*- does the skew adjusted 97.5% threshold exceed the maximum (98.75%)?    -!*/
/*- If so, set it to the maximum.  Determine the acceptable skew adjustment-!*/
/*- and add it to the original 82.5% threshold to get a new threshold.     -!*/
/*--------------------------------------------------------------------------!*/
   if ( b6thresh1 > b6thresh_max ) 
     {
     b6thresh1 = b6thresh_max;
     b6thresh2 = (b6thresh1 - b6thresh1_save) + b6thresh2_save;
     sprintf(pb,"Threshold exceeds maximum allowed (1.25%%): %6.2f ",
           b6thresh_max);      pr(pst);
     sprintf(pb," 1.25%% threshold shifted: %6.2f",b6thresh1);   pr(pst);
     sprintf(pb,"17.50%% threshold shifted: %6.2f",b6thresh2);   pr(pst);
     }
   else 
     {
     prb(pst,"shifted threshold employed.. maximum not exceeded ");
     sprintf(pb,"1.25%% threshold shifted: %6.2f ",b6thresh1);      pr(pst);
     sprintf(pb,"17.50%% threshold shifted: %6.2f ",b6thresh2);     pr(pst);
     }
   }
  else
   {
   prb(pst,"Threshold Skew Disabled (Negative)");
   sprintf(pb," 2.5%% thermal threshold (shift): %6.2f ",b6thresh1); pr(pst);
   sprintf(pb,"17.5%% thermal threshold (shift): %6.2f ",b6thresh2); pr(pst);
   }
 cenb(pst," "); cenb(pst,"-"); 
/*--------------------------------------------------------------------------!*/
/*-                               begin loop                               -!*/
/*--------------------------------------------------------------------------!*/
  if ( !virinit( &b6clouds200, "b6clouds200.temp.file",bfact ) )
   ERROR("Error opening virtual buffer b6clouds2","csm" );

  count9=  0;
  count10= 0;
  vir_reinit( &b6clouds );
  cloudvalues=  0;
  cloudvalues2= 0;

  for (iy=0; iy<nls ; iy++ )
    {
    if ( !get_line( &b6tempeture_File, (unsigned char*)b6tempeture, iy ) )
        ERROR("getline from b6tempeture","csm" );
    if ( !get_line( &clmaskb_File,   (unsigned char*)clmaskb,   iy ) )
        ERROR("getline from clmaskb","csm" );

    for ( ix=0; ix<nps; ix++)
      {
      if ( clmaskb[ix] == b0 ) 
        {
        if ( b6tempeture[ix] <= b6thresh1 ) 
          {
          clmaskb[ix]= b125;
          count9++;
          if ( b6tempeture[ix] <= b6thresh2 ) 
            {
            clmaskb[ix]= b200;
            count10++;
            virput( &b6clouds200, b6tempeture[ix] );
            cloudvalues2++;
            }
          virput( &b6clouds, b6tempeture[ix] );
          cloudvalues++;                  /* b6clouds histogram */
          }
        }
      }
    if ( !put_line( &clmaskb_File, (char*)clmaskb, iy ) )
      ERROR("putline clmaskb file","csm" );
    if ( ( iy==0 || iy ==(nls-1) || iy%100==0 ) && odometer_flag )
      {
      printf("--- main pass2 loop, line %d of %d --- \r",iy,nls-1); 
      fflush(stdout); 
      }
    }
  if ( odometer_flag )printf("\n");
  virflush( &b6clouds );
  virflush( &b6clouds200 );
/*--------------------------------------------------------------------------!*/
/*-                            second histogram                            -!*/
/*--------------------------------------------------------------------------!*/
   index= cloudvalues;
   if ( index > 0 )
     {
     histo( &b6clouds, index, b6histogram, HMIN,  HMAX, &b6min, &b6max1 ) ;
     bins= nint( b6max1 ) - nint( b6min ) + 1; /* number of bins used*/

     if ( HIST_LOG_FLAG ) 
       {
       cenb(pst," ");
       sprintf(pb,"second histogram (2) nbins=%d ",bins);        pr(pst);
       sprintf(pb,"histogram min=%f max=%f ",b6min, b6max1);     pr(pst);
       cenb(pst," ");
       for(ihi=nint(b6min)-1; ihi<nint(b6max1)+1; ihi++)
         {
         sprintf(pb,"hist[%3.3i]=%10d",ihi,b6histogram[ihi]);    pr(pst);
         }
       }
     }
 cenb(pst,"-"); 
 cenb(pst,"2.5% Threshold Statistics");                cenb(pst,"-");
 sprintf(pb," cloud minimum: %6.2f",b6min);          pr(pst);
 sprintf(pb," cloud maximum: %6.2f",b6max1);         pr(pst);

 if ( cloudvalues >  0)
   {
   moment( &b6clouds, cloudvalues, result ); 
   b6mean1 = result[0];
   sprintf(pb," cloud mean   : %6.2f",b6mean1); pr(pst);
   cenb(pst,"-"); cenb(pst," "); 
   }
 else
   {
   b6mean1 = b6mean;
   }
/*--------------------------------------------------------------------------!*/
/*-                             third histogram                            -!*/
/*--------------------------------------------------------------------------!*/
   index= cloudvalues2;
   if ( index > 0 )
     {
     histo( &b6clouds200, index, b6histogram, HMIN,  HMAX, &b6min, &b6max2 ) ;
     bins= nint( b6max2 ) - nint( b6min ) + 1; /* number of bins used*/
     if ( HIST_LOG_FLAG ) 
       {
       sprintf(pb,"third histogram (3) nbins=%d",bins);           pr(pst);
       sprintf(pb,"histogram min=%f max=%f",b6min,b6max2);        pr(pst);
       for(ihi=nint(b6min)-1; ihi<nint(b6max2)+1; ihi++)
         {
         sprintf(pb,"hist[%3.3i]=%10d",ihi,b6histogram[ihi]);     pr(pst);
         }
       }
     }
/*--------------------------------------------------------------------------!*/
/*--------------------------------------------------------------------------!*/
   if ( cloudvalues2 > 1 ) 
     {
     cenb(pst," ");   cenb(pst,"-");
     cenb(pst,"17.5% Threshold Statistics");   cenb(pst,"-");   cenb(pst," "); 
     sprintf(pb," cloud minimum: %6.2f ",b6min);   pr(pst);
     sprintf(pb," cloud maximum: %6.2f ",b6max2);  pr(pst);     
     moment( &b6clouds200, cloudvalues2, result );
     b6mean2 = result[0];
     sprintf(pb," cloud mean   : %6.2f ",b6mean2);   pr(pst); cenb(pst," ");
     }
   else
     {
     b6mean2 = b6mean;
     }

/*--------------------------------------------------------------------------!*/
   sprintf(pb,"pixel count through  2.5%% thermal threshold: %d ",count9); 
   pr(pst); 
   sprintf(pb,"pixel count through 17.5%% thermal threshold: %d ",count10);
   pr(pst);  cenb(pst," "); cenb(pst,"-"); 
/*--------------------------------------------------------------------------!*/
/* determine cloud percentage for 5% threshold and thermal effect          -!*/
/*--------------------------------------------------------------------------!*/
   thermeffecta =  (float)count9 / (float)totpix * 100.0;
   thermeffectb = (float)count10 / (float)totpix * 100.0;

   if(snpercent > 1.0 || desert_index <= di_threshold )
     clpercentpass1 = clpercent254;
   else
     clpercentpass1 = clpercentorig;

   clpercenta = clpercentpass1 + thermeffecta;
   clpercentb = clpercentpass1 + thermeffectb;

   sprintf(pb,"baselined pass 1 cloud percentage: %6.2f %%",clpercentpass1);
   pr(pst); 
   sprintf(pb,"net thermal effect for 2.5%% filter %6.2f%%",thermeffecta);
   pr(pst); 
   sprintf(pb,"net thermal effect for 17.5%% filter %6.2f%%",thermeffectb);  
   sprintf(pb,"cloud %% after 2.5%% thermal filter:   %6.2f%%",clpercenta);
   pr(pst); 
   sprintf(pb,"cloud %% after 17.5%% thermal filter: %6.2f%%",clpercentb );
   pr(pst); 

/*--------------------------------------------------------------------------!*/
/*- if a radical thermal effect took place or if snow is present step down -!*/
/*- to a more conservative filter. change to 40 from 30 due to 164 20.     -!*/
/*--------------------------------------------------------------------------!*/

   b6max_thresh1_diff = b6max - b6thresh1;
   if ( thermeffecta >  thermal_effect_high || snpercent > snow_thresh  ||
        b6mean1      >= cloud_mean_temp     || b6max_thresh1_diff <= 2.0 ) 
     {
     if(thermeffectb > thermal_effect_low || b6mean2 > cloud_mean_temp )
       {
        filter_used = 3; /*  clmask(index) = 0 */
       }
     else
       {
       filter_used = 2; /*  clmask(index) = 255 */
       }
     }
   else
     {
     filter_used = 1;
     }
/*--------------------------------------------------------------------------!*/
/*                                 40% to JK                               -!*/
/*--------------------------------------------------------------------------!*/
   cenb(pst,"-");  sprintf(pb,"filter used=%d",filter_used); cen(pst);
   cenb(pst,"-"); cenb(pst," ");
/*--------------------------------------------------------------------------!*/
/*-              index = where(clmask eq 200, cloudvalues)                 -!*/
/*--------------------------------------------------------------------------!*/
   cloudvalues=  0;
   for (iy=0; iy<nls ; iy++ )
     {
     if ( !get_line( &clmaskb_File,   (unsigned char*)clmaskb,   iy ) )
       ERROR("getline from clmaskb","csm" );

     for ( ix=0; ix<nps; ix++)
       {
       if ( clmaskb[ix] == b200 || clmaskb[ix]  == b125 ) 
         {
         if ( filter_used == 1 )
           {
           clmaskb[ix]= b255;
           }
         else if ( filter_used == 2 )
           {
           if ( clmaskb[ix] == b200 )clmaskb[ix]= b255;
           }
         else if ( filter_used == 3 )
           {
           if ( clmaskb[ix] == b200 )clmaskb[ix]= b0;
           }
         }
/*--------------------------------------------------------------------------!*/
/*- now zero out everything but clouds and compute cloud percentage before -!*/
/*- filtering                                                              -!*/
/*--------------------------------------------------------------------------!*/
       if ( clmaskb[ix] == b255 )
         cloudvalues++;
       else
         clmaskb[ix] == b0;

       }
     if ( !put_line( &clmaskb_File, (char*)clmaskb, iy ) )
       ERROR("putline clmaskb file","csm" );

     if ( ( iy==0 || iy ==(nls-1) || iy%100==0 ) && odometer_flag )
       {
       printf("--- pass2 filter, line %d of %d --- \r",iy,nls-1); 
       fflush(stdout); 
       }
     }
   if ( odometer_flag )printf("\n");
   clpercent =  (float)cloudvalues / (float)totpix * 100.0;

/*--------------------------------------------------------------------------!*/
/*- determine its effect... if still radical forget the filter and goback  -!*/
/*- to the original cloud determination result.                            -!*/
/*--------------------------------------------------------------------------!*/
  
   if ( filter_used == 1 )
     {
     sprintf(pb,"2.5 percent filter used ");  pr(pst);
     }
   if ( filter_used == 2 )
     {
     sprintf(pb,"17.5 percent filter used due to: ");
     if ( snpercent > 1.0 )
       {
       sprintf(pb,"snow percentage %6.2f %% is greater than 1%%",snpercent);
       pr(pst);
       }
    if( thermeffecta > thermal_effect_high ) 
      {
      sprintf(pb," 2.5%% thermal effect was greater than 40%%");     pr(pst);
      }
   if( b6max_thresh1_diff <= 2.0 )
      {
      sprintf(pb," 2.5%% threshold is within 2 bins of the band 6 maximum");
      pr(pst);
      }
    }
   if ( filter_used == 3 )
     {
     sprintf(pb,"questionable results, thermal pass discarded  ");   pr(pst);
     }
   sprintf(pb,"cloud %% before filtering: %6.2f",clpercent);         pr(pst);
   sprintf(pb,"before filtering count: %d", cloudvalues);            pr(pst);
   }
 else
   {  
/*--------------------------------------------------------------------------!*/
/*- endif else begin  ;                                                    -!*/
/*- only if clouds were found endelse (i.e > 4%) and desert index          -!*/
/*- use the conservative cloud count because desert features are present.  -!*/
/*--------------------------------------------------------------------------!*/
   if ( b6mean < cloud_mean_temp )
     {
     clpercent = clpercent254;
     sprintf(pb,"thermal analysis bypassed with a valid b6 mean!");  pr(pst);
     sprintf(pb,"band 6 mean: %f (<threshold=%f) cloud percent=%f",
     b6mean,cloud_mean_temp,clpercent);                              pr(pst);
     }
   else
     {
     cenb(pst," ");
     cenb(pst,"thermal analysis bypassed!");   cenb(pst," ");   cenb(pst,"-"); 
     sprintf(pb,"invalid results, band 6 mean: %f zero out all clouds!",
             b6mean);                                                pr(pst); 

     clpercent = 0.0;
/*--------------------------------------------------------------------------!*/
/*-                            clmask(index) = 0                           -!*/
/*--------------------------------------------------------------------------!*/
     for (iy=0; iy<nls; iy++ )
       {
       if ( !get_line(&clmaskb_File,  (unsigned char*)clmaskb,   iy ) )
         ERROR("getline from clmaskb","csm" );
       for ( ix=0; ix<nps; ix++)
         {
         if ( clmaskb[ix] == b255 )
           {
           clmaskb[ix]= b0;
           }
         }
       if ( !put_line(&clmaskb_File, (char*)clmaskb, iy) )
         ERROR("putline clmaskb file","csm" );
       if ( ( iy==0 || iy ==(nls-1) || iy%100==0 ) && odometer_flag )
         {
         printf("--- conservative clear clouds pass, line %d of %d --- \r",
         iy,nls-1); 
         fflush(stdout); 
         }
      }
     if ( odometer_flag )printf("\n");
     }/*  endelse */

   cenb(pst," "); 
   cenb(pst,"thermal analysis bypassed!");        cenb(pst," "); 
   sprintf(pb,"cloud cover percentage: %f",clpercent); pr(pst); 
   sprintf(pb,"desert index: %f",desert_index);        pr(pst); 
   }
 cenb(pst," "); cenb(pst,"-"); 
/*--------------------------------------------------------------------------!*/
/*- endelse    ; only if clouds were found endif                           -!*/
/*--------------------------------------------------------------------------!*/
/*- set up for spatial filter                                              -!*/
/*- spatial filter modified to work on both snow and cloud masks           -!*/
/*- clmaskb image set to 1 for cloud 16 for snow                           -!*/
/*--------------------------------------------------------------------------!*/
 cloudvalues=0;
 for (iy=0; iy<nls; iy++ )
   {
   if ( !get_line(&clmaskb_File,  (unsigned char*)clmaskb,   iy ) )
     ERROR("getline from clmaskb","csm" );
   if ( !get_line( &clmask_File,   (unsigned char*)clmask,   iy ) )
      ERROR("getline from clmask","csm" );

   for ( ix=0; ix<nps; ix++)
     {
     clmaskb2[ix]= ( clmaskb[ix] == b255 || clmaskb[ix]== b254 ) ? b1 : b0;
     if ( clmask[ix] == CLSTAT_S ) clmaskb2[ix]+= 16;
     if ( clmaskb2[ix]==b1 ) cloudvalues++;
     }

   if ( !put_line(&clmaskb2_File, (char*)clmaskb2, iy) )
     ERROR("putline clmaskb2 file","csm" );

   if ( ( iy==0 || iy ==(nls-1) || iy%100==0 ) && odometer_flag )
     {
     printf("--- pre spatial filter pass line %d of %d --- \r",iy,nls-1); 
     fflush(stdout); 
     }
   }
 if ( odometer_flag )printf("\n");
 sprintf(pb,"cloudvalues prior to final filter=%d ",cloudvalues);   pr(pst); 
/*--------------------------------------------------------------------------!*/
/*-                             spatial filter                             -!*/
/*-                  fill in cloud holes if clouds exist                   -!*/
/*--------------------------------------------------------------------------!*/
 memset(clmaskb2, 0, nps*4*sizeof(unsigned char));
 filt_new= 0;
 new_snow= 0;
 memset(cloud,0,sizeof(int)*2*2);
 cloudvalues= 0;
 if ( !get_line( &clmaskb2_File, (unsigned char*)&clmaskb2[0*nps], 0 ) )
  ERROR("getline from clmaskb2","csm" );

 for (iy=0; iy<nls; iy++ )
   {
   if ( !get_line( &clmaskb_File,  (unsigned char*)clmaskb,  iy ) )
    ERROR("getline from clmaskb","csm" );
   if ( !get_line( &clmask_File,   (unsigned char*)clmask,   iy ) )
    ERROR("getline from clmask","csm" );

   p= (iy+1)%3;
   if ( (iy+1)<nls )
     if ( !get_line( &clmaskb2_File, (unsigned char*)&clmaskb2[p*nps], iy+1 ) )
      ERROR("getline from clmaskb2","csm" );

   clmaskb2p= &clmaskb2[p*nps];
   clmaskb2m= &clmaskb2[((p+1)%3)*nps];
   clmaskb2c= &clmaskb2[((p+2)%3)*nps];

   if ( iy==0       )clmaskb2m= clmaskb2c;
   if ( iy==(nls-1) )clmaskb2p= clmaskb2c;

   iyp= iy + 1;
   iym= iy - 1;
   if ( iym < 0       )iym= iym + 1;
   if ( iyp > (nls-1) )iyp= iyp - 1;

   for ( ix=0; ix<nps; ix++)
     {
     if ( clmaskb[ix] < b254 )
       {
       ixp= ix + 1;
       ixm= ix - 1;
       if ( ixm < 0       )ixm= ixm + 1;
       if ( ixp > (nps-1) )ixp= ixp - 1;

       isum_cloud= (int)( clmaskb2m[  ixm  ] & 0x0f ) +
                   (int)( clmaskb2m[  ix   ] & 0x0f ) +
                   (int)( clmaskb2m[  ixp  ] & 0x0f ) +
                   (int)( clmaskb2c[  ixm  ] & 0x0f ) +
                   (int)( clmaskb2c[  ixp  ] & 0x0f ) +
                   (int)( clmaskb2p[  ixm  ] & 0x0f ) +
                   (int)( clmaskb2p[  ix   ] & 0x0f ) +
                   (int)( clmaskb2p[  ixp  ] & 0x0f );

       isum_snow=  (int)( clmaskb2m[  ixm  ] >> 4 ) +
                   (int)( clmaskb2m[  ix   ] >> 4 ) +
                   (int)( clmaskb2m[  ixp  ] >> 4 ) +
                   (int)( clmaskb2c[  ixm  ] >> 4 ) +
                   (int)( clmaskb2c[  ixp  ] >> 4 ) +
                   (int)( clmaskb2p[  ixm  ] >> 4 ) +
                   (int)( clmaskb2p[  ix   ] >> 4 ) +
                   (int)( clmaskb2p[  ixp  ] >> 4 );

       if ( isum_cloud >= 5 ) /* more than 5 cloudy neighbors= cloud */
         {
         clmaskb[ix]= b255 ;
         filt_new++;
         }
                              /* more than 5 snow/icd neighbors= snow */
       if ( isum_snow >= 5 && clmask[ix] != CLSTAT_S )  
         {
         clmask[ix]=  CLSTAT_S;
         new_snow++;
         }
       }
     if ( clmaskb[ix] == b255 || clmaskb[ix] == b254   )
       {
       clmask[ix]= CLSTAT_C;
       cloudvalues++;
       cloud[(iy/yhalf)*2+(ix/xhalf)]++;
       }
     /* eliminate the snow mask */
     if ( clmask[ix] == CLSTAT_S )clmask[ix]= CLSTAT_L;
     }
   if ( param->sieve_thresh>0 || param->apply_kernel>0 )
     memcpy(&imask_img[nps*iy],clmask,nps);

   if ( !put_line(&clmaskb_File, (char*)clmaskb, iy ) )
     ERROR("putline clmaskb file","csm" );

   if ( !put_line(&clmask_File, (char*)clmask, iy ) )
     ERROR("putline clmask file","csm" );

   if ( ( iy==0 || iy ==(nls-1) || iy%100==0 ) && odometer_flag )
      {
      printf("--- 5 neighbor filter, line %d of %d --- \r",iy,nls-1);
      fflush(stdout); 
      }
   }
 if ( odometer_flag )printf("\n");
 
/*--------------------------------------------------------------------------!*/
/*-                                sieve filter                            -!*/
/*--------------------------------------------------------------------------!*/
 if ( param->sieve_thresh>0 )
   {
  sieve( nps        , nls    ,
         imask_img           ,
         sive_flag           ,
         sive_cflag          ,
         param->sieve_thresh , 
         odometer_flag
       );
   }
/*--------------------------------------------------------------------------!*/
/*-                                write sieved                            -!*/
/*--------------------------------------------------------------------------!*/
   if ( WRITE_SIEVE )
     {
     for (iy=0; iy<nls; iy++ )
       if ( !put_line(&sive_File, (char*)&imask_img[iy*nps], iy ) )
         ERROR("putline clmask file","csm" );
     img_close(&sive_File);
     }
/*--------------------------------------------------------------------------!*/
/*-                             spatial filter (2)                         -!*/
/*-                   apply kernel to expand cloud mask                    -!*/
/*--------------------------------------------------------------------------!*/
 sprintf(pb,"filtered mask, ksize= %d, apply=%d sieve_thresh=%d",
   param->ksize,     param->apply_kernel, param->sieve_thresh); pr(pst); 

 if ( param->apply_kernel>0 )
   {
   ksize= param->ksize;
   for (ik=0; ik<param->apply_kernel; ik++ )
     {
     for (iy=0; iy<nls; iy++ )
       {
       for ( ix=0; ix<nps; ix++)
         {
         oval= imask_img[iy*nps+ix];
         for (yy=0; yy<ksize; yy++)
           {
           yyy= iy-(ksize/2)+yy;
           for (xx=0; xx<ksize; xx++)
             {
             xxx= ix-(ksize/2)+xx;
             if ( xxx>=0 && xxx<nps && yyy>=0 && yyy<nls )
               {
               ival= imask_img[yyy*nps+xxx];
               if ( ival == CLSTAT_C )oval= CLSTAT_C;
               }
             }
           }
         omask_img[iy*nps+ix]= oval;
         }
       }
     if ( ik < (param->apply_kernel-1) )
       {
       imask_img= & bmask_img[ ((bufptr++)%2)*nps*nls ];
       omask_img= & bmask_img[ ((bufptr  )%2)*nps*nls ]; 
       }
    else
       {
       sprintf(pb,"last ik=%d",ik); pr(pst); 
       }
     }
/*--------------------------------------------------------------------------!*/
/*-                               write mask                               -!*/
/*--------------------------------------------------------------------------!*/
     for (iy=0; iy<nls; iy++ )
       {
       if ( !put_line(&clmask_File, (char*)&omask_img[iy*nps], iy ) )
         ERROR("putline clmask file","csm" );

       if ( ( iy==0 || iy ==(nls-1) || iy%100==0 ) && odometer_flag )
          {
          printf("--- expand with kernel, line %d of %d --- \r",iy,nls-1);
          fflush(stdout); 
          }
       }
     if ( odometer_flag )printf("\n");
   }
/*--------------------------------------------------------------------------!*/
/*-                            write final mask                            -!*/
/*--------------------------------------------------------------------------!*/
 for (iy=0; iy<nls; iy++ )
   {
   if ( !get_line( &clmask_File, (unsigned char*)clmask, iy ) )
    ERROR("getline from clmask","csm" );
    if (!PutOutputLine(output, CLMASK, iy, clmask))
      ERROR("writing output data for a line (CLMASK)", "csm");
   if ( ( iy==0 || iy ==(nls-1) || iy%100==0 ) && odometer_flag )
      {
      printf("--- write final mask, line %d of %d --- \r",iy,nls-1);
      fflush(stdout); 
      }
   }
 if ( odometer_flag )printf("\n");


/*--------------------------------------------------------------------------!*/
/*-     compute final quadrant cloud cover percent for log file output     -!*/
/*--------------------------------------------------------------------------!*/
 clpercent =  ( (float)cloudvalues / (float)totpix ) * 100.0;
 c_cover_tot= 0.0;
 for (iy=0; iy<2; iy++)
   for (ix=0; ix<2; ix++)
     {
     c_cover[iy*2+ix]= (float)cloud[iy*2+ix] * 100.0/(float)( xhalf*yhalf );
     c_cover_tot+= c_cover[iy*2+ix]/4.0;
     }
 sprintf(pb,"filtered total: %d",filt_new);           pr(pst); 
 sprintf(pb,"filtered snow:  %d",new_snow);           pr(pst); 
 sprintf(pb,"final cloud count: (%d,%f,%f)",
              cloudvalues,clpercent,c_cover_tot);     pr(pst); 
 sprintf(pb,"cloud cover scores=(%f,%f,%f,%f)",   
       c_cover[0],c_cover[1],c_cover[2],c_cover[3]);  pr(pst); cenb(pst," "); 
 cenb(pst,"processing completed");     cenb(pst," "); prb(pst,"*"); 
/*--------------------------------------------------------------------------!*/
/*-                         close temporary images                         -!*/
/*--------------------------------------------------------------------------!*/
 img_close_rm( &b6tempeture_File );
 img_close_rm( &clmaskb_File   );
 img_close_rm( &clmaskb2_File  );
 img_close_rm( &clmask_File    );
/*--------------------------------------------------------------------------!*/
/*-                         close virutal buffers                          -!*/
/*--------------------------------------------------------------------------!*/
 virclose( &b6clouds    );
/* virclose( &b6clouds200 ); */
 virclose( &b6clouds254 );
 if ( LOG_FLAG )fclose( mout );
/*--------------------------------------------------------------------------!*/
/*-                           free temp buffers                            -!*/
/*--------------------------------------------------------------------------!*/
 FREE_IT:
 free( b6tempeture );
 free( clmask    );
 free( clmaski2  );
 free( clmaskb   );
 free( clmaskb2  );
 free( ib2buf    );
 free( ib3buf    );
 free( ib4buf    );
 free( ib5buf    );
 free( ib6buf    );
 free(line_in_buf);
 return true;
}
/*--------------------------------------------------------------------------!*/
/*--------------------------------------------------------------------------!*/
/*-                                                                        -!*/
/*-                          histogram function                            -!*/
/*-                                                                        -!*/
/*--------------------------------------------------------------------------!*/
/*--------------------------------------------------------------------------!*/
void histo( vbuf_t* vb, int num, int* outhist, const int HMIN, const int HMAX, 
            float* omin, float* omax ) 
{
int i,j,num_above,num_below;
/*--------------------------------------------------------------------------!*/
/*                                initalize                                -!*/
/*--------------------------------------------------------------------------!*/
      num_above= 0;
      num_below= 0;

      *omin= 1024;
      *omax=-1024;

      for (i=0; i<1024; i++)outhist[i]= 0;
/*--------------------------------------------------------------------------!*/
/*                             fill histogram                              -!*/
/*--------------------------------------------------------------------------!*/
      for (i=0; i<num; i++)
	{
        float data= virget( vb,i );
        j= (int)data ;
        if ( j == 0 )
            j= HMIN;
         else if ( j < HMIN )
	   {
           sprintf(pb,"bad histogram value=%d ",j);  pr(pst); 
           j= HMIN;
           num_below= num_below + 1;
           }
         else if ( j > HMAX )
           {
           sprintf(pb,"bad histogram value=%d ",j);   pr(pst); 
           j= HMAX;
           num_above= num_above + 1;
           }
         outhist[ j-1 ]= outhist[ j-1 ] + 1;
         if ( data < *omin )*omin= data;
         if ( data > *omax )*omax= data;
	}
/*--------------------------------------------------------------------------!*/
/*                    check for histogram extent error                     -!*/
/*--------------------------------------------------------------------------!*/
      if ( num_below > 0 )
        {
        sprintf(pb,"warning histogram min=%d number below=%d "
             ,HMIN,num_below);    pr(pst); 
        }
      if ( num_above > 0 )
        {
        sprintf(pb,"warning histogram max=%d number above=%d "
            ,num_above,HMAX);    pr(pst); 
        }
      if ( num_above > 0 || num_below > 0 )
	{
        cenb(pst," histogram limit error"); 
	}
}
/*--------------------------------------------------------------------------!*/
/*--------------------------------------------------------------------------!*/
/*-                                                                        -!*/
/*-                            moment function                             -!*/
/*-                                                                        -!*/
/*--------------------------------------------------------------------------!*/
/*--------------------------------------------------------------------------!*/
void moment( vbuf_t* vb, const int num, float* output ) 
{
int  i;
double  mean
      , variance
      , skewness
      , kurtosis
      , meanabsdev
      , stdv
      , sum    
      , sumsq  
      , sumdsq 
      , sumadif
      , sumskew
      , sumkert;
/*--------------------------------------------------------------------------!*/
/*                                                                         -!*/
/*  the purpose of this routine is to compute some statistics              -!*/
/*  the output is in the array output                                      -!*/
/*  the contents are as follows:                                           -!*/
/*  output(1)= mean                                                        -!*/
/*  output(2)= variance                                                    -!*/
/*  output(3)= skewness                                                    -!*/
/*  output(4)= kurtosis                                                    -!*/
/*  output(5)= mean absolute deviation                                     -!*/
/*  output(5)= standard deviation                                          -!*/
/*  warning elimination                                                    -!*/
/*                                                                         -!*/
/*--------------------------------------------------------------------------!*/
 sum=     0.0;
 sumsq=   0.0;
 sumdsq=  0.0;
 sumadif= 0.0;
 sumskew= 0.0;
 sumkert= 0.0;


 for (i=0; i<num; i++)
   {
   double data= (double)virget( vb,i );
   sum= sum + data;
   sumsq= sumsq + ( data * data );
   }
 mean= sum/ (double)num;

 for (i=0; i<num; i++)
   {
   double data= (double)virget( vb,i );
   sumdsq= sumdsq + pow( data - mean,2);
   sumadif= sumadif + fabs( data - mean );
   }

 variance= sumdsq / (float)max(num-1,1);
 stdv= max( sqrt( variance ) , 0.01 );
 meanabsdev= sumadif/(float)num;

 for (i=0; i<num; i++)
   {
   double data= (double)virget( vb,i );
   sumskew= sumskew + pow( ( data - mean )/stdv ,3);
   sumkert= sumkert + pow( ( data - mean )/stdv ,4);
   }

 skewness= sumskew / (double)num;
 kurtosis= sumkert / (double)num;

 output[1-1]= (double)mean;
 output[2-1]= (double)variance;
 output[3-1]= (double)skewness;
 output[4-1]= (double)kurtosis;
 output[5-1]= (double)meanabsdev;
 output[6-1]= (double)stdv;
}
