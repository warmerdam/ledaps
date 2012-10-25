#ifndef CSM_H
#define CSM_H

#include "lndcsm.h"
#include "bool.h"
#include "lut.h"
#include "param.h"
#include "input.h"
#include "output.h"
#include "tiff.h"     /* for tiff package*/
#include "virbuf.h"   /* for virtual buffer*/

bool CloudMask(  
                Input_t *input
             ,  Input_t *input_th
             ,  Output_t *output
             ,  Lut_t *lut
             ,  Param_t *param
             ,  int odometer_flag
                );
/*****************************************************************************/
void moment( vbuf_t* vb,       const int num, float* output );
void histo( vbuf_t* vb, int num, int* outhist, const int hmin, const int hmax, 
            float* omin, float* omax ) ;
void compute_std(vbuf_t* vb, const int num, float* mean_std);
/*****************************************************************************/
#endif
/*****************************************************************************/
/*****************************************************************************/
/***                                                                       ***/
/*** the following IDL source code is the program for the original         ***/
/*** ACCA algorithm on which this C program was based.                     ***/
/***                                                                       ***/
/***                                                                       ***/
/*** I made two modifications to Richards IDL programs                     ***/
/***                                                                       ***/
/*** 1)   This change keeps floating point precision of temperature        ***/
/***    ;b6_temp(0:samples-1,l) = double(b6 / 100) + 273                   ***/
/***     b6_temp(0:samples-1,l) = (double(b6) / 100.0) + 273.0             ***/
/***                                                                       ***/
/*** 2) This change to filtmaj keeps the kernel symetric by not modifying  ***/
/***    the input image as the kernel is applied.                          ***/
/***                                                                       ***/
/***      image(x-1,y-1) = fill_val                                        ***/
/***    ; imf(x,y) = fill_val                                              ***/
/***                                                                       ***/
/*** ;return, imf( radius:(Lx+radius), radius:(Ly+radius) ) ; fixed(jk)    ***/
/*** return, image                                                         ***/
/***                                                                       ***/
/*****************************************************************************/
				
/*
;-----------------------------------------------------------------------------------
;
;            		Automatic Cloud Cover Assessment Algorithm
;		         for application to Landsat 5 TM 1G imagery
;                           Algorithm ID:  ACCA_10May03_LPSO.ALG
;-----------------------------------------------------------------------------------
pro accal5_1g

true = 1
false = 0

;-----------------------------------------------------------------------------------
openr,2,'b2.raw'
openr,3,'b3.raw'
openr,4,'b4.raw'
openr,5,'b5.raw'
openr,6,'b6.raw'
openw,7,'p1.mask'
openw,8,'p2.mask'	
openw,9,'snow.mask'
openw,10,'acca_results'
;-----------------------------------------------------------------------------------		

samples = 9349L
lines = 8524L
path = 38
row = 22
sun_el = 33.
deg2rad =  atan(1.0) / 45.
pi = double(3.1415926536)

;-----------------------------------------------------------------------------------
; Initialize counters.
;-----------------------------------------------------------------------------------

count1 = 0L
count2 = 0L
count3 = 0L
count4 = 0L
count5 = 0L
count6 = 0L
count7 = 0L
count9 = 0L
count10= 0L
count42 = 0L
cloudclount = 0L
cloudvalues = 0L
snowtab = 0L
hotpixels = 0L
landtab = 0L
B3tab = 0L
NDSItab = 0L
index = 0L
index2 = 0L
b6delete = 0L
snpercent = 0

b2 = intarr(samples)
b3 = intarr(samples)
b4 = intarr(samples)
b5 = intarr(samples)
b6 = intarr(samples)

b6_temp = fltarr(samples,lines)

clmask = bytarr(samples,lines)
snowmask = bytarr(samples,lines)
;desert_mask = bytarr(samples,lines)
printf,10,'************************************************************************'
printf,10,''

print,'samples: ',samples,' lines: ',lines
printf,10, format = '("* Path: ",I4," Row: ",I4,49x,"*")',path,row
printf,10, format = '("* Sun Elevation Angle: ",F8.3,40x,"*")',sun_el

totpix = samples * lines

sum = 0L
s = 0L

di_threshold = .5
b3thresh = .08
thresh_b56_high = 225.
thresh_b56_low = 210.
b45_thresh = 1.00
b42_thresh = 2.00
b43_thresh = 2.00
b6_thresh = 300.
thresh_NDSI_max = .7
thresh_NDSI_min = -.25
thresh_b3_lower = .07
thermal_effect_high = 35.0
thermal_effect_low = 25.
cloud_mean_temp = 295.
cloud_254_thresh = .2
thresh_NDSI_snow = .8
snow_thresh = 1.
;b43_thresh = 2.35
;b42_thresh = 2.16

targ_samp = 6910
targ_line = 2539

for l = 0L, (lines - 1) do begin
        readu, 2, b2
        readu, 3, b3
        readu, 4, b4
        readu, 5, b5
        readu, 6, b6
        fill_index = where(b2 eq -9999, fill_values)
        b2rflect = double(b2 / 10000.)
	b3rflect = double(b3 / 10000.)
	b4rflect = double(b4 / 10000.)
	b5rflect = double(b5 / 10000.)
;	b6_temp(0:samples-1,l) = double(b6 / 100) + 273   ; fixed(jk) below
  	b6_temp(0:samples-1,l) = (double(b6) / 100.0) + 273.0
	if(fill_values gt 0)then begin
	  b2rflect(fill_index) = 0
	  b3rflect(fill_index) = 0
	  b4rflect(fill_index) = 0
	  b5rflect(fill_index) = 0
	  b6_temp(fill_index,l) = 0
	endif  
;-----------------------------------------------------------------------------------
; band 3 reflectance threshold - screens out bright targets.
;-----------------------------------------------------------------------------------
	for s = 0, samples-1 do begin	
	  if b3rflect(s) gt b3thresh then begin 
          if(s eq targ_samp and l eq targ_line)then begin
            print,' target pixel b3rflect: ',b3rflect(s)
          endif
		ndsi = (b2rflect(s) - b5rflect(s)) / (b2rflect(s) + b5rflect(s))
                if(s eq targ_samp and l eq targ_line)then begin
                  print,' target pixel b5rflect: ',b5rflect(s)
                  print,' target pixel: ndsi: ',ndsi
                endif
		count1 = count1 + 1
;-----------------------------------------------------------------------------------
; normalized difference snow index. Less that .4? Snow eliminated, but not ice.
;-----------------------------------------------------------------------------------
		if (ndsi lt thresh_NDSI_max and ndsi gt thresh_NDSI_min) then begin
		  count2 = count2 + 1	
                  if(s eq targ_samp and l eq targ_line)then begin
                    print,' target pixel: b3rflect: ',b3rflect(s)
                  endif	
		  if(b6_temp(s,l) lt b6_thresh) then begin
                     if(s eq targ_samp and l eq targ_line)then begin
                       print,' target pixel: b6 temp: ', b6_temp(s,l)
                     endif
		     b6delete = b6delete + 1
		     b56thresh = ( 1. - b5rflect(s) ) * b6_temp(s,l)
                     if(s eq targ_samp and l eq targ_line)then begin
                        print,' target pixel: b56thresh: ',b56thresh
                     endif
;-----------------------------------------------------------------------------------		     
; band 6/5 composite. Eliminate ice
; placed liberally at 225....secondary filter at 210
;-----------------------------------------------------------------------------------
		     if(b56thresh lt thresh_b56_high) then begin 
			   count6 = count6 + 1
			   b43ratio = b4rflect(s) / b3rflect(s)
;-----------------------------------------------------------------------------------
; band 4/3 ratio (equals about 1.0 for clouds), elimintes bright veg and soil.
;-----------------------------------------------------------------------------------
			   if(b43ratio lt b43_thresh)then begin			
            			count3 = count3 + 1
            			b42ratio = b4rflect(s) / b2rflect(s)
;-----------------------------------------------------------------------------------
; 4/2 ratio for senescing vegetation
;-----------------------------------------------------------------------------------
            			if(b42ratio lt b42_thresh)then begin
            			  count42 = count42 + 1
				  b45ratio = b4rflect(s) / b5rflect(s)
;----------------------------------------------------------------------------------- 
; critical for desert regions, 174-47 - Sudan
; band 4/5 ratio. Eliminates rocks and desert    
;-----------------------------------------------------------------------------------      				  
				  if(b45ratio gt b45_thresh)then begin				    
				    count4 = count4 + 1			
				    clmask(s,l) = 255
				    if(b56thresh lt thresh_b56_low)then begin
				      clmask(s,l) = 254
;				      desert_mask(i) = 255	
				    endif  					
			   	 endif		 ; b45 ratio endif
			        endif		 ; b42 endif
			       endif		 ; b43 endif
			    endif else begin 	 ; b56 endif	
			      if(b5rflect(s) lt .08)then clmask(s,l) = 55			    
			    endelse	
			endif else begin	 ; b6 endif
			   clmask(s,l) = 55
			   hotpixels = hotpixels + 1
			endelse
		      endif else begin		 ; ndsi endif else
		        if(ndsi gt thresh_NDSI_snow)then begin
		          snowtab = snowtab + 1
		          snowmask(s,l) = 255
		        endif
		        clmask(s,l) = 55
		        NDSItab = NDSItab + 1		  
  		      endelse			 ; ndsi endelse
  		    endif else begin		 ; b3thresh endif
  		      if(b3rflect(s) lt thresh_b3_lower)then begin
  		        clmask(s,l) = 55
		        B3tab = B3tab + 1
  		      endif 	      
  		    endelse			 ; b3thresh endelse
    endfor ; samples endfor
    if(fill_values gt 0)then begin
	  clmask(fill_index,l) = 1
    endif  	
    line100 = (l+1) mod 10
    if(line100 eq 0)then print,'lines processed: ',l+1    
endfor ; lines endfor

print,'NDSI minimum:      ',thresh_NDSI_min
print,'NDSI maximum:      ',thresh_NDSI_max
print,'b45 threshold:     ',b45_thresh
print,'b42 threshold:     ',b42_thresh
print,'b43 threshold:     ',b43_thresh
print,'thermal effect high: ',thermal_effect_high
print,'thermal effect low:  ',thermal_effect_low
print,'cloud_mean_temp:     ',cloud_mean_temp
print,'cloud 254 threshold: ',cloud_254_thresh

printf,10, format = '("* Desert Threshold: ",F8.3,43X,"*")',di_threshold
printf,10, format = '("* Band 3 Threshold: ",F8.3,43x,"*")',b3thresh
printf,10, format = '("* NDSI Minimum:     ",F8.3,43x,"*")',thresh_NDSI_min
printf,10, format = '("* NDSI Maximum:     ",F8.3,43x,"*")',thresh_NDSI_max
printf,10, format = '("* B6 Threshold:     ",F8.3,43x,"*")',b6_thresh
printf,10, format = '("* B43 Threshold:    ",F8.3,43x,"*")',b43_thresh
printf,10, format = '("* B42 Threshold:    ",F8.3,43x,"*")',b42_thresh
printf,10, format = '("* B45 Threshold:    ",F8.3,43x,"*")',b45_thresh
printf,10, format = '("* B45 Threshold:    ",F8.3,43x,"*")',thresh_b56_high
printf,10, format = '("* B45 Threshold:    ",F8.3,43x,"*")',thresh_b56_low
printf,10, format = '("* B3 Lower Threshold: ",F8.3,41x,"*")',thresh_b3_lower
printf,10, format = '("* NDSI Snow Minimum: ",F8.3,42x,"*")',thresh_NDSI_snow
printf,10, format = '("*",70x,"*")'
printf,10, format = '("* Pass Two Thresholds ",49x,"*")'
printf,10, format = '("*",70x,"*")'
printf,10, format = '("* Snow Threshold: ",F8.3,45x,"*")',snow_thresh
printf,10, format = '("* Cloud 254 Threshold:   ",F8.3,38x,"*")',cloud_254_thresh
printf,10, format = '("* Allowable Cloud Mean:  ",F8.3,38x,"*")',cloud_mean_temp
printf,10, format = '("* Thermal Effect High:   ",F8.3,38x,"*")',thermal_effect_high
printf,10, format = '("* Thermal Effect Low:    ",F8.3,38x,"*")',thermal_effect_low
;----------------------------------------------------------------------------------- 
; write pass 2 mask
;-----------------------------------------------------------------------------------                 

writeu,7,clmask

landtab = B3tab + NDSItab + hotpixels

print,'band 3 threshold:         ',b3thresh
print,'************************************************************************'
print,'*                          REAL LANDSAT 7 DATA                         *'
print,'*                            ----------------                          *'
print,'*                          PIXEL FILTER TALLIES                        *'
print,'*                                                                      *'
print, format = '("* -> total image pixels:      ",i8,"                                 *")', $
totpix
print, format = '("* -> band 3 reflectance count:",i8,"                                 *")', $
count1
print, format = '("* -> band 5/2 NDSI count:     ",i8,"                                 *")', $
count2
print, format = '("* -> band 6 < 300K count:     ",i8,"                                 *")', $
b6delete
print, format = '("* -> band 5/6 threshold:      ",i8,"                                 *")', $
count6
print, format = '("* -> band 43 ratio count:     ",i8,"                                 *")', $
count3
print, format = '("* -> band 42 ratio count:     ",i8,"                                 *")', $
count42
print, format = '("* -> band 4/5 threshold count:",i8,"                                 *")', $
count4
print, format = '("* -> NDSI pixels:             ",i8,"                                 *")', $
NDSItab
print, format = '("* -> band 3 pixels eliminated:",i8,"                                 *")', $
B3tab
print, format = '("* -> hot pixels eliminated:   ",i8,"                                 *")', $
hotpixels
print, format = '("* -> snow pixels eliminated:  ",i8,"                                 *")', $
snowtab
print, format = '("* -> unambiguous land total:  ",i8,"                                 *")', $
landtab

printf,10,'band 3 threshold:         ',b3thresh
printf,10,'************************************************************************'
printf,10,'*                          REAL LANDSAT 7 DATA                         *'
printf,10,'*                            ----------------                          *'
printf,10,'*                          PIXEL FILTER TALLIES                        *'
printf,10,'*                                                                      *'
printf,10, format = '("* -> total image pixels:      ",i8,"                                 *")', $
totpix
printf,10, format = '("* -> band 3 reflectance count:",i8,"                                 *")', $
count1
printf,10, format = '("* -> band 5/2 NDSI count:     ",i8,"                                 *")', $
count2
printf,10, format = '("* -> band 6 < 300K count:     ",i8,"                                 *")', $
b6delete
printf,10, format = '("* -> band 5/6 threshold:      ",i8,"                                 *")', $
count6
printf,10, format = '("* -> band 43 ratio count:     ",i8,"                                 *")', $
count3
printf,10, format = '("* -> band 42 ratio count:     ",i8,"                                 *")', $
count42
printf,10, format = '("* -> band 4/5 threshold count:",i8,"                                 *")', $
count4
printf,10, format = '("* -> NDSI pixels:             ",i8,"                                 *")', $
NDSItab
printf,10, format = '("* -> band 3 pixels eliminated:",i8,"                                 *")', $
B3tab
printf,10, format = '("* -> hot pixels eliminated:   ",i8,"                                 *")', $
hotpixels
printf,10, format = '("* -> snow pixels eliminated:  ",i8,"                                 *")', $
snowtab
printf,10, format = '("* -> unambiguous land total:  ",i8,"                                 *")', $
landtab

clczech = 0L
cloudvalues = 0L

index = where(clmask eq 255 or clmask eq 254, cloudvalues)
clpercentorig =  cloudvalues / float(totpix) * 100
index = where(clmask eq 254, cloudcount)
clpercent254 = cloudcount / float(totpix) * 100
print,'************************************************************************'
print,'*              CLOUD ASSESSMENT PRIOR TO THERMAL FILTER                *'
print,'*                                                                      *'
print, format = '("* -> moderate cloud count    : ",i8,"  ",F5.2,"%                        *")', $
cloudvalues, clpercentorig
print, format = '("* -> conservative cloud count: ",i8,"  ",F5.2,"%                        *")', $
cloudcount, clpercent254

printf,10,'************************************************************************'
printf,10,'*              CLOUD ASSESSMENT PRIOR TO THERMAL FILTER                *'
printf,10,'*                                                                      *'
printf,10, format = '("* -> moderate cloud count    : ",i8,"  ",F5.2,"%                        *")', $
cloudvalues, clpercentorig
printf,10, format = '("* -> conservative cloud count: ",i8,"  ",F5.2,"%                        *")', $
cloudcount, clpercent254

landvalues = 0L
;-----------------------------------------------------------------------------------
; compute cloud mean values and standard deviation but only if clouds were 
; at the half of percent level.
;-----------------------------------------------------------------------------------
desert_index = float(count4) / float(count42)
snpercent =  snowtab / float(totpix) * 100
print,format = '("* -> desert index:   ",F5.3,"                                             *")', $
desert_index
printf,10,format = '("* -> desert index:   ",F5.3,"                                             *")', $
desert_index
;-----------------------------------------------------------------------------------
;If snow is present in the image then ice is too. Use the more conservative 
; cloud mask (254) for developing the training statistics or for determining
; the cloud count if we're over desert.
;-----------------------------------------------------------------------------------

if(snpercent gt 1. or desert_index le di_threshold)then begin
  index = where(clmask eq 255, cloudvalues)
  if(cloudvalues gt 0)then clmask(index) = 0
  index = where(clmask eq 254, cloudvalues)
  print,'************************************************************************'
  print,'*                  Conservative statistics employed!                   *'
  print,format = '("* Snow (Ice) present ",F6.3,"%                                           *")', $
  snpercent
  print,format = '("* Desert Region - index: ",F6.3,"                                        *")', $
  desert_index
  print,'************************************************************************'

  printf,10,'************************************************************************'
  printf,10,'*                  Conservative statistics employed!                   *'
  printf,10,format = '("* Snow (Ice) present ",F6.3,"%                                           *")', $
  snpercent
  printf,10,format = '("* Desert Region - index: ",F6.3,"                                        *")', $
  desert_index
  printf,10,'************************************************************************'
endif else begin
    index = where(clmask eq 255 or clmask eq 254, cloudvalues)
    print,'************************************************************************'
    print,'* Snow (Ice) absent ',snpercent,'% moderate statistics employed!       *'
    print,'************************************************************************'

    printf,10,'************************************************************************'
    printf,10,'* Snow (Ice) absent ',snpercent,'% moderate statistics employed!       *'
    printf,10,'************************************************************************'
endelse
  

;this if was added 3/6/99 and is not in the LPS code.... arctic images made this necessary

if(cloudvalues gt 0)then begin
  clmask(index) = 255
  b6clouds = b6_temp(index)
  result = moment(b6clouds)
  b6mean = result(0)
endif else begin
  print,'b6 mean to be zeroed because there is no band 6 mean'
  print,'cloudcount:',cloudvalues

  printf,10,'b6 mean to be zeroed because there is no band 6 mean'
  printf,10,'cloudcount:',cloudvalues
  b6mean = 0
endelse  

if(clpercent254 gt cloud_254_thresh and desert_index gt di_threshold and b6mean lt cloud_mean_temp)then begin

  print, format = '("* cloud thermal mean:      ",F6.2,38x,"*")',result(0)
  print, format = '("* cloud thermal variance : ",F6.2,38x,"*")',result(1)
  print, format = '("* cloud thermal skew:      ",F6.2,38x,"*")',result(2)
  print, format = '("* cloud kurtosis:          ",F6.2,38x,"*")',result(3)

  printf,10, format = '("* cloud thermal mean:      ",F6.2,38x,"*")',result(0)
  printf,10, format = '("* cloud thermal variance : ",F6.2,38x,"*")',result(1)
  printf,10, format = '("* cloud thermal skew:      ",F6.2,38x,"*")',result(2)
  printf,10, format = '("* cloud kurtosis:          ",F6.2,38x,"*")',result(3)
;--------------------------------------------------------------------------------------------------------------
; this .95 thermal rule was born form 190 44 (JK)
;--------------------------------------------------------------------------------------------------------------
; compute cloud signature standard deviation and histogram the cloud temperatures.
;--------------------------------------------------------------------------------------------------------------
     b6sd = stdev(b6clouds,result(0))
     b6histogram = histogram(b6clouds,  OMAX = b6max,  OMIN = b6min)

     bins = ceil(b6max - b6min)
     print, format = '("* cloud maximum:           ",F6.2,38x,"*")',b6max
     print, format = '("* cloud minimum:           ",F6.2,38x,"*")',b6min
     print, format = '("* standard deviation:      ",F6.2,38x,"*")',b6sd
     print, format = '("* number of histogram bins: ",I6,37x,"*")',bins

     printf,10, format = '("* cloud maximum:           ",F6.2,38x,"*")',b6max
     printf,10, format = '("* cloud minimum:           ",F6.2,38x,"*")',b6min
     printf,10, format = '("* standard deviation:      ",F6.2,38x,"*")',b6sd
     printf,10, format = '("* number of histogram bins: ",I6,37x,"*")',bins

     for i = 0L, (bins-1) do begin
       print,'histogram bin: ',i+1,' ',b6histogram(i)
       printf,10,'histogram bin: ',i+1,' ',b6histogram(i)
     endfor
;--------------------------------------------------------------------------------------------------------------
; examine the thermal skew and set the skew shift factor accordingly.  The skew shift
; factor is limited to between -1 and 1.
;--------------------------------------------------------------------------------------------------------------
     shiftfact = result(2)
     print, format = '("* thermal threshold shift factor: ",F6.2,31x,"*")',shiftfact
     printf,10, format = '("* thermal threshold shift factor: ",F6.2,31x,"*")',shiftfact

     if(result(2) gt 1.)then shiftfact = 1.
     if(result(2) lt -1.)then shiftfact = -1.
     print, format = '("* thermal threshold new shift factor: ",F6.2,27x,"*")',shiftfact
     printf,10, format = '("* thermal threshold new shift factor: ",F6.2,27x,"*")',shiftfact
     shift = b6sd * shiftfact
     print, format = '("* threshold skew shift (shiftfact * b6sd): ",F6.2,22x,"*")',shift
     printf,10, format = '("* threshold skew shift (shiftfact * b6sd): ",F6.2,22x,"*")',shift
;--------------------------------------------------------------------------------------------------------------
; calculate 3 cloud temperature thresholds, 1.25%, 2.5%, and 17.5% below the max.
; first determine the total pixel count at each percentage levelÉ..
;--------------------------------------------------------------------------------------------------------------
print,'************************************************************************'
     print,'*                      Threshold Calculations                          *'
     print,'*                      ----------------------                          *'
     cloudcut0125 = long(cloudvalues * .0125)
     print, format = '("* 1.25% of cloud pixels: ",I8,38x,"*")',cloudcut0125
     cloudcut025 = long(cloudvalues * .025) 
     print, format = '("* 2.50% of cloud pixels: ",I8,38x,"*")',cloudcut025
     cloudcut175 = long(cloudvalues * .175)
     print, format = '("* 17.5% of cloud pixels: ",I8,38x,"*")',cloudcut175

printf,10,'************************************************************************'
     printf,10,'*                      Threshold Calculations                          *'
     printf,10,'*                      ----------------------                          *'
     printf,10, format = '("* 1.25% of cloud pixels: ",I8,38x,"*")',cloudcut0125
     printf,10, format = '("* 2.50% of cloud pixels: ",I8,38x,"*")',cloudcut025
     printf,10, format = '("* 17.5% of cloud pixels: ",I8,38x,"*")',cloudcut175
;     print,b6histogram
;--------------------------------------------------------------------------------------------------------------
; now determine the actual histogram bin number represented by each of the 3 %'s.
;--------------------------------------------------------------------------------------------------------------
     cloudsum = 0
     for i = 1, bins do begin
        cloudsum = cloudsum + b6histogram(bins - i)
        if (cloudsum ge cloudcut0125) then goto, gotit
     endfor
     gotit:
     b6thresh_max = b6max - i
     print, format = '("* 1.25% bin number: ",I3," cloud pixels excluded:",i8,17x,"*")',i,cloudsum
     printf,10, format = '("* 1.25% bin number: ",I3," cloud pixels excluded:",i8,17x,"*")',i,cloudsum     
     cloudsum = 0
     for i = 1, bins do begin
       cloudsum = cloudsum + b6histogram(bins - i)
       if (cloudsum ge cloudcut025) then goto, gotit2
     endfor
     gotit2:
     b6thresh1 = b6max - i 
     print, format = '("* 2.50% bin number: ",I3," cloud pixels excluded:",i8,17x,"*")',i,cloudsum
     printf,10, format = '("* 2.50% bin number: ",I3," cloud pixels excluded:",i8,17x,"*")',i,cloudsum
     cloudsum = 0
     for i = 1, bins do begin
         cloudsum = cloudsum + b6histogram(bins - i)
         if (cloudsum ge cloudcut175) then goto, gotit3
     endfor
     gotit3:   
     b6thresh2 = b6max - i
     print, format = '("* 17.5% bin number: ",I3," cloud pixels excluded:",i8,17x,"*")',i,cloudsum
     print, format = '("* band 6 threshold 98.75%: ",I8,36x,"*")',b6thresh_max
     print, format = '("* band 6 threshold 97.50%: ",I8,36x,"*")',b6thresh1
     print, format = '("* band 6 threshold 82.50%: ",I8,36x,"*")',b6thresh2

     print,'************************************************************************'
     print,'*                    SECOND PASS RESULTS                               *'

     printf,10, format = '("* 17.5% bin number: ",I3," cloud pixels excluded:",i8,17x,"*")',i,cloudsum
     printf,10, format = '("* band 6 threshold 98.75%: ",I8,36x,"*")',b6thresh_max
     printf,10, format = '("* band 6 threshold 97.50%: ",I8,36x,"*")',b6thresh1
     printf,10, format = '("* band 6 threshold 82.50%: ",I8,36x,"*")',b6thresh2

     printf,10,'************************************************************************'
     printf,10,'*                    SECOND PASS RESULTS                               *'

;--------------------------------------------------------------------------------------------------------------
; the b6skew threshold was borne from the 169/28 image. The old algorithm had
; the mean + 2SD as the threshold. This yielded an enormous threshold which
; classified the entire image as cloud (large neg skew) 
;--------------------------------------------------------------------------------------------------------------
; exploit the cloud distribution skew characteristics but only if the skew is positive.
; change the two working thresholds (2.5%, 17.5%) by the computed skew shift.
;--------------------------------------------------------------------------------------------------------------
    if(shiftfact gt 0)then begin
      print,'* Threshold Skew Enabled                                              *'
      printf,10,'* Threshold Skew Enabled  
      b6thresh1_save = b6thresh1
      b6thresh2_save = b6thresh2
      
      b6thresh1 = b6thresh1 + shift
      b6thresh2 = b6thresh2 + shift 
;--------------------------------------------------------------------------------------------------------------
; does the skew adjusted 97.5% threshold exceed the maximum (98.75%)?  
; If so, set it to the maximum.  Determine the acceptable skew adjustment and
; add it to the original 82.5% threshold to get a new threshold.
;--------------------------------------------------------------------------------------------------------------

      if(b6thresh1 gt b6thresh_max)then begin        
        b6thresh1 = b6thresh_max
        b6thresh2 = (b6thresh1 - b6thresh1_save) + b6thresh2_save
        print, format = '("* Threshold exceeds maximum allowed (1.25%): ",F6.2,19x,"*")',b6thresh_max
        print, format = '("*  1.25% threshold shifted: ",F6.2,36x,"*")',b6thresh1
        print, format = '("* 17.50% threshold shifted: ",F6.2,36x,"*")',b6thresh2

        printf,10, format = '("* Threshold exceeds maximum allowed (1.25%): ",F6.2,19x,"*")',b6thresh_max
        printf,10, format = '("*  1.25% threshold shifted: ",F6.2,36x,"*")',b6thresh1
        printf,10, format = '("* 17.50% threshold shifted: ",F6.2,36x,"*")',b6thresh2

      endif else begin      
        print, format = '("*  shifted threshold employed.. maximum not exceeded",18x,"*")'    
        print, format = '("*  1.25% threshold shifted: ",F6.2,36x,"*")',b6thresh1
        print, format = '("* 17.50% threshold shifted: ",F6.2,36x,"*")',b6thresh2   

        printf,10, format = '("*  shifted threshold employed.. maximum not exceeded",18x,"*")'    
        printf,10, format = '("*  1.25% threshold shifted: ",F6.2,36x,"*")',b6thresh1
        printf,10, format = '("* 17.50% threshold shifted: ",F6.2,36x,"*")',b6thresh2        
     
      endelse
    endif else begin
       print,'* Threshold Skew Disabled (Negative)                                   *'
       printf,10,'* Threshold Skew Disabled (Negative) 
       print, format = '("*  2.5% thermal threshold (shift): ",F6.2,30x,"*")',b6thresh1
       print, format = '("* 17.5% thermal threshold (shift): ",F6.2,30x,"*")',b6thresh2
       printf,10, format = '("*  2.5% thermal threshold (shift): ",F6.2,30x,"*")',b6thresh1
       printf,10, format = '("* 17.5% thermal threshold (shift): ",F6.2,30x,"*")',b6thresh2
    endelse

    for i = 0L, lines - 1 do begin
      for j = 0L, samples - 1 do begin
        if(clmask(j,i) eq 0)then begin	
	  if(b6_temp(j,i) le b6thresh1)then begin
	    clmask(j,i) = 125
	    count9 = count9 + 1
	    if(b6_temp(j,i) le b6thresh2)then begin	  	
		clmask(j,i) = 200
	    	count10 = count10 + 1
	    endif
	  endif
        endif
      endfor
    endfor
    index = where(clmask eq 125 or clmask eq 200, cloudvalues)
    b6clouds = b6_temp(index)
    
    b6histogram = histogram(b6clouds,  OMAX = b6max1,  OMIN = b6min)
    print,'*----------------------------------------------------------------------*'
    print,'*                     2.5% Threshold Statistics                        *'  
    print,'*                                                                      *' 
    print, format = '("*  cloud minimum: ",F6.2,47x,"*")',b6min
    print, format = '("*  cloud maximum: ",F6.2,47x,"*")',b6max1

    printf,10,'*----------------------------------------------------------------------*'
    printf,10,'*                     2.5% Threshold Statistics                        *'  
    printf,10,'*                                                                      *' 
    printf,10, format = '("*  cloud minimum: ",F6.2,47x,"*")',b6min
    printf,10, format = '("*  cloud maximum: ",F6.2,47x,"*")',b6max1

    result = moment(b6clouds)
    b6mean1 = result(0)
    print, format = '("*  cloud mean   : ",F6.2,47x,"*")',b6mean1
    print,'*----------------------------------------------------------------------*'

    printf,10, format = '("*  cloud mean   : ",F6.2,47x,"*")',b6mean1
    printf,10,'*----------------------------------------------------------------------*'

    index = where(clmask eq 200, cloudvalues)
    if(cloudvalues gt 1)then b6clouds = b6_temp(index)
    b6histogram = histogram(b6clouds,  OMAX = b6max2,  OMIN = b6min)

    if(cloudvalues gt 1)then begin
        print,'*                    17.5% Threshold Statistics                       *' 
        print,'*                                                                     *' 
        print, format = '("*  cloud minimum: ",F6.2,46x,"*")',b6min
        print, format = '("*  cloud maximum: ",F6.2,46x,"*")',b6max2
    	result = moment(b6clouds)
   	b6mean2 = result(0)
   	print, format = '("*  cloud mean   : ",F6.2,46x,"*")',b6mean2

        printf,10,'*                    17.5% Threshold Statistics                       *' 
        printf,10,'*                                                                     *' 
        printf,10, format = '("*  cloud minimum: ",F6.2,46x,"*")',b6min
        printf,10, format = '("*  cloud maximum: ",F6.2,46x,"*")',b6max2
   	printf,10, format = '("*  cloud mean   : ",F6.2,46x,"*")',b6mean2

    endif else begin
    	b6mean2 = b6mean
    endelse

    
  print,'***********************************************************************'
  print, format = '("* pixel count through  2.5% thermal threshold: ",I8,15x,"*")',count9
  print, format = '("* pixel count through 17.5% thermal threshold: ",I8,15x,"*")',count10

  printf,10,'***********************************************************************'
  printf,10, format = '("* pixel count through  2.5% thermal threshold: ",I8,15x,"*")',count9
  printf,10, format = '("* pixel count through 17.5% thermal threshold: ",I8,15x,"*")',count10

;  tv,clmask,channel=0
;----------------------------------------------------------------------------------
; determine cloud percentage for 5% threshold and thermal effect
;----------------------------------------------------------------------------------
  thermeffecta =  count9 / float(totpix) * 100
  thermeffectb =  count10 / float(totpix) * 100
  if(snpercent gt 1. or desert_index le di_threshold)then begin
    clpercentpass1 = clpercent254
  endif else begin
    clpercentpass1 = clpercentorig
  endelse
  clpercenta = clpercentpass1 + thermeffecta
  clpercentb = clpercentpass1 + thermeffectb
  print, format = '("* baselined pass 1 cloud percentage:      ",F6.2,"%",21x,"*")',clpercentpass1
  print, format = '("* net thermal effect for the  2.5% filter ",F6.2,"%",21x,"*")',thermeffecta
  print, format = '("* net thermal effect for the 17.5% filter ",F6.2,"%",21x,"*")',thermeffectb
  print, format = '("* cloud % after  2.5% thermal filter:   ",F6.2,"%",23x,"*")',clpercenta
  print, format = '("* cloud % after 17.5% thermal filter:   ",F6.2,"%",23x,"*")',clpercentb

  printf,10, format = '("* baselined pass 1 cloud percentage:      ",F6.2,"%",21x,"*")',clpercentpass1
  printf,10, format = '("* net thermal effect for the  2.5% filter ",F6.2,"%",21x,"*")',thermeffecta
  printf,10, format = '("* net thermal effect for the 17.5% filter ",F6.2,"%",21x,"*")',thermeffectb
  printf,10, format = '("* cloud % after  2.5% thermal filter:   ",F6.2,"%",23x,"*")',clpercenta
  printf,10, format = '("* cloud % after 17.5% thermal filter:   ",F6.2,"%",23x,"*")',clpercentb

;----------------------------------------------------------------------------------
; if a radical thermal effect took place or if snow is present step down
; to a more conservative filter. change to 40 from 30 due to 164 20.
;----------------------------------------------------------------------------------
  b6max_thresh1_diff = b6max - b6thresh1
  if(thermeffecta gt thermal_effect_high or snpercent gt snow_thresh or b6mean1 ge cloud_mean_temp or b6max_thresh1_diff le 2.)then begin         ; 40% to JK
      index = where(clmask eq 200, cloudvalues)
      if(thermeffectb gt thermal_effect_low or b6mean2 gt cloud_mean_temp)then begin
        clmask(index) = 0
        filter_used = 3
      endif else begin
        clmask(index) = 255
        filter_used = 2
      endelse
  endif else begin    
      index = where(clmask eq 200 or clmask eq 125, cloudvalues)
      clmask(index) = 255
      filter_used = 1
  endelse
;----------------------------------------------------------------------------------
; now zero out everything but clouds and compute cloud percentage before
; filtering
;----------------------------------------------------------------------------------
  index = where(clmask ne 255, cloudvalues)
  clmask(index) = 0
  index = where(clmask eq 255, cloudvalues)
  clpercent =  cloudvalues / float(totpix) * 100
;----------------------------------------------------------------------------------
; determine its effect... if still radical forget the filter and go back
; to the original cloud determination result.
;----------------------------------------------------------------------------------
  if(filter_used eq 1)then begin
    print,'*  2.5 percent filter used                                            *'
    printf,10,'*  2.5 percent filter used                                            *'
  endif
  if(filter_used eq 2)then begin
    print,'* 17.5 percent filter used due to:                                    *'
    printf,10,'* 17.5 percent filter used due to:                                    *'
    if(snpercent gt 1.)then begin
      print,format = '("*     the snow percentage",F6.2,"% is greater than 1%",19x,"*")',$
      snpercent
      printf,10,format = '("*     the snow percentage",F6.2,"% is greater than 1%",19x,"*")',$
      snpercent
    endif
    if(thermeffecta gt thermal_effect_high)then begin
      print,format = '("*     the 2.5% thermal effect was greater than 40%",20x,"*")'
      printf,10,format = '("*     the 2.5% thermal effect was greater than 40%",20x,"*")'
    endif
    if(b6max_thresh1_diff le 2.)then begin
      print,format = '("*     the 2.5% threshold is within 2 bins of the band 6 maximum",7x,"*")'
      printf,10,format = '("*     the 2.5% threshold is within 2 bins of the band 6 maximum",7x,"*")'
    endif
  endif
  if(filter_used eq 3)then begin
    print,'* questionable results, thermal pass discarded                        *'
    printf,10,'* questionable results, thermal pass discarded                        *'
  endif
  print, format = '("* cloud % before filtering:  ",F6.2,"                                   *")', $
  clpercent
  print, format = '("* before filtering count:  ",i8,"                                   *")', $
  cloudvalues
  printf,10, format = '("* cloud % before filtering:  ",F6.2,"                                   *")', $
  clpercent
  printf,10, format = '("* before filtering count:  ",i8,"                                   *")', $
  cloudvalues

endif else begin	 ; only if clouds were found endelse (i.e > 4%) and desert index
;----------------------------------------------------------------------------------
; use the conservative cloud count because desert features are present.
;----------------------------------------------------------------------------------
  if(b6mean lt cloud_mean_temp)then begin
    clpercent = clpercent254
    print,'thermal analysis bypassed with a valid b6 mean!'
    print,'band 6 mean:',b6mean
    printf,10,'thermal analysis bypassed with a valid b6 mean!'
    printf,10,'band 6 mean:',b6mean
  endif else begin
    print,'*              thermal analysis bypassed!               *'
    print,' invalid results, band 6 mean: ',b6mean,' zero out all clouds!'
    printf,10,'*              thermal analysis bypassed!               *'
    printf,10,' invalid results, band 6 mean: ',b6mean,' zero out all clouds!'
    clmask(index) = 0
    clpercent = 0.
  endelse
  print,'thermal analysis bypassed!'
  print,'cloud cover percentage: ',clpercent
  print,'desert index:          ',desert_index

  printf,10,'thermal analysis bypassed!'
  printf,10,'cloud cover percentage: ',clpercent
  printf,10,'desert index:          ',desert_index
  
endelse			 ; only if clouds were found endif
;writeu, 7, clmask
;clmaskf = median(clmask,2)

;----------------------------------------------------------------------------------
; fill in cloud holes if clouds exist
;----------------------------------------------------------------------------------
if(clpercent gt 0)then begin
  lineshalf = lines / 2
  lines2 = lineshalf * 2
  if(lines2 ne lines)then begin
    lines_quad12 = lineshalf
    lines_quad34 = lineshalf + 1
  endif else begin
    lines_quad12 = lineshalf
    lines_quad34 = lineshalf 
  endelse
  clmaskfilt = filtmaj(clmask, WIDTH = 3, THRESH = 1275, CELLFLAG = 0, FILL=255,/ALL)
		
  index = where(clmaskfilt eq 255, cloudvalues)
  clpercent =  cloudvalues / float(totpix) * 100
  print, format = '("* final cloud pixel count: ",i8,"                                   *")', $
  cloudvalues
  printf,10, format = '("* final cloud pixel count: ",i8,"                                   *")', $
  cloudvalues


  print, format = '("* -> final cloud cover %: ",F6.2,"                                      *")', $
  clpercent
  printf,10, format = '("* -> final cloud cover %: ",F6.2,"                                      *")', $
  clpercent

  cloud_quad = extrac(clmaskfilt, 15, 0, 192, lines_quad12)
  index = where(cloud_quad eq 255, quadvalues)
  quad1percent = quadvalues / float(totpix/4) * 100
  cloud_quad = extrac(clmaskfilt, 207, 0, 192, lines_quad12)
  index = where(cloud_quad eq 255, quadvalues)
  quad2percent = quadvalues / float(totpix/4) * 100
  cloud_quad = extrac(clmaskfilt, 15, lines/2, 192, lines_quad34)
  index = where(cloud_quad eq 255, quadvalues)
  quad3percent = quadvalues / float(totpix/4) * 100
  cloud_quad = extrac(clmaskfilt, 207, lines/2, 192, lines_quad34)
  index = where(cloud_quad eq 255, quadvalues)
  quad4percent = quadvalues / float(totpix/4) * 100
  print, format = '("* -> quad 1 cloud count:  ",F7.3,"                                     *")', $
quad1percent
  print, format = '("* -> quad 2 cloud count:  ",F7.3,"                                     *")', $
quad2percent
  print, format = '("* -> quad 3 cloud count:  ",F7.3,"                                     *")', $
quad3percent
  print, format = '("* -> quad 4 cloud count:  ",F7.3,"                                     *")', $
quad4percent
  printf,10, format = '("* -> quad 1 cloud count:  ",F7.3,"                                     *")', $
quad1percent
  printf,10, format = '("* -> quad 2 cloud count:  ",F7.3,"                                     *")', $
quad2percent
  printf,10, format = '("* -> quad 3 cloud count:  ",F7.3,"                                     *")', $
quad3percent
  printf,10, format = '("* -> quad 4 cloud count:  ",F7.3,"                                     *")', $
quad4percent
;  tv,clmaskfilt
;  tv,snowmask
;tv,desert_mask
endif

;clmaskfilt2 = filtmaj(clmaskfilt, WIDTH=3, THRESH = 255, CELLFLAG = 255, FILL=0, /ALL)
;tv,clmaskfilt2,channel=0

;index = where(clmaskfilt2 eq 255, cloudvalues)
;clpercent =  cloudvalues / float(totpix) * 100

;print, format = '(" cloud noise removed: ",i8,"  ",F6.2,"%")', $
;cloudvalues, clpercent


;zz = rebin(clmask,512,512,/sample) & tv,zz,channel=0
print,'*                                                                     *'
print,'*                     Processing Completed                            *'
print,'***********************************************************************'

printf,10,'*                                                                     *'
printf,10,'*                     Processing Completed                            *'
printf,10,'***********************************************************************'

writeu,8,clmask
writeu,9,snowmask

close, 2,3,4,5,6, 7, 8, 9, 10

end




function filtmaj,image,WIDTH=box_wid,THRESH=thrsh_hold,CELLFLAG = pix_val,FILL= fill_val,ALL_PIXELS=all_pixels
		
;
; NAME:
;	filter_image
;
; PURPOSE:
;	Computes the average and/or median of pixels in moving box,
;	replacing center pixel with the computed average and/or median,
;		(using the IDL smooth or median functions).
;	The main reason for using this function is the options to
;	also process the pixels at edges and corners of image, and,
;	to apply iterative smoothing simulating convolution with Gaussian,
;	and/or to convolve image with a Gaussian kernel.
;
; CALLING:
;	Result = filter_image( image,/ALL )
;
; INPUT:
;	image = 2-D array (matrix)
;
; KEYWORDS:
;
;	/ALL_PIXELS causes the edges of image to be filtered as well,
;		accomplished by reflecting pixels adjacent to edges outward.
;
;
; RESULT:
;	Function returns the majority filtered image.
;
; EXAMPLES:
;
; EXTERNAL CALLS:
;	function convolve
;	pro factor
;	function prime		;all these called only if FWHM is specified.
;
; PROCEDURE:
;	If /ALL_PIXELS or  then
;	create a larger image by reflecting the edges outward,
;	then call the IDL median and/or smooth function on the larger image,
;	and just return the central part (the orginal size image).
; HISTORY:
;	Written, 1991, Frank Varosi, NASA/GSFC.
;-
	sim = size( image )
	Lx = sim(1)-1
	Ly = sim(2)-1
	
	xsize = sim(1)
	ysize = sim(2)

	if (sim(0) NE 2) OR (sim(4) LE 4) then begin
		message,"input must be an image (a matrix)",/INFO
		return,image
	   endif


	if keyword_set( all_pixels ) then begin
		
		box_wid = fix( box_wid )
		radius = (box_wid/2) > 1
		Lxr = Lx+radius
		Lyr = Ly+radius
		rr = 2*radius
		imf = bytarr( sim(1)+rr, sim(2)+rr )
		imf(radius,radius) = image		; reflect edges outward, 
							; to make larger image.
		imf(  0,0) = rotate( imf(radius:rr,*), 5 )	;Left -   Transpose: True
		imf(Lxr,0) = rotate( imf(Lx:Lxr,*), 5 )		;right -  Transpose: True
		imf(0,  0) = rotate( imf(*,radius:rr), 7 )	;bottom - Transpose: True
		imf(0,Lyr) = rotate( imf(*,Ly:Lyr), 7 )		;top -    Transpose: True
	  endif else begin

		radius=0
		imf = image
	   endelse
	   filtab = 0L
;-----------------------------------------------------------------------------------
; salt
;-----------------------------------------------------------------------------------
	   if(thrsh_hold eq 255) then begin
	   	for y = 1, ysize do begin
	   		for x = 1, xsize do begin
	   			if(imf(x,y) eq pix_val) then begin
					sum = total(imf(x - 1: x + 1, y - 1: y + 1))
					if(sum eq thrsh_hold) then begin
					  image(x,y) = fill_val
;					  imf(x-1,y-1) = fill_val
					  filtab = filtab + 1
					endif
				endif
			endfor
	   	endfor
	   endif else begin

;print,image(310:319,0)
;print,image(320:329,0)
;print,image(330:339,0)
;print,image(340:349,0)
;print,image(350:359,0)
;print,image(360:369,0)
;print,image(370:381,0)

tab = 0
	   for y = 1, ysize do begin
	   		for x = 1, xsize do begin
	   			if(imf(x,y) eq pix_val) then begin
					sum = total(imf(x - 1: x + 1, y - 1: y + 1))
					if(sum ge thrsh_hold) then begin
					  tab = tab + 1
;					  if(tab lt 25)then print,'pixel ',tab+1,' x: ',x,' y:',y
					  image(x-1,y-1) = fill_val
;					  imf(x,y) = fill_val
					  filtab = filtab + 1
					endif
				endif
			endfor
	   	endfor
	   endelse
printf,10, format = '("* filtered pixel total:    ",i8,"                                   *")',filtab
;return, imf( radius:(Lx+radius), radius:(Ly+radius) ) ; fixed(jk) 
return, image  
end

 */
