#include "util.h"
#include "error.h"
/*
!C****************************************************************************

!File: util.c
  
!Description: Utility Functions (these read radiometric tables)

!Revision History:
 Revision 1.0 2004/06/16
 Jonathan Kutler
 Original Version.


!Team Unique Header:
  This software was written for the LEDAPS project and is based on the software 
  written by the MODIS Land Science Team Support Group for the Laboratory for 
  Terrestrial Physics (Code 922) at the   National Aeronautics and Space Administration,
  Goddard Space Flight Center

 ! References and Credits:

      Robert E. Wolfe (Code 922)
      MODIS Land Team Support Group     Raytheon ITSS
      robert.e.wolfe.1@gsfc.nasa.gov    4400 Forbes Blvd.
      phone: 301-614-5508               Lanham, MD 20770  
  
      J. Kutler 
      LEDAPS,                           Raytheon ITSS
      jonathan_l_kutler@raytheon.com    7501 Forbes Blvd, ste 103
      phone: 301-352-2152               Seabrook MD, 20706
  
 ! functions included : 
 
   mygetline   - C line read
   get_wo    - reads the Landsat Work Order file to find the "g-old" gains for each band
   get_gnew  - reads the G-new system table to find current "g-new" gains for each band

!END****************************************************************************
*/
/********************************************************************************/
/************************************ get_wo ************************************/
/********************************************************************************/
Param_wo_t *get_wo(const char* fname)
{
 Param_wo_t *this;
 FILE *iFile;

  char  /*                  char                       */
      line[1024]
    , msgbuf[1024]
    , algorithm[9]
    , completion_date[11]
   ;
  int   /*                  int                        */
      ib
    , jb
    , id
    , band
    , refdet
    , det
      ;

  float   /*                float                        */
        gain
      , bias
      , roff
      , rgain
      , foff
      , fgain
      , gain_sum
      , offset_sum
     ;

 this = (Param_wo_t*)malloc(sizeof(Param_wo_t));
 iFile= fopen(fname,"r");
 if ( !iFile ) 
   { 
   sprintf(msgbuf,"unable to open file %s",fname); 
   ERROR(msgbuf, "get_wo");
   }

 if ( !iFile ) return (Param_wo_t *)NULL;
  while ( mygetline(line,1024,iFile)>=0 )
   if ( (int)strstr(line,"Completion date")>0 )
     {
     strncpy(completion_date,&line[strlen("Completion date:    ")],10); 
     completion_date[10]='\0';
     sscanf(completion_date,"%d %d %d",
      &this->completion_year,&this->completion_month,&this->completion_day);
     break;
     }
  while ( mygetline(line,1024,iFile)>=0 )
   if ( (int)strstr(line,"RADIOMETRIC CORRECTION")>0 ) break;
  while ( mygetline(line,1024,iFile)>=0 )
   if ( (int)strstr(line,"Algorithm:")>0 ) 
     {
     strncpy(algorithm,&line[strlen("Algorithm: ")],8); algorithm[8]='\0';
     break;
     }

  this->nasa_flag= strncmp(algorithm,"NASA",4) ? false : true;

  if ( !strncmp(algorithm,"NASA CPF",8) )
    this->algorithm= ALG_NASA_CPF;
  else if ( !strncmp(algorithm,"NASA",4) )
    this->algorithm= ALG_NASA;
  else
    this->algorithm= ALG_CCRS;

  while( mygetline(line,1024,iFile)>=0 )
    if ( (int)strstr(line,"Detector")>0 && (int)strstr(line,"gain")>0 )break;

  skipline(iFile,1);
  for (ib=0; ib<7; ib++)
    {
    mygetline(line,1024,iFile);
    sscanf(line," %d    |   %d        %f   %f       %5s",
            &band,&refdet,&gain,&bias,msgbuf);
    this->DN_to_Radiance_gain[ib]=   gain;
    this->DN_to_Radiance_offset[ib]= bias;
    this->reference_detector[ib]=    refdet;
    }

  for (ib=0; ib<7; ib++)
    {
    jb= ib==6 ? 5 : ib;
    gain_sum=   0.0;
    offset_sum= 0.0;
    while( mygetline(line,1024,iFile)>=0 )if ( (int)strstr(line,"Band")>0 )break;
     
    sscanf(line,"Band %d",&band);
    if ( band != (ib+1) )
      {
      sprintf(msgbuf,"in file %s expecting band %d found band %d",
 	fname,ib+1,band);
      ERROR(msgbuf, "get_wo");
      }
    skipline(iFile,4);
 
    for (id=0; id<(ib==5?4:16); id++)
      {
      mygetline(line,1024,iFile);
      sscanf(line," %d      | %f   %f    %f   %f",
                   &det,    &fgain,&foff,&rgain,&roff);
      if ( det != (id+1) )
        {
        sprintf(msgbuf,
         "in file %s band %d expecting detector %d found detector %d",
         fname,ib+1,id+1,det);
         ERROR(msgbuf, "get_wo");
        }

    if ( ib==5 )
      {
      this->thermal_forward_gains[id]=       fgain;
      this->thermal_forward_offsets[id]=     foff;
      this->thermal_reverse_gains[id]=       rgain;
      this->thermal_reverse_offsets[id]=     roff;
      }
    else
      {
      this->reflective_forward_gains[jb][id]=   fgain;
      this->reflective_forward_offsets[jb][id]= foff;
      this->reflective_reverse_gains[jb][id]=   rgain;
      this->reflective_reverse_offsets[jb][id]= roff;
      }
    gain_sum+=   ( fgain + rgain );
    offset_sum+= ( foff + roff );
    }
  if ( this->nasa_flag )
    {
    this->final_gain[ib]=   gain_sum   / (ib==5 ? 8.0 : 32.0);
    this->final_offset[ib]= offset_sum / (ib==5 ? 8.0 : 32.0);
    }
  else
    {
    refdet= this->reference_detector[ib];
    this->final_gain[ib]= ( ib==5 ) ? 
           ( this->thermal_forward_gains[refdet-1] +
             this->thermal_reverse_gains[refdet-1] )  / 2.0 :
           ( this->reflective_forward_gains[jb][refdet-1] +
             this->reflective_reverse_gains[jb][refdet-1] ) / 2.0 ;

     
    this->final_offset[ib]=  ( ib==5 ) ? 
           ( this->thermal_forward_offsets[refdet-1] +
             this->thermal_reverse_offsets[refdet-1] )  / 2.0 :
           ( this->reflective_forward_offsets[jb][refdet-1] +
             this->reflective_reverse_offsets[jb][refdet-1] ) / 2.0;
    }
  }

 return this;

}
/********************************************************************************/
/*********************************** get_gnew ***********************************/
/********************************************************************************/
Gains_t* get_gnew(const char* fname, Date_t* adate)
{
 Gains_t* this;
 FILE *iFile;

 char  /*                  char                       */
      line[1024]
    , msgbuf[1024]
   ;

 float   /*                float                        */
        gain1
      , gain2
      , gain3
      , gain4
      , gain5
      , gain7
      , gaino5
      , gaino7
   ;


 double  /*                double                       */
         dyear
   ;

 int   /*                  int                        */
      dsl
    , doy
   ;
 bool equation=false;

 this = (Gains_t*)malloc(sizeof(Gains_t));
 iFile= fopen(fname,"r");
 if ( !iFile ) 
   { 
   sprintf(msgbuf,"unable to open file %s",fname); 
   ERROR(msgbuf, "get_gnew");
   }

  mygetline(line,1024,iFile);
  if ( (int)strstr(line,"Gnew[ib]=")>0 )equation=true;

  skipline(iFile,(equation?2:1));
  this->valid_flag= false;

  while ( mygetline(line,1024,iFile)>=0 )
   {
   if ( equation )
     {
     sscanf(line,
         "%d	%lg	%d	%f	%f	%f	%f	%f	%f",
     &dsl,&dyear,&doy,&gain1,&gain2,&gain3,&gain4,&gain5,&gain7);
     }
   else
     {
     sscanf(line,
  "%d	%lg	%d	%f	%f	%f	%f	%f	%f	%f	%f",
     &dsl,&dyear,&doy,&gain1,&gain2,&gain3,&gain4,&gain5,&gain7,&gaino5,&gaino7);
     }


    if ( (int)dyear == adate->year && doy == adate->doy )
      {
      this->gains[0]= gain1;
      this->gains[1]= gain2;
      this->gains[2]= gain3;
      this->gains[3]= gain4;
      this->gains[4]= equation ? gain5 : gaino5;
      this->gains[6]= equation ? gain7 : gaino7;
      this->valid_flag= true; 
      break;
      }
   }
 return this;
}
/********************************************************************************/
/*********************************** get_gold ***********************************/
/********************************************************************************/
Gains_t* get_gold(const char* fname, Date_t* adate)
{
 Gains_t *this;
 FILE *iFile;

 char  /*                  char                        */
      line[1024]
    , msgbuf[1024]
   ;

 double  /*                double                      */
        gain1
      , gain2
      , gain3
      , gain4
      , gain5
      , gain7
      , ysl
   ;

 int   /*                  int                        */
      dsl
    , doy
    , year
   ;

 this = (Gains_t*)malloc(sizeof(Gains_t));
 iFile= fopen(fname,"r");
 if ( !iFile ) 
   { 
   sprintf(msgbuf,"unable to open file %s",fname); 
   ERROR(msgbuf, "get_gold");
   }

  skipline(iFile,2);
  this->valid_flag= false;

  while ( mygetline(line,1024,iFile)>=0 )
   {
   sscanf(line,"%d %d %d %lg %lg %lg %lg %lg %lg %lg",
          &year,&doy,&dsl,&ysl,&gain1,&gain2,&gain3,&gain4,&gain5,&gain7);

    if ( year == adate->year && doy == adate->doy )
      {
      this->gains[0]= (float)gain1;
      this->gains[1]= (float)gain2;
      this->gains[2]= (float)gain3;
      this->gains[3]= (float)gain4;
      this->gains[4]= (float)gain5;
      this->gains[6]= (float)gain7;
      this->valid_flag= true; 
      break;
      }
   }
 return this;
}
/********************************************************************************/
/*********************************** mygetline  ***********************************/
/********************************************************************************/
int mygetline(char* l,int m,FILE* f)
{
if ( fgets(l,m,f) == NULL )return -1;
else { if ( l[strlen(l)-1]=='\n' )l[strlen(l)-1]='\0'; return strlen(l);  }
}
/********************************************************************************/
/*********************************** skipline ***********************************/
/********************************************************************************/
int skipline(FILE* f, int n)
{
char l[1024];
int i;
for ( i=0; i<n; i++)if ( fgets(l,1024,f) == NULL )break;
return i+1;
}
/*************************************************************************
 *** this program zooms an integer array                               ***
 *************************************************************************/
void zoomIt(int*o,int*i,int n,int z){int j,c=n*z;for(j=0;j<c;j++)o[j]=i[j/z];}
/********************************************************************************/
void zoomIt8(char*o,char*i,int n,int z){int j,c=n*z;for(j=0;j<c;j++)o[j]=i[j/z];}
/************************************* eof **************************************/
/********************************************************************************/
int big_endian(){long int i=1,j=0; char b=1; memcpy(&j,&b,1); return i!=j;}
