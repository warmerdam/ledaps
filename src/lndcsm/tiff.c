#include <stdio.h>       /* for true, false */
#include <string.h>     /* for strlen      */
#include <math.h>
#include <stdlib.h>
#include <ctype.h> 
#include <string.h>
#include <time.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "tiff.h"
#include "error.h"
char msgbuf[1024];
bool tiffHdr (char* hdrdata, imgfile_t* iFile) 
{
const short int BITS_PER_SAMPLE_ADDR=174; /*-> for color move to BitsPerSamp*/
int  NB       = iFile->nb;

char* buf= hdrdata;
int  color=0;
int  NumPixels= iFile->np;
int  NumLines=  iFile->nl;
int  SAMPLE_COUNT=1;
int  SampPer=  iFile->ps;
int  Photo;

const short int  fortytwo= 42;
const       int  FIRST_OFFSET=8;
      short int  NUM_DIRECT= 11;

const short int IMG_WIDTH_TAG=256;         /* IDF ENTRY 2*/
const short int IMG_WIDTH_TYPE=4;
const      int IMG_WIDTH_COUNT=1;
const      int IMG_WIDTH= NumPixels;

const short int IMG_LENGTH_TAG=257;        /* IDF ENTRY 3*/
const short int IMG_LENGTH_TYPE=4;
const      int IMG_LENGTH_COUNT=1;
const      int IMG_LENGTH=NumLines;

const short int BITS_PER_SAMPLE_TAG=258;   /* IDF ENTRY 4*/
const short int BITS_PER_SAMPLE_TYPE=3;
            int BITS_PER_SAMPLE_COUNT;
            int BITS_PER_SAMPLE_c;
const short int BITS_PER_SAMPLE=8*iFile->ps;
const short int BITS_PER_SAMPLE_FILL=0;

const short int COMPRESSION_TAG=259;       /* IDF ENTRY 5*/
const short int COMPRESSION_TYPE=3;
const      int COMPRESSION_COUNT=1;
const short int COMPRESSION=1;
const short int COMPRESSION_FILL=0;

const short int PHOTOINTERP_TAG=262;       /* IDF ENTRY 6*/
const short int PHOTOINTERP_TYPE=3;
const      int PHOTOINTERP_COUNT=1;
      short int PHOTOINTERP= Photo;
const short int PHOTOINTERP_FILL=0;

const short int STRIP_OFFSETS_TAG=273;     /* IDF ENTRY 7*/
const short int STRIP_OFFSETS_TYPE=4;
const      int STRIP_OFFSETS_COUNT=1;
const      int STRIP_OFFSETS=512;

const short int SAMP_PIX_TAG=277;          /* IDF ENTRY 8*/
const short int SAMP_PIX_TYPE=3;
const      int SAMP_PIX_COUNT=1;
      short int SAMP_PIX= SampPer;
const short int SAMP_PIX_FILL=0;

const short int ROWS_PER_STRIP_TAG=278;    /* IDF ENTRY 9*/
const short int ROWS_PER_STRIP_TYPE=4;
const      int ROWS_PER_STRIP_COUNT=1;
const      int ROWS_PER_STRIP= NumLines ;

const short int STRIP_BYTE_COUNT_TAG=279;  /* IDF ENTRY 10*/
const short int STRIP_BYTE_COUNT_TYPE=4;
const      int STRIP_BYTE_COUNT_COUNT=1;
           int STRIP_BYTE_COUNT= ( NumPixels * NumLines * SampPer  );

const short int X_RESOLUTION_TAG=282;      /* IDF ENTRY 11*/
const short int X_RESOLUTION_TYPE=5;
const      int X_RESOLUTION_COUNT=1;
const      int X_RESOLUTION_ADDR=158;

const short int Y_RESOLUTION_TAG=283;      /* IDF ENTRY 12*/
const short int Y_RESOLUTION_TYPE=5;
const      int Y_RESOLUTION_COUNT=1;
const      int Y_RESOLUTION_ADDR=166;

const short int RESOLUTION_UNIT_TAG=296;   /* IDF ENTRY 13*/
const short int RESOLUTION_UNIT_TYPE=3;
const      int RESOLUTION_UNIT_COUNT=1;
const short int RESOLUTION_UNIT=1;
const short int RESOLUTION_UNIT_FILL=0;

const      int NEXT_IDF_OFFSET=0;
const      int X_RESOLUTION[2]={1,1};
const      int Y_RESOLUTION[2]={1,1};

const short int BITS_PER_SAMPLE_COLOR[3]={8,8,8};
int L=4;
int i=0;

if ( NB > 1 )color= 1;
NumPixels= iFile->np;
NumLines=  iFile->nl;
SAMPLE_COUNT= 1;
SampPer=   iFile->ps;
Photo=        1;

if ( color ) 
   {
   SAMPLE_COUNT= 3;
   SampPer=      3;
   Photo=        2;
   }

PHOTOINTERP= Photo;
if ( NB>1 )NUM_DIRECT++;

STRIP_BYTE_COUNT= ( NumPixels * NumLines * SampPer  );
SAMP_PIX= SampPer;
BITS_PER_SAMPLE_COUNT= SAMPLE_COUNT;
BITS_PER_SAMPLE_c=BITS_PER_SAMPLE_ADDR;

for (i=0; i<512; i++) buf[i]=0;

i=0;
memcpy( &buf[i+=0], (big_endian()? "MM" : "II"),      2);  /* 0*/
memcpy( &buf[i+=2], &fortytwo,     2);  /* 2*/
memcpy( &buf[i+=2], &FIRST_OFFSET, 4);  /* 4*/
memcpy( &buf[i+=4], &NUM_DIRECT,   2);  /* 8*/

memcpy( &buf[i+=2], &IMG_WIDTH_TAG,2);          /* IDF ENTRY 2*/
memcpy( &buf[i+=2], &IMG_WIDTH_TYPE,2);  
memcpy( &buf[i+=2], &IMG_WIDTH_COUNT,4);   
memcpy( &buf[i+=4], &IMG_WIDTH,4);

memcpy( &buf[i+=4], &IMG_LENGTH_TAG,2);        /* IDF ENTRY 3*/
memcpy( &buf[i+=2], &IMG_LENGTH_TYPE,2);   
memcpy( &buf[i+=2], &IMG_LENGTH_COUNT,4);   
memcpy( &buf[i+=4], &IMG_LENGTH,4);

memcpy( &buf[i+=4], &BITS_PER_SAMPLE_TAG,2);      /* IDF ENTRY 4*/
memcpy( &buf[i+=2], &BITS_PER_SAMPLE_TYPE,2);   
memcpy( &buf[i+=2], &BITS_PER_SAMPLE_COUNT,4);  
if ( color ) 
   { 
   memcpy( &buf[i+=4], &BITS_PER_SAMPLE_c,4);
   i+=2;
   }
else
   {
   memcpy( &buf[i+=4], &BITS_PER_SAMPLE,2);
   memcpy( &buf[i+=2], &BITS_PER_SAMPLE_FILL,2);
   }

memcpy( &buf[i+=2], &COMPRESSION_TAG,2);        /* IDF ENTRY 5*/
memcpy( &buf[i+=2], &COMPRESSION_TYPE,2);   
memcpy( &buf[i+=2], &COMPRESSION_COUNT,4);  
memcpy( &buf[i+=4], &COMPRESSION,2);   
memcpy( &buf[i+=2], &COMPRESSION_FILL,2);   

memcpy( &buf[i+=2], &PHOTOINTERP_TAG,2);         /* IDF ENTRY 6*/
memcpy( &buf[i+=2], &PHOTOINTERP_TYPE,2);  
memcpy( &buf[i+=2], &PHOTOINTERP_COUNT,4);   
memcpy( &buf[i+=4], &PHOTOINTERP,2);  
memcpy( &buf[i+=2], &PHOTOINTERP_FILL,2);   

memcpy( &buf[i+=2], &STRIP_OFFSETS_TAG,2);        /* IDF ENTRY 7*/
memcpy( &buf[i+=2], &STRIP_OFFSETS_TYPE,2);   
memcpy( &buf[i+=2], &STRIP_OFFSETS_COUNT,4);   
memcpy( &buf[i+=4], &STRIP_OFFSETS,4);  
L=4;

if ( NB>1 ) 
  {
  memcpy( &buf[i+=4], &SAMP_PIX_TAG,2);             /* IDF ENTRY 8*/
  memcpy( &buf[i+=2], &SAMP_PIX_TYPE,2);   
  memcpy( &buf[i+=2], &SAMP_PIX_COUNT,4);    
  memcpy( &buf[i+=4], &SAMP_PIX,2);   
  memcpy( &buf[i+=2], &SAMP_PIX_FILL,2);   
  L=2;
  }

memcpy( &buf[i+=L], &ROWS_PER_STRIP_TAG,2);       /* IDF ENTRY 9*/
memcpy( &buf[i+=2], &ROWS_PER_STRIP_TYPE,2);  
memcpy( &buf[i+=2], &ROWS_PER_STRIP_COUNT,4);   
memcpy( &buf[i+=4], &ROWS_PER_STRIP,4); 

memcpy( &buf[i+=4], &STRIP_BYTE_COUNT_TAG,2);     /* IDF ENTRY 10*/
memcpy( &buf[i+=2], &STRIP_BYTE_COUNT_TYPE,2);   
memcpy( &buf[i+=2], &STRIP_BYTE_COUNT_COUNT,4);  
memcpy( &buf[i+=4], &STRIP_BYTE_COUNT,4);   
 
memcpy( &buf[i+=4], &X_RESOLUTION_TAG,2);         /* IDF ENTRY 11*/
memcpy( &buf[i+=2], &X_RESOLUTION_TYPE,2);  
memcpy( &buf[i+=2], &X_RESOLUTION_COUNT,4);   
memcpy( &buf[i+=4], &X_RESOLUTION_ADDR,4);   
 
memcpy( &buf[i+=4], &Y_RESOLUTION_TAG,2);         /* IDF ENTRY 12*/
memcpy( &buf[i+=2], &Y_RESOLUTION_TYPE,2);   
memcpy( &buf[i+=2], &Y_RESOLUTION_COUNT,4);  
memcpy( &buf[i+=4], &Y_RESOLUTION_ADDR,4);   

memcpy( &buf[i+=4], &RESOLUTION_UNIT_TAG,2);      /* IDF ENTRY 13*/
memcpy( &buf[i+=2], &RESOLUTION_UNIT_TYPE,2);   
memcpy( &buf[i+=2], &RESOLUTION_UNIT_COUNT,4);   
memcpy( &buf[i+=4], &RESOLUTION_UNIT,2);   
memcpy( &buf[i+=2], &RESOLUTION_UNIT_FILL,2);  

if ( NB==1 ) i+=12;
memcpy( &buf[i+=2], &NEXT_IDF_OFFSET,4);   
memcpy( &buf[i+=4], &X_RESOLUTION[0],8);   
memcpy( &buf[i+=8], &Y_RESOLUTION[0],8);  

if ( color )memcpy( &buf[i+=8], &BITS_PER_SAMPLE_COLOR[0],6);  
return true;
}
bool getTiff(imgfile_t* iFile, const char* buf)
{
int i,j;
int first_offset;
int bits_per;
short int forty_two,num_direct,num_bands,bps[3];
memcpy(&forty_two,&buf[2],2);
memcpy(&first_offset,&buf[4],4);
if ( forty_two == 42 )
  {
  iFile->tiff_flag= 1;
  num_bands= 1;
  for (j=0; j<10; j++)
  {
  memcpy(&num_direct,&buf[first_offset],2);
  }
  for ( i=0; i<num_direct; i++)
    {
    unsigned    short int idf[6];
    memcpy(idf,&buf[first_offset+2+i*12],12);
    if ( idf[0]==256 )memcpy( &iFile->np,  &idf[4], 4);
    if ( idf[0]==257 )memcpy( &iFile->nl,  &idf[4], 4);
    if ( idf[0]==277 )memcpy( &num_bands,  &idf[4], 2);
    if ( idf[0]==258 )
      {
      bits_per= idf[4];
      if ( bits_per!=8 && bits_per!=16 )
        {
        memcpy(&bps,&buf[idf[4]+0],6);
        bits_per= bps[0];
        }
      }
    }
   iFile->ps= bits_per/8;
   iFile->nb= num_bands;
   iFile->sz= iFile->np*iFile->nl*num_bands*iFile->ps+512;
  }
 else
   { 
   iFile->tiff_flag= 0;
   sprintf(msgbuf,"*** no tiff header on file \"%s\"",iFile->fname); 
   RETURN_ERROR(msgbuf,"getTiff", false);
   }
 return true;
}
bool img_openw( imgfile_t* iFile, const char* fname, 
               const int np, const int nl, const int nb, const int ps )
{
  char hdr[512];
  int io_size;
  iFile->fp= open( fname, (O_CREAT|O_RDWR), (S_IRUSR|S_IWUSR|S_IRGRP|S_IWGRP) );
  if ( iFile->fp==-1 )
    { 
    sprintf(msgbuf,"*** error opening file \"%s\"",fname); 
    RETURN_ERROR(msgbuf,"img_openw", false);
    }
  iFile->np= np;
  iFile->nl= nl;
  iFile->nb= nb;
  iFile->ps= ps;
  tiffHdr( hdr,iFile );          /* generate tiff hdr*/
  iFile->tiff_flag= 1;
  io_size= write(iFile->fp,hdr,512);
  if ( io_size==-1 )
    {
    sprintf(msgbuf,"*** error, can not write to file \"%s\" ***",fname);
    RETURN_ERROR(msgbuf,"img_openw", false);
    }
 iFile->fname= (char*)malloc( strlen(fname)+1 );
 strcpy(iFile->fname,fname);
 return true;
}
bool img_openr( imgfile_t* iFile, const char* fname )
{
  int n_read;
  char buf[512];
 
  iFile->fp= open( fname, O_RDONLY );
  if ( iFile->fp==-1 )
    {
    sprintf(msgbuf,"*** error opening file \"%s\"",fname);
    RETURN_ERROR(msgbuf,"img_openr", false);
    }
  n_read= read( iFile->fp, buf, 512 );      /* read hdr*/
  if ( n_read==-1 )
    { 
    sprintf(msgbuf,"*** error reading header from file \"%s\"",fname);
    RETURN_ERROR(msgbuf,"img_openr", false);
    }
  getTiff( iFile, buf );
 iFile->fname= (char*)malloc( strlen(fname)+1 );
 strcpy(iFile->fname,fname);

  return true;
}
bool img_close( imgfile_t* iFile )
{
  free( iFile->fname );
  close( iFile->fp );
  return true;
}
bool img_close_rm( imgfile_t* iFile )
{
  remove( iFile->fname );
  img_close( iFile );
  return true;
}
bool put_line(imgfile_t* iFile, char* buf, int line )
{
int stat=0;
int np= iFile->np;
int ps= iFile->ps;
int nb= iFile->nb;
stat=lseek(iFile->fp, 512*iFile->tiff_flag + line*np*ps*nb, SEEK_SET );
if ( stat==-1 )
  {
  sprintf(msgbuf,"error seeking line number %d to file %s",line,iFile->fname);
  RETURN_ERROR(msgbuf,"put_line", false);
  }
stat= write(iFile->fp,(char*)buf, np*ps*nb);
if ( stat==-1 )
  {
  sprintf(msgbuf,"error writing line number %d to file %s",line,iFile->fname);
  RETURN_ERROR(msgbuf,"put_line", false);
  }
return true;
}

bool get_line(imgfile_t* iFile, unsigned char* buf, int line )
{
int stat=0;
int np= iFile->np;
int ps= iFile->ps;
int nb= iFile->nb;
stat=lseek(iFile->fp, 512*iFile->tiff_flag + line*np*ps*nb, SEEK_SET );
if ( stat==-1 )
 {
 sprintf(msgbuf,"error seeking line number %d to file %s",line,iFile->fname);
 RETURN_ERROR(msgbuf,"get_line", false);
 }
stat= read(iFile->fp,(unsigned char*)buf, np*ps*nb);
if ( stat==-1 )
 {
 sprintf(msgbuf,"error reading line number %d \n to file %s",
         line,iFile->fname);
 RETURN_ERROR(msgbuf,"get_line", false);
 }

return true;
}
int big_endian(){long int i=1,j=0; char b=1; memcpy(&j,&b,1); return i!=j;}

