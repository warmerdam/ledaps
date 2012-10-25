#ifndef UTIL_HPP
#define UTIL_HPP
/****************************************************************************/
#include <stdio.h>      /* for true, false                                  */
#include <string.h>     /* for strlen                                       */
#include <math.h>
#include <stdlib.h>
#include <ctype.h> 
#include <string.h>
#include <time.h>
#define max(A,B) (A>B?A:B)
#define min(A,B) (A>B?B:A)
#define nint(A)(A<0?(int)(A-0.5):(int)(A+0.5))
#define bounded(A,B,C)(A>B?(A<C?A:C):B)
#define FABS(A)(A>0?A:-1.0*A)
#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "csm.h"
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
void scale_buf(unsigned char* out, unsigned short int* in, int n, double scale, double rf1);
char* format_deg2dms(const double deg, char* dms, const char delim_in);
char* deg2dms(const double deg, char* dms);
char* deg2dms0(const double deg, char* dms);
double dms2deg(const char* dms);
void compute_std(vbuf_t* vb, const int num, float* mean_std);
void printit(FILE* mfout, char* buf,int wid, char edge);
void center(FILE* mfout, char* buf,int wid, char edge);
typedef struct {
  char *buffer;             /* print buffer */
  int wid;                  /* print buffer width */
  FILE* ofile;              /* output file */
  char edge;                /* edge character */
} Pbuf_t;
Pbuf_t *InitPr(FILE* mfout, char* buf, int wid, const char edge);
void pr(Pbuf_t* pb);
void cen(Pbuf_t* pb);
void prb(Pbuf_t* pb,char* buf);
void cenb(Pbuf_t* pb,char* buf);
                                                                       /*****************************************************************************/
#define CLOUD_LAND   (0)
#define CLOUD_CLOUD  (1)    /* bit 0 */
#define CLOUD_SHADOW (2)    /* bit 1 */
#define CLOUD_WATER  (4)    /* bit 2 */
#define CLOUD_SNOW   (8)    /* bit 3 */
#define CLOUD_FILL   (255)
typedef enum {
  DM_NULL = -1,
  DM_LANDV = 0,
  DM_CLOUD, 
  DM_SNOWV,
  DM_FILLV
} indx_t;
#define FILLV (255)
#define LANDV (0)
#define CLOUD (1)
#define SNOWV (8)
int incr_nr[3];
int incr_nd[3];
int index_lut[256];
/*****************************************************************************/
bool sieve( int            xsize
          , int            ysize
          , unsigned char* sive_buf
          , unsigned int*  sive_flag
          , unsigned char* sive_cflag
          , int            sive_size
          , int odometer_flag
		);
int touching( unsigned char* img_buf
            , int nps
            , int nls
            , int ix
            , int iy
            , int ib
            , int val
            , int* ntouch
            , int maxp
            , unsigned int* flag
            , int fv
            , int* nr
            , int* nd
            );
/*****************************************************************************/
#endif /* UTIL_HPP */
