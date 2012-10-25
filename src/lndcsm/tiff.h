#ifndef TIFF_HPP
#define TIFF_HPP
#include "error.h"
#include "bool.h"
typedef struct 
{
int fp;
int ps;
int np;
int nl;
int nb;
int sz;
int tiff_flag;
char* fname;
} imgfile_t;
bool getTiff(imgfile_t* iFile, const char* buf);
bool tiffHdr (char* hdrdata, imgfile_t* iFile) ;
bool img_openr( imgfile_t* iFile, const char* fname );
bool img_openw( imgfile_t* iFile, const char* fname, 
               const int np, const int nl, const int nb, const int ps );
bool put_line(imgfile_t* iFile, char* buf, int line );
bool get_line(imgfile_t* iFile, unsigned char* buf, int line );
bool img_close( imgfile_t* iFile );
bool img_close_rm( imgfile_t* iFile );
int big_endian();
#endif /* TIFF_HPP */
