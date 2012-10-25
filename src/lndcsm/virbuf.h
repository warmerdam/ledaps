#ifndef VIRB_HPP
#define VIRB_HPP
#include "bool.h"
#include "error.h"
typedef struct 
{
int size;
int all;
int blocking;
int current_line;
int file;
float* buffer;
float* data_buffer;
char* fname;
} vbuf_t ;

bool virinit( vbuf_t *virbuf, char* fname, const int blocking );
bool virclose( vbuf_t* virbuf );
bool vir_reinit( vbuf_t* virbuf );
bool virput( vbuf_t* virbuf, float value );
bool virflush( vbuf_t* virbuf );
float virget( vbuf_t* virbuf, int index );

#endif /* VIRB_HPP                                                   */
