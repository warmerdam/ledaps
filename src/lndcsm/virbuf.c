#include <stdio.h>       /* for true, false */
#include <string.h>     /* for strlen       */
#include <math.h>
#include <stdlib.h>
#include <ctype.h> 
#include <time.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "virbuf.h"
char msgbuf[1024];
bool virinit( vbuf_t* virbuf, char* fname, const int blocking )
{
virbuf->file= open(fname,(O_CREAT|O_RDWR),(S_IRUSR|S_IWUSR|S_IRGRP|S_IWGRP));
if ( virbuf->file==-1 )
  { 
  sprintf(msgbuf,"*** error opening file \"%s\"",fname); 
  RETURN_ERROR(msgbuf,"virnit", false);
  }
virbuf->size=0;
virbuf->current_line=0;
virbuf->blocking= blocking;
virbuf->buffer = (float*)malloc ( blocking* sizeof(float) );
virbuf->fname = (char*)malloc ( (strlen(fname)+1)* sizeof(char) );
strcpy(virbuf->fname,fname);
return true;
}
bool virclose( vbuf_t* virbuf )
{
int status;
virbuf->size=0;
virbuf->current_line=0;
status= close( virbuf->file );
free( virbuf->buffer );
remove( virbuf->fname );
free( virbuf->fname );
return true;
}
bool vir_reinit( vbuf_t* virbuf )
{
int seek_status, seek_addr;
virbuf->size=0;
virbuf->current_line=0;
seek_addr= virbuf->blocking * virbuf->current_line * sizeof(float) ;
seek_status= lseek(virbuf->file, seek_addr, SEEK_SET );
return true;
}
bool virput( vbuf_t* virbuf, float value )
{
virbuf->buffer[ virbuf->size ]= value;
virbuf->size++;
if ( virbuf->size == virbuf->blocking )
  {
  write(virbuf->file, (char*)virbuf->buffer, virbuf->blocking*sizeof(float) );
  if ( virbuf->file==-1 )
    { 
    sprintf(msgbuf,"*** error writing to file \"%s\"",virbuf->fname); 
    RETURN_ERROR(msgbuf,"virput", false);
    }
  virbuf->current_line++;
  virbuf->size=0;
  }
return true;
}
bool virflush( vbuf_t* virbuf )
{
write(virbuf->file, (char*)virbuf->buffer,virbuf->blocking*sizeof(float) ); 
return true;
}
float virget( vbuf_t* virbuf, int index )
{
int line= index/virbuf->blocking;
if ( line != virbuf->current_line )
  {
  int seek_status,seek_addr,n_read;
  virbuf->current_line= line;
  seek_addr= virbuf->blocking * virbuf->current_line * sizeof(float) ;
  seek_status= lseek(virbuf->file, seek_addr, SEEK_SET );
  n_read=read(virbuf->file,(char*)virbuf->buffer,virbuf->blocking*sizeof(float));
  if ( virbuf->file==-1 )
    { 
    sprintf(msgbuf,"*** error reading  file \"%s\"",virbuf->fname); 
    RETURN_ERROR(msgbuf,"virget", false);
    }
  }
return virbuf->buffer[ index%virbuf->blocking ];
}
