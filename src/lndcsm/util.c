#include "util.h"

void scale_buf(unsigned char* out, unsigned short int* in, int n, double scale, double rf1)
{
int i;
for (i=0; i<n; i++)out[i]= bounded( nint( ((float)in[i]-rf1)*scale ),0,255 );
}
char* format_deg2dms(const double deg, char* dms, const char delim_in)
{
  double fdeg, fmin, fsec, frac;
  int    ideg, imin, isec, irac, sign;
  char* delim=":";
  memcpy(delim,&delim_in,1);
  sign= deg<0 ? -1 : 1;
  fdeg= deg*(double)sign;
  ideg= (int)fdeg;
  fmin= ( fdeg-(double)ideg )*60.0;
  imin= (int)fmin;
  fsec= ( fmin-(double)imin )*60.0;
  isec= (int)fsec;
  frac= ( fsec-(double)isec );
  irac= (int)(frac*10000.0);
  sprintf(dms,"%4d%s%2.2d%s%2.2d.%4.4d",
          sign*ideg,delim,imin,delim,isec,irac);
  return dms;
}
char* deg2dms(const double deg, char* dms){return format_deg2dms(deg,dms,':');}
char* deg2dms0(const double deg, char* dms){return format_deg2dms(deg,dms,'0');}
double dms2deg(const char* dms)
{
char temp[1024];
const double sixty=60.0;
int start=0, itok=0;
double answer= 0.0,q=1.0,sign=1.0,v,av;
char* token= (char*)dms;
while ( token!=(char*)NULL && itok<3 )
  {
  strcpy(temp, &dms[start]);
  token= strtok(temp," :|");
  start+= strlen(temp);
  v= atof( token );
  if ( itok==0 && v<0 )sign=-1.0;
  av= FABS(v);
  answer+= (av*q) ;
  q/= sixty;
  itok++;
  }
return sign*answer;
}

void compute_std(vbuf_t* vb, const int num, float* mean_std)
{
int i;
double data;
double sum=0.0;
double sumsq=0.0;
for (i=0; i<num; i++)
  {
  data= (double)virget( vb,i );
  sum    += data;
  sumsq  += data*data;
  }
mean_std[0]= (float)(sum/(double)num);
mean_std[1]= (float)(sqrt((sumsq- sum*sum/(double)num)/(double)(num-1)));
}


Pbuf_t *InitPr(FILE* mfout, char* buf, int wid, const char edge)
{
  Pbuf_t *xthis;
  xthis = (Pbuf_t *)malloc(sizeof(Pbuf_t));
  xthis->ofile = mfout;
  xthis->wid=    wid;
  xthis->buffer= buf;
  xthis->edge= (char)edge;
  return xthis;
}

void pr(Pbuf_t* pb){printit(pb->ofile, pb->buffer, pb->wid, pb->edge);}
void cen(Pbuf_t* pb){center(pb->ofile, pb->buffer, pb->wid, pb->edge);}
void prb(Pbuf_t* pb,char* buf){printit(pb->ofile, buf, pb->wid, pb->edge);}
void cenb(Pbuf_t* pb,char* buf){center(pb->ofile, buf, pb->wid, pb->edge);}

void printit(FILE* mfout, char* buf,int wid, char edge)
{
int s,i,same=1;
char st[2]= " ", e[2], end[3], ebl[4];
e[0]=edge; e[1]=(char)0;
end[0]=edge; end[1]='\n'; end[2]=(char)0;
ebl[0]=edge; ebl[1]=' '; ebl[2]=' ';  ebl[3]=(char)0;
for(i=1; i<strlen(buf); i++)if ( buf[i]!=buf[0] )same=0;
if ( same ) 
  {
  st[0]=buf[0];
  fprintf(mfout,e);
  for(i=0; i<wid-2; i++)fprintf(mfout,st);
  fprintf(mfout,end);
  }
else
  {
  fprintf(mfout,ebl);
  fprintf(mfout,"%s",buf);
  s= wid - ( strlen(buf)+4 );
  if ( s<0 ) s=0;
  for (i=0; i<s; i++)fprintf(mfout," ");
  fprintf(mfout,end);
  }
fflush(mfout);
}
void center(FILE* mfout, char* buf,int wid, char edge)
{
int s,s1,s2,i,same=1;
char st[2]= " ", e[2], end[3], ebl[3];
e[0]=edge; e[1]=(char)0;
end[0]=edge; end[1]='\n'; end[2]=(char)0;
ebl[0]=edge; ebl[1]=' ';  ebl[2]=(char)0;
for(i=1; i<strlen(buf); i++)if ( buf[i]!=buf[0] )same=0;
if ( same && buf[0]!=' ' ) 
  {
  st[0]=buf[0];
  fprintf(mfout,e);
  for(i=0; i<wid-2; i++)fprintf(mfout,st);
  fprintf(mfout,end);
  }
else
  {
  fprintf(mfout,ebl);
  s= wid - ( strlen(buf)+3 );
  if ( s<0 ) s=0;
  s1= s/2;
  s2= s-s1;
  for (i=0; i<s1; i++)fprintf(mfout," ");
  fprintf(mfout,"%s",buf);
  for (i=0; i<s2; i++)fprintf(mfout," ");
  fprintf(mfout,end);
  }
fflush(mfout);
}
bool sieve( int            xsize
          , int            ysize
          , unsigned char* sive_buf
          , unsigned int*  sive_flag
          , unsigned char* sive_cflag
          , int            sive_size
          , int            odometer_flag
          )
{  

  int 
      max_particle_size
    , i
    , maxint=2147483647
    , jtime=0
    , value
    , ix
    , iy
    , ib
    , x
    , y
    , fv
    , ntouch
    , num_same
    , nr=0
    , nd=0
    ;
  const int fillv= FILLV
          , landv= LANDV
          , cloud= CLOUD
          , snowv= SNOWV
    ;
  unsigned char 
                mask_lut[4]
              , cflag
    ;
  int old_val;
  bool candidate_flag= true;
  /***********************************************************************/
  /***************************** begin sieve  ****************************/ 
  /***********************************************************************/
  max_particle_size  = 5;
  max_particle_size  = sive_size;
  mask_lut[DM_LANDV+1]= landv;
  mask_lut[DM_CLOUD+1]= cloud;
  mask_lut[DM_SNOWV+1]= snowv;

  for (i=0; i<256;i++)index_lut[i]= (int)DM_NULL;
  index_lut[landv]= (int)DM_LANDV;
  index_lut[cloud]= (int)DM_CLOUD;
  index_lut[snowv]= (int)DM_SNOWV;
  index_lut[fillv]= (int)DM_FILLV;

  incr_nr[DM_LANDV]= cloud;
  incr_nr[DM_CLOUD]= landv;
  incr_nr[DM_SNOWV]= cloud;

  incr_nd[DM_LANDV]= snowv;
  incr_nd[DM_CLOUD]= snowv;
  incr_nd[DM_SNOWV]= landv;

  fv=1;
  nr= nd= 0;

  memset(sive_cflag, 0,xsize*ysize*sizeof(unsigned char));
  memset(sive_flag,  0,xsize*ysize*sizeof(unsigned int));

  for (iy=0; iy<ysize; iy++ )
    {
    if ( ( iy==0 || iy%25==0 || iy==(ysize-1) ) && odometer_flag )
      {  
      printf("line: %4d \r",iy);  fflush(stdout); 
      }
    for (ix=0; ix<xsize; ix++ )
      { 
      candidate_flag= true;
      if ( max_particle_size<9 )
        {
        num_same=0;
        for (y= iy-1; y<iy+2; y++)
          for (x= ix-1; x<ix+2; x++)
            if ( (ix!=x || iy!=y) && x>=0 && y>=0 && x<xsize && y<ysize )
              if ( sive_buf[xsize*iy+ix]==sive_buf[xsize*y+x] )num_same++;
        if ( num_same>max_particle_size )candidate_flag= false;
        }
      value= (int)sive_buf[xsize*iy+ix];
      if( candidate_flag )
        {
        ntouch=0;
        sive_flag[xsize*iy+ix]= fv;
        nr=nd=0;

        if ( !touching( sive_buf
                      , xsize
                      , ysize
                      , ix
                      , iy 
                      , ib
                      , value
                      , &ntouch
                      , max_particle_size-1          
                      , sive_flag
                      , fv 
                      , &nr
                      , &nd
                      ) 
           ) 
          {
          switch (value) 
            {
            case LANDV: cflag= nr>nd ? DM_CLOUD+1 : DM_SNOWV+1; break;
            case CLOUD: cflag= nr>nd ? DM_LANDV+1 : DM_SNOWV+1; break;
            case SNOWV: cflag= nr>nd ? DM_CLOUD+1 : DM_LANDV+1; break;
            default:    cflag= 0;
            }
          sive_cflag[xsize*iy+ix]= cflag;
          }
        fv= fv<(maxint-1) ? fv+1 : 1;
        if ( fv==1 )memset(sive_flag,0,xsize*ysize*sizeof(int));
        if ( fv==1 )printf("*** mem_setting time=%d ***\n",jtime++);
        }
      }
    }
  for (iy=0; iy<ysize; iy++ )
    for (ix=0; ix<xsize; ix++ )
      {
      old_val= sive_buf[xsize*iy+ix];
      if ( sive_cflag[xsize*iy+ix]>0 )
        sive_buf[xsize*iy+ix]=mask_lut[ (int)sive_cflag[xsize*iy+ix] ];
      }
  if ( odometer_flag )printf("\n");
  return true;
}
int touching( unsigned char* ibuf
            , int xsize
            , int ysize
            , int ix
            , int iy
            , int ib
            , int ival
            , int* ntouch
            , int maxp
            , unsigned int* flag
            , int fv
            , int* nr
            , int* nd
            )
{
  int x,y,value;

  for (y= iy-1; y<iy+2; y++)
    for (x= ix-1; x<ix+2; x++)
      {
      if ( (ix!=x || iy!=y) &&(ix==x ||  iy==y) &&
            x>=0 && y>=0 && x<xsize && y<ysize )
        {
        value= (int)ibuf[xsize*y+x];
        if ( value!=ival && flag[xsize*y+x]!=fv )
	  {
          if ( value == incr_nr[ index_lut[ival] ] )(*nr)++;
          if ( value == incr_nd[ index_lut[ival] ] )(*nd)++;
	  }
        if ( value==ival && flag[xsize*y+x]!=fv )
	  {
          if ( ++(*ntouch)>maxp )return 1;
          flag[xsize*y+x]= fv;
          if (touching(ibuf,xsize,ysize,x,y,ib,ival,ntouch,maxp,flag,fv,nr,nd) 
             )return 1;
	  }
	}
      }
  return 0;
}
