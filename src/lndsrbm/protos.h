int decode2scans(uchar *buffer,short *ch[5][2],float *flat[],float *flon[],
                 float *ts[]);
void intergeom(float *lat_51,float *lon_51,float *ts_51,float *lat,float *lon,
               float *ts,short *type);
int getscantime(uchar *buffer,short *year,short *day,float *gmt);
int IsaTape(char *input);
int initializepar(char *file,t_project *proj);
int rdparamsline(FILE *fd,int *par,char *value);
int rdtbmheader(void);
int rddataheader(void);
int rdrecord(void);
short getsatid(unsigned char c);
int openorbitfile(FILE *fdlist,FILE **fdinput);
int ProjectInit(t_project *proj);
int ProjSetup(t_project *proj);
int ProjFOR(t_project *project,int *y,int *x,double lat,double lon);
int ProjINV(t_project *project,int i,int j,double *lat,double *lon);
int getavhrrparams(short scid);
int getsatname(short scid,char *satname);
int getprojname(long proj,char *projname);
int gettypename(short type,char *typename);
int prescan(t_project *proj,double taos,
            int *xmin,int *xmax,int *ymin,int *ymax,int *scanmin,int *scanmax,
            int *nbportions,float *scanAN);
int domap(char *file,short color,short *layer,t_project *project);
int doframe(short color,short *layer,t_project *project);
int dogrid(short color,short *layer,t_project *project);

