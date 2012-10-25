#include <stdio.h> 
#include <stdlib.h> 
#include "mfhdf.h"
#define D2R     1.745329251994328e-2

/*
 * cc -DINV -O0 -o xy_to_geo LS_geoloc.o LS_geoloc_driver.c -I$HDFINC -L$HDFLIB -lmfhdf -ldf -ljpeg -lz -lm
      OR 
 * cc       -O0 -o geo_to_xy LS_geoloc.o LS_geoloc_driver.c -I$HDFINC -L$HDFLIB -lmfhdf -ldf -ljpeg -lz -lm
 */

/* Had been using these values from a sample LS file: */
/*  *zonecode = 18;*/
/*  *sphercode = 8;*/
/*  *orientationangle = 0.0; */
/*  *pixelsize = 28.5; */
/*  *upperleftx = 395038.500000;*/
/*  *upperlefty = 4728748.500000;*/

int main(int argc, char **argv)
{

float coordinates[8];
double parm[13];
double radius, lat, lon, dl, ds;
double corner[2];
int i, l, s, ret;
int zonecode, sphercode, rows, cols;
float orientationangle, pixelsize, upperleftx, upperlefty;
double arg2, arg3;

int LSsphdz(float coordinates[8], double *parm, double *radius, double corner[2]);
int LSutminv(double s, double l, double *lon, double *lat);
int LSutmfor(double *s, double *l, double lon, double lat);
int get_data(char *filename, int *zonecode, int *sphercode, float *orientationangle, float *pixelsize, 
             float *upperleftx, float *upperlefty, int *rows, int *cols);

if (argc < 4) {
#ifdef INV
   printf("usage: %s <file> <sample> <line>\n", argv[0]);
#else
   printf("usage: %s <file> <longitude> <latitude>\n", argv[0]);
#endif
   printf("Jim Ray, SSAI, %s\n\n", __DATE__);
   exit(0);
   }

if ( (ret = get_data(argv[1], &zonecode, &sphercode, &orientationangle, &pixelsize, 
         &upperleftx, &upperlefty, &rows, &cols)) != 0) {
   printf("Error reading file %s, cannot continue\n", argv[1]);
   exit(1);
      } 

arg2 = atof(argv[2]);
arg3 = atof(argv[3]);

coordinates[4] = (double)zonecode;
coordinates[5] = (double)sphercode;
coordinates[6] = (double)orientationangle;
coordinates[7] = (double)pixelsize;
for(i=0;i<13;i++) parm[i] = 0.0;
corner[0] = (double)upperleftx;
corner[1] = (double)upperlefty;

LSsphdz(coordinates, parm, &radius, corner);


#ifdef INV
ds = arg2;
dl = arg3;

if (ds > (double)cols) {
   printf("Sample argument (%s) exceeds number of columns in file (%d): will use %d\n", 
               argv[2], cols, cols);
   ds = (double)cols;
   }

if (dl > (double)rows) {
   printf("Sample argument (%s) exceeds number of rows in file (%d): will use %d\n", 
               argv[3], rows, rows);
   dl = (double)rows;   
   }

ret = LSutminv(ds, dl, &lon, &lat);
printf("line   %5.1f  samp   %5.1f  => long %f lat %f\n", dl, ds, lon, lat);
#else
lon = arg2;
lat = arg3;

/* We need sanity checks on these as well */

ret = LSutmfor(&ds, &dl, lon, lat);
printf("long %f lat %f => line   %f  samp   %f  \n",  lon, lat, dl, ds);
#endif
}

	     
/**********************************************************************************
 **********************************************************************************/

int get_data(char *filename, int *zonecode, int *sphercode, float *orientationangle, float *pixelsize, 
             float *upperleftx, float *upperlefty, int *rows, int *cols)
{
int32 i, j, n, sd, n_sets, n_gattr, count, structmetadata_exists, number_type;
float floatattr;
double doubleattr;
char projection[256];
char coordinate[256];
char attrib[256];
char attribu[256];
void metareader(int32 sd_id, char *type_of_meta, char *metastring, int32 *count, char *data);


*zonecode = *sphercode = *rows = *cols = -1;
*orientationangle = *pixelsize = -999.0;

/* Make sure input file is HDF, is 
 * readable, has SDSs, et cetera
 ***********************************/
if ((sd = SDstart(filename, DFACC_RDONLY)) == -1) {
     printf("Error: file '%s' can't be opened with SDstart(): cannot continue...\n", filename);
     exit(-2);
    }
       
n_sets = 0;
SDfileinfo(sd, &n_sets, &n_gattr);
if (n_sets == 0) {
    printf("Error: file %s doesn't seem to have any SDSs: cannot continue...\n", filename);
    SDend(sd);
    exit(-4);
   } 
if (n_gattr == 0) {
    printf("Error: file %s doesn't seem to have any global attributes: cannot continue...\n", filename);
    SDend(sd);
    exit(-4);
   } 
   
structmetadata_exists = 0;
for (j=0;j<n_gattr;j++) {
     SDattrinfo(sd, j, attrib, &number_type, &count);
     n = strlen(attrib);
     for (i=0;i<n;i++) attribu[i] = toupper(attrib[i]);
     
     if (strstr(attribu, "ORIENTATIONANGLE")) {
        SDreadattr(sd, j, &doubleattr);
	*orientationangle = (float)doubleattr;
	}
     if (strstr(attribu, "PIXELSIZE")) {
        SDreadattr(sd, j, &doubleattr);
	*pixelsize = (float)doubleattr;
	}
	
     if (strstr(attribu, "STRUCTMETADATA")) {
        structmetadata_exists = 1;
	}
    }

if (structmetadata_exists == 0) {
    printf("Error: file %s doesn't seem to have any StructMetadata: cannot continue...\n", filename);
    SDend(sd);
    exit(-4);
   } 

metareader(sd, "STRUCTMETADATA\0", "Projection\0", &count, projection);

if (  !strstr(projection, "GCTP_UTM\0") ) {
    printf("Error: file %s has projection %s, not UTM: cannot continue...\n", filename, projection);
    SDend(sd);
    exit(-4);
   } 

metareader(sd, "STRUCTMETADATA", "UpperLeftPointMtrs\0", &count, coordinate);
parse_comma_sep(coordinate, upperleftx, upperlefty);
    
coordinate[0] = '\0';    
metareader(sd, "STRUCTMETADATA", "XDim\0", &count, coordinate);
*cols = atoi(coordinate);

coordinate[0] = '\0';    
metareader(sd, "STRUCTMETADATA", "YDim\0", &count, coordinate);
*rows = atoi(coordinate);

coordinate[0] = '\0';    
metareader(sd, "STRUCTMETADATA", "ZoneCode\0", &count, coordinate);
*zonecode = atoi(coordinate);

coordinate[0] = '\0';    
metareader(sd, "STRUCTMETADATA", "SphereCode\0", &count, coordinate);
*sphercode = atoi(coordinate);

SDend(sd);

if ( *zonecode == -1) {
    printf("Error reading zone code, cannot continue...\n");
    exit(-5);
  } 
if ( *sphercode == -1) {
    printf("Error reading sphere code, cannot continue...\n");
    exit(-5);
  }
if ( *rows == -1) {
    printf("Error reading number of rows, cannot continue...\n");
    exit(-5);
   }
if ( *cols == -1) {
    printf("Error reading number of columns, cannot continue...\n");
    exit(-5);
  }
if ( *orientationangle < -998.0 )  {
    printf("Error reading orientation angle, cannot continue...\n");
    exit(-5);
  }
if ( *pixelsize < 0.0 )  {
    printf("Error reading pixel size, cannot continue...\n");
    exit(-5);
  }

/*printf("As read from file %s (%d rows  by %d columns)\n", filename, *rows, *cols);
printf("Zonecode %d: sphere code %d; Orientation angle %f; Pixelsize %f; UL: %f, %f\n", 
        *zonecode, *sphercode, *orientationangle, *pixelsize, *upperleftx, *upperlefty); */
	
return(0);
}
	     
/**********************************************************************************
 **********************************************************************************/

int parse_comma_sep(char *argument, float *left, float *right)
{  /* designed to parse StructMetadata's UpperLeft/LowerRight metadata in the GCTP_UTM projection */
int i, j, n, n_non, where;
double tmpdb;
char c, x, tmpstr[256];
char tmparg[256];
int in_paren;

n = strlen(argument);
if (n < 3) return(-1);

/* first strip off the parentheses */
j = in_paren = 0;
for(i=0;i<n;i++) {
   x = argument[i];
   if ((x == '(') && (in_paren == 0)) in_paren = 1; 
   if ((x == ')') && (in_paren == 1)) in_paren = 0; 
   if (x == ')') break;
   if (x == '(') continue;
   if (in_paren == 1) tmparg[j++] = x;   
   }
/*tmparg[j++] = ',';*/
tmparg[j] = '\0';
/*printf("%s\n", tmparg);*/

n_non = where = 0;
n = strlen(tmparg);
for(i=0;i<n;i++) {
   x = tmparg[i];
   if (( !isdigit(x)) && ( x != '.') && ( x != '-')) {
      n_non++;
      c = x;
      where = i;
      }
   }

if ((n_non != 1) || (c != ',') || (where ==  0)) return(-1);

for(i=0;i<where;i++) {
   x = tmparg[i];
   tmpstr[i] = x;
   }
tmpstr[i] = '\0';

tmpdb = atof(tmpstr);
*left = (float)tmpdb;

j=0;
i++;
for(;i<n;i++) {
   x = tmparg[i];
   tmpstr[j++] = x;
   }
tmpstr[j] = '\0';

tmpdb = atof(tmpstr);
*right = (float)tmpdb;

return(0);
}



/*****************************************************************************************************************
 **************  start of metareader() ***************************************************************************
 *****************************************************************************************************************/
 
void metareader(int32 sd_id, char *type_of_meta, char *metastring, int32 *count, char *data)
{
#define XMAXLENGTH 1000
int32 i, ii, j, n_attr1, n_sets1, count1, number_type1, n_val;
int start,within,wasjustobject;
int obj_offset[10],n_obj;
char attrib[XMAXLENGTH], attrib1[XMAXLENGTH];
char *charattr;
char line[XMAXLENGTH], objs[XMAXLENGTH];
char lhs[XMAXLENGTH], rhs[XMAXLENGTH];
char uppers[XMAXLENGTH];
void get_a_line(char *text, int lengthoftext, int *start, char line[XMAXLENGTH]);

SDfileinfo(sd_id, &n_sets1, &n_attr1);
for (j=0;j<n_attr1;j++) {
     SDattrinfo(sd_id, j, attrib1, &number_type1, &count1);
     attrib[0] = '\0';  
     for (i=0;i<strlen(attrib1);i++) attrib[i] = toupper(attrib1[i]);  
     attrib[i] = '\0';  
     start = 0;
     if (strstr(attrib,type_of_meta)) {
         if ((charattr = (char *)malloc((count1+1)*sizeof(char))) == NULL) {
	         printf("Out of memory, array 'charattr'\n");
		 return;
	                                                                   }
         SDreadattr(sd_id, j, charattr);
         count1 = strlen(charattr);
         if (!strcmp(type_of_meta, metastring)) {  /* signal for no need to parse metadata */
	    strcpy(data, charattr);
	    return;
	    }


	 n_obj=0;
	 wasjustobject=0;
	 objs[0] = '\0';
	 do { 
	      line[0] = '\0';
	      get_a_line(charattr,count1+1,&start,line);

	      /* Get rid of whitespace, if not inside "..." characters. */
	      within=0;
              for (i=0,ii=0;i<strlen(line);i++) {
	          if (line[i] == '"') {
		     if (!within) within=1;
		     else if (within) within=0;
		     }
	          if (within) line[ii++] = line[i]; 
		  else {
		     if (line[i] != ' ') line[ii++] = line[i];  
		     }
		  }
	      line[ii] = '\0';

	      /* Get left-hand-side, rhs of equation */
	      lhs[0] = rhs[0] = '\0';
	      if (strchr(line,'=')) {
                 for (i=0,ii=0;i<(strchr(line,'=')-line);i++) {
		        lhs[ii++]=line[i];
	             }
	         lhs[ii] = '\0';
		 i++;
                 for (ii=0;i<strlen(line)-1;i++,ii++) {
		        rhs[ii]=line[i];
	             }
	         rhs[ii] = '\0';
		 }
	
	      if (!strcmp(lhs,"OBJECT")) {
	           wasjustobject=1;
	           obj_offset[n_obj++]=strlen(objs);
	           strcat(objs,rhs);
	          }
	      if (!strcmp(lhs,"GROUP")) wasjustobject=0;
	      if ((!strcmp(lhs,"CLASS"))&&(wasjustobject)) {
		   for (i=0;i<strlen(rhs)-1;i++) rhs[i]=rhs[i+1];
	           rhs[i-1] = '\0';
	           strcat(objs,rhs);
	          }
	      if (!strcmp(lhs,"END_OBJECT")) {
	           if (n_obj > 0)
	              objs[obj_offset[--n_obj]] = '\0';
	          }
	      if (!strcmp(lhs,"NUM_VAL")) {
	           n_val=atoi(rhs);
	          }
	      
	      if (!strcmp(lhs,"VALUE")) {
		     /*printf("!%s! %d *%s*%s*\n", objs,n_obj,lhs,rhs);*/
	         if (!strcmp(objs, metastring)) {
	             strcpy(data,rhs);
		     *count=n_val;
		     /*printf("!%s! %d *%s*%s*\n", objs,n_obj,lhs,rhs);*/
	                }
	            }
	      if ((!strcmp(lhs,"\0"))&&(wasjustobject)) {
		     /*printf("!%s! %d *%s*%s*\n", objs,n_obj,lhs,rhs);*/
	         if (!strcmp(objs, metastring)) {
	             strcat(data,line);
		     *count+=n_val;
		     /*printf("!%s! %d *%s*%s*\n", objs,n_obj,lhs,rhs);*/
	                }
	            }
		    
	      /* special case for StructMetadata */
	      if  ( (strstr(type_of_meta, "STRUCTMETADATA")) && ( strstr(lhs, metastring))  ) strcat(data, rhs);
	      	    
	    } while ( line[0] != '\0') ;
	 
	 
	 free(charattr);  
                                      }


                        }  /* for (j=0;j<n_attr1;j++) */

}

/*****************************************************************************************************************
 **************  start of get_a_line() ***************************************************************************
 *****************************************************************************************************************/
 
void get_a_line(char *text, int lengthoftext, int *start, char line[XMAXLENGTH])
{
#define XMAXLENGTH 1000
int i=0;
int where;
int getout = 0;

where = *start;
line[i] = '\0';
if (where >= lengthoftext) return;
if (text[where] == '\0') return;
else {
   while (getout == 0) {  
	 if ((text[where] == '\0')||(text[where] == '\n')||(i>=XMAXLENGTH)||(where>=lengthoftext)) getout=1;
         line[i++] = text[where++];
      }
      
    *start = where;
    line[i] = '\0';
    return;
     }
}



