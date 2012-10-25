/**
 * !Description
 *  This program appends thermal brightness temperature to vnir reflectance file 
 *
 * !Input  
 *  surface reflectance file (lndsr*) in HDF format (from "lndsr")
 *  thermal band temperature file (lndth*) in HDF format (from "lndcal")
 *
 * !Output
 *  updated lndsr* file in HDF format
 *
 * !Developer: 
 *  Feng Gao (Feng.Gao@nasa.gov)
 *
 * ! History 
 *  first on 8/28/2006  - produce a separate combined HDF file 
 *  revised on 7/2/2010 - update original lndsr file to prepare input 
 *          for stand-alone Ledaps package that include SR based mask 
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "hdf.h"
#include "mfhdf.h"
#include "HdfEosDef.h"

#define MAX_NCOLS   10000 
#define Max_StrLen  1000
#define FAILURE     -1
#define SUCCESS     1

typedef struct {

  char  fileName[Max_StrLen];
  int32 nrows;           
  int32 ncols;
  float ulx;
  float uly;

  int32 SD_ID;            
  int32 band_id; 

  char    sdsName[Max_StrLen];
  char    longName[Max_StrLen];
  char    units[Max_StrLen];
  int16   range[2];
  int16   fillValue;
  float64 scale;

} GRID_SR;    

void usage(char *);
int  getMetaInfo(GRID_SR *);
int  openForWrite(GRID_SR *sr, GRID_SR *th);
int  cleanUp(GRID_SR *sr);

/*#define DEBUG*/

main(int argc, char *argv[])
{
  int i, irow, icol, ret;
  int32 start[2];
  int32 length[2];
  int16 *buf;
  GRID_SR sr, th;

  if(argc != 5) {
    usage(argv[0]);
    exit(1);
  }
  
  /* parse command line */
  for(i=1; i<argc; i++){
    if(strcmp(argv[i],"-sr")==0)
      strcpy(sr.fileName, argv[++i]);
    else
      if(strcmp(argv[i],"-th")==0)
	strcpy(th.fileName, argv[++i]);
      else{
	printf("\nWrong option:%s\n",argv[i]);
	usage(argv[0]);
	exit(1);
      }    
  }

  strcpy(sr.sdsName, "band1");
  strcpy(th.sdsName, "band6");

  /* get metadata from lndsr* file */
  if((ret=getMetaInfo(&sr))==FAILURE) {
    printf("Retrieve %s metadata error\n", sr.fileName);
    exit(1);
  }

  /* get metadata from lndth* file */
  if((ret=getMetaInfo(&th))==FAILURE) {
    printf("Retrieve %s metadata error\n", th.fileName);
    exit(1);
  }
  
  if(sr.nrows!=th.nrows || sr.ncols!=th.ncols || sr.ulx!=th.ulx || sr.uly!=th.uly)
    printf("Warning: lndsr and lndth images does not match!\n");

  /* allocate memory */
  buf = calloc(th.ncols, sizeof(int16));
  if(!buf) {
    printf("Memory allocation error\n");
    exit(1);
  }
  
  /* create output and write metadata */
  printf("Updating lndsr file ...\n");
  if((ret=openForWrite(&sr, &th))==FAILURE) {
    printf("Open for write error\n");
    exit(1);
  }

  if(ret == 0)
    printf("Appending band6 for line --- ");
  else
    printf("Updating band6 for line --- ");

  for(irow=0; irow<sr.nrows; irow++) {
    printf("%4d\b\b\b\b",irow);
   
    start[0] = irow;
    start[1] = 0;
    length[0] = 1;
    length[1] = th.ncols;

    if((SDreaddata(th.band_id, start, NULL, length, buf)) == FAILURE) {
      printf("Error: reading band6 error irow: %d\n", irow);
      exit(1);
    }
    
    if((SDwritedata(sr.band_id, start, NULL, length, buf)) == FAILURE) {
      printf("Error: writing to combined file error irow: %d\n", irow);
      exit(1);
    }
  } /* end of irow */

  /* clean up files */
  cleanUp(&sr);
  cleanUp(&th);
  free(buf);
  printf("\n");
  return EXIT_SUCCESS;      /* success */
}


/* display usage */
void usage(char *command)
{
  printf("Usage: %s [-sr][-th] \n", command);
  printf("   -sr   <lndsr_HDF_File>     surface reflectance file\n");
  printf("   -th   <lndth_HDF_File>     thermal band temperature file\n");
}


/* get metadata and open specific sds */
int getMetaInfo(GRID_SR *sr) {

  char name[100];
  int i, index, ret;

  char GD_gridlist[100];
  int32 gfid=0, ngrid=0, gid=0;
  int32 bufsize=100;
  float64 GD_upleft[2];
  float64 GD_lowright[2];

  int32 rank,data_type,attributes,dim_sizes[2];
  int32 att_id;

  /* open a hdf file */
  gfid = GDopen(sr->fileName, DFACC_READ);
  if(gfid==FAILURE){
      printf ("Not successful in retrieving grid file ID/open\n");
      return FAILURE;
  }

  /* find out about grid type */
  ngrid=GDinqgrid(sr->fileName, GD_gridlist,&bufsize);
  if(ngrid==FAILURE){
      printf ("Not successful in retrieving grid name list\n");
      return FAILURE;
  }

  /* attach grid */
  gid = GDattach(gfid, GD_gridlist);
  if(gid==FAILURE){
      printf ("Not successful in attaching grid.\n");
      return FAILURE;
  }

  /* get grid info */
  ret = GDgridinfo(gid, &sr->ncols, &sr->nrows, GD_upleft, GD_lowright);
  if(ret==FAILURE){
      printf ("Failed to read grid info.\n");
      return FAILURE;
  }
  sr->ulx = GD_upleft[0];
  sr->uly = GD_upleft[1];

  /* detach grid */
  ret = GDdetach(gid);
  if(ret==FAILURE){
      printf ("Failed to detach grid.\n");
      return FAILURE;
  }

  /* close for grid access */
  ret = GDclose(gfid);
  if(ret==FAILURE){
      printf ("GD-file close failed.\n");
      return FAILURE;
  }

  /* open hdf file and get sds_id from given sds_name */  
  if ((sr->SD_ID = SDstart(sr->fileName, DFACC_READ))<0) {
    printf("Can't open file %s\n",sr->fileName);
    return FAILURE;
  }

  /* open band1 */
  if ((index=SDnametoindex(sr->SD_ID, sr->sdsName))<0) {
    printf("Not successful in convert sdsName %s to index\n", sr->sdsName);
    return FAILURE;
  }
  sr->band_id = SDselect(sr->SD_ID, index);

  /* retrieve SR Fill Value */
  if ((att_id = SDfindattr(sr->band_id, "_FillValue")) == FAILURE) {
    printf ("Can't retrieve fill value from SR SDS attr.\n");
    return FAILURE;
  }
  if (SDreadattr(sr->band_id, att_id, &sr->fillValue) == FAILURE) {
    printf ("Can't retrieve fill value from SR SDS attr.\n");
    return FAILURE;
  }

  if ((att_id = SDfindattr(sr->band_id, "long_name")) == FAILURE) {
    printf ("Can't retrieve long name from SR SDS attr.\n");
    return FAILURE;
  }
  if (SDreadattr(sr->band_id, att_id, sr->longName) == FAILURE) {
    printf ("Can't retrieve long name from SR SDS attr.\n");
    return FAILURE;
  }

  if ((att_id = SDfindattr(sr->band_id, "units")) == FAILURE) {
    printf ("Can't retrieve units from SR SDS attr.\n");
    return FAILURE;
  }
  if (SDreadattr(sr->band_id, att_id, sr->units) == FAILURE) {
    printf ("Can't retrieve units from SR SDS attr.\n");
    return FAILURE;
  }

  if ((att_id = SDfindattr(sr->band_id, "valid_range")) == FAILURE) {
    printf ("Can't retrieve valid range from SR SDS attr.\n");
    return FAILURE;
  }
  if (SDreadattr(sr->band_id, att_id, sr->range) == FAILURE) {
    printf ("Can't retrieve valid range from SR SDS attr.\n");
    return FAILURE;
  }

  if ((att_id = SDfindattr(sr->band_id, "scale_factor")) == FAILURE) {
    printf ("Can't retrieve scale factor from SR SDS attr.\n");
    return FAILURE;
  }
  if (SDreadattr(sr->band_id, att_id, &sr->scale) == FAILURE) {
    printf ("Can't retrieve scale factor from SR SDS attr.\n");
    return FAILURE;
  }

  
#ifdef DEBUG
  printf("File: %s\n",sr->fileName);
  printf("nrows:%d, ncols:%d\t",sr->nrows, sr->ncols);
  printf("ulx:%8.1f, uly:%8.1f\t", sr->ulx, sr->uly);
  printf("fillV:%d\n", sr->fillValue); 
#endif

  cleanUp(sr);
  return SUCCESS;
}


/* create output HDF file and open sds_id for write */
int openForWrite(GRID_SR *sr, GRID_SR *th)
{
  char str[Max_StrLen];
  int ret, index, exist;
  int32 GDfid, GDid; 
  
  /* open hdf file and check if thermal band has already been appened */  
  if ((sr->SD_ID = SDstart(sr->fileName, DFACC_READ))<0) {
    printf("Can't open file %s",sr->fileName);
    return FAILURE;
  }
  
  /* check if band 6 already existed */
  if ((index=SDnametoindex(sr->SD_ID, th->sdsName))<0) {
    exist = 0;
  }
  else exist = 1;

  if((SDend(sr->SD_ID)) == FAILURE) {
    printf("SDend for %s error!\n", sr->fileName);
    return FAILURE;
  }

  /* create output use HDF-EOS functions */
  GDfid = GDopen(sr->fileName, DFACC_RDWR);
  if(GDfid == FAILURE){
    printf("Not successful in creating grid file %s", sr->fileName);
    return FAILURE;
  }
  
  GDid = GDattach(GDfid, "Grid");
  if(GDid == FAILURE) {
    printf("Not successful in attaching grid ID/ create for %s", sr->fileName);
    return FAILURE;
  }

  if(exist == 0) {
    ret = GDdeffield(GDid, th->sdsName, "YDim,XDim", DFNT_INT16, HDFE_NOMERGE);
    if(ret==FAILURE){
      printf ("Not successful in defining %s SDS", th->sdsName);
      return FAILURE;
    }
  }
  
  /* detach grid */
  ret = GDdetach(GDid);
  if(ret==FAILURE){
    printf ("Failed to detach grid.");
    return FAILURE;
  }

  /* close for grid access */
  ret = GDclose(GDfid);
  if(ret==FAILURE){
    printf ("GD-file close failed.");
    return FAILURE;
  }
 
  /* open hdf file and get sds_id from given sds_name and then write metedata */  
  if ((sr->SD_ID = SDstart(sr->fileName, DFACC_RDWR))<0) {
    printf("Can't open file %s",sr->fileName);
    return FAILURE;
  }
  
  if ((index = SDnametoindex(sr->SD_ID, th->sdsName))<0) {
    printf("Not successful in convert SR sdsName band6 to index");
    return FAILURE;
  }
  sr->band_id = SDselect(sr->SD_ID, index);

  /* write SDS metadata */
  ret = SDsetattr(sr->band_id, "long_name", DFNT_CHAR8, strlen(th->longName), th->longName);
  if (ret == FAILURE) {
    printf("Can't write SR long_name for SDS %s\n", th->sdsName);
    return FAILURE;
  } 
    
  ret = SDsetattr(sr->band_id, "units", DFNT_CHAR8, strlen(th->units), th->units);
  if (ret == FAILURE) {
    printf("Can't write SR units for SDS %s", th->sdsName);
    return FAILURE;
  } 

  ret = SDsetattr(sr->band_id, "valid_range", DFNT_INT16, 2, th->range);
  if (ret == FAILURE) {
    printf ("Can't write SR valid_range for SDS %s", th->sdsName);
    return FAILURE;
  } 
  
  ret = SDsetattr(sr->band_id, "_FillValue", DFNT_INT16, 1, &(th->fillValue));
  if (ret == FAILURE) {
    printf ("Can't write SR _FillValue for SDS %s", th->sdsName);
    return FAILURE;
  } 
  
  ret = SDsetattr(sr->band_id, "scale_factor", DFNT_FLOAT64, 1, &(th->scale));
  if (ret == FAILURE) {
    printf ("Can't write SR scale_factor for SDS %s", th->sdsName);
    return FAILURE;
  } 

  /* open thermal hdf file and get sds_id from given sds_name */  
  if ((th->SD_ID = SDstart(th->fileName, DFACC_READ))<0) {
    printf("Can't open file %s\n",th->fileName);
    return FAILURE;
  }

  if ((index=SDnametoindex(th->SD_ID, th->sdsName))<0) {
    printf("Not successful in convert sdsName %s to index\n", th->sdsName);
    return FAILURE;
  }
  th->band_id = SDselect(th->SD_ID, index);

  return exist;
}


int cleanUp(GRID_SR *sr) 
{

  if((SDendaccess(sr->band_id)) == FAILURE) {
    printf("SDendaccess for %s aot error!\n", sr->fileName);
    return FAILURE;
  }

  if((SDend(sr->SD_ID)) == FAILURE) {
    printf("SDend for %s error!\n", sr->fileName);
    return FAILURE;
  }
  return SUCCESS;

}


