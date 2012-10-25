#ifndef uchar
# define uchar unsigned char
#endif

#ifndef M_PI
# define M_PI 3.141592653589793
#endif
#define ROVRDG (M_PI/180.)

/****************************************************************************
***** General purpose MACROS
****************************************************************************/
#define max(a,b)	((a)<(b) ? (b) : (a))
#define min(a,b)	((a)>(b) ? (b) : (a))
#define ABS(F) (((F) < 0)? -(F):(F))

#define   SAT_PROJ		99
#define	  LEN_TBM_HEADER	122
#define	  LEN_GAC_RECORD	6440
#define	  LEN_LAC_RECORD	7400

/*****************************************************************/
/* Error Codes */
/*****************************************************************/
#define 	OK 			0

#define		SOURCE_NOT_OPENED	10
#define		READ_ERROR		11
#define		WRITE_ERROR		12
#define		MEMORY_ERROR		13

#define		BAD_TBM_FORMAT		40
#define		NO_TBM_HEADER           41
#define		UNKOWN_DATA_TYPE	42

#define		BAD_SIGNATURE		90

#define   NO_ENTRY -9999

#define   SPHERE_RAD 6370997. /* radius of reference sphere in m */

/*****************************************************************/
/* data types */
/*****************************************************************/
#define TYPE_LAC	1
#define TYPE_GAC	2
#define TYPE_HRPT	3

/*****************************************************************/
/* satellites */
/*****************************************************************/
#define NOAA_06		2
#define NOAA_07		4
#define NOAA_08		6
#define NOAA_09		7
#define NOAA_10		8
#define NOAA_11		1
#define NOAA_12		5
#define NOAA_14		3
#define UNKNOWN		-99

/*****************************************************************/
/* parameters */
/*****************************************************************/
#define PAR_UNKNOWN	-99
#define PAR_DEFAULTSAT	1
#define PAR_DEFAULTTYPE	2
#define PAR_DIGITALCNTS	3
#define PAR_ASCENDING	4
#define PAR_PROJECTION	5
#define PAR_LAT0	6
#define PAR_LATMIN	7
#define PAR_LATMAX	8
#define PAR_DELTALAT	9
#define PAR_LON0	10
#define PAR_LONMIN	11
#define PAR_LONMAX	12
#define PAR_DELTALON	13
#define PAR_XSIZE	14
#define PAR_YSIZE	15
#define PAR_AUTOMATIC	16
#define PAR_NBOFTAPES	17
#define PAR_DARKCH1	18
#define PAR_DARKCH2	19
#define PAR_CALCOEF1	20
#define PAR_CALCOEF2	21
#define PAR_NIGHT	22
#define PAR_OUTPUTPREF	23
#define PAR_SAVE	24
#define PAR_NBOFDAYS	25
#define PAR_MAPFILE	26
#define PAR_PATCHMAP	27
#define PAR_PRESCAN	28
#define PAR_CELLWIDTH	29
#define PAR_CELLHEIGHT	30
#define PAR_GRIDSIZE	31
#define PAR_GRIDCOLOR	32
#define PAR_OUTPUT	33
#define PAR_GRIDLATMIN		34
#define PAR_GRIDLATMAX		35
#define PAR_GRIDLONMIN		36
#define PAR_GRIDLONMAX		37
#define PAR_GRIDLATSTEP		38
#define PAR_GRIDLONSTEP		39
#define PAR_GSHHS		40
#define PAR_LAND		41
#define PAR_ELEVATION		42
#define PAR_LATLON		43
#define PAR_HDF		44

/****
***** parameters related to projection
****/
#define PAR_STDPAR1		70
#define PAR_STDPAR2		71
#define PAR_CTRMERID		72
#define PAR_ORIGINLAT		73
#define PAR_SCALEFACTOR		74
#define PAR_TRUESCALELAT	75
#define PAR_LONBELOWPOLE	76
#define PAR_PIXELWIDTH		77
#define PAR_PIXELHEIGHT		78
#define PAR_CENTRALLON		79
#define PAR_CENTRALLAT		80

/****
***** units
****/
#define UNIT_UNKNOWN		-1
#define UNIT_M			1
#define UNIT_KM			2
#define UNIT_DEG		3
#define UNIT_RD			4

/*****************************************************************/
/* projections */
/*****************************************************************/
#define PROJ_GEO	0
#define PROJ_UTM	1
#define PROJ_SPCS	2
#define PROJ_ALBERS	3
#define PROJ_LAMCC	4
#define PROJ_MERCAT	5
#define PROJ_PS		6
#define PROJ_POLYC	7
#define PROJ_EQUIDC	8
#define PROJ_TM		9
#define PROJ_STEREO	10
#define PROJ_LAMAZ	11
#define PROJ_AZMEQD	12
#define PROJ_GNOMON	13
#define PROJ_ORTHO	14
#define PROJ_GVNSP	15
#define PROJ_SNSOID	16
#define PROJ_EQRECT	17
#define PROJ_MILLER	18
#define PROJ_VGRINT	19
#define PROJ_HOM	20
#define PROJ_ROBIN	21
#define PROJ_SOM	22
#define PROJ_ALASKA	23
#define PROJ_GOOD	24
#define PROJ_MOLL	25
#define PROJ_IMOLL	26
#define PROJ_HAMMER	27
#define PROJ_WAGIV	28
#define PROJ_WAGVII	29
#define PROJ_OBLEQA	30
#define PROJ_SAT	99
#define PROJ_UNKNOWN	-1
