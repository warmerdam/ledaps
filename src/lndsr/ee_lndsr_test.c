/**
 * lndsr_ee_test.c - Main simulating EarthEngine calls to lndsr_ee.c.
 *
 * Author: Frank Warmerdam <warmerdam@google.com>
 */

#include "ee_lndsr.h"

/************************************************************************/
/*                               Usage()                                */
/************************************************************************/
static void Usage() {
    printf( "Usage: ee_lndsr [-check is il] <source_directory>\n" );
    exit(1);
}

/************************************************************************/
/*                           ReadWholeFile()                            */
/************************************************************************/

static void *ReadWholeFile( const char *filename, long expected_size ) {
    FILE *fp = fopen(filename, "rb");
    long size;
    void *data;

    if( fp == NULL ) {
        perror(filename);
        exit(2);
    }

    fseek(fp, 0, SEEK_END);
    size = ftell(fp);
    fseek(fp, 0, SEEK_SET);

    if (expected_size != 0 && size != expected_size) {
        fprintf( stderr, "file %s not expected %ld bytes long, got %ld bytes.\n", 
                 filename, expected_size, size );
        exit(2);
    }

    data = calloc(size+1,1);
    if (fread(data, size, 1, fp) != 1) {
        perror("fread");
        exit(2);
    }

    fclose(fp);
    
    return data;
}

/************************************************************************/
/*                          WriteRasterFile()                           */
/************************************************************************/

static void WriteRasterFile( const char *filename, short *data,
                             int size_s, int size_l ) {

    FILE *fp = fopen(filename, "wb");
    long size = size_s * size_l * sizeof(short);

    if( fp == NULL ) {
        perror(filename);
        exit(2);
    }

    if (fwrite(data, size, 1, fp) != 1) {
        perror("fread");
        exit(2);
    }

    fclose(fp);
}

/************************************************************************/
/*                                main()                                */
/************************************************************************/
int main(int argc, char **argv) {

    const char *src_dir = NULL;

    char filename[2048], line[128], *r;
    int i;

    char *metadata, *sixs_filename;

    int full_input_size_s=0, full_input_size_l=0;
    int input_offset_s=0, input_offset_l=0;
    int input_size_s=0, input_size_l=0, nband;
    short *band1, *band2, *band3, *band4, *band5, *band7;
    short *thermal;
    int8 *lndcal_QA;
    char *csm_mask;

    int utm_zone;
    double ul_x, ul_y, pixel_size;

    int ar_size_s, ar_size_l;
    short *aerosol;
    float *anc_SP, *anc_WV, *anc_ATEMP, *anc_O3, *anc_dem;

    FILE *fpParms;

/* -------------------------------------------------------------------- */
/*      Parse arguments.                                                */
/* -------------------------------------------------------------------- */
    for( i=1; i < argc; i++ ) 
    {
        if (strcmp(argv[i],"-check") == 0 && i < argc-2)
        {
            add_check_pixel( atoi(argv[i+1]), atoi(argv[i+2]) );
            i += 2;

        } else if( src_dir == NULL ) {
            src_dir = argv[i];

        } else { 
            printf( "Unrecognised argument: %s\n", argv[i] );
            Usage();
        }
    }

    if( src_dir == NULL )
        Usage();

/* -------------------------------------------------------------------- */
/*      Read the parms file.                                            */
/* -------------------------------------------------------------------- */
    sprintf( filename, "%s/ee_lndsr_parms.txt", src_dir );
    fpParms = fopen(filename, "r");

    r = fgets(line, sizeof(line), fpParms);
    sixs_filename = strdup(line);
    sixs_filename[strlen(sixs_filename)-1] = '\0'; /* trim newline */
    
    r = fgets(line, sizeof(line), fpParms);
    full_input_size_s = atoi(line);

    r = fgets(line, sizeof(line), fpParms);
    full_input_size_l = atoi(line);

    r = fgets(line, sizeof(line), fpParms);
    input_offset_s = atoi(line);

    r = fgets(line, sizeof(line), fpParms);
    input_offset_l = atoi(line);

    r = fgets(line, sizeof(line), fpParms);
    input_size_s = atoi(line);

    r = fgets(line, sizeof(line), fpParms);
    input_size_l = atoi(line);

    r = fgets(line, sizeof(line), fpParms);
    ar_size_s = atoi(line);

    r = fgets(line, sizeof(line), fpParms);
    ar_size_l = atoi(line);

    r = fgets(line, sizeof(line), fpParms);
    utm_zone = atoi(line);

    r = fgets(line, sizeof(line), fpParms);
    ul_x = atof(line);

    r = fgets(line, sizeof(line), fpParms);
    ul_y = atof(line);

    r = fgets(line, sizeof(line), fpParms);
    pixel_size = atof(line);

    fclose(fpParms);

    adjust_check_pixels( input_offset_s, input_offset_l );

/* -------------------------------------------------------------------- */
/*      Read the metadata file.                                         */
/* -------------------------------------------------------------------- */
    sprintf( filename, "%s/metadata.txt", src_dir );
    metadata = (char *) ReadWholeFile(filename, 0);

/* -------------------------------------------------------------------- */
/*      Read the source landsat data.                                   */
/* -------------------------------------------------------------------- */
    sprintf( filename, "%s/lndcal_band1", src_dir );
    band1 = (short *) ReadWholeFile(filename, 2*input_size_s*input_size_l);

    sprintf( filename, "%s/lndcal_band2", src_dir );
    band2 = (short *) ReadWholeFile(filename, 2*input_size_s*input_size_l);

    sprintf( filename, "%s/lndcal_band3", src_dir );
    band3 = (short *) ReadWholeFile(filename, 2*input_size_s*input_size_l);

    sprintf( filename, "%s/lndcal_band4", src_dir );
    band4 = (short *) ReadWholeFile(filename, 2*input_size_s*input_size_l);

    sprintf( filename, "%s/lndcal_band5", src_dir );
    band5 = (short *) ReadWholeFile(filename, 2*input_size_s*input_size_l);

    sprintf( filename, "%s/lndcal_band7", src_dir );
    band7 = (short *) ReadWholeFile(filename, 2*input_size_s*input_size_l);

    sprintf( filename, "%s/thermal", src_dir );
    thermal = (short *) ReadWholeFile(filename, 2*input_size_s*input_size_l);

    sprintf( filename, "%s/lndcal_QA", src_dir );
    lndcal_QA = (int8 *) ReadWholeFile(filename, input_size_s*input_size_l);

    sprintf( filename, "%s/csm", src_dir );
    csm_mask = (int8 *) ReadWholeFile(filename, input_size_s*input_size_l);

/* -------------------------------------------------------------------- */
/*      Read the ancillary data.                                        */
/* -------------------------------------------------------------------- */
    sprintf( filename, "%s/ar", src_dir );
    aerosol = (short *) ReadWholeFile(filename, 2*ar_size_s*ar_size_l);

    sprintf( filename, "%s/slp", src_dir );
    anc_SP = (float *) ReadWholeFile(filename, 4*4*ar_size_s*ar_size_l);

    sprintf( filename, "%s/pr_wtr", src_dir );
    anc_WV = (float *) ReadWholeFile(filename, 4*4*ar_size_s*ar_size_l);

    sprintf( filename, "%s/air", src_dir );
    anc_ATEMP = (float *) ReadWholeFile(filename,4*4*input_size_s*input_size_l);

    sprintf( filename, "%s/ozone", src_dir );
    anc_O3 = (float *) ReadWholeFile(filename, 4*ar_size_s*ar_size_l);

    sprintf( filename, "%s/dem", src_dir );
    anc_dem = (float *) ReadWholeFile(filename, 4*ar_size_s*ar_size_l);

/* -------------------------------------------------------------------- */
/*      Allocate additional output arrays.                              */
/* -------------------------------------------------------------------- */
    short *atmos_opacity = (short *) calloc(input_size_s*input_size_l, 2);
    short *lndsr_QA = (short *) calloc(input_size_s*input_size_l, 2);

/* -------------------------------------------------------------------- */
/*      Invoke the mainline.                                            */
/* -------------------------------------------------------------------- */
    int error;

    error = ee_lndsr_main( metadata, sixs_filename,
                           full_input_size_s, full_input_size_l,
                           input_offset_s, input_offset_l,
                           input_size_s, input_size_l, 6 /* nband */,
                           band1, band2, band3, band4, band5, band7,
                           thermal, lndcal_QA, csm_mask,
                           atmos_opacity, lndsr_QA,
                           utm_zone, ul_x, ul_y, pixel_size,
                           ar_size_s, ar_size_l,
                           aerosol, anc_SP, anc_WV, anc_ATEMP,
                           anc_O3, anc_dem);

    if (error != 0) {
        printf( "Error in ee_lndsr_main()!\n" );
    }

/* -------------------------------------------------------------------- */
/*      Write out the results.                                          */
/* -------------------------------------------------------------------- */
    sprintf( filename, "%s/lndsr_band1", src_dir );
    WriteRasterFile(filename, band1, input_size_s, input_size_l);
    sprintf( filename, "%s/lndsr_band2", src_dir );
    WriteRasterFile(filename, band2, input_size_s, input_size_l);
    sprintf( filename, "%s/lndsr_band3", src_dir );
    WriteRasterFile(filename, band3, input_size_s, input_size_l);
    sprintf( filename, "%s/lndsr_band4", src_dir );
    WriteRasterFile(filename, band4, input_size_s, input_size_l);
    sprintf( filename, "%s/lndsr_band5", src_dir );
    WriteRasterFile(filename, band5, input_size_s, input_size_l);
    sprintf( filename, "%s/lndsr_atmos_opacity", src_dir );
    WriteRasterFile(filename, atmos_opacity,  input_size_s, input_size_l);
    sprintf( filename, "%s/lndsr_QA", src_dir );
    WriteRasterFile(filename, lndsr_QA, input_size_s, input_size_l);

    exit(0);
}
