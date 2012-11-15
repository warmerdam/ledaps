#!/usr/bin/env python

import os
import sys

from osgeo import gdal
from osgeo import osr

def extract_subwin(srcfile, dstfile, srcwin, pixel_type):

    cmd = "gdal_translate -srcwin %d %d %d %d -ot %s '%s' %s -of PAux" \
          % (srcwin[0], srcwin[1], srcwin[2], srcwin[3],
             pixel_type, srcfile, dstfile)
    print cmd
    os.system(cmd)

def extract_cal_bands(srcwin, src_hdf_file):
    src = 'HDF4_EOS:EOS_GRID:"'+src_hdf_file+'":Grid:band'
    bands = ['1', '2', '3', '4', '5', '7']
    for band in bands:
        extract_subwin(src + band, 'lndcal_band'+band, srcwin, 'Int16')

def extract_lndcal_QA(srcwin, src_hdf_file, srs, gt):
    srcfile = 'HDF4_SDS:UNKNOWN:"'+src_hdf_file+'":6'
    cmd = "gdal_translate -srcwin %d %d %d %d -ot %s '%s' %s -of PAux -a_srs %s -a_ullr %.2f %.2f %.2f %.2f" \
          % (srcwin[0], srcwin[1], srcwin[2], srcwin[3],
             'Byte', srcfile, 'lndcal_QA', srs,
             gt[0], gt[3],
             gt[0] + srcwin[2] * gt[1],
             gt[3] + srcwin[3] * gt[5])
    print cmd
    os.system(cmd)

    

def extract_cal_info(src_hdf_file):
    src = 'HDF4_EOS:EOS_GRID:"'+src_hdf_file+'":Grid:band1'
    ds = gdal.Open(src)

    return ((ds.RasterXSize,ds.RasterYSize),
            'EPSG:32615',
            ds.GetGeoTransform())

def extract_anc( srs, geotransform, size, srcfile, dstfile):
    xmin = geotransform[0]
    xmax = xmin + geotransform[1] * size[0]
    ymax = geotransform[3]
    ymin = ymax + geotransform[5] * size[1]
    
    cmd = "gdalwarp -overwrite -r bilinear -t_srs %s -ts %d %d -te %f %f %f %f -ot Float32 -of PAux '%s' %s" \
          % (srs, size[0], size[1], xmin, ymin, xmax, ymax, srcfile, dstfile)
    print cmd
    os.system(cmd)


def extract_metadata(src_file, dst_file):
    ds = gdal.Open(src_file)
    md = ds.GetMetadata()
    ds = None

    of = open(dst_file,"w")
    for key in md.keys():
        of.write(key + '\n')
        of.write(md[key]+'\n')

def write_parms(filename, parm_list):
    of = open(filename,'w')
    for parm in parm_list:
        of.write(str(parm) + '\n')
       
# Mainline

ar_subsample = 40
src_cal='../lndcal.L71023039_03919990805.hdf'
src_th='../lndth.L71023039_03919990805.hdf'
src_csm='../lndcsm.L71023039_03919990805.hdf'
srcwin=(1200,3600,800,400)

(full_input_size, srs, full_gt) = extract_cal_info(src_cal)

win_gt = (full_gt[0] + srcwin[0] * full_gt[1],
          full_gt[1],
          full_gt[2],
          full_gt[3] + srcwin[1] * full_gt[5],
          full_gt[4],
          full_gt[5])

utm_zone = int(srs[-2:])
ul_x = win_gt[0]
ul_y = win_gt[3]
pixel_size = win_gt[1]

aerosol_geotransform = (win_gt[0],
                        win_gt[1]*ar_subsample,
                        win_gt[2]*ar_subsample,
                        win_gt[3],
                        win_gt[4]*ar_subsample,
                        win_gt[5]*ar_subsample)

extract_cal_bands(srcwin, src_cal)
extract_lndcal_QA(srcwin, src_cal, srs, win_gt)
extract_subwin(src_csm, 'csm', srcwin, 'Byte')
extract_subwin(src_th, 'thermal', srcwin, 'Int16')

aerosol_size = (srcwin[2]/ar_subsample, srcwin[3]/ar_subsample)

extract_anc( srs, aerosol_geotransform, aerosol_size,
             '../anc/REANALYSIS_1999217_pr_wtr.tif',
             'pr_wtr')
extract_anc( srs, aerosol_geotransform, aerosol_size,
             '../anc/REANALYSIS_1999217_slp.tif',
             'slp')
extract_anc( srs, aerosol_geotransform, aerosol_size,
             '../anc/TOMS_1999217.tif',
             'ozone')
extract_anc( srs, aerosol_geotransform, aerosol_size,
             '../anc/CMGDEM.tif',
             'dem')
extract_anc( srs, win_gt, (srcwin[2], srcwin[3]),
             '../anc/REANALYSIS_1999217_air.tif',
             'air')

extract_subwin('../AEROSOL_RUN_RESULTS.hdf', 'ar',
               (srcwin[0]/ar_subsample,
                srcwin[1]/ar_subsample,
                srcwin[2]/ar_subsample,
                srcwin[3]/ar_subsample),
               'Int16')

extract_metadata(src_cal,'metadata.txt')

write_parms('ee_lndsr_parms.txt',
            [ 'SIXS_RUN_RESULTS.TXT',
              full_input_size[0], full_input_size[1],
              srcwin[0], srcwin[1], srcwin[2], srcwin[3],
              aerosol_size[0], aerosol_size[1],
              utm_zone, ul_x, ul_y, pixel_size ] )
