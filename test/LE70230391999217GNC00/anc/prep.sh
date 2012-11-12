#!/bin/bash

ANC_BOUND="-181.25 91.25 178.75 -91.25"

gdal_translate -a_srs WGS84 -a_ullr $ANC_BOUND \
  -unscale -ot Float32 \
  HDF4_SDS:UNKNOWN:"REANALYSIS_1999217.hdf":2 \
  REANALYSIS_1999217_slp.tif

gdal_translate -a_srs WGS84 -a_ullr $ANC_BOUND \
  -unscale -ot Float32 \
  HDF4_SDS:UNKNOWN:"REANALYSIS_1999217.hdf":3 \
  REANALYSIS_1999217_pr_wtr.tif

gdal_translate -a_srs WGS84 -a_ullr $ANC_BOUND \
  -unscale -ot Float32 \
  HDF4_SDS:UNKNOWN:"REANALYSIS_1999217.hdf":4 \
  REANALYSIS_1999217_air.tif

gdal_translate -a_srs WGS84 -a_ullr -180 90 180 -90 -a_nodata -9999 \
  HDF4_SDS:UNKNOWN:"CMGDEM.hdf":0 CMGDEM.tif

gdal_translate -a_srs WGS84 -a_ullr -180 90 180 -90 \
  TOMS_1999217.hdf TOMS_1999217.tif

rm REANALYSIS*tif.aux.xml 

