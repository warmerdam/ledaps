#!/usr/bin/env python

import os
import sys
import filecmp

from osgeo import gdal, osr

#######################################################
def diff_band( output_filename, template_ds, golden_band, new_band):

    # Create output file.
    drv = gdal.GetDriverByName('GTiff')
    out_db = drv.Create(output_filename,
                        golden_band.XSize,
                        golden_band.YSize,
                        1, gdal.GDT_Float32)

    out_db.SetGeoTransform(template_ds.GetGeoTransform())
    out_db.SetProjection(template_ds.GetProjection())

    out_band = out_db.GetRasterBand(1)

    for line in range(golden_band.YSize):
        golden_line = golden_band.ReadAsArray(0, line, golden_band.XSize, 1)
        new_line = new_band.ReadAsArray(0, line, golden_band.XSize, 1)
        diff_line = golden_line - new_line
        out_band.WriteArray(diff_line, 0, line)

#######################################################
#
# Mainline
#

if __name__ == '__main__':

  # Default GDAL argument parsing.
  argv = gdal.GeneralCmdLineProcessor( sys.argv )
  if argv is None:
    sys.exit( 0 )

  if len(argv) == 1:
    Usage()

  # Script argument parsing.
  golden_file = None
  new_file = None

  i = 1
  while  i < len(argv):

    if golden_file is None:
      golden_file = argv[i]

    elif new_file is None:
      new_file = argv[i]

    else:
      print('Urecognised argument: ' + argv[i])
      Usage()

    i = i + 1
    # next argument


  golden_ds = gdal.Open(golden_file)
  new_ds = gdal.Open(new_file)

  diff_band('diff.tif', golden_ds,
            golden_ds.GetRasterBand(1), new_ds.GetRasterBand(1))
  
  
