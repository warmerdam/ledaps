export TOP=`pwd`
export BIN=`echo $TOP| sed 's/src/bin/'`

export GEOTIFF_INC=$TOP/include
export GEOTIFF_LIB=$TOP/lib
export HDFINC=$TOP/include
export HDFEOS_INC=$TOP/include
export TIFFINC=$TOP/include
export TIFFLIB=$TOP/lib
export HDFEOS_LIB=$TOP/lib
export HDFLIB=$TOP/lib

export PATH=$BIN:$PATH
