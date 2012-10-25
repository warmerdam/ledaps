#! /bin/sh
# shell script file to update the cloud mask
#
# where is this executable?
exe_dir=`echo $0 | sed -e "s|/[^/]*$||"`

lndsr_inp=$1
if test -z "$lndsr_inp"
then
  echo "FAIL  : no input filename"
  echo "USAGE : lndsrbm.ksh lndsr.*.txt"
  exit
fi
if test ! -f "$lndsr_inp"
then
  echo "FAIL  : '$lndsr_inp' not a valid file"
  echo "USAGE : lndsrbm.ksh lndsr.*.txt"
  exit
fi

# find lndsr*.hdf file
filein=`grep SR_FILE $lndsr_inp | awk '{print $3}'`
if test ! -f "$filein"
then
  echo "FAIL  : can't find '$filein'"
  exit
fi

#fine lndth*.hdf file
fileth=`grep TEMP_FILE $lndsr_inp | awk '{print $3}'`
if test ! -f "$fileth"
then
  echo "FAIL  : can't find '$fileth'"
  exit
fi

# find REANALYSIS ancillary file 
fileanc=`grep PRWV_FIL $lndsr_inp | awk '{print $3}'`
if test ! -f "$fileanc"
  then
    echo "FAIL  : no ancillary file '$fileanc'"
    exit
fi

echo "using ancillary data '$fileanc'"

case=`basename $filein .hdf`
case="$case.txt"

# append thermal band to lndsr*.hdf file
echo "appending thermal file $fileth to lndsr file $filein"
lndapp -sr $filein -th $fileth

# work with metadata
ncdump -h $filein > tmp.meta

# compute lat,long of the center of the scene
cat tmp.meta | grep Coordinate >tmp.updatecm
lon1=`grep WestBoundingCoordinate tmp.updatecm | awk '{print $3}'`
lon2=`grep EastBoundingCoordinate tmp.updatecm | awk '{print $3}'`
lat1=`grep NorthBoundingCoordinate tmp.updatecm | awk '{print $3}'`
lat2=`grep SouthBoundingCoordinate tmp.updatecm | awk '{print $3}'`
echo $lon1 $lon2 $lat1 $lat2
lonc=`echo $lon1 $lon2 | awk '{print ($1+$2)/2.}'`
latc=`echo $lat1 $lat2 | awk '{print ($1+$2)/2.}'`
echo $lonc $latc
ygrib=`echo $latc | awk '{print int((90.-$1)*73/180.)}'` 
xgrib=`echo $lonc | awk '{print int((180.+$1)*144/360.)}'` 
echo $ygrib $xgrib
$exe_dir/SDSreader3.0 -f $fileanc -w "$ygrib $xgrib 1 1" -v >tmp.dumpfileanc
grep SDS tmp.dumpfileanc | grep air | awk '{print $8}' | tr -d "," | awk '{print ($1*0.01)+512.81}' >tmp.airtemp
sctest=`grep  AcquisitionDate tmp.meta | awk '{print $3}' | awk -F "T" '{print $2}' | tr -d '"'`
if [[ $sctest != "00:00:00.000000Z" ]] 
then
scenetime=`grep  AcquisitionDate tmp.meta | awk '{print $3}' | awk -F "T" '{print $2}' | awk -F : '{print $1+$2/60.}'`
else
scenetime=`echo 10.5 $lonc | awk '{print $1-$2/15.}'`
fi
scenetimet=`echo $scenetime | awk '{printf "%d\n", int($1*1000000)}'`
scenetime=`echo $scenetime | awk '{printf "%f\n", int($1*100000)/100000}'`
if [[ $scenetimet -lt 0 ]]
then
scenetime=`echo $scenetime | awk '{print $1+24}'`
echo "WARNING WE ASSUME THE DATE IS GMT IS IT?"
fi
echo $scenetime
cat tmp.airtemp
tclear=`$exe_dir/comptemp $scenetime <tmp.airtemp`
echo $tclear

# compute scene orientation
# get row and col of image center
ccol=`grep  "XDim_Grid =" tmp.meta | tail -1 | awk '{print $3/2}'`
crow=`grep  "YDim_Grid =" tmp.meta | tail -1 | awk '{print $3/2}'`
echo $ccol $crow
# compute lat,long of the center
str=`$exe_dir/xy2geo $filein $ccol $crow`
clat=`echo $str | awk '{print $9}'`
clon=`echo $str | awk '{print $7}'`
echo $clat $clon
# compute lat,lon of the point 100 pixels north from the center 
cpcol=$ccol
cprow=`echo $crow | awk '{print $1-100}'`
str=`$exe_dir/xy2geo $filein $cpcol $cprow`
cplat=`echo $str | awk '{print $9}'`
cplon=`echo $str | awk '{print $7}'`
# now move to the longitude of the center of the image
cplon=$clon
# and compute the deviation in pixel
str=`$exe_dir/geo2xy $filein $cplon $cplat`
cscol=`echo $str | awk '{print $9}'`
csrow=`echo $str | awk '{print $7}'`
echo $cscol $csrow
deltay=`echo $csrow $crow | awk '{print $2-$1}'`
deltax=`echo $cscol $ccol | awk '{print $1-$2}'`
echo $deltax $deltay
adjnorth=`$exe_dir/compadjn $deltax $deltay`
echo $adjnorth
ts=`grep SolarZenith tmp.meta | awk '{print $3}' | sed -e "s/f//"`
fs=`grep SolarAzimuth tmp.meta | awk '{print $3}' | sed -e "s/f//"`
echo $tclear $ts 0.0 $fs $adjnorth
echo $tclear $ts 0.0 $fs $adjnorth >anc-$case
echo "Updating cloud mask"
$exe_dir/cmrbv1.0 $filein <anc-$case

rm -f tmp.airtemp tmp.dumpfileanc tmp.meta tmp.updatecm anc-$case