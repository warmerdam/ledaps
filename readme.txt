Installation Guide 
(tested in Linux bash, Feb. 2012)

1. Install ancillary files

tar -xvzf ledapsAnc.20120210.tar.gz

2. Install source files

tar -xvzf ledapsSrc.20120210.tar.gz
cd ledapsSrc/src
source env.sh
make
make install
make clean

This will create a list of executable files under bin/
lndapp lndcal  lndcsm  lndpm  lndsr  lndsrbm.ksh sixsV1.0B
cmrbv1.0  compadjn  comptemp  geo2xy SDSreader3.0 xy2geo
(tested in gcc and gfortran compiler)

* Note that new cloud mask program (lndsrbm) calls command 
"ncdump" which is a standard HDF command. This release includes 
a linux version ncdump command in bin directory. You need replace 
ncdump program if it runs from different OS systems. If you don't 
want to run optional lndsr-based cloud mask (lndsrbm), you will 
not need a ncdump program)

3. Setup environment

export ANC_PATH="directory_saved_ancillary_files"
(or in c shell use 
setenv ANC_PATH "directory_saved_ancillary_files")

(note to replace "directory_saved_ancillary_files" with 
the actual path that stores the ancillary files)

4. Test 

download Landsat files and then run following commands separately
lndpm <Landsat_meta_file>
lndcal <lndcal_input_text>
lndcsm <lndcsm_input_text>
lndsr <lndsr_input_text>
lndsrbm.ksh <lndsr_input_text> (optional)

OR simply run the combined commands to generate all (includes lndsrbm)
do_ledaps.csh <Landsat_meta_file>

(verify result with the provided example)
mkdir test
cp example.tar.gz test/
cd test
tar -xvzf example.tar.gz
do_ledaps.csh L5014033_03320060801_MTL.txt

(then compare result lndsr.L5014033_03320060801.hdf with 
L5TSR.a2006213.w2p014r033.020.2011045152828.hdf from LEDAPS 
operational system)


5. Output

lndcal.*.hdf: top-of-atmosphere (TOA) reflectance in HDF format
lndth.*.hdf: themal brightness temperature in HDF format
lndcsm.*.hdf: cloud mask from ACCA algorithm
lndsr.*.hdf: surface reflectance in HDF format
(each file is also associated with a ENVI header file that contains geographic information)

REFERENCES

Jeffrey G. Masek, Eric F. Vermote, Nazmi E. Saleous, Robert Wolfe,
Forrest G. Hall, Karl F. Huemmrich, Feng Gao, Jonathan Kutler, and 
Teng-Kui Lim, 2006, A Landsat Surface Reflectance Dataset for North
America, 1990-2000, IEEE Geoscience and Remote Sensing Letters, 
Vol. 3, No. 1, pp. 68-72.

