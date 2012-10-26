#!/bin/sh

set errexit

meta_file=L71012011_01120010829_MTL.txt
meta=`echo $meta_file | sed -e 's/.txt//' -e 's/_MTL//' -e 's/.met//'`

export ANC_PATH=../../ANC
LEDAPS_BIN=../../bin
export PATH=$LEDAPS_BIN:$PATH

# run LEDAPS modules
set -x

lndpm $meta_file
lndcal lndcal.$meta.txt
lndcsm lndcsm.$meta.txt
lndsr lndsr.$meta.txt

# Updates lndsr_QA and adds band6 in lndsr.*.hdf file
#lndsrbm.ksh lndsr.$meta.txt

