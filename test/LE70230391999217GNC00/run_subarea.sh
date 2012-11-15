#!/bin/bash

lndsr lndsr.L71023039_03919990805.txt \
	-debug \
        -srcwin 1200 3600 800 400 \
	-check 1210 3610  \
	-check 1260 3860  \
	-check 1600 3800 \
        -check 1850 3970  \
	|& tee subarea/lndsr.log

ee_lndsr subarea \
	-check 1210 3610  \
	-check 1260 3860  \
	-check 1600 3800 \
        -check 1850 3970  \
	|& tee subarea/ee_lndsr.log

