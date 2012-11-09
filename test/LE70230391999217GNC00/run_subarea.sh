#!/bin/bash

lndsr lndsr.L71023039_03919990805.txt \
	-srcwin 2400 800 800 400 \
	-check 2410 805  \
	-check 2460 860  \
	-check 2800 1000 \
        -check 3050 970  \
	|& tee subarea/lndsr.log

ee_lndsr subarea \
	-check 2410 805  \
	-check 2460 860  \
	-check 2800 1000 \
        -check 3050 970  \
	|& tee subarea/ee_lndsr.log

