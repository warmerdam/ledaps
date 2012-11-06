#!/bin/sh

lndsr lndsr.L71023039_03919990805.txt \
	-srcwin 2400 800 800 400 \
	-check 2450 850  \
	-check 2460 860  \
	-check 2800 1000 \
        -check 3050 970  \
        -checkll 30.991165 -92.0439054 \
	|& tee subarea/lndsr.log
