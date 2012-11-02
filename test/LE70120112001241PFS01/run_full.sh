#!/bin/sh

lndsr lndsr.L71012011_01120010829.txt \
	-check 2450 850  \
	-check 2460 860  \
	-check 2800 1000 \
        -check 3050 970  \
        -checkll 70.46175 -54.99569 \
	|& tee full/lndsr.log

