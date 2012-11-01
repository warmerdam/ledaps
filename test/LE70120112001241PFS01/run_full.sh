#!/bin/sh

lndsr lndsr.L71012011_01120010829.txt \
	-check 2450 850 \
	-check 2460 860 \
	-check 2800 1000 \
	|& tee full/lndsr.log

