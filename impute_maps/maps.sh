#!/bin/sh

for m in $(seq 1 30);\
	do python map_maker.py ../chrsplit/170112_merged.chr$m.bim imputemap.chr$m.map; done
