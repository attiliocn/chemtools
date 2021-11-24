#!/bin/bash

FOLDER_PREFIX='*phenol*/'

mkdir results
mkdir results/tbp results/complex results/ts

for i in $FOLDER_PREFIX; do cp $i/tbp/*xtb_tbp.xyz results/tbp/; done
for i in $FOLDER_PREFIX; do cp $i/complex/*-const.xyz results/complex/; done
 
for i in $FOLDER_PREFIX; do cp $i/ts_forward/*ts-fw.xyz results/ts/; done
for i in $FOLDER_PREFIX; do cp $i/ts_reverse/*ts-rv.xyz results/ts/; done

for i in $( ls -d $FOLDER_PREFIX | sort -V); do echo $i $(cat $i/ts_forward/*forward.output | grep "forward  barrier (kcal)") >> results/barrier_ts-fw.txt; done
for i in $( ls -d $FOLDER_PREFIX | sort -V); do echo $i $(cat $i/ts_forward/*forward.output | grep "reaction energy  (kcal)") >> results/reaction-energy_ts-fw.txt; done


