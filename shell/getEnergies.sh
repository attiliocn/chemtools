#!/bin/bash

MAIN_DIR=$PWD
FOLDER_PREFIX='*phenol*/'
#FOLDER_PREFIX='*num*/'
TS_CTR_FOLDER='results/ts_constr/'

GENERIC_RUN=false
GENERIC_PREFIX='*traj*/'

if [ $GENERIC_RUN = false ]
then

	# TBP #
	#############################
	for i in $(ls -d $FOLDER_PREFIX | sort -V)
	do
		cd $i/tbp
		#echo $PWD
		echo ${i%/}_tbp,$(grep "TOTAL ENERGY" *.output | awk '{print $4}') >> $MAIN_DIR/results/energies_tbp.csv
		cd $MAIN_DIR	
	done

	# TS FW PATHFINDER #
	#############################
	for i in $(ls -d $FOLDER_PREFIX | sort -V)
	do
		cd $i/ts_forward
		#echo $PWD
		echo ${i%/}_ts-fw,$(head -n 2 xtbpath_ts.xyz | tail -1 | awk '{print $2}') >> $MAIN_DIR/results/energies_ts_semopt.csv
		cd $MAIN_DIR
	done

	# TS CONSTRAINED #
	#############################
	cd $TS_CTR_FOLDER
	for i in $(ls *.output | sort -V)
	do
		echo ${i%_*}_xtb-ts-fw,$(grep "TOTAL ENERGY" $i | awk '{print $4}') >> $MAIN_DIR/results/energies_ts_optconstrained.csv
	done
	cd $MAIN_DIR

	# COMPLEX  #
	#############################
	for i in $(ls -d $FOLDER_PREFIX | sort -V)
	do
		cd $i/complex
		#echo $PWD
		echo ${i%/}_complex,$(grep "TOTAL ENERGY" *.output | awk '{print $4}') >> $MAIN_DIR/results/energies_complex.csv
		cd $MAIN_DIR
	done
	exit
fi

# GENERIC RUN  #
#############################
if [ $GENERIC_RUN == 'true' ]; then

	for i in $(ls -d $GENERIC_PREFIX | sort -V)
	do
		cd $i/
		echo ${i%/},$(grep "TOTAL ENERGY" *.output | awk '{print $4}') >> $MAIN_DIR/results/energies_generic.csv
		cd $MAIN_DIR
	done

fi
