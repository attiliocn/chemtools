#!/bin/bash

element=$1
auto=$2

PERIODIC_TABLE=(H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt)


if [ $auto == 'true' ]; then
  
	for k in "${PERIODIC_TABLE[@]}"; do
		counter=1
		echo "==================="
		echo $k

		for i in $(ls *.xyz|sort -V); do
		    if [ ! -z $(tail -n +3 $i | grep $k |tail -1|awk '{print $1}') ];then 
			echo $counter $i
			counter=$((counter+1)) 
		    fi
		done
		echo "==================="
	done
fi

#####
#unitary run
#####
if [ $auto == 'false' ]; then
	for i in $(ls *.com|sort -V)
	do 

	    if [ ! -z $(grep "$element" $i|tail -1|awk '{print $1}') ] 
	    then
		echo $i 
	    fi

	done
fi
