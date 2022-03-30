#!/bin/bash

if [ -z $1 ];
then
	files=$(ls *.log | sort -V)
else
	files=$@
fi


echo 'filename,total_energy,zpe,thermal,enthalpy,free_energy,imaginary_freq_qtd,lowest_freq' # > $output 

for file in $files;
do

	total_energy=$(grep 'Recovered' $file | tail -1 | awk '{print $3}')
	zero_point=$(grep 'and zero-point Energies' $file | awk '{print $7}')
	thermal_energy=$(grep 'thermal Energies' $file | awk '{print $7}')
	enthalpy=$(grep 'thermal Enthalpies' $file | awk '{print $7}')
	free_energy=$(grep 'thermal Free Energies' $file | awk '{print $8}')
	imaginary_freq_qtd=$(grep 'imaginary frequencies (negative Signs)' $file | awk '{print $2}')
	if [ -z $imaginary_freq_qtd ]; then imaginary_freq_qtd=0; fi 
	lowest_freq=$(grep 'Frequencies' $file | head -n 1 | awk '{print $3}')
	
	#write energies of $file to output
	echo $file,$total_energy,$zero_point,$thermal_energy,$enthalpy,$free_energy,$imaginary_freq_qtd,$lowest_freq  #>> $output
done

exit
