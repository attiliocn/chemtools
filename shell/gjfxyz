#!/bin/bash

for i in $@; do
	echo "GJF ---> XYZ $i"

	tr -d '\15\32' < $i > temp
	mv temp $i
	coordinatesLine1=$(grep -n "^$" $i | head -n 2 | tail -1 | sed 's/://')
	coordinatesLine1=$((coordinatesLine1 + 2))
	#coordinatesLine1=$(grep -nE "^\s+" $i | head -n 1 | sed 's/:/\t/'| awk '{print $1}')
	tail -n +"$coordinatesLine1" $i | sed '/^$/d'| sed '${/^ $/d}' > temp
	cat temp | wc -l > ${i%.*}.xyz
	echo ${i%.*} >> ${i%.*}.xyz
	cat temp >> ${i%.*}.xyz
	echo "-----------"
done
rm temp

