#!/bin/bash

for i in $@; do
	echo "GJF ---> XYZ $i"

    # convert DOS carriage return to UNIX CR
	tr -d '\15\32' < $i > temp
	mv temp $i

    # remove Gaussian Conectivity Table
    head -n $(grep -n "^$" "$i" | head -n 3 | tail -1 | sed 's/://') "$i" > temp
	mv temp $i

	coordinates_first_line=$(grep -n "^$" $i | head -n 2 | tail -1 | sed 's/://')
	coordinates_first_line=$((coordinates_first_line + 2))
	tail -n +"$coordinates_first_line" $i | sed '/^$/d'| sed '${/^ $/d}' > temp
	cat temp | wc -l > ${i%.*}.xyz
	echo ${i%.*} >> ${i%.*}.xyz
	cat temp >> ${i%.*}.xyz
    rm temp
	echo "-----------"
done

