#!/bin/bash

esiFile=$@

splitCriteria="^M06" #match the xyz name
#splitCriteria="^[0-9'’a-zA-Z]*-" #match the xyz name
#splitCriteria="^[a-zA-Z]\{1,2\}[(]" #used for monfort2014


#split
csplit -zs $esiFile "/$splitCriteria/" {*}
sed -i '/^[A-Z][ ]*=/d' xx* #delete lines containing E and G energies (example match G = or E = )
sed -i '/^S[0-9]/d' xx* #delete lines containing page numbers (example match S123)
sed -i '1 s/ //g' xx* #remove any space in the first line of the file 
sed -i '/^$/d' xx* #delete any blank line
sed -i '${/^ $/d}' xx* #delete last line if is blank but contain space
find . -empty -type f -delete

exit #DELETE THIS EXIT!!!
 
#insert atom numbers
for i in xx*
do
    numberOfLines=$(cat $i|wc -l)
    numberOfAtoms=$((numberOfLines-1))
    echo $numberOfAtoms > temp
    cat $i >> temp
    mv temp $i
done
#rename and run babel to standardize
c=1
for i in xx*
do
    filename=$(head -2 $i|tail -1)_$c.xyz
    echo $i $filename >> err
    obabel -ixyz $i -oxyz > $filename 2>> err
    #mv $i $filename
    c=$((c+1))
done

rm xx*
find . -empty -type f -delete
