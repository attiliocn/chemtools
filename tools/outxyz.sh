#!/bin/bash

for outputFile in $@
do

    xMolFile="${outputFile%.*}.xyz"

    startline=$(grep -n "CARTESIAN COORDINATES (ANGSTROEM)" $outputFile | tail -1 | sed 's/\:/ /'| awk '{print $1}')
    endline=$(grep -n "CARTESIAN COORDINATES (A.U.)" $outputFile | tail -1 | sed 's/\:/ /'| awk '{print $1}')
    
    startLine=$(($startline+2))
    endLine=$(($endline-3))

    sedCommand="sed -n $startLine,"$endLine"p $outputFile"

    $sedCommand | wc -l > $xMolFile
    echo $outputFile >> $xMolFile
    $sedCommand >> $xMolFile

done

