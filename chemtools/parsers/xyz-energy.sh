#!/bin/bash

echo "file,header" > energy.csv
for f in $@
do
    energies=$(head -n 2 $f | tail -1 | grep -oP "\-[0-9]+\.[0-9]+")
    i=1
    for h in $energies
    do
        echo "${f/.*/}"_"$i".xyz,$h >> energy.csv
        i=$((i+1))
    done
done