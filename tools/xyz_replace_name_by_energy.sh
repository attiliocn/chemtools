#!/bin/bash

# requires energies.csv where the first column
# is the filename (ending with .log) and the second column
# is the energy for the entry

for i in $@; do
    basename=${i%.*}
    energy=$(awk -F "," -v output="$basename".log '$1==output {print $2}'  energies.csv)
    sed -i "2s/.*/$energy/" $i
done
