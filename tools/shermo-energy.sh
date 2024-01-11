#!/bin/bash

output_file="shermo.csv"

echo filename,electronic energy,enthalpy,free energy > $output_file
for f in $(find . -type f -name "*.shermo" | sort -V); do
    electronic_energy=$(grep "Electronic energy:" $f | tail -1 | awk '{print $3}')
    enthalpy=$(grep "Sum of electronic energy and thermal correction to H:" $f | tail -1 | awk '{print $((NF - 1))}')
    free_energy=$(grep "Sum of electronic energy and thermal correction to G:" $f | tail -1 | awk '{print $((NF - 1))}')
    echo "$f","$electronic_energy","$enthalpy","$free_energy" >> $output_file
done
