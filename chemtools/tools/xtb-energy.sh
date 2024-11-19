#!/bin/bash

echo "filename,electronic energy,enthalphy,free energy" > energies.csv

for xtb_output in $(find . -maxdepth 1 -type f -name "*.xtb" | sort -V); do
    electronic_energy=$(grep "TOTAL ENERGY" "$xtb_output" | awk '{print $4}')
    enthalpy_energy=$(grep "TOTAL ENTHALPY" "$xtb_output" | awk '{print $4}')
    free_energy=$(grep "TOTAL FREE ENERGY" "$xtb_output" | awk '{print $5}')

    echo "$xtb_output","$electronic_energy","$enthalpy_energy","$free_energy" >> energies.csv
done

cat energies.csv

