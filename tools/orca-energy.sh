#!/bin/bash

echo "filename,electronic energy,enthalphy,free energy" > energies.csv

for orca_output in $(find . -maxdepth 1 -type f -name "*.out" | sort -V); do
    electronic_energy=$(grep "FINAL SINGLE POINT ENERGY" "$orca_output" | tail -1 | awk '{print $NF}')
    enthalpy_energy=$(grep "Total enthalpy" "$orca_output" | tail -1 | awk '{print $(NF-1)}')
    free_energy=$(grep "Final Gibbs free energy" "$orca_output" | tail -1 | awk '{print $(NF-1)}')

    echo "$orca_output","$electronic_energy","$enthalpy_energy","$free_energy" >> energies.csv
done

cat energies.csv

