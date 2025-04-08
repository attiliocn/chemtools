#!/bin/bash

echo "filename,E(el),H,G,ifreq quantity,ifreq" > energies.csv

for orca_output in $(find . -maxdepth 1 -type f -name "*.out" | sort -V); do
    electronic_energy=$(grep "FINAL SINGLE POINT ENERGY" "$orca_output" | tail -1 | awk '{print $NF}')
    enthalpy_energy=$(grep "Total enthalpy" "$orca_output" | tail -1 | awk '{print $(NF-1)}')
    free_energy=$(grep "Final Gibbs free energy" "$orca_output" | tail -1 | awk '{print $(NF-1)}')

    last_spectra_line=$(grep -n "VIBRATIONAL FREQUENCIES" "$orca_output" | tail -1 | awk -F ":" '{print $1}')
    if [[ -n $last_spectra_line ]]; then
        tail -n +"$last_spectra_line" "$orca_output" > spectra.tmp
        grep "imaginary mode" spectra.tmp > imaginary.tmp
        imaginary_frequencies=$(awk '{print $2}' imaginary.tmp | tr '\n' ';')
        imaginary_frequencies_qty=$(wc -l imaginary.tmp | awk '{print $1}')
        rm spectra.tmp imaginary.tmp
    else
        imaginary_frequencies=''
        imaginary_frequencies_qty=''
    fi

    echo "$orca_output","$electronic_energy","$enthalpy_energy","$free_energy","$imaginary_frequencies_qty","$imaginary_frequencies" >> energies.csv
done

cat energies.csv

