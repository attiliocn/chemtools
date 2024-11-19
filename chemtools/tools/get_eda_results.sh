#!/bin/bash

output_files=$@

echo "output,charge_transfer,electrostatic,polarization,exchange,fra1_def,frag2_def,frag1_se,frag2_se" > eda_energies.csv

for output in $output_files; do
    charge_transfer=$(grep "CT =" $output | awk '{print $NF}')
    electrostatic=$(grep "ES =" $output | awk '{print $NF}')
    polarization=$(grep "POL =" $output | awk '{print $NF}')
    exchange=$(grep "XC =" $output | awk '{print $NF}')
    frag1_def=$(grep "^ 1\. " $output | sed 's/[()]/ /g' | awk '{print $((NF-1))}')
    frag1_se=$(grep "^ 1\. " $output | sed 's/[()]/ /g' | awk '{print $NF}')
    frag2_def=$(grep "^ 2\. " $output | sed 's/[()]/ /g' | awk '{print $((NF-1))}')
    frag2_se=$(grep "^ 2\. " $output | sed 's/[()]/ /g' | awk '{print $NF}')

    grep "^ 1\. " $output | sed 's/[()]/ /g' | awk '{print $((NF-1))}' |sed 's/[()]/ /g'

    echo $output,$charge_transfer,$electrostatic,$polarization,$exchange,$frag1_def,$frag2_def,$frag1_se,$frag2_se >> eda_energies.csv
done