#!/bin/bash

for xyz_file in $@; do
    bn="${xyz_file%.*}"
    smi=$(obabel -ixyz $xyz_file -ocan 2> /dev/null | awk '{print $1}' )

    printf ""$bn"\t"$smi"\n"

done
