#!/bin/bash

for xyz_file in $@; do
    bn="${xyz_file%.*}"
    input_file="$bn".mol
    obabel -ixyz $xyz_file -omol > $input_file 2> /dev/null
cat $input_file
done
