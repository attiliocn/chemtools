#!/bin/bash

for xyz_file in $@; do
    bn="${xyz_file%.*}"
    input_file="$bn".inp

    obabel -ixyz $xyz_file -oorcainp > $input_file 2> /dev/null
    sed -i  '/^#/d' $input_file

    sed -i 's/^! insert inline commands here/! insert_route\n\n%insert_section\n/' $input_file
    sed -i 's/! insert_route/! pal8\n! insert_route/' $input_file
    sed -i 's/%insert_section/%maxcore 5000\n%insert_section/' $input_file

cat $input_file

done
