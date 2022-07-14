#!/bin/bash

for xyz_file in $@; do
    bn="${xyz_file%.*}"
    input_file="$bn".com
    chk_file="$bn".chk

    echo -e "%mem=16gb\n%nprocs=8\n%chk="$chk_file"" > $input_file
    obabel -ixyz $xyz_file -ocom >> $input_file 2> /dev/null
    sed -i  '/!Put/d' $input_file

    sed -i 's/^#/# insert_route/' $input_file
    sed -i -E "s/^0 \+1/0 1/" $input_file

done
