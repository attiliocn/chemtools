#!/bin/bash

for xyz_file in $@; do
    bn="${xyz_file%.*}"

    cp -v $xyz_file "$bn"_rv.xyz
    cp -v $xyz_file "$bn"_fw.xyz
    rm -v $xyz_file
done
