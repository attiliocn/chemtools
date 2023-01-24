#!/bin/bash

reference_file="$1"
shift

for file in $@; do
    if [[ "$file" != "$reference_file" ]]; then
        temp_file="temp_"$(echo $(date "+%M%S%N") | sha1sum | cut -c1-5)".xyz"
        obrms -m "$reference_file" "$file" -o $temp_file &>> rmsd.txt
        echo "" >> rmsd.txt
        mv $temp_file $file
    fi
done