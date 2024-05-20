#!/bin/bash

smarts=$1
reference_file=$2
shift
shift

for file in $@; do
    if [[ "$file" != "$reference_file" ]]; then
        obfit $smarts $reference_file $file 2> /dev/null 1> /tmp/geometry.tmp
        mv /tmp/geometry.tmp $file
    fi
done