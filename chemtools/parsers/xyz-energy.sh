#!/bin/bash

echo "file,header" > headers.csv
for f in $@
do
    headers=$(grep -P "^-[0-9]+\.[0-9]" $f)
    i=1
    for h in $headers
    do
        echo "${f/.*/}"_"$i".xyz,$h >> headers.csv
        i=$((i+1))
    done
done