#!/bin/bash

for i in $@; do
    basename=${i%.*}
    convert $i -trim $i
done
