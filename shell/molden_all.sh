#!/bin/bash

num_entries=$#

j=$num_entries
for i in $@; do
    echo "Current molecule:"
    echo $i
    echo "Remaining molecules: "$j""
    molden $i &> /dev/null
    echo "Done"
    j=$(("$j"-1))
    echo ''
done

