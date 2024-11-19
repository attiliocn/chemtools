#!/bin/bash

# set search variables based on arg1 and arg2
# arg1 -> the software
# arg2 -> file extension of the log file
if [[ $1 == "xtb" ]]; then
    search_string="normal termination of xtb"
    extension=".xtb"
elif [[ $1 == "crest" ]]; then
    search_string="CREST terminated normally"
    extension=".crest"
else
    printf "Compatible outputs are xtb,crest\n"
    printf "Exiting\n"
    exit
fi

# set extension to arg2 if provided
if [ ! -z $2 ]; then
    extension="$2"
fi

counts=$(grep "$search_string" -c *"$extension")

for count in $counts; do
    file=$(echo $count | awk -F ":" '{print $1}')
    count=$(echo $count | awk -F ":" '{print $2}')

    if [ $count -eq 0 ]; then
        printf "$file $count error\n"
    else
        printf "$file $count normal\n"
    fi
done