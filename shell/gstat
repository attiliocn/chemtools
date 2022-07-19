#!/bin/bash

function grep_gstat() {
for i in $(cat "$1"); do
    if [ -z $3 ] || [ $3 = "fast" ]; then
        tail -n 10 "$i" | grep -m 1 "$2" | xargs -I{} echo $i: {}
    elif [ $3 = "full" ]; then
        grep -H "$2" "$i"
    fi
done
}

required_stat=$1
search_in=$2
temp_file="/tmp/$(timestamp)_tempfile.txt"

find . -type f -name "*.log" | sort -V > $temp_file
sed -i '/slurm/d' $temp_file

if [ "$required_stat" = 'error'  ]; then
    grep_string='Error termination'
    grep_gstat "$temp_file" "$grep_string" "$search_in"
elif [ "$required_stat" = 'normal' ]; then
    grep_string='Normal termination'
    grep_gstat "$temp_file" "$grep_string" "$search_in"
elif [ "$required_stat" = 'time' ]; then
    grep_string='Elapsed time'
    grep_gstat "$temp_file" "$grep_string" "$search_in"
elif [ "$required_stat" = 'steps' ]; then
    grep_string='Step number'
    for i in $(cat "$temp_file"); do
        echo $i $(grep "$grep_string" "$i" | wc -l)
    done
else
    grep_string=$required_stat
    grep_gstat $temp_file $grep_string $search_in
fi

rm $temp_file

exit
