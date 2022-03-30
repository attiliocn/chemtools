#!/bin/bash

grep "Normal termination" *.log -m 1 | sed 's/:/ /' | awk '{print $1}' | sort -V > normal.temp
ls *.log | sort -V  > log.temp

diff normal.temp log.temp > non_finished.temp

grep "^>" non_finished.temp | sed 's/> //g' > temp.temp
mv temp.temp non_finished.temp

for i in $(cat non_finished.temp); do
	last_mod=$(stat -c %y "$i")
	spacer=$'\t'
	echo "$i""$spacer""$spacer""$last_mod"

done

rm normal.temp log.temp non_finished.temp
