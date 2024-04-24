#!/bin/bash

files=$@

for i in $files; do 

	tr -d '\15\32' < $i > temp
	mv temp $i

	head -n $(grep -n "^$" $i | head -n 3 | tail -1 | sed 's/://') $i > temp
	mv temp $i
done
