#!/bin/bash

# files passed through argv are set to files variable
files=$@

for i in $files; do
  echo "convertendo $i"
  obabel -ig09 $i -oxyz >> ${i%.*}.xyz 2> /dev/null
  echo ""
done 
