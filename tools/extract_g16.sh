#! /bin/bash

if [[ -d results_extracted ]]; then
    echo "The folder results_extracted already exists."
    echo "Delete this folder and try again"
    exit
else
    mkdir results_extracted 
    find . -mindepth 2 -type f \( -name "*.com" -o -name "*.log" \) | xargs -I{} cp -v {} results_extracted
fi




