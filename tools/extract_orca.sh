#!/bin/bash

if [[ -d results_extracted ]]; then
    mkdir results_extracted
    find . -mindepth 2 -type f \( -name "*.out" -o -name "*.inp" -o -name "*.hess" \) -a -not \( -name "*smd*" -o -name "*scf*" -o -name "*atom*" \) -exec cp -v {} results_extracted/ \;
else
    exit

