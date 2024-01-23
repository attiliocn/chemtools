#!/bin/bash

if [[ ! -d results_extracted ]]; then
    mkdir results_extracted

    find . -mindepth 2 -type f \
        \( \
            -name "*.out" \
            -o -name "*.inp" \
            -o -name "*.hess" \
            -o -name "*_trj.xyz" \
        \) \
        -a -not \
        \( \
            -name "*smd*" \
            -o -name "*scf*" \
            -o -name "*atom*" \
            -o -name "*proc*" \
            -o -name "*appr*" \
        \) \
        -exec cp -v {} results_extracted/ \;
else
    echo "Error: The results_extracted directory already exists"
fi