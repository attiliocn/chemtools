#!/bin/bash
mkdir results_extracted && find . -mindepth 2 -type f \( -name "*.out" -o -name "*.inp" \) -exec cp -v {} results_extracted/ \;