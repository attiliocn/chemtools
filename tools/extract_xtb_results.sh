#!/bin/bash

if [ -d results_extracted/ ]; then
    echo 'WARNING: Impossible to extract the results. The results_extracted folder already exists!'
    exit
else
    mkdir results_extracted
    mkdir results_extracted/logs
fi

for i in $(ls -d */); do
    i=${i/\//}
    if test -f "$i"/crest_clustered.xyz; then
        echo "Found CREGEN clustered ensemble!"
        cp "$i"/crest_clustered.xyz results_extracted/"$i"_confs.xyz
        cp "$i"/cregen.output results_extracted/logs/"$i".cregen
        cp "$i"/crest.output results_extracted/logs/"$i".crest
        cp "$i"/cluster.order results_extracted/logs/"$i".clusters
    
    elif test -f "$i"/crest.output; then
        echo "Standard CREST output!"
        cp "$i"/crest_conformers.xyz results_extracted/"$i"_confs.xyz
        cp "$i"/crest.output results_extracted/logs/"$i".crest
        cp "$i"/xtbopt.log results_extracted/logs/"$i".log

    elif test -f "$i"/xtbscan.log; then
        echo "XTB output (scan)"
        cp "$i"/xtbopt.xyz results_extracted/"$i".xyz
        cp "$i"/xtbscan.log results_extracted/logs/"$i".scan
        cp "$i"/xtbopt.log results_extracted/logs/"$i".log
        cp "$i"/xtb.output results_extracted/logs/"$i".xtb

    elif test -f "$i"/xtblast.xyz; then
        echo "XTB output (non terminated)"
        cp "$i"/xtblast.xyz results_extracted/"$i"_error.xyz
        cp "$i"/xtbopt.log results_extracted/logs/"$i"_error.log
        cp "$i"/xtb.output results_extracted/logs/"$i"_error.xtb

    elif test -f "$i"/xtb.output; then
        echo "Standard XTB output"
        cp "$i"/xtbopt.xyz results_extracted/"$i".xyz
        cp "$i"/xtbopt.log results_extracted/logs/"$i".log
        cp "$i"/xtb.output results_extracted/logs/"$i".xtb
  
    fi
done

echo 'Done!'
