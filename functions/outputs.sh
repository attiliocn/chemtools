#!/bin/bash
function orca.opt(){
    grep -li 'hurray' *.out | while read -r f; do echo mv ${f/.out/.*} converged/; done
}

function gaussian.opt(){
    grep -l 'Optimization completed.' *.log | while read -r f; do echo mv ${f/.log/.*} converged/; done
}
