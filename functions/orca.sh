#!/bin/bash
function orca.opt(){
    grep -li 'hurray' *.out | while read -r f; do echo mv $f converged/; done
}
