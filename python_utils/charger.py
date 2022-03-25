#!/usr/bin/env python3

import sys
from chemparser.tools.tools import electron_counter

#requires chemparser tool as simlink on the same folder

print(electron_counter(sys.argv[1:]))


#shell one liner
# for i in $(ls *.xyz | sort -V); do echo $i $(cat $i | awk '{print $1}' | tail -n +3 | xargs ./charger); done
