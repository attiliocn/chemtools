#!/bin/bash

smarts=$1
mol_fix=$2
mol_mov=$3

obfit $smarts $mol_fix $mol_mov 2>&1 >/dev/null | head -n 1 | awk '{print $2}'
