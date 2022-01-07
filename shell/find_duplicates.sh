#!/bin/bash

# https://stackoverflow.com/questions/16276595/how-to-find-duplicate-filenames-recursively-in-a-given-directory-bash #

find $1 -type f -printf '%p/ %f\n' | sort -k2 | uniq -f1 --all-repeated=separate
