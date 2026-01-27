#!/bin/bash

for f in *.xyz; do
cat << EOF > tmpfile
! pal8
! insert_route

%maxcore 3000
%insert_section

* xyz 0 1
EOF
    tail -n +3 $f >> tmpfile
    echo "*" >> tmpfile

    mv tmpfile ${f/xyz/inp}

done
