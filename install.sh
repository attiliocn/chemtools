#!/bin/bash

for f in "$CHEMTOOLS_DIR"/functions/*.sh; do
    source $f
done

for f in "$CHEMTOOLS_DIR"/run-softwares/*.sh; do
    source $f
done

export PATH=""$PATH":"$CHEMTOOLS_DIR"/tools"
export PATH=""$PATH":"$CHEMTOOLS_DIR"/converters"
export PATH=""$PATH":"$CHEMTOOLS_DIR"/shell"

