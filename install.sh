#!/bin/bash

for f in "$CHEMTOOLS_DIR"/chemtools/functions/*.sh; do
    source $f
done

for f in "$CHEMTOOLS_DIR"/chemtools/run-softwares/*.sh; do
    source $f
done

export PATH=""$PATH":"$CHEMTOOLS_DIR"/chemtools"
export PATH=""$PATH":"$CHEMTOOLS_DIR"/chemtools/tools"
export PATH=""$PATH":"$CHEMTOOLS_DIR"/chemtools/parsers"
export PATH=""$PATH":"$CHEMTOOLS_DIR"/chemtools/converters"
export PATH=""$PATH":"$CHEMTOOLS_DIR"/chemtools/shell"

