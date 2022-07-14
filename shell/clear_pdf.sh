#!/bin/bash
echo $1
pdftk \"$1\" output - uncompress | sed '/^\/Annots/d' | pdftk - output out.pdf compress
