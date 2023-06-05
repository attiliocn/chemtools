#!/bin/bash
function gvv(){
    temp_directory="/tmp/geometry"
    
    if [[ -d "$temp_directory" ]]; then
        rm -rf "$temp_directory"
        mkdir -p "$temp_directory"
    else
        mkdir -p "$temp_directory"
    fi

    if [[ $1 == 'conv' ]]; then
        force_conversion=true
        shift
    else
        force_conversion=false
    fi

    for filepath in $@; do
        filepath=$(realpath "$filepath")
        filename=$(basename "$filepath")
        file_extension=$(echo "$filename" | awk -F . '{print $NF}')
        converted_filename=""$temp_directory"/${filename/"$file_extension"/mol}"
        if [[ $file_extension == 'xyz' || $force_conversion == true ]]; then
            obabel -ixyz $filepath -omol > $converted_filename
        else
            cp $filepath "$temp_directory"
        fi
    done
    gv "$temp_directory"/*
  }

#gvv $@
