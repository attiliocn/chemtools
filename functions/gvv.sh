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

    for file in $@; do
        file_extension=$(echo "$file" | awk -F . '{print $NF}')
        file_basename=${file/."$file_extension"/}
        
        if [[ $file_extension == 'xyz' || $force_conversion == true ]]; then
            obabel -ixyz $file -omol > "$temp_directory"/"$file_basename".mol
        elif [[ $file_extension == 'out' ]]; then
            echo "ORCA output. The output will be converted using OFakeG"
            OfakeG $file
            mv "$file_basename"_fake.out "$temp_directory"/"$file_basename".log
        else
            cp $file "$temp_directory"
        fi
    done
    gv "$temp_directory"/* &
  }