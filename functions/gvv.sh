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
        # elif [[ $file_extension == 'out' && $force_conversion == false ]]; then
        #     echo "ORCA output. The optimized geometry will be extracted manually"
        #     outxyz.sh $filepath
        #     obabel -ixyz "${filepath/."$file_extension"/.xyz}" -omol > $converted_filename
        elif [[ $file_extension == 'out' ]]; then
            echo "ORCA output. The output will be converted using OFakeG"
            OfakeG "$filepath"
            mv ${filepath/."$file_extension"/_fake."$file_extension"} $temp_directory/${filename/"$file_extension"/log}
        else
            cp $filepath "$temp_directory"
        fi
    done
    gv "$temp_directory"/* &
  }

#gvv $@
