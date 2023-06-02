function gvv(){
    if [[ $1 == 'conv' ]]; then
        force_conversion=true
        shift
    else
        force_conversion=false
    fi

    for file in $@; do
        file_extension="${file#*.}"
        converted_filename="/tmp/${file/"$file_extension"/.mol}"
        if [[ $file_extension == 'xyz' || $force_conversion ]]; then
            obabel -ixyz $file -omol > $converted_filename
            gv $converted_filename
            rm $converted_filename
        else
          gv $file
        fi
    done
  }
