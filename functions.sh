function gvv(){
    for file in $@; do
        file_extension="${file#*.}"
        if [[ $file_extension == 'xyz' ]]; then
            obabel -ixyz $file -omol > /tmp/"${file/.xyz/.mol}"
            gv "/tmp/${file/.xyz/.mol}"
            rm "/tmp/${file/.xyz/.mol}"
        else
          gv $file
        fi
    done
  }
