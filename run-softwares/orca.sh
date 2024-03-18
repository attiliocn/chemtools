function orca.run() {
    input_file=$1
    echo Running "$input_file"
    $ORCA/orca $input_file 1> ${input_file/inp/out} &
}