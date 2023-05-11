#!/bin/bash

if [ $# -ne 2 ]; then
    echo "Usage: $0 <tar_file_name> <folder_name>"
    exit 1
fi

tar_file_name=$1
folder_name=$2
folder_name=${folder_name/\//}

# Check if the specified folder exists
if [ ! -d "$folder_name" ]; then
    echo "Error: Folder '$folder_name' does not exist."
    exit 1
fi

# Check if the specified tar file already exists
if [ -e "$tar_file_name".tar.gz ]; then
    echo "Error: Tar file '$tar_file_name' already exists."
    exit 1
fi

# Compress the folder and rename it inside the tar file
tar -czvf "$tar_file_name".tar.gz "$folder_name" --transform "s/^"$folder_name"/"$tar_file_name"/" --show-transformed-names

echo "Compression complete. Tar file '$tar_file_name' created with the renamed folder."
