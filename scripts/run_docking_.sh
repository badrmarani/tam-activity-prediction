#!/bin/bash

pmp="CC1=NC=C(C(=C1O)CN)COP(=O)(O)O"
enz_dir="enzymes"
sub_dir="enzymes"

num_cores=$(grep -c ^processor /proc/cpuinfo)

# grep ATOM 3ERK.pdb > rec.pdb
# obabel rec.pdb -O rec.pdb
for file in "$enz_dir"/*;
do
    # Skip directories
    if [[ -d "$file" ]]; then
        continue
    fi

    if grep -q ".pdb" <<< $file; then
        filename=$(echo $file | sed 's/\.[^.]*$//')
        grep ATOM "$(readlink -e $file)" > "$filename".temp.pdb
        obabel -i $(readlink -e "$filename".temp.pdb) -O $(readlink -e "$filename".temp.pdb)
    fi
done
