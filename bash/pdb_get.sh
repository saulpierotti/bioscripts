#!/bin/bash
#
# Author: Saul Pierotti
# Mail: saulpierotti.bioinfo@gmail.com
# Last updated: 16/05/2020
#
# Takes as first argument a file containing a list of PDB IDs, or a single pdb id directly
# Downloads all of them in the current folder in pdb format

ext="cif"

case $2 in
"--pdb")
    ext="pdb"
    ;;
"*")
    echo "Extension not given. Defaulting to mmcif."
    ;;
esac

if [ -f "$1" ]; then
    while read -r line; do
        wget -nc "https://files.rcsb.org/download/$line.$ext"
    done <"$1"
else
    echo "No input file given. Trying to interpret the argument as a PDB ID"
    wget -nc "https://files.rcsb.org/download/$1.$ext"
fi
