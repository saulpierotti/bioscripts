#!/bin/bash
#
# Author: Saul Pierotti
# Mail: saulpierotti.bioinfo@gmail.com
# Last updated: 16/05/2020
#
# Takes as first argument a file containing a list of PDB IDs, or a single pdb id directly
# Downloads all of them in the current folder in pdb format

if [ -f "$1" ]; then
	for i in `cat $1`; do
		wget https://files.rcsb.org/download/$i.cif
	done
else
	echo "No input file given. Trying to interpret the argument as a PDB ID"
	wget https://files.rcsb.org/download/$1.cif
fi
