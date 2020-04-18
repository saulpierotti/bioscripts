#!/bin/bash
#
# Author: Saul Pierotti
# Mail: saulpierotti.bioinfo@gmail.com
# Last updated: 17/04/2020
#
# Takes as first argument a file containing a list of PDB IDs
# Downloads all of them in the current folder in pdb format

if [ -f "$1" ]; then
	for i in `cat $1`; do
		wget https://files.rcsb.org/download/$i.pdb
	done
else
	echo "No input file given. This script requires a file containing a list of PDB IDs as first argument"
fi
