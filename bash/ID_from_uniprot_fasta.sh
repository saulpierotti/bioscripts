#!/bin/bash
#
# Author: Saul Pierotti
# Mail: saulpierotti.bioinfo@gmail.com
# Last updated: 17/04/2020
#
# Extracts a list of UniProt IDs from a Uniprot fasta file

fasta=$1

if [ -f "$fasta" ]; then
	grep ">" $fasta|cut -d "|" -f 2|sort
else
	echo "No input file given. This script requires an uniprot fasta file as first argument."
fi
