#!/bin/bash
#
# Author: Saul Pierotti
# Mail: saulpierotti.bioinfo@gmail.com
# Last updated: 17/04/2020
#
# Extracts a list of UniProt IDs from a Uniprot fasta file

fasta=$1

grep ">" $fasta|cut -d "|" -f 2|sort
