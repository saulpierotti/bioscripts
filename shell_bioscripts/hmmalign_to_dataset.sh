#!/bin/bash
#
# Author: Saul Pierotti
# Mail: saulpierotti.bioinfo@gmail.com
# Last updated: 17/04/2020
#
# Takes a hmmsearch output file in tablular format as first argument and prints
# to STDOUT a simplified version of it.
# It also keeps only the lowest E-value line for each entry.
#
# The output is a simplified version containing only the query UniProt ID
# and the associated E-value.
# 
# It is possible to give an optional second argument that will be appended
# at the end of each line. This is useful for instance for differentiating
# classes for machine learning methods.
# I am taking the local E-value for now (field $8 in awk). This can be easily changed to the global one (field $5).

hmmfile="$1"
class=$2

if [ -f "$hmmfile" ]; then
	cat $hmmfile|grep -v "^#"|awk -v class=$class '{eval=$8;split($1,v,"|"); print v[2],eval,class}'\
	|awk '{if (!($1 in out) ||$2 < evalue[$1]) {out[$1]=$0; evalue[$1]=$2}} END{for(key in out){print out[key]}}'|sort -k 1

else
	echo "No input file given. This script requires an hmmsearch output file in tabular format as first argument."
fi
