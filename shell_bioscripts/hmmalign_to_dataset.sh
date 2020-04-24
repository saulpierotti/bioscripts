#!/bin/bash
#
# Author: Saul Pierotti
# Mail: saulpierotti.bioinfo@gmail.com
# Last updated: 23/04/2020
#
# Takes a hmmsearch output file in tablular format as first argument and prints
# to STDOUT a simplified version of it.
# The hmmsearch output seemsto be uniquem but just for safety this script also
# keeps only the highest bitscore line for each entry.
#
# The output is a simplified version containing only the query UniProt ID
# and the associated bitscore.
# 
# It is possible to give an optional second argument that will be appended
# at the end of each line. This is useful for instance for differentiating
# classes for machine learning methods.
# By default it outputs the global E-value in the second field.
# If I want the local corealue put 1 as fourth argument.

hmmfile="$1"
class="$2"
is_local="$3"

if [ -f "$hmmfile" ]; then
	cat $hmmfile|grep -v "^#"|awk -v class=$class -v is_local=$is_local '{bitscore_loc=$9;bitscore_glob=$6;split($1,v,"|");\
		if (is_local)\
			{bitscore=bitscore_loc}\
		else\
			{bitscore=bitscore_glob};\
		print v[2],bitscore,class}'\
	|awk '{if (!($1 in out) ||$2 < bitscore[$1]) {out[$1]=$0; bitscore[$1]=$2}} END{for(key in out){print out[key]}}'|sort -k 1

else
	echo "No input file given. This script requires an hmmsearch output file in tabular format as first argument."
fi
