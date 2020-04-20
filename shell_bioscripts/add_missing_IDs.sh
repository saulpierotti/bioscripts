#!/bin/bash
#
# Author: Saul Pierotti
# Mail: saulpierotti.bioinfo@gmail.com
# Last updated: 17/04/2020
#
# It adds missing IDs to a dataset to be used for classification.
# Input must be sorted
# Takes in input as first argument a dataset for in the form
#
# ID score class
#
# Takes in input as second argument a list of IDs
# Takes in input as third argument a score to be used for the missing IDs
# Takes in input as fourth a class to be used for the missing IDs

dataset=$1
id_list=$2
score=$3
class=$4

cat $dataset
for id in $(cut -d " " -f 1 $dataset|comm -3 -1 - $id_list); do
	echo "$id $score $class"
done
