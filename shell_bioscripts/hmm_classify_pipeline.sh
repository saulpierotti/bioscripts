#!/bin/bash
#
# Author: Saul Pierotti
# Mail: saulpierotti.bioinfo@gmail.com
# Last updated: 01/05/2020
#
# Just a convenient wrapper for all the computation done in building an HMM from a seed alignment
# and classifying a positive and negative dataset with it.
# If sequences of the seed alignment are present in the dataset they will be removed.
# Just execute this file without any argument.
# It will source data from the folder `input`, put intermediate data in `data` and the final results in `output`.
#
# The seed alignment should be in ./input/seed.aln
# The positive dataset in ./input/positives.fasta.gz (also uncompressed is ok)
# The negative dataset in ./input/negatives.fasta.gz (also uncompressed is ok)

echo "Unpacking the data..."
cat ./input/seed.aln | awk '$1 ~ /^>/ {print $0} $1 !~ /^>/ {gsub(/-/,"");print $0}' > ./data/seed.fasta
gunzip ./input/positives.fasta.gz
gunzip ./input/negatives.fasta.gz

echo "Blasting the datasets against the seed..."
makeblastdb -in ./data/seed.fasta -dbtype prot
cat ./input/positives.fasta|blastp -db ./data/seed.fasta -out ./data/positives_against_seed.bltab -outfmt 6 -evalue 1e-03
cat ./input/negatives.fasta|blastp -db ./data/seed.fasta -out ./data/negatives_against_seed.bltab -outfmt 6 -evalue 1e-03
cat ./data/positives_against_seed.bltab|awk '{if ($3==100){print $1}}' > ./data/to_be_removed_from_positives.list
cat ./data/negatives_against_seed.bltab|awk '{if ($3==100){print $1}}' > ./data/to_be_removed_from_negatives.list

echo "Stripping seqeunces that are in the seed from the dataset..."
./script/strip_from_fasta.py ./input/positives.fasta ./data/to_be_removed_from_positives.list > ./data/positives_cleaned.fasta
./script/strip_from_fasta.py ./input/negatives.fasta ./data/to_be_removed_from_negatives.list > ./data/negatives_cleaned.fasta

echo "Creating the HMM..."
hmmbuild ./data/seed.hmm ./input/seed.aln

echo "Screening the positive dataset with hmmsearch..."
hmmsearch -Z 1 --max --tblout ./data/positives.hmmsearch.tbl ./data/seed.hmm ./data/positives_cleaned.fasta
echo "Screening the negative dataset with hmmsearch..."
hmmsearch -Z 1 --max --tblout ./data/negatives.hmmsearch.tbl ./data/seed.hmm ./data/negatives_cleaned.fasta

echo "Extracting useful data from the hmmsearch output..."
hmmalign_to_dataset.sh ./data/positives.hmmsearch.tbl 1 1 > ./data/positives_local.dat
hmmalign_to_dataset.sh ./data/negatives.hmmsearch.tbl 0 1 > ./data/negatives_local.dat
hmmalign_to_dataset.sh ./data/positives.hmmsearch.tbl 1 0 > ./data/positives_global.dat
hmmalign_to_dataset.sh ./data/negatives.hmmsearch.tbl 0 0 > ./data/negatives_global.dat

echo "Extracting IDs list from the datasets..."
ID_from_uniprot_fasta.sh ./data/positives_cleaned.fasta > ./data/positives.list
ID_from_uniprot_fasta.sh ./data/negatives_cleaned.fasta > ./data/negatives.list

echo "Adding missing IDs..."
add_missing_IDs.sh ./data/positives_local.dat  ./data/positives.list 10 1 > ./data/positives_local_complete.dat
add_missing_IDs.sh ./data/negatives_local.dat  ./data/negatives.list 10 0 > ./data/negatives_local_complete.dat
add_missing_IDs.sh ./data/positives_global.dat ./data/positives.list 10 1 > ./data/positives_global_complete.dat
add_missing_IDs.sh ./data/negatives_global.dat ./data/negatives.list 10 0 > ./data/negatives_global_complete.dat

echo "Merging positives and negatives..."
cat ./data/positives_local_complete.dat  ./data/negatives_local_complete.dat  > ./data/final_dataset_local.dat
cat ./data/positives_global_complete.dat ./data/negatives_global_complete.dat > ./data/final_dataset_global.dat

echo "Just testing that everything worked well before doing cross validation..."
model_stats.py ./data/final_dataset_local.dat 1e-03
model_stats.py ./data/final_dataset_global.dat 1e-03

echo "Doing cross validation in python..."
cp -l ./data/final_dataset_global.dat Global # I use the filename as class name for the plots in cross_val.py, so I make an hard link
cp -l ./data/final_dataset_local.dat Best\ domain
cross_val.py Global Best\ domain

echo "Final cleanup.."
rm Global # remove the hard links
rm Best\ domain
gzip ./input/positives.fasta
gzip ./input/negatives.fasta

echo "End of the script"
