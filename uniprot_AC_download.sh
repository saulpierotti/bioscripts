for i in `cat $1`; do
	wget https://www.uniprot.org/uniprot/$i.fasta
done
