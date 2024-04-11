#!/usr/bin/bash

#DOC https://github.com/seqan/raptor

# assuming this folders tree

#   competitors/
#       install_softs.sh
#       data/
#           salmonella_genomes.list
#           queries.fasta
#           genomes/
#               genome1.fasta
#               genome2.fasta
#       indexes/
#       cobs/
#       PAC/
#       RAMBO/
#       competitors_env/


list_of_files="../data/salmonella_genomes.list"
queries_file="../data/queries.fasta"

k=31
fpr=0.05

#INDEX
#from run_scripts/ folder
raptor build --input $list_of_files --output ../indexes/salmonella_index.raptor --kmer $k --fpr $fpr

#QUERIES
raptor search --index ../indexes/salmonella_index.raptor --query $queries_file --output ../indexes/raptor_output.out
