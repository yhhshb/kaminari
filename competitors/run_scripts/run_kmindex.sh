#!/usr/bin/bash

#DOC https://tlemane.github.io/kmindex/construction/

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

#create good format list of files
ls -1 ../data/Salmonella_enterica/Genomes/*  | sort -n -t/ -k 2 |awk '{print ++count" : "$1}' > ../data/salmonella_genomes.txt

list_of_files="../data/salmonella_genomes.txt"
queries_file="../data/queries.fasta"

k=31
bloom_size=1000000 #power of 2 above this number 

#INDEX
#from run_scripts/ folder
kmindex build --fof $list_of_files --run-dir ../indexes/kmindex_index --index ../indexes --register-as salmonella_index.kmindex --hard-min 1 --kmer-size $k --bloom-size $bloom_size

#QUERIES
kmindex query --index ../indexes --fastx $queries_file --names salmonella_index.kmindex --zvalue 3 --threshold 0
