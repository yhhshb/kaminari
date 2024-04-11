#!/usr/bin/bash

#DOC https://github.com/kalininalab/metaprofi

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
fpr=0.05

#INDEX
#from run_scripts/ folder
metaprofi build $list_of_files /path/to/config.yml

#QUERIES
metaprofi search_index -f $queries_file -i nucleotide -t 50
