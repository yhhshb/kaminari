#!/usr/bin/bash

#DOC : https://github.com/RUSH-LAB/RAMBO

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


#INDEX
#from run_scripts/ folder
../RAMBO/bin/rambo build $list_of_files -o ../indexes/rambo_index

#QUERIES
../RAMBO/bin/rambo query $queries_file --database ../indexes/rambo_index