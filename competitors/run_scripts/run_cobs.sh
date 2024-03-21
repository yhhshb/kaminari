#!/usr/bin/bash

#DOC : https://github.com/bingmann/cobs

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
../cobs/build/src/cobs compact-construct $list_of_files ../indexes/salmonella_index.cobs

#QUERIES
../cobs/build/src/cobs query -i ../indexes/salmonella_index.cobs -f $queries_file