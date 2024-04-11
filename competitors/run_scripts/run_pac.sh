#!/usr/bin/bash

#DOC :  https://github.com/Malfoy/PAC

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
bloom_size=1000000 #power of 2 above this number 

#INDEX
#from run_scripts/ folder
../PAC/PAC -f $list_of_files -d ../indexes/pac_index -k $k -b $bloom_size

#QUERIES
../PAC/PAC -l ../indexes/pac_index -q $queries_file