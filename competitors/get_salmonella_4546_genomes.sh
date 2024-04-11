#!/usr/bin/bash

mkdir data
cd data
wget https://zenodo.org/records/1323684/files/Salmonella_enterica.zip\?download\=1 -O Salmonella_enterica.zip
unzip Salmonella_enterica.zip
find Salmonella_enterica/Genomes/ -type f > salmonella_genomes.list