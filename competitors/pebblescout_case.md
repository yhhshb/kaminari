
${toc}

* Open date: 18/07/2023
* Close date: 21/07/2023 (au moins 2 mois de calcul sur les 50 jeux de données)
* What: Tester PebbleScout sur les 50 jeux Tara QQSS
* Where:  cl1n051:/WORKS/expes_indexations/expe_PebbleScout
* Tool versions:


## Environnement
```
tmux
srun -p seqdigger --mem 512G --threads=16 --cpus-per-task=16 -w cl1n051 --pty bash

export TIME="\t%E real,\t%U user,\t%S sys,\t%K amem,\t%M mmem"


cd /WORKS/expe_PebbleScout
```

## Install
wget https://ftp.ncbi.nlm.nih.gov/pub/agarwala/pebblescout/v2.25_release.zip
(lien donné dans le papier)

Environnement:
```
mamba create -n 'env_pebblescout' python=3.9 sqlite=3 gcc=11 cmake gxx=11
mamba activate env_pebblescout
```

Modification de "software/pebblescout/Lookup/CmakeLists.txt"

```
target_link_libraries ( pebblesearch /home/symbiose/ppeterlo/.conda/envs/env_pebblescout/lib/libsqlite3.so )
#target_link_libraries ( pebblesearch /usr/local/sqlite/3.35.5/lib/libsqlite3.a )
target_link_libraries ( pebblesearch ${CMAKE_DL_LIBS})
target_include_directories(  pebblesearch PUBLIC /home/symbiose/ppeterlo/.conda/envs/env_pebblescout/include/ )
#target_include_directories(  pebblesearch PUBLIC /usr/local/sqlite/3.35.5/include )
```
 
Puis suivi le README:
```
Judy array code is freely available at https://sourceforge.net/projects/judy/
We downloaded file Judy-1.0.5.tar.gz in this directory. Do following to compile this code:
$ gunzip Judy-1.0.5.tar.gz
$ tar xvf Judy-1.0.5.tar
$ cd judy-1.0.5
$ ./configure
$ make
$ make check
$ cd ..

Software for Pebblescout is in pebblescout.zip file in this directory.
Unzippping this file creates a pebblescout directory.
Access to sqlite3 is expected.
Edit pebblescout/Lookup/CMakeLists.txt to point to the location of sqlite headers and library.
After that, do:
$ cd pebblescout
$ cp ../judy-1.0.5/src/obj/.libs/libJudy.a .
$ g++ -I ../judy-1.0.5/src/ -O2 -Wall -DSUBSAMPLESIZE=18 -o harvester25 -L . harvester25.cpp -l Judy
$ g++ -Wall -O2 -o indexbuilder indexbuilder.cpp
$ cmake Lookup -D CMAKE_BUILD_TYPE=RELEASE
$ make
```


## Build index tiny example
```
cd /WORKS/expes_indexations/expe_PebbleScout
for file in 11SUR1QQSS11 142SUR1QQSS11 143SUR1QQSS11
do
    zcat /WORKS/expes_indexations/expe_kmindex/data_per_station/${file}.fastq.gz | head -n 40000 | gzip > head_${file}.fastq.gz
done
ls -d $PWD/head* > tiny_fof.txt

#create the tiny_list.txt (cf in the directory)

#create the index
 /usr/bin/time sh run_test_one_tara.sh
#edit the the json files: adjust paths in db.json and db.with_suppressed.json to the directory where you built the database and check that file sizes for *.tr.bin and .vocab are correct in the corresponding json files under tranTablesSZ and vocabSZ, respectively. Entry for vocabFLDS should reflect metadata fields in the list.txt file.
```



## Build index full example
```

#create the list.txt :
cnt=0
while read file; do
        base_file=$(basename $file .fastq.gz)
    cnt=$((cnt + 1))
        echo  -e "${cnt}\t ${base_file}\t mandatory_metadata"
done < fof.txt > full_list.txt
ln -s full_list.txt list.txt

#create the index
 /usr/bin/time sh run_test_one_tara.sh
 
#edit the the json files: adjust paths in db.json and db.with_suppressed.json to the directory where you built the database and check that file sizes for *.tr.bin and .vocab are correct in the corresponding json files under tranTablesSZ and vocabSZ, respectively. Entry for vocabFLDS should reflect metadata fields in the list.txt file.
```


[run_test_one_tara.sh](:/42958a44d61741d8b931e666bf359b78)

21/07/2023 16:13 : killed - entre 20 et 22h pour harves le premier jeu de données. Il y en a 50...


## SEARCH
```
#create the query
zcat head_11SUR1QQSS11.fastq.gz | head -n 2 | tr "@" ">"  > query.fa

#make the search
/usr/bin/time ./pebblescout_v2.25/software/pebblescout/pebblesearch -f query.fa -m 2 -F "QueryID,SubjectID,%coverage,PBSscore,BioSample,Sample,Host" -i db.json -o score.out 2> score.log
```
-- 

Pierre Peterlongo
Inria Research Director
Head of the GenScale team

Webpage: http://people.rennes.inria.fr/Pierre.Peterlongo/
Visio: https://inria.webex.com/join/pierre.peterlongo
Tel : +33 (0)6 64 71 15 95

