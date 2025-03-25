# kaminari
雷 - kaminari (thunder/lightning)

Kaminari is a *minimizer-based* approximate index for large-scale matching and color queries, in the same vein as [Fulgor](https://github.com/jermp/fulgor), but powered by [PTHash](https://github.com/jermp/pthash).

**A pre-print describing how the index works can be found [here (TODO)]().**

### Table of contents
* [Dependencies](#dependencies)
* [Compiling the code](#compiling-the-code)
* [Tools](#tools)
* [Example](#Demo)

Dependencies
------------

#### zlib

If you do not have `zlib` installed, you can do

    sudo apt-get install zlib1g

if you are on Linux/Ubuntu, or

    brew install zlib

if you are using macOS.

Compiling the code
------------------

The code is tested on Linux with `gcc` and on MacOS with `clang`.
To build the code, [`CMake`](https://cmake.org/) is required.

First clone the repository with

    git clone --recursive https://github.com/yhhshb/kaminari.git

for HTTPS or 

    git clone --recursive git@github.com:yhhshb/kaminari.git

for SSH.

If you forgot `--recursive` when cloning, do

    git submodule update --init --recursive

before compiling.

To compile the code for a release environment (see file `CMakeLists.txt` for the used compilation flags), it is sufficient to do the following, within the parent `kaminari` directory:

    mkdir build
    cd build
    cmake ..
    make -j

For a testing environment, use the following instead:

    mkdir debug_build
    cd debug_build
    cmake .. -D CMAKE_BUILD_TYPE=Debug -D KAMINARI_USE_SANITIZERS=On
    make -j

Tools
-----

There is one executable called `kaminari` after the compilation, which can be used to run a tool.
Run `./kaminari` to see a list of available tools.

	Usage: ./kaminari [--help] [--version] {build, query}

    Optional arguments:
    -h, --help     shows help message and exits 
    -v, --version  prints version information and exits 

    Subcommands:
    build         Build a kaminari index from a list of datasets 
    query         Query a kaminari index 



### Build 
`./kaminari build` allows you to build an index. `./kaminari build -h` displays a reminder of differents parameters that are mentionned below.

 + `-i` or `--input-list` is required. The list of files to index. You can provide multiple fasta/fastq (compressed or not) or one single file of files (which name extension should be ".list").
 + `-o` or `--output-filename`. The path and filename of the created index, base filename is ".". We encourage file extension ".kaminari". Default: index.kaminari
 + `-k`. The length of the k-mers (words of length k) that are indexed and queried. Reminder : only the shortest (according to murmurhash) minimizer (m-mer) is indexed for each k-mer. k=31 is a classic value. Default: 31
 + `-m`. The length of the minimizers (see k parameter). m=19 is a classic value. Default: 19
 + `-a` or `--canonical`. Consider canonical minimizers or not. DNA has two strands with the same information in reverse-complementary forms. Using canonical minimizers means choosing the same representation for both strands, effectively reducing redundancy and factoring the information. Default: OFF
 + `-b` or `--bit-check`. The number of bits used as an additionnal fingerprint for each ColorSetID. Reduce the number of alien k-mers in a negative query but add `b` bits for each unique minimizer indexed in kaminari. Default: 1 
 + `--metagenome`. Modifies the strategy of encoding colors during merging step. Recommanded in case of sparse colors, i.e. very complex data (metagenomics) combined with a large collection (>1000 metagenomic docs) Default: OFF
 + `-d` or `--tmp-dir`. The directory where temporary files are stored. Files are not deleted in case of manual interruption of the execution. Default: "."
 + `-g` or `--max-ram`. A maximum amount of memory allocated for kaminari in GB. Kaminari will try to use the maximum amount given by the user for better performances. Default: 4
 + `-t` or `--threads`. The number of threads allocated for kaminari. Kaminari will try to use them all for better performances. Default: 1
 + `-s` or `--seed`. The seed for murmurhash hash function used for determining minimizers. Two indexes with different seeds will have different results but this parameter should not be changed unless you know what you are doing. Default: 42
 + `-c` or `--pthash-constant`. The PTHash (MPHF) build constant, see pthash paper for details. Higher `c` leads to slower construction but more space efficient MPHF. Default: 4
 + `-v` or `--verbose`: The output verbosity level. Goes from 1 (global steps) to 5 (debug). Default: 0


 ### Query 
`./kaminari query` allows you to query an index. `./kaminari query -h` displays a reminder of differents parameters that are mentionned below.

 + `-x` or `--index` is required. The filepath to the index to be loaded and queried. Usually ends with ".kaminari".
 + `-i` or `--input-list` is required. List of fasta filenames to be queried. If only one file ending with ".list" is provided, it is assumed to be a file of files. Each sequence in the files is going to be queried and the header will be the name of the answer in the output.
 + `-o`. The filename of the output of the queries. Each line is the answer to a query. Answers are in the form of `name \t number_of_docs \t doc1 doc2 doc3` if `ranking` is OFF, else `name \t number_of_docs \t (doc1, count) (doc2, count) (doc3, count)`. Doc *n* corresponds to the *n*th doc provided for building the index. Default: "kaminari_results.txt"
 + `-r` or `--ratio`. The ratio of k-mers needed for a doc to be selected in the answer. Example: if my doc 5 contains (according to the index) 85 k-mers of my query of size 100 k-mers, it is selected. If `ranking` is ON, (5, 85) will be reported. (Reminder: with `k`=31, a query Q has |Q|-`k`+1 k-mers). Default: 1.0 (100% of the kmers)
 + `--no-ranking`. Kaminari returns the ranking by default (i.e. (doc, count) in the result). This can be turned OFF using this flag to have the set of answers unsorted.
 + `-d` or `--tmp-dir`. The directory where temporary files are stored. Files are not deleted in case of manual interruption of the execution. Default: "."
 + `-g` or `--max-ram`. A maximum amount of memory allocated for kaminari in GB. Kaminari will try to use the maximum amount given by the user for better performances. Default: 4
 + `-t` or `--threads`. The number of threads allocated for kaminari. Kaminari will try to use them all for better performances. Default: 1
 + `-v` or `--verbose`: The output verbosity level. Goes from 1 (global steps) to 5 (debug). Default: 0 




Demo
----

This short demo shows how to index the 10-genome collection in the folder `example`. The collection can be downloaded by running the script `download_test_datasets.sh` which generates `example/data/salmonella10/`.

From `kaminari/example`, run

    ./download_test_datasets.sh
    find $(pwd)/data/salmonella10/* > salmonella_10_filenames.list

We will use the standard values k = 31 and m = 19.

`kaminari` can take a single text file listing the inputs (one input file per line) as its `-i` option, for instance from the `build` directory: 

    ./kaminari build -i ../example/salmonella_10_filenames.list -o ../example/salmonella10.kaminari -k 31 -m 19 -a -d /tmp -g 4 -t 1 -v 1
    
in order to build an index that will be serialized to the file `example/salmonella10.kaminari`.

Now that the index is generated, we can query the presence of a sequence in the documents.

From `kaminari/build`, run

    ./kaminari query -x ../example/salmonella10.kaminari -i ../example/two_queries.fasta -o ../example/two_queries_result.txt -r 0.8 -d /tmp -g 4 -t 4 -v 1
    cat ../example/two_queries_result.txt

You can check the result in `example/two_queries_result.txt` and see that the first sequence is mainly present in docs 0, 3, 7 and 8 which corresponds to the 1st, 4th, 8th and 9th lines in `example/salmonella_10_filenames.list` (start from 0). The second sequence is reported for all 10 documents, but note how only documents 0,4,7,8 have the maximal count. The other ones share less k-mers with the sequence, but still over 80% (`-r 0.8`), so we report them.

