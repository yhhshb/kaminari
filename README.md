# kaminari
雷 - kaminari (thunder/lightning)

Kaminari is a *minimizer-based* approximate index for large-scale matching and color queries, in the same vein as [Fulgor](https://github.com/jermp/fulgor), but powered by [PTHash](https://github.com/jermp/pthash).

**A pre-print describing how the index works can be found [here](https://www.biorxiv.org/content/10.1101/2025.05.16.654317v1).**

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
`./kaminari build` allows you to build an index. `./kaminari build -h` displays a reminder of the different parameters mentioned below.

* **`-i` or `--input-list`** (Required): The list of files to index. You can provide multiple FASTA/FASTQ files (compressed or not) or a single file of files (where the file extension must be ".list").
* **`-o` or `--output-dirname`**: The output directory name where the index will be created. Default: `kaminari_index`
* **`-k`**: The length of the k-mers (words of length k) that are indexed and queried. Reminder: only the shortest (according to murmurhash) minimizer (m-mer) is indexed for each k-mer. k=31 is a classic value. Default: `31`
* **`-m`**: The length of the minimizers (must be < k). m=19 is a classic value. For highly redundant data (e.g. a human pangenome), a higher value is advised (21 to 23). Default: `19`
* **`-n` or `--non-canonical`**: Turn OFF canonical minimizers. By default, Kaminari considers canonical minimizers (using the same representation for a DNA strand and its reverse-complement to reduce redundancy). Use this flag to disable this behavior and index sequences exactly as they appear. 
* **`-b` or `--bit-check`**: The number of bits used as an additional fingerprint for each ColorSetID. Reduces the number of alien k-mers in a negative query but adds `b` bits for each unique minimizer indexed in Kaminari. Default: `1` 
* **`-g` or `--max-ram`**: The maximum amount of memory allocated for Kaminari, expressed in **MB**. Kaminari will try to use the maximum amount given by the user for better performance. Default: `8192` (8 GB)
* **`-t` or `--threads`**: The number of threads allocated for Kaminari. Default: `1`
* **`-s` or `--seed`**: The random seed for the murmurhash function used for determining minimizers. Two indexes with different seeds will have different results, so this parameter should not be changed unless you know what you are doing. Default: `42`
* **`-c` or `--pthash-constant`**: The PTHash (MPHF) build constant. A higher `c` leads to slower construction but a more space-efficient MPHF. Default: `4.0`
* **`-K` or `--keep-tmp-files`**: Keep temporary files after execution. *Warning: Might cause a crash if there is not enough space on disk.* Default: `OFF`
* **`-v` or `--verbose`**: The tool verbosity level. Options are `0` (silent), `1` (times), `2` (times & stats), or `3` (debug). Default: `0`

### Query 
`./kaminari query` allows you to query an index. `./kaminari query -h` displays a reminder of the different parameters mentioned below.

* **`-x` or `--index`** (Required): The Kaminari index **directory** to use (e.g., if your index is in `my_index/index.kaminari`, use `-x my_index`).
* **`-i` or `--input-list`** (Required): List of FASTA files to query. If exactly 1 file ending with ".list" is provided, it is assumed to be a file of filenames. Each sequence in the files will be queried, and the header will be the name of the answer in the output.
* **`-o` or `--output-filename`**: The filename for the query output. For each sequence of each input file, it outputs its color set. Default: `kaminari_results.txt`
* **`-r` or `--ratio`**: The ratio of k-mers needed to select a color. For example, `r=0.3` means at least 30% of the k-mers must belong to color `c1` to select `c1`. Default: `1.0` (100% of the k-mers)
* **`-t` or `--threads`**: The number of threads allocated for Kaminari. Default: `1`
* **`-v` or `--verbose`**: Increases output verbosity. Default: `0`



## Demo

This short demo shows how to index the 10-genome collection in the folder `example`. The collection can be downloaded by running the script `download_test_datasets.sh` which generates `example/data/salmonella10/`.

From `kaminari/example`, run:

    ./download_test_datasets.sh
    find "$(pwd)/data/salmonella10" -type f > salmonella_10_filenames.list

We will use the standard values k = 31 and m = 19.

`kaminari` can take a single text file listing the inputs (one input file per line) as its `-i` option. For instance, from the `build` directory: 

    ./kaminari build -i ../example/salmonella_10_filenames.list -o ../example/salmonella10_index -k 31 -m 19 -g 4096 -t 4 -v 1
    
This will build an index that will be saved into the directory `example/salmonella10_index`. Canonical minimizers are enabled by default.

Now that the index is generated, we can query the presence of a sequence in the documents.

Still from `kaminari/build`, run:

    ./kaminari query -x ../example/salmonella10_index -i ../example/two_queries.fasta -o ../example/two_queries_result.txt -r 0.8 -t 4 -v 1
    cat ../example/two_queries_result.txt

You can check the result in `example/two_queries_result.txt` and see that the first sequence is mainly present in docs 0, 5, 7, and 9, which corresponds to the 1st, 6th, 8th, and 10th lines in `example/salmonella_10_filenames.list` (document IDs start from 0). 

The second sequence is reported for all 10 documents, but note how only documents 0, 5, 7, and 9 have the maximal count. The other ones share fewer k-mers with the sequence, but still over 80% (`-r 0.8`), so Kaminari reports them.
