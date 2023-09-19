# kaminari
雷 - kaminari (thunder/lightning)

Kaminari is a *colored compacted de Bruijn graph* approximate index for large-scale matching and color queries, in the same vein as [Fulgor](https://github.com/jermp/fulgor), but powered by [LPHash](https://github.com/jermp/lphash).
Colors are still pre-computed [GGCAT](https://github.com/algbio/GGCAT). 

**A pre-print describing how the index works can be found [here (TODO)]().**

### Table of contents
* [Dependencies](#dependencies)
* [Compiling the code](#compiling-the-code)
* [Tools](#tools)
* [Demo](#Demo)
* [Indexing an example Salmonella pan-genome]

Dependencies
------------

#### GGCAT

The code uses the [GGCAT](https://github.com/algbio/GGCAT) Rust library,
so make sure you have Rust installed. 
If not, Rust can be installed as recommended [here](https://www.rust-lang.org/tools/install), with

	curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

Then, switch to the [Nightly](https://doc.rust-lang.org/book/appendix-07-nightly-rust.html#rustup-and-the-role-of-rust-nightly) release channel.

	rustup toolchain install nightly

#### zlib

If you do not have `zlib` installed, you can do

    sudo apt-get install zlib1g

if you are on Linux/Ubuntu, or

    brew install zlib

if you are using macOS.

Compiling the code
------------------

The code is tested on ~~Linux with `gcc`~~ and on MacOS with `clang`.
To build the code, [`CMake`](https://cmake.org/) is required.

First clone the repository with

    git clone --recursive https://github.com/yhhshb/kaminari.git

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

	Usage: ./kaminari [--help] [--version] {build}

    Optional arguments:
    -h, --help     shows help message and exits 
    -v, --version  prints version information and exits 

    Subcommands:
    build         Build a kaminari index from a colored compacted de Bruijn Graph

Demo
----

This short demo shows how to index the 10-genome collection in the folder `data/salmonella_10` which can be downloaded by running the script `download_test_datasets.sh` that can be found in the `scripts` subdirectory of `kaminari`.

We will use the standard value k = 31.

`kaminari` can take multiple filenames as its -i option

    ./kaminari build -i ../data/salmonella10/SAL_*.fasta.gz -o ../data/salmonella10.kmn.index -k 31 -m 19 -d tmp_dir -t 1 --verbose --check

Or, alternatively, a single text file listing the inputs (one input file per line).
First create a list of filenames (with absolute paths) for the files in `data/salmonella10`.
From `kaminari/data`, do

    find $(pwd)/salmonella10/* > salmonella_10_filenames.txt

Then, from `kaminari/build`, run

    ./kaminari build -i ../data/salmonella_10_filenames.txt -o ../data/salmonella10.kmn.index -k 31 -m 19 -d tmp_dir -t 1 --verbose --check

in order to build an index that will be serialized to the file `data/salmonella10.kmn.index`.