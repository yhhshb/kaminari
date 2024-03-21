#!/usr/bin/bash

# Script to install COBS, PAC, RAMBO, kmindex, RAPTOR, metaprofi
# assuming all requirements are met, it installs everything but it might 
# be better to go one by one and check for errors

##############################################################
# MANUAL INSTALL SOFTS
##############################################################

#COBS https://github.com/bingmann/cobs  ========================================

# test with older cmake 
git clone --recursive https://github.com/bingmann/cobs.git;
mkdir cobs/build;
cd cobs/build;
cmake ..;
make -j;

cd ../..;

#RUN with ./cobs/build/src/cobs ARGS


#PAC https://github.com/Malfoy/PAC  ============================================

#had to install zlib from https://zlib.net/
#then add `-I../zlib/include` at the end of CFLAGS=
#and `-L../zlib/lib` at the end of LDFLAGS=

git clone --depth 1 --recursive https://github.com/Malfoy/PAC.git
cd PAC
make

cd ..

#RUN with ./PAC/PAC ARGS


#RAMBO https://github.com/RUSH-LAB/RAMBO  ======================================

#had to run these commands after error on ./configure
#libtoolize --force
#aclocal
#autoheader
#automake --force-missing --add-missing
#autoconf
#./configure

git clone --recursive https://github.com/RUSH-LAB/RAMBO.git
cd RAMBO
./configure
make 
make install

cd ..

#RUN with ./RAMBO/bin/rambo ARGS


##############################################################
# CONDA SOFTS
##############################################################

conda create -p competitors_env python==3.8 pigz
conda activate ./competitors_env

#KMINDEX https://github.com/tlemane/kmindex  ===================================

conda install -c conda-forge -c tlemane kmindex

#git clone --recursive https://github.com/tlemane/kmindex 
#cd kmindex && ./install.sh && cd .. 

#RUN with kmindex ARGS


#RAPTOR https://github.com/seqan/raptor  =======================================

conda install -c bioconda -c conda-forge raptor

#git clone --recursive https://github.com/seqan/raptor
#cd raptor && mkdir build && cd build && cmake .. && make && cd ../..

#RUN with raptor ARGS


#METAPROFI https://github.com/kalininalab/metaprofi  ===========================

git clone https://github.com/kalininalab/metaprofi.git
#pip install Cython==0.29.28 numpy==1.22.3 if error with SharedArray
pip install ./metaprofi

#RUN with metaprofi ARGS