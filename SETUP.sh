#!/bin/bash
echo ========================================
echo Setting up Genetic Data analysis Library and compiling utility tools
echo ========================================
cd util
sh INSTALL.sh
cd ..
echo ========================================
echo Downloading and installing FragGeneScan
echo ========================================
cd util/ext
wget -O FragGeneScan.tar.gz https://sourceforge.net/projects/fraggenescan/files/latest/download
tar -zxvf FragGeneScan.tar.gz
mv FragGeneScan*/ FragGeneScan
cd FragGeneScan
make clean
make fgs
cd ../../..
echo ========================================
echo Downloading and extracting gold standard datasets
echo ========================================
mkdir data
cd data
wget http://giant.csail.mit.edu/carnelian/EC-2010-DB.tar.gz
tar -zxvf EC-2010-DB.tar.gz
rm EC-2010-DB.tar.gz
cd ..
