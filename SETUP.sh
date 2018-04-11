#!/bin/bash
echo ========================================
echo Setting up Genetic Data analysis Library and compiling utility tools
echo ========================================
cd util
sh INSTALL.sh
cd ..
echo ========================================
echo Downloading and extracting gold standard datasets
echo ========================================
mkdir data
cd data
wget http://giant.csail.mit.edu/carnelian/EC-2192-DB.tar.gz
tar -zxvf EC-2192-DB.tar.gz
rm EC-2192-DB.tar.gz
cd ..
