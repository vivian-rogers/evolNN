#!/bin/bash


wget https://gitlab.com/libeigen/eigen/-/archive/3.3.9/eigen-3.3.9.zip 
unzip eigen-3.3.9.zip 
#rm eigen-3.3.9.zip
mv ./eigen-3.3.9 ./eigen
#rm -r eigen-3.9.9
cd mc
make
