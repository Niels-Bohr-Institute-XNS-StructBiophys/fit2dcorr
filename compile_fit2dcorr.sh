#!/bin/bash
# ./compile_fit2dcorr.sh g++ 3 NONE
make clean
echo "make CC=$1 OLEVEL=$2 OOPTIONS=$3"
make CC=$1 OLEVEL=$2 OOPTIONS=$3

