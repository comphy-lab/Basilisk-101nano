#!/bin/bash

mkdir -p $1

cp $1.c $1/
cd $1
qcc -O2 -Wall -disable-dimensions $1.c -o $1 -lm
./$1