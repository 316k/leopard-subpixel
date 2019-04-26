#!/bin/bash

for i in *unwrapping/
do
    cd $i
    ./all.sh
    cd ..
done
