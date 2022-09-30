#!/bin/bash

declare -a xyzs=(12 18 24 28 30 34 36 40 44 46 6 14 22 27 29 31 35 38 42 45 47 8)

for xyz in "${xyzs[@]}"
do 
    cp "thf_file${xyz}.xyz" ../utils/xyz_dump/
done