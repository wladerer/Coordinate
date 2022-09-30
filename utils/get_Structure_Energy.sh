#!/bin/bash

#Claire notes: 
#####Run in directory where all your mass-directories were created
#####Greps molecule # and energy output from .log file and puts it in UFF-Energies.txt file

for file in /p/work1/wladerer/uffopt/ltnf/xyz/*; do
    mkdir final_g16_structs 
    name=$(basename "$file" .xyz)
    #Make name for directory path
    cd UFF-$name
    export name
    echo $(grep "Energy=" $name.out | perl -ne 'print "$ENV{name} | $_"') >> ../UFF-Energies.txt
    grep "Standard orientation" -A 75 *.out | tail -76 > ${name}.log
    obabel -ig16 ${name}.log -oxyz
    cp ${name}.xyz ../final_g16_structs/
    cd ..
done

tar -cvzf final_ltnf_uff.tar.gz final_g16_structs/