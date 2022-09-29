#!/bin/bash

###This script is used to setup OpenBabel + include python bindings 
install_dir=$1
obabeltar="openbabel-2.3.1.tar.gz"
echo "Make sure you know that you're in the install directory bc I can't be bothered to check that for you :)"
echo "Unpacking OpenBabel archive: ${obabeltar}"
tar -xzf "$obabeltar"
cd "openbabel-2.3.1" || echo "Failed to cd into openbabel" exit; mkdir ob-build; cd ob-build || exit
cmake -DCMAKE_INSTALL_PREFIX=${install_dir} ../
make -j install