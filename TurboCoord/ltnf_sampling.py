from TurboCoord import *
# import required module
import os
import numpy as np
from TurboCoord import generateStructures
from TurboCoord import filterStructures
# assign directory
directory = '../utils/test_folder'
 
file = "/home/wladerer/github/Coordinate/utils/s-bridge-thf-def2-TZVPP.xyz"
thf_file = "thf.xyz"

structures = generateStructures(file, thf_file)

#save the structures to a folder
for i, structure in enumerate(structures):
    structure.save(f"{directory}/sbridge_{i}.xyz", freeze_ligand=False)