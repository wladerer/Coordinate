from TurboCoord import *
# import required module
import os

from TurboCoord import generateStructures
# assign directory
directory = '../utils/test_folder'
 
file = "/home/wladerer/github/Coordinate/utils/s-bridge-thf-def2-TZVPP.xyz"
thf_file = "thf.xyz"

generateStructures(file, thf_file)