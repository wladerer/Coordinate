from randomSpheres import *
# import required module
import os
# assign directory
directory = '../utils/xyz_dump'
 
file = "/home/wladerer/github/Coordinate/less-terts-noferr-thf-def2-TZVPP.xyz"
thf_file = "thf.xyz"

find_all_sites(file, thf_file, 2.5)

