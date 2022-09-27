from cgi import test
from randomSpheres import *
# import required module
import os
# assign directory
directory = '../utils/xyz_dump'
 
file = "/home/wladerer/github/Coordinate/utils/noTHF-molecule.xyz"
testAlgorithm(file, "/home/wladerer/github/Coordinate/utils/thf.xyz", 2.5)


# # iterate over files in
# # that directory
# good_count = 0
# good_files = []
# for filename in os.listdir(directory):
#     f = os.path.join(directory, filename)
    
#     if checkBadness(f, 1):
#         good_count += 1
#         good_files.append(f)

# file_string = []
# for file in good_files:
#     file_string.append(file)
#     os.system(f"cp {file} ../")