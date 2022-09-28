from randomSpheres import *
# import required module
import os
# assign directory
directory = '../utils/'
 
file = ### put Yb xyz file here
thf_file = "thf.xyz"
testAlgorithm(file, thf_file, 2.5)
find_all_sites(file, thf_file, 2.5)

os.system(f"tar -cvzf xyzs.tar.gz {directory}")


#you can use this to see if the structures are good on first pass
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