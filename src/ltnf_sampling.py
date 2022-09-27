from randomSpheres import *
# import required module
import os
# assign directory
directory = '../utils/xyz_dump'
 
# iterate over files in
# that directory
good_count = 0
good_files = []
for filename in os.listdir(directory):
    f = os.path.join(directory, filename)
    
    if checkBadness(f, 1):
        good_count += 1
        good_files.append(f)


for file in good_files:
    print(file)