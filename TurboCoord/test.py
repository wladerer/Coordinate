from better_turbo import *

 
file = "/home/wladerer/Downloads/s-bridge-thf-def2-TZVPP.xyz"
thf_file = "/home/wladerer/github/Coordinate/utils/thf.xyz"
 
test = CoordinationComplex(file, thf_file)
test.save("test1.xyz")