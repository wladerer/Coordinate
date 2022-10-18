from better_turbo import *

file = "/Users/wladerer/github/Coordinate/utils/s-bridge-thf-def2-TZVPP.xyz"
thf_file = "/Users/wladerer/github/Coordinate/utils/thf.xyz"
test = CoordinationComplex(file, thf_file, [0,1,1])

test.save("test1.xyz")