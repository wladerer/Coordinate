from better_turbo import *

file = "/home/wladerer/github/Coordinate/utils/s-bridge-thf-def2-TZVPP.xyz"
thf_file = "/home/wladerer/github/Coordinate/utils/thf.xyz"
test = CoordinationComplex(file, thf_file)
sphere = Sphere(file, 100, 1.5)

point = sphere.valid_points[0]
test.orient_ligand(point)
test.plot_CoordinationComplex(point)

