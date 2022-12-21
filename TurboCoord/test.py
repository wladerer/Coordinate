from TurboCoord import *
import spatialTools as st

file = "/home/wladerer/github/Coordinate/utils/s-bridge-thf-def2-TZVPP.xyz"
thf_file = "/home/wladerer/github/Coordinate/utils/thf.xyz"

complex = Complex.from_xyz(file, thf_file)
points = Sphere(file).valid_points
ligand_axis = complex.ligand.coords[0] - ( complex.ligand.coords[2] + complex.ligand.coords[3] ) / 2

for i,point in enumerate(points):
   complex.orient_ligand(ligand_axis)
   complex.to_xyz(f"s-bridge_{i}.xyz", freeze = False)
