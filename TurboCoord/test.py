from TurboCoord import *
import spatialTools as st

file = "/home/wladerer/github/Coordinate/utils/s-bridge-thf-def2-TZVPP.xyz"
thf_file = "/home/wladerer/github/Coordinate/utils/thf.xyz"

complex = Complex.from_xyz(file, thf_file)
points = Sphere(file).valid_points
ligand_axis = complex.ligand.coords[0] - ( complex.ligand.coords[2] + complex.ligand.coords[3] ) / 2


# X = "X"

# filename = "aaaaaaahhhh.xyz"
# n_atoms = len(complex.substrate.atoms) + len(points)
# with open(filename, 'w') as file:
#       file.write(f"{n_atoms}\n\n")

#       for point in points:
#          line = f"{X:<2}{point[0]:>15.5f}{point[1]:>15.5f}{point[2]:>15.5f}\n"
#          file.write(line)

#       for atom in complex.atoms:
#          line = f"{atom.element:<2}{atom.x:>15.5f}{atom.y:>15.5f}{atom.z:>15.5f}\n"
#          file.write(line)


for i,point in enumerate(points):
   complex.orient_ligand(-1*ligand_axis, point)
   complex.to_xyz(f"s-bridge_{i}.xyz", freeze = False)
print(points[35])
# complexes = [complex.orient_ligand(-1*ligand_axis, point) for point in points]



