from TurboCoord import *
import spatialTools as st

file = "/Users/wladerer/github/Coordinate/utils/s-bridge-thf-def2-TZVPP.xyz"
thf_file = "/Users/wladerer/github/Coordinate/utils/thf.xyz"
molecule_info = st.from_xyz(file)
ligand_info = st.from_xyz(file)

ligand = Ligand(*ligand_info)
test = CoordinationComplex(*molecule_info, ligand)
sphere = Sphere(file, 100, 1.5)

for i,point in enumerate(sphere.valid_points):
    new = test.orient_ligand(point)
    new.save(filename=f"s-bridge_{i}")





