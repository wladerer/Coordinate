from gcutil import distance_matrix, readxyz
import numpy as np


def rotation_matrix_from_vectors(vec1, vec2):
    """ Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector -- this would be our thf axis
    :param vec2: A 3d "destination" vector --- this would be the negative of the Yb - O bond
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rotation_matrix


def fibonacci_sphere(samples=100):

    points = []
    phi = np.pi * (3. - np.sqrt(5.))  # golden angle in radians

    for i in range(samples):
        y = 1 - (i / float(samples - 1)) * 2  # y goes from 1 to -1
        radius = np.sqrt(1 - y * y)  # radius at y

        theta = phi * i  # golden angle increment

        x = np.cos(theta) * radius
        z = np.sin(theta) * radius

        points.append((x, y, z))

    return points

def writeLines(lines, filename="xyzs.xyz"):
    n_lines = len(lines)

    with open(filename, 'w') as f:
        
        f.write(f"{n_lines}\n\n") #create header for xyz file

        for line in lines:
            f.write(f"{line}\n")


def validPoints(sphere, xyzarr, cutoff=1.5):
    points = []
    for point in sphere:
        dists = []
        for atom in xyzarr:
            dist = np.linalg.norm(point -atom)
            dists.append(dist)

        list_check = True
        for dist in dists:
            if dist < cutoff:
               list_check = False
    
        if list_check == True:
            points.append(point)
    return points

def set_origin(xyzarr):

    xyzarr = xyzarr - xyzarr[0]
    return xyzarr


def getInfo(file):
    xyzarr, atoms = readxyz(file) 
    xyzarr = set_origin(xyzarr)

    return xyzarr, atoms

def sphericalSample(radius, xyzarr):
    sphere = radius*np.array(fibonacci_sphere())
    points = validPoints(sphere, xyzarr)

    return points

def makeLigand(file, pointsInPlane=[0,2,3]):
    file = "/home/wladerer/github/minis/thf.xyz"
    ligand_arr, ligand_atoms = getInfo(file)
    ligand_axis = pointsInPlane[0] - ((pointsInPlane[2] + pointsInPlane[3])/2)

    return ligand_arr, ligand_atoms, ligand_axis

def find_all_sites(molecule_file, ligand_file, radius):

    points = sphericalSample(radius, ligand_file)
    xyzarr, atoms = getInfo(molecule_file)
    ligand_arr, ligand_atoms, ligand_axis = makeLigand(ligand_file)


    lines = []
    for xyz, atom in zip(xyzarr,atoms):
        lines.append(f"{atom:<2}{xyz[0]:>15.5f}{xyz[1]:>15.5f}{xyz[2]:>15.5f}")

    for i,point in enumerate(points):
        ligand_rot_mat = rotation_matrix_from_vectors(-1*ligand_axis, point)
        ligand_lines = []
        for xyz, atom in zip(ligand_arr,ligand_atoms):
            xyz = ligand_rot_mat @ xyz
            xyz = xyz + point
            ligand_lines.append(f"{atom:<2}{xyz[0]:>15.5f}{xyz[1]:>15.5f}{xyz[2]:>15.5f}")

            writeLines(lines + ligand_lines, filename=f"../utils/xyz_dump/thf_file{i}.xyz")


def testAlgorithm(molecule_file, ligand_file, radius):

    points = sphericalSample(radius, ligand_file)
    xyzarr, atoms = getInfo(molecule_file)

    lines = []
    for xyz, atom in zip(xyzarr,atoms):
        lines.append(f"{atom:<2}{xyz[0]:>15.5f}{xyz[1]:>15.5f}{xyz[2]:>15.5f}")

    dummy = "XX"
    for xyz in points:
        lines.append(f"{dummy:<2}{xyz[0]:>15.5f}{xyz[1]:>15.5f}{xyz[2]:>15.5f}")

    writeLines(lines, "../utils/all_points.xyz")

    