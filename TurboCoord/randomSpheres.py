import json

import numpy as np
import pandas as pd




def rotation_matrix_from_vectors(vec1, vec2):
    """ Find the rotation matrix that aligns vec1 to vec2
    vec1 is the ligand axis
    vec2 is the negative of the Yb - O bond
    returns a transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rotation_matrix


def fibonacci_sphere(samples=100):
    """
    Generates a sphere around the coordination complex
    Number of points is adjustable
    """
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
    """
    Writes generated strings to an xyz file
    """
    n_lines = len(lines)
    with open(filename, 'w') as f:
        
        f.write(f"{n_lines}\n\n") #create header for xyz file

        for line in lines:
            f.write(f"{line}\n")


def validPoints(sphere, xyzarr, cutoff=1.5):
    """
    Tests if a point within the generated sphere is too close to a atom within a designated cutoff distance
    """
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
    """
    Sets the first atom in the xyz array as the origin
    It is recommended that you format your xyz file such that the coordination complex is the first atom
    """
    xyzarr = xyzarr - xyzarr[0]
    return xyzarr


def getInfo(file):
    """
    Primes the xyzarr and sets the origin for other vector caclulus needs
    """
    xyzarr, atoms = readxyz(file) 
    xyzarr = set_origin(xyzarr)

    return xyzarr, atoms

def sphericalSample(radius, xyzarr):
    """
    Returns the points on the sample sphere of a given radius
    """
    sphere = radius*np.array(fibonacci_sphere())
    points = validPoints(sphere, xyzarr)

    return points

def makeLigand(file = "/home/wladerer/github/minis/thf.xyz"):
    """
    Currently only works for THF - missing a way to define the axis of the ligand needed for rotation
    """
    ligand_arr, ligand_atoms = getInfo(file)
    ligand_axis = ligand_arr[0] - ((ligand_arr[1] + ligand_arr[2])/2)

    return ligand_arr, ligand_atoms, ligand_axis

def runTurboCoord(molecule_file, ligand_file, radius):

    xyzarr, atoms = getInfo(molecule_file)
    points = sphericalSample(radius, xyzarr)
    ligand_arr, ligand_atoms, ligand_axis = makeLigand(ligand_file)
    

    lines = []
    for xyz, atom in zip(xyzarr,atoms):
        lines.append(f"{atom:<2}{-1:^4}{xyz[0]:>15.5f}{xyz[1]:>15.5f}{xyz[2]:>15.5f}")

    for i,point in enumerate(points):
        ligand_rot_mat = rotation_matrix_from_vectors(-1*ligand_axis, point)
        ligand_lines = []
        for xyz, atom in zip(ligand_arr,ligand_atoms):
            xyz = (ligand_rot_mat @ xyz) + point

            ligand_lines.append(f"{atom:<2}{0:^4}{xyz[0]:>15.5f}{xyz[1]:>15.5f}{xyz[2]:>15.5f}")

            writeLines(lines + ligand_lines, filename=f"thf_file{i}.xyz")

        


def testAlgorithm(molecule_file, ligand_file, radius):

    xyzarr, atoms = getInfo(molecule_file)
    points = sphericalSample(radius, xyzarr)

    lines = []
    for xyz, atom in zip(xyzarr,atoms):
        lines.append(f"{atom:<2}{xyz[0]:>15.5f}{xyz[1]:>15.5f}{xyz[2]:>15.5f}")

    dummy = "XX"
    for xyz in points:
        lines.append(f"{dummy:<2}{xyz[0]:>15.5f}{xyz[1]:>15.5f}{xyz[2]:>15.5f}")

    writeLines(lines, "all_points.xyz")


def distance_matrix_from_points(xyzs):
    """
    Creates a distance matrix for the xyz array
    """
    n_atoms = len(xyzs)
    dist_mat = np.zeros((n_atoms, n_atoms))
    for i in range(n_atoms):
        for j in range(n_atoms):
            dist_mat[i,j] = np.linalg.norm(xyzs[i] - xyzs[j])
    return dist_mat



def checkBadness(xyzs, filter):
    import numpy as np

    dists = distance_matrix_from_points(xyzs)
    sorted_dists = np.sort(dists)

    averages = []
    
    for i,dist in enumerate(sorted_dists):
        avg = np.average(dist[0:5])
        averages.append(avg)

    averages = np.sort(averages) 

    shortest_bond_avg = np.average(averages[1:10])
    if shortest_bond_avg > filter:
        return True
    else:
        return False


def dipole(ligand_xyz):
    import pandas as pd
    r, ligand_atoms = getInfo(ligand_xyz)

    ptable = pd.read_csv("/home/wladerer/github/Coordinate/utils/ptable.csv")

    q = [ptable.loc[ptable['Symbol'] == ligand_atom, 'Electronegativity'].values[0] for ligand_atom in ligand_atoms] #overly complicated way to match an e.neg to a symbol
    
    mu_x = []
    mu_y = []
    mu_z = []
    dipole = []
    for i,q_i in enumerate(q): #standard dipole calculation
        mu = (q_i * r[i])
        mu_x.append(mu[0])
        mu_y.append(mu[1])
        mu_z.append(mu[2])
        dipole += mu

    ligand_axis = dipole

    return mu_x, mu_y, mu_z, ligand_atoms, ligand_axis


def checkfiles(dir):
    import os
    files = os.listdir(dir)
    files = [file for file in files if file.endswith(".xyz")]

    return files

