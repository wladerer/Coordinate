import numpy as np
from scipy.spatial.distance import cdist

def set_origin(coords):
    """
    Sets the first atom in the xyz array as the origin
    It is recommended that you format your xyz file such that the coordination complex is the first atom
    """
    coords = coords - coords[0]
    return coords


def from_xyz(xyzfile: str):
    """
    Reads in an xyz file and returns an a numpy array of xyz coordinates and atom names
    """

    #Read the first line of the file to get the number of atoms
    with open(xyzfile, 'r') as f:
        n_atoms = int(f.readline())
        #skip the second line
        f.readline()
        #Initialize the Nx3 xyz array with the number of atoms
        coords = np.zeros((n_atoms, 3))
        #Initialize the list of atom names
        atoms = []
        indices = []
        #Loop over the lines in the file
        for i, line in enumerate(f):
            #Split the line into a list of strings
            line = line.split()
            #Add the atom name to the list
            atoms.append(line[0])
            indices.append(i)
            #Add the xyz coordinates to the array
            coords[i, :] = line[1:]
    
    #set the origin to the first atom
    coords = set_origin(coords)
        
    #Finally, find the distance array
    dist_mat = cdist(coords, coords) #cdist is a function from scipy

    return coords, atoms, indices, dist_mat


def rotation_matrix(vec1, vec2):
    """ Create a rotation matrix that aligns two vectors """
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
    phi = np.pi *(3.0 - np.sqrt(5.0))  # golden angle in radians
    for i in range(samples):
        y = 1 - (i / float(samples - 1)) * 2  # y goes from 1 to -1
        radius = np.sqrt(1 - y * y)  # radius at y
        theta = phi * i  # golden angle increment
        x = np.cos(theta) * radius
        z = np.sin(theta) * radius
        points.append((x, y, z))

    return points


# def validPoints(points, coords, cutoff=0.3):
#     """
#     Tests if one point in sphere is too close to a point in coords
#     """


    
    

   



def generateSphere(xyzfile: str, samples: int, cutoff: float):
    """
    Generates a sphere around the coordination complex
    """

    if xyzfile == None:
        raise Exception("No xyz file provided")

    
    coords, *_ = from_xyz(xyzfile)
    points = fibonacci_sphere(samples)
    # points = validPoints(sphere, coords, cutoff)

    return points


