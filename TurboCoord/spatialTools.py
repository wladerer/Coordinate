import numpy as np
from scipy.spatial.distance import cdist

def set_origin(coords):
    """
    Sets the first atom in the xyz array as the origin.
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


def rotation_matrix(vec1, vec2) -> np.ndarray:
    """ Find the rotation matrix that aligns vec1 to vec2
    vec1 is the ligand axis
    vec2 is the negative of the Yb - O bond
    returns a transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c: np.ndarray = np.dot(a, b)
    s: float = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rotation_matrix

def fibonacci_sphere(samples=100) -> list:
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

def generateSphere(xyzfile: str, samples: int, cutoff: float=1.26, radius: float = 2.5):
    """
    Generates a sphere around the coordination complex
    """

    if xyzfile == None:
        raise Exception("No xyz file provided")

    sphere = radius*np.array(fibonacci_sphere(samples))
    coords, *_ = from_xyz(xyzfile)
    points = validPoints(sphere, coords, cutoff)

    return points

def validPoints(points, coords, cutoff=1.3) -> tuple[list, list]:
    """
    Tests if one point in sphere is too close to a point in coords, if it is too close, remove the point
    """
    valid = []
    invalid = []
    for point in points:
        if np.all(cdist([point], coords) > cutoff):
            valid.append(point)
        else:
            invalid.append(point)
    return valid, invalid


def plot_sphere(xyzfile: str, samples: int, cutoff: float, radius: float = 2.5):
    """
    Generates a 3d plot of the sphere around the coordination complex using plotly 3d scatter plot
    """
    import plotly.graph_objects as go

    validpoints, invalidpoints  = generateSphere(xyzfile, samples, cutoff, radius)
    x, y, z = zip(*validpoints)


    fig = go.Figure(data=[go.Scatter3d(x=x, y=y, z=z, mode='markers', name='Ligand' ,marker=dict(size=4, color='blue', opacity=0.8))])
    
    # now plot the coordination complex as a new trace
    
    coords, *_ = from_xyz(xyzfile)
    cx, cy, cz = zip(*coords)
    fig.add_trace(go.Scatter3d(x=cx, y=cy, z=cz, mode='markers', name='Complex', marker=dict(size=4, color='red')))


    fig.show()


def plot_complex(xyzfile: str, ligand_xyzfile: str):
    """
    Generates a 3d plot of the coordination complex using plotly 3d scatter plot
    """
    import plotly.graph_objects as go

    cx, cy, cz = vector_decomposition(xyzfile)
    lx, ly, lz = vector_decomposition(ligand_xyzfile)

    fig = go.Figure(data=[go.Scatter3d(x=cx, y=cy, z=cz, mode='markers', marker=dict(size=4, color='red'))]) #plot the coordination complex
    fig.add_trace(go.Scatter3d(x=lx, y=ly, z=lz, mode='markers', marker=dict(size=4, color='blue'))) # now plot the ligand as a new trace

    fig.show()

def vector_decomposition(ligand_xyzfile: str):
    """
    Decomposes each point of an xyz file into x, y, and z components
    """
    coords, *_ = from_xyz(ligand_xyzfile)
    x, y, z = zip(*coords)
    return x,y,z


def plot_3d_vectors(*vectors):
    """
    Plots multiple 3d vectors using a plotly image

    vectors:: list of tuples of 3d vectors
    """
    import plotly.graph_objects as go
    fig = go.Figure()
    for i, vec in enumerate(vectors):
        fig.add_trace(go.Scatter3d(x=[0, vec[0]], y=[0, vec[1]], z=[0, vec[2]], mode='lines', name=f'vec{i}'))
    fig.show()


def align_vectors(vec1, vec2, translation: float =0.0):
    """
    Aligns vec1 with vec2 using a rotation matrix

    vec1:: 3d vector
    vec2:: 3d vector
    """

    mat = rotation_matrix(vec1, vec2)
    rvec2 = vec2 @ mat + translation

    return rvec2


def plot_cdist(xyzfile: str, cutoff: float):
    """
    Plots the distance between each point in the xyz file using a plotly heat map
    """
    import plotly.graph_objects as go
    import plotly.express as px
    coords, *_ = from_xyz(xyzfile)
    dist = cdist(coords, coords)
    fig = px.imshow(dist)
    fig.show()


def plot_cdist_sphere(xyzfile: str, samples: int, cutoff: float, radius: float = 2.5):
    """
    Plots the distance between each point in the xyz file and the sphere using a plotly heat map
    """
    import plotly.graph_objects as go
    import plotly.express as px
    coords, *_ = from_xyz(xyzfile)
    sphere = radius*np.array(fibonacci_sphere(samples))
    dist = cdist(coords, sphere)
    fig = px.imshow(dist)
    fig.show()

def rotate(coords, theta):
    """
    Rotate a set of coordinates around an axis by an angle theta
    """
    
    #create an arbitrary 3d rotation matrix
    R_z = np.array([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
    R_y = np.array([[np.cos(theta), 0, np.sin(theta)], [0, 1, 0], [-np.sin(theta), 0, np.cos(theta)]])
    R_x = np.array([[1, 0, 0], [0, np.cos(theta), -np.sin(theta)], [0, np.sin(theta), np.cos(theta)]])
    
    #create a dictionary representing the rotation matrices
    rotations = {'x': R_x, 'y': R_y, 'z': R_z}
    R = rotations[np.random.choice(['x', 'y', 'z'])] #choose a random rotation matrix
    
    #rotate the coordinates
    rotated_coords = coords @ R
    
    return rotated_coords


def rotate_ligand_away(xyzfile: str, ligand_xyzfile: str, theta: float):
    """
    Rotates the ligand away from the coordination complex by an angle theta
    """
    coords, *_ = from_xyz(xyzfile)
    ligand_coords, *_ = from_xyz(ligand_xyzfile)
    
    #calculate the distance between the ligand atoms and the coordination complex atoms
    dist = cdist(ligand_coords, coords)

    #rotate the ligand away from the coordination complex until the average cdist is minimized
    min_dist = np.inf
    for angle in np.linspace(0, theta, 100):
        rotated_coords = rotate(ligand_coords, angle)
        new_dist = cdist(rotated_coords, coords)
        if new_dist.mean() < min_dist:
            min_dist = new_dist.mean()
            min_coords = rotated_coords
    
    return min_coords
    

def plot_rotation(xyzfile: str, ligand_xyzfile: str, theta: float):
    """
    Plots the rotation of the ligand away from the coordination complex
    """
    import plotly.graph_objects as go
    coords, *_ = from_xyz(xyzfile)
    ligand_coords, *_ = from_xyz(ligand_xyzfile)
    
    min_coords = rotate_ligand_away(xyzfile, ligand_xyzfile, theta)
    
    #plot the rotation
    fig = go.Figure(data=[go.Scatter3d(x=coords[:,0], y=coords[:,1], z=coords[:,2], mode='markers', marker=dict(size=4, color='red'))])
    fig.add_trace(go.Scatter3d(x=ligand_coords[:,0], y=ligand_coords[:,1], z=ligand_coords[:,2], mode='markers', marker=dict(size=4, color='blue')))
    fig.add_trace(go.Scatter3d(x=min_coords[:,0], y=min_coords[:,1], z=min_coords[:,2], mode='markers', marker=dict(size=4, color='green')))
    fig.show()

