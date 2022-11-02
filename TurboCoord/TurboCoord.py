import numpy as np
import spatialTools as st
import pandas as pd
from scipy.spatial.distance import cdist


def writeLines(lines, filename="xyzs.xyz"):
    """
    Writes generated strings to an xyz file
    """
    n_lines = len(lines)
    with open(filename, 'w') as f:
        
        f.write(f"{n_lines}\n\n") #create header for xyz file

        for line in lines:
            f.write(f"{line}\n")

class Atom:
    """
    A class that represents an atom in a molecule

    atom_type: str
        The element of the atom
    xyz: np.ndarray
        The xyz coordinates of the atom
    index: int
        The index of the atom in the xyz file
    """

    def __init__(self, xyz, atom_type, index):
        self.point = xyz
        self.coords = xyz
        self.x = xyz[0]
        self.y = xyz[1]
        self.z = xyz[2]
        self.index = index
        self.element: str = atom_type
        self.frozen_line = f"{self.element:<2}{-1:^4}{xyz[0]:>15.5f}{xyz[1]:>15.5f}{xyz[2]:>15.5f}"
        self.labile_line = f"{self.element:<2}{0:^4}{xyz[0]:>15.5f}{xyz[1]:>15.5f}{xyz[2]:>15.5f}"
        self.line = f"{self.element:<2}{xyz[0]:>15.5f}{xyz[1]:>15.5f}{xyz[2]:>15.5f}"
    
    def __str__(self) -> str:
        return self.frozen_line

class Sphere:
    """
    A class that represents a sphere of points around a central atom
    Mainly used for sampling
    """

    def __init__(self, xyzfile, n_points=100, cutoff=1.5):

        self.samples = n_points
        self.cutoff = cutoff
        self.xyzfile = xyzfile
        self.points = st.generateSphere(xyzfile, n_points, cutoff)
        self.valid_points = self.points[0] #points that are within the cutoff
        self.invalid_points = self.points[1] #points that are outside the cutoff
    
    def __str__(self) -> str:
        return f"A sphere with {self.samples} points and a cutoff of {self.cutoff}"

class Complex:
    """
    A generic class that encodes atomic coordinates and other useful information
    """

    def __init__(self, xyzfile) -> None:
        
        #Unpack the xyz file
        self.coords, self.atoms, self.indices, self.dists = st.from_xyz(xyzfile)
        #Populate the complex with atom objects
        self.Atoms = [Atom(self.coords[i], self.atoms[i], self.indices[i]) for i in range(len(self.atoms))]

    def as_dataframe(self):
        #combine all the data into a pandas dataframe
        return pd.DataFrame({'atom': self.atoms, 'x': self.coords[:, 0], 'y': self.coords[:, 1], 'z': self.coords[:, 2]})

    def __str__(self): #honestly, this is the nicest way to print the complex
        return self.as_dataframe().__str__()

class Ligand(Complex):
    """
    Represents a ligand
    """

    #currently only works for thf
    def __init__(self, xyzfile) -> None:
        super().__init__(xyzfile)
        self.ligand_axis = self.coords[0] - ((self.coords[1] + self.coords[2])/2) #this needs to be changed later


class CoordinationComplex(Complex):

    def __init__(self, xyzfile: str, ligand_xyzfile: str) -> None:
        super().__init__(xyzfile)
        self.ligand = Ligand(ligand_xyzfile)
        self.ligand_axis = self.ligand.ligand_axis
        self.ligand_coords = self.ligand.coords
        self.ligand_atoms = self.ligand.atoms
        self.ligand_indices = self.ligand.indices
        self.ligand_Atoms = self.ligand.Atoms
        
        
        self.complex_atoms = self.atoms + self.ligand_atoms
        self.complex_coords = np.concatenate((self.coords, self.ligand_coords))
        self.distance_matrix = cdist(self.complex_coords, self.complex_coords)

    def orient_ligand(self, point):
        """
        Orients the ligand axis so that it is pointing towards the central atom
        """
        ligand_axis = self.ligand_axis #this is the axis of the ligand
        ligand_coords = self.ligand_coords #these are the coordinates we must rotate

        mat = st.rotation_matrix(ligand_axis, point) #get the rotation matrix
        new_vecs = [ mat @ coord + point for coord in ligand_coords]
        ligand_coords = np.reshape(new_vecs, (len(new_vecs),3)) #rotate the ligand around the central atom and translate

        setattr(self, 'ligand_coords', ligand_coords) #set the new coordinates

        return self


    def plot_CoordinationComplex(self, point=None):
        import plotly.graph_objects as go
        if point is None:

            #plot the complex and ligand atoms separately in a plotly scatter plot
            fig = go.Figure(data=[go.Scatter3d(x=self.coords[:, 0], y=self.coords[:, 1], z=self.coords[:, 2], mode='markers', marker=dict(size=2, color='blue')),
                                go.Scatter3d(x=self.ligand_coords[:, 0], y=self.ligand_coords[:, 1], z=self.ligand_coords[:, 2], mode='markers', marker=dict(size=2, color='red'))])
            
            fig.show()
        
        if point is not None:

            #rotate coordinates with respect to the point specified
            self.orient_ligand(point)
            #plot the complex and ligand atoms separately in a plotly scatter plot
            fig = go.Figure(data=[go.Scatter3d(x=self.coords[:, 0], y=self.coords[:, 1], z=self.coords[:, 2], mode='markers', marker=dict(size=2, color='blue')),
                                go.Scatter3d(x=self.ligand_coords[:, 0], y=self.ligand_coords[:, 1], z=self.ligand_coords[:, 2], mode='markers', marker=dict(size=2, color='red'))])
            
            fig.show()                    
        

    def rotate_ligand(self, theta, axis):
        """
        Rotates the ligand around an axis by an angle theta
        """
        self.ligand_coords = st.rotate(self.ligand_coords, theta, axis)
        self.ligand.axis = st.rotate(self.ligand.axis, theta, axis)
        self.ligand.Atoms = [Atom(self.ligand_coords[i], self.ligand.atoms[i], self.ligand.indices[i]) for i in range(len(self.ligand.atoms))]

    def rotate(self, theta, axis):
        """
        Rotates the complex around an axis by an angle theta
        """
        self.coords = st.rotate(self.coords, theta, axis)
        self.ligand_coords = st.rotate(self.ligand_coords, theta, axis)
        self.ligand.axis = st.rotate(self.ligand.axis, theta, axis)
        self.axis = st.rotate(self.axis, theta, axis)
        self.Atoms = [Atom(self.coords[i], self.atoms[i], self.indices[i]) for i in range(len(self.atoms))]
        self.ligand.Atoms = [Atom(self.ligand_coords[i], self.ligand.atoms[i], self.ligand.indices[i]) for i in range(len(self.ligand.atoms))]

    def translate(self, vector):
        """
        Translates the complex by a vector
        """
        self.coords = st.translate(self.coords, vector)
        self.ligand_coords = st.translate(self.ligand_coords, vector)
        self.ligand.axis = st.translate(self.ligand.axis, vector)
        self.axis = st.translate(self.axis, vector)
        self.Atoms = [Atom(self.coords[i], self.atoms[i], self.indices[i]) for i in range(len(self.atoms))]
        self.ligand.Atoms = [Atom(self.ligand_coords[i], self.ligand.atoms[i], self.ligand.indices[i]) for i in range(len(self.ligand.atoms))]

    def __str__(self) -> str:
        return f"Coordination Complex with {len(self.atoms)} atoms and {len(self.ligand.atoms)} ligand atoms"

    def save(self, filename=f"coordination_complex.xyz", freeze_ligand=False):
        """
        Saves the complex to an xyz file
        """
        
        if freeze_ligand:

            lines = [atom.labile_line for atom in self.Atoms]
            #replace the ligand line with the frozen ligand line
            lines.extend([self.Atoms[i].frozen_line if i in self.ligand.indices else self.Atoms[i].line for i in range(len(self.Atoms))])
            writeLines(lines, filename)
                
        else:
            
            lines = [atom.line for atom in self.Atoms]
            lines.extend([atom.line for atom in self.ligand.Atoms])
            writeLines(lines, filename)


    def as_dataframe(self):
        """
        combine all the coordination complex attributes into a pandas dataframe
        """
        import pandas as pd 
        return pd.DataFrame({'atom': self.complex_atoms, 'index': self.complex_indices, 'x': self.complex_coords[:, 0], 'y': self.complex_coords[:, 1], 'z': self.complex_coords[:, 2]})


def generateStructures(xyzfile: str, ligand_xyzfile: str) -> list[CoordinationComplex]:
    """
    Generates a list of CoordinationComplex objects from an xyz file
    Returns a list of CoordinationComplex objects
    """
    complex = CoordinationComplex(xyzfile, ligand_xyzfile) #generate a sampling sphere around the central atom
    points = Sphere(xyzfile, 100, 1.5).valid_points #collect only the valid points
    structures = [complex.orient_ligand(point) for point in points] #create a list of structures with the ligand oriented around the central atom

    return structures


def filterStructures(structures: list[CoordinationComplex], cutoff_distance: float) -> list[CoordinationComplex]:
    '''
    Filters structures that have atoms too close to other atoms
    '''
    #get distances between all atoms in the complex
    distances = [distance(structures[i].Atoms[j], structures[i].Atoms[k]) for i in range(len(structures)) for j in range(len(structures[i].Atoms)) for k in range(len(structures[i].Atoms)) if j != k]
    #get the indices of the structures with distances less than the cutoff distance
    indices = [i for i in range(len(structures)) if any([distances[i*len(structures[i].Atoms) + j*len(structures[i].Atoms) + k] < cutoff_distance for j in range(len(structures[i].Atoms)) for k in range(len(structures[i].Atoms)) if j != k])]
    #filter the structures
    structures = [structures[i] for i in range(len(structures)) if i not in indices]
    return structures

def distance(atom1: Atom, atom2: Atom) -> float:
    """
    Calculates the distance between two atoms
    """
    return np.linalg.norm(atom1.coords - atom2.coords)