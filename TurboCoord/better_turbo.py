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

    def __init__(self, xyz, atom_type, index):
        self.point = xyz
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

    def __init__(self, xyzfile, n_points=100, cutoff=1.5):

        self.samples = n_points
        self.cutoff = cutoff
        self.xyzfile = xyzfile
        self.points = st.generateSphere(xyzfile, n_points, cutoff)
    
    def __str__(self) -> str:
        return f"A sphere with {self.samples} points and a cutoff of {self.cutoff}"

class Complex:

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

    #currently only works for thf
    def __init__(self, xyzfile) -> None:
        super().__init__(xyzfile)
        self.ligand_axis = self.coords[0] - ((self.coords[1] + self.coords[2])/2) #this needs to be changed later


class CoordinationComplex(Complex):

    def __init__(self, xyzfile: str, ligand_xyzfile: str, point=None) -> None:
        super().__init__(xyzfile)
        self.ligand = Ligand(ligand_xyzfile)
        self.ligand_coords = self.ligand.coords
        self.ligand_axis = self.ligand.ligand_axis
        self.ligand_atoms = self.ligand.atoms
        self.ligand_indices = self.ligand.indices
        self.ligand_Atoms = self.ligand.Atoms
        self.complex_atoms = self.atoms + self.ligand_atoms
        self.complex_coords = np.concatenate((self.coords, self.ligand_coords))
        self.complex_dists = cdist(self.complex_coords, self.complex_coords)

        def orient_ligand(ligand_xyzfile, point=None):
            """
            Orient the ligand to point towards the central atom
            Points must be chosen from the Sphere object
            """
            if point == None:
                #if no point is given, orient the ligand to point towards the first point on the sphere
                points = Sphere(ligand_xyzfile, n_points=100, cutoff=1.5).points
                point = points[0]

            #create a rotation matrix to orient the ligand
            rotation_matrix = st.rotation_matrix(-1*self.ligand_axis, point)
            #rotate the ligand and displace the ligand from the central atom
            oriented_ligand_coords = (rotation_matrix @ self.ligand_coords) + point


            return oriented_ligand_coords

        self.ligand_coords = orient_ligand(ligand_xyzfile, point) 

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

