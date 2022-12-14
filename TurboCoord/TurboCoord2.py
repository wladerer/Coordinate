import numpy as np
import spatialTools as st
import pandas as pd
from scipy.spatial.distance import cdist

# a molecule is made from atoms - a complex is made from two molecules - therefore we need to 
# faithfully represent atoms, molecules, and complexes

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


class Atom:

    def __init__(self, element: str, coordinates: np.ndarray) -> None:
        self.element = element
        self.coordinates = coordinates
        self.x = coordinates[0]
        self.y = coordinates[1]
        self.z = coordinates[2]
    
    def update_position(self, new_coordinates: np.ndarray):
        self.coordinates = new_coordinates
        self.x = new_coordinates[0]
        self.y = new_coordinates[1]
        self.z = new_coordinates[2]

    def __repr__(self) -> str:
        return "Atom of type {self.element} located at {self.x},{self.y},{self.z}"

class Molecule:

    def __init__(self, atoms: list[Atom]) -> None:
        self.atoms = atoms
        self.elements = [atom.element for atom in atoms]
        self.coords = [atom.coordinates for atom in atoms]


class Complex:

    def __init__(self, substrate: Molecule, ligand: Molecule) -> None:
        self.substrate = substrate
        self.ligand = ligand
    
    
    def orient_ligand(self, point, ligand_axis: np.ndarray):
        """
        Orients the ligand axis so that it is pointing towards the central atom
        """

        mat = st.rotation_matrix(ligand_axis, point) #get the rotation matrix
        new_vecs = [ mat @ coord + point for coord in ligand_coords]
        ligand_coords = np.reshape(new_vecs, (len(new_vecs),3)) #rotate the ligand around the central atom and translate

        #create a new ligand object and re-assign it to the complex
        self.ligand = Molecule()
        
        return self
    
    def to_xyz(self, filename: str = f"complex.xyz", freeze = False):

        substrate_atoms = self.substrate.atoms
        ligand_atoms = self.ligand.atoms
        n_atoms = len(substrate_atoms.append(ligand_atoms))
        
        with open(filename, 'w') as file:
            file.write(f"{n_atoms}\n\n")
            

        if freeze == False:

            atoms = substrate_atoms.append(ligand_atoms) 
            for atom in atoms:
                line = f"{atom.element:<2}{atom.x:>15.5f}{atom.y:>15.5f}{atom.z:>15.5f}\n"
                file.write(line)

        else:

            for atom in substrate_atoms:
                line = f"{self.element:<2}{0:^4}{atom.x:>15.5f}{atom.y:>15.5f}{atom.z:>15.5f}\n"
                file.write(line)

            for atom in ligand_atoms:
                line = f"{self.element:<2}{-1:^4}{atom.x:>15.5f}{atom.y:>15.5f}{atom.z:>15.5f}\n"
                file.write(line)
            


