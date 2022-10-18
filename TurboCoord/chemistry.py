import numpy as np
import spatialTools as st
from TurboCoord import Atom, Ligand, CoordinationComplex



def dipole(Ligand):
    """Calculate the dipole moment of a ligand"""
    dipole = np.zeros(3)
    for atom in Ligand.Atoms:
        dipole += atom.charge * (atom.coords)
    return dipole


def molecular_mechanics(Ligand):
    """
    Slightly displace or rotate the ligand if it is too close to an atom in the complex
    """
    for atom in Ligand.Atoms:
        for atom2 in Ligand.Atoms:
            if atom != atom2:
                if np.linalg.norm(atom.coords - atom2.coords) < 0.5:
                    atom.coords += np.random.rand(3) * 0.1
                    atom.coords = st.rotate(atom.coords, np.random.rand() * 0.1, np.random.rand(3))
    return Ligand



