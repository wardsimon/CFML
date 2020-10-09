import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))

import matplotlib.pyplot as plt
import numpy as np

import CFML_api

filename = "../Examples/Data/SrTiO3.cif"
cif_file = CFML_api.CIFFile(filename )

cell = cif_file.cell
space_group = cif_file.space_group
atom_list = cif_file.atom_list
print("========\nInitial atom list")
atom_list.print_description()

# Switch atoms
atom1 = atom_list[0]
atom2 = atom_list[1]
atom_list[0]=atom2
atom_list[1]=atom1
print("========\nSwitched atom list")
atom_list.print_description()

# Modify atom position
atom = atom_list[0]
atom.xyz=np.asarray([0.2, 0.3, 0.4])
atom_list[0] = atom
print("========\nModify 1st atom position")
atom_list.print_description()

# Modify atom with string
atom = atom_list[0]
atom.from_string('ATOM C C 0.0 0.0 0.0 0.2 0.5')
atom_list[0] = atom
print("========\nModified 1st atom")
atom_list.print_description()

# Create atom from string
atom2 = CFML_api.Atom('ATOM C C 0.0 0.0 0.0 0.2 0.5')
print("========\nAtom created from string")
print(atom2)

# Replace atom
print("========\nReplace 4th atom in atom list")
atom_list[3] = atom2
atom_list.print_description()