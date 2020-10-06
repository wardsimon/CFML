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
atom_list.print_description()

atom = atom_list[0]

atom.from_string('ATOM C C 0.0 0.0 0.0 0.2 0.5')
atom_list[0] = atom 
atom_list.print_description()

atom.xyz=np.asarray([0.2, 0.3, 0.4])
atom_list[0] = atom
atom_list.print_description()