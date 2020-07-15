import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))

import CFML_api

filename = sys.argv[1]

cif_file = CFML_api.CIFFile(filename)

Cell = cif_file.cell
SpG = cif_file.space_group
A = cif_file.atom_list

Cell.print_description()
SpG.print_description()
A.print_description()

#PowPat_conditions = CFML_api.Read_PowPatConditions(filename)

#hkl = CFML_api.hkl_Uni(Cell, SpG, lFriedel, 0,  PowPat_conditions.get_stlmax())

#hkl.StructureFactors(SpG, A)

#Pat = CFML_api.Calc_powder_pattern(PowPat_conditions, hkl)
