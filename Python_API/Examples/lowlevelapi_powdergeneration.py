import sys
import os

import CFML_api

filename = sys.argv[1]

cif_file = CFML_api.Readn_set_Xtal_Structure(filename)
Cell = cif_file.cell
SpG = cif_file.space_group
A = cif_file.atom_list

PowPat_conditions = CFML_api.Read_PowPatConditions(filename)

hkl = CFML_api.hkl_Uni(Cell, SpG, lFriedel, 0,  PowPat_conditions.get_stlmax())

hkl.StructureFactors(SpG, A)

Pat = CFML_api.Calc_powder_pattern(PowPat_conditions, hkl)
