import sys
import os

import CFML_api

filename = sys.argv[1]

(Cell, SpG, A) = CFML_api.Readn_set_Xtal_Structure(filename)

PowPat_conditions = CFML_api.Read_PowPatConditions(filename)

hkl = CFML_api.hkl_Uni(Cell, SpG, lFriedel, 0,  PowPat_conditions.get_stlmax())

hkl.StructureFacors(SpG, A)

Pat = CFML_api.Calc_powder_pattern(PowPat_conditions, hkl)
