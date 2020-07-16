# **************************************************************************
#
# CrysFML API Atoms
#
# @file      Src/API_Atom_TypeDef.py
# @brief     Atoms properties based on CrysFML
#
# @homepage  https://code.ill.fr/scientific-software/crysfml
# @license   GNU LGPL (see LICENSE)
# @copyright Institut Laue Langevin 2020-now
# @authors   Scientific Computing Group at ILL (see AUTHORS)
#
# **************************************************************************

import CFML_api.crysfml_api
import CFML_api.FortranBindedClass

class AtomList(CFML_api.FortranBindedClass):   
    def print_description(self):
        CFML_api.crysfml_api.atom_typedef_write_atom_list(self.get_fortran_address())