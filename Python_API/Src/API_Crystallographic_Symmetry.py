# **************************************************************************
#
# CrysFML API Crystallographic Symmetry
#
# @file      Src/API_Crystallographic_Symmetry.py
# @brief     Symmetry groups properties based on CrysFML
#
# @homepage  https://code.ill.fr/scientific-software/crysfml
# @license   GNU LGPL (see LICENSE)
# @copyright Institut Laue Langevin 2020-now
# @authors   Scientific Computing Group at ILL (see AUTHORS)
#
# **************************************************************************

import CFML_api.crysfml_api as crysfml_api

class SpaceGroup():
    def __init__(self, group_id=None, address=None):
        if address is not None:
            self.__address = address
        elif group_id is not None:
            self.__address = crysfml_api.crystallographic_symmetry_set_spacegroup(group_id)["address"]
    
    def printDescription(self):
        crysfml_api.crystallographic_symmetry_write_spacegroup(self.__address)
        
    def getLatticeTranslation(self):
        return crysfml_api.crystallographic_symmetry_get_latt_trans(self.__address)["Array"]

    lattice_translation = property(getLatticeTranslation)