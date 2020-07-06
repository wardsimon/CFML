# **************************************************************************
#
# CrysFML API
#
# @file      Src/SymmetryGroups.py
# @brief     Symmetry groups properties based on CrysFML
#
# @homepage  https://code.ill.fr/scientific-software/crysfml
# @license   GNU LGPL (see LICENSE)
# @copyright Institut Laue Langevin 2020-now
# @authors   Scientific Computing Group at ILL (see AUTHORS)
#
# **************************************************************************

import CFML_api.crysfml_symmetry as crysfml_symmetry

class SymmetryGroups():
    def __init__(self, group_id):
        self.__index = crysfml_symmetry.create_space_group(group_id)["Index"]
    
    def printDescription(self):
        crysfml_symmetry.get_description(self.__index)
        
    def getLatticeTranslation(self):
        return crysfml_symmetry.get_latt_trans(self.__index)["Array"]
    
    def setLatticeTranslation(self, array):
        print("Unmutable")
        
    lattice_translation = property(getLatticeTranslation, setLatticeTranslation)