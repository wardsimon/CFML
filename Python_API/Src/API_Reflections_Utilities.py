# **************************************************************************
#
# CrysFML API Reflection Utilities
#
# @file      Src/API_Reflections_Utilities.py
# @brief     Reflections utilities based on CrysFML
#
# @homepage  https://code.ill.fr/scientific-software/crysfml
# @license   GNU LGPL (see LICENSE)
# @copyright Institut Laue Langevin 2020-now
# @authors   Scientific Computing Group at ILL (see AUTHORS)
#
# **************************************************************************

import CFML_api.crysfml_api as crysfml_api

class ReflectionList():
    def __init__(self, cell, spg, lfriedel, value1, value2):
        self.__address = crysfml_api.reflections_utilities_hkl_uni_reflist(
            cell.as_fortran_object(), spg.as_fortran_object(), lfriedel, value1, value2)["address"]
    
    def compute_structure_factors(self, space_group, atom_list):
        crysfml_api.structure_factors_structure_factors(
            atom_list.as_fortran_object(), space_group.as_fortran_object(), self.__address)
        crysfml_api.structure_factors_write_structure_factors(self.__address)
