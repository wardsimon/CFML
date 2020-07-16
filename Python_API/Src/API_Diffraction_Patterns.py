# **************************************************************************
#
# CrysFML API Diffraction Patters
#
# @file      Src/API_Diffraction_Patterns.py
# @brief     Diffraction Patterns utilities based on CrysFML
#
# @homepage  https://code.ill.fr/scientific-software/crysfml
# @license   GNU LGPL (see LICENSE)
# @copyright Institut Laue Langevin 2020-now
# @authors   Scientific Computing Group at ILL (see AUTHORS)
#
# **************************************************************************

import CFML_api.crysfml_api
import CFML_api.FortranBindedClass

class DiffractionPattern(CFML_api.FortranBindedClass):
    def __init__(self, simulation_conditions, reflection_list):
        self._set_fortran_address(CFML_api.crysfml_api.diffraction_patterns_compute_powder_pattern(
            simulation_conditions, reflection_list))
        
