# **************************************************************************
#
# CrysFML API IO Format
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
        self.__address = crysfml_api.reflections_utilities_hkl_uni_reflist(cell, spg, lfriedel, value1, value2)["address"]
