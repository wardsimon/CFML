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
    def __init(self, address=None):
        self.__address = None
        if address is not None:
            self.__address = address
