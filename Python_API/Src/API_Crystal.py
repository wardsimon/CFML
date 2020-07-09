# **************************************************************************
#
# CrysFML API Crystal
#
# @file      Src/API_Crystal.py
# @brief     Crystal properties based on CrysFML
#
# @homepage  https://code.ill.fr/scientific-software/crysfml
# @license   GNU LGPL (see LICENSE)
# @copyright Institut Laue Langevin 2020-now
# @authors   Scientific Computing Group at ILL (see AUTHORS)
#
# **************************************************************************

import CFML_api.crysfml_api as crysfml_api

class Cell():
    def __init__(self, address):
        self.__address = None
        if address is not None:
            self.__address = address
