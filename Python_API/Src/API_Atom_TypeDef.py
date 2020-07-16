# **************************************************************************
#
# CrysFML API Atoms
#
# @file      Src/API_Atom_Typedef.py
# @brief     Atoms properties based on CrysFML
#
# @homepage  https://code.ill.fr/scientific-software/crysfml
# @license   GNU LGPL (see LICENSE)
# @copyright Institut Laue Langevin 2020-now
# @authors   Scientific Computing Group at ILL (see AUTHORS)
#
# **************************************************************************

import CFML_api.crysfml_api as crysfml_api

class AtomList():
    def __init__(self, address=None):
        self.__address = None
        if address is not None:
            self.__address = address

    def print_description(self):
        crysfml_api.atom_typedef_write_atom_list(self.__address)
