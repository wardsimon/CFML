# **************************************************************************
#
# CrysFML API IO Format
#
# @file      Src/API_IO_Format.py
# @brief     Symmetry groups properties based on CrysFML
#
# @homepage  https://code.ill.fr/scientific-software/crysfml
# @license   GNU LGPL (see LICENSE)
# @copyright Institut Laue Langevin 2020-now
# @authors   Scientific Computing Group at ILL (see AUTHORS)
#
# **************************************************************************

import CFML_api.crysfml_api as crysfml_api
import CFML_api.API_Atoms_TypeDef
import CFML_api.API_Crystal_Metrics
import CFML_api.API_Crystallographic_Symmetry

class CIFFile():
    def __init__(self, filename):
        dict = crysfml_api.IO_formats_readn_set_xtal_structure(filename)
        self.cell = Cell(address=dict["Cell"])
        self.space_group = SpaceGroup(address=dict["SpG"])
        self.atom_list = AtomList(address=dict["A"])
