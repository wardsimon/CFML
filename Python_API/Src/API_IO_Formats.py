# **************************************************************************
#
# CrysFML API IO Format
#
# @file      Src/API_IO_Formats.py
# @brief     IO Formats based on CrysFML
#
# @homepage  https://code.ill.fr/scientific-software/crysfml
# @license   GNU LGPL (see LICENSE)
# @copyright Institut Laue Langevin 2020-now
# @authors   Scientific Computing Group at ILL (see AUTHORS)
#
# **************************************************************************

import CFML_api.crysfml_api
import CFML_api.API_Crystal_Metrics
import CFML_api.API_Crystallographic_Symmetry
import CFML_api.API_Atom_TypeDef

class CIFFile():
    def __init__(self, filename):
        dict = CFML_api.crysfml_api.IO_formats_readn_set_xtal_structure(filename)
        self.cell = CFML_api.API_Crystal_Metrics.Cell.from_fortran_address(dict["Cell"])
        self.space_group = CFML_api.API_Crystallographic_Symmetry.SpaceGroup.from_fortran_address(dict["SpG"])
        self.atom_list = CFML_api.API_Atom_TypeDef.AtomList.from_fortran_address(dict["A"])
