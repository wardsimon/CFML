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
import CFML_api.API_Atoms
import CFML_api.API_Crystal
import CFML_API.API_Crystallographic_Symmetry

def parse_cif_file(filename):
    dict = crysfml_api.IO_formats_readn_set_xtal_structure(filename)
    return Cell(address=dict["Cell"]), SpaceGroup(address=dict["SpG"]), AtomList(address=dict["A"])