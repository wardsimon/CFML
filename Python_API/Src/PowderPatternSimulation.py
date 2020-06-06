# **************************************************************************
#
# CrysFML API
#
# @file      Src/PowderPatternSimulation.py
# @brief     Powder pattern simulation based on CrysFML
#
# @homepage  https://code.ill.fr/scientific-software/crysfml
# @license   GNU LGPL (see LICENSE)
# @copyright Institut Laue Langevin 2020-now
# @authors   Scientific Computing Group at ILL (see AUTHORS)
#
# **************************************************************************

import CFML_api.powder_generation as powder_generation

import os
import numpy as np

class PowderPatternSimulator():       
    def __init__(self):
        self.y = None
        self.x = None

    def compute(self, file):
        # Removing suffix (if any)
        if file[:-4] == ".cif" or file[:-4] == ".cfl":
            mfile = file[:-4]

        # Check if file exists
        cfl_exists = False
        cif_exists = False
        if os.path.exists(file+".cif"):
            cif_exists = True
        if os.path.exists(file+".cfl"):
            cfl_exists = True
        if not (cfl_exists or cif_exists):
            raise IOError("No cfl or cif file with name: " + file)
        else:
            # Compute powder pattern
            ret = powder_generation.crysfml_powder_compute_powder_pattern(file)
            # Get back data 
            self.y = ret["Powder pattern"]#.copy()
            self.x = np.arange(ret["xmin"], ret["xmax"]+ret["xstep"]/2, ret["xstep"])