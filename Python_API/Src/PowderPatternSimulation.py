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
import logging

class PowderPatternSimulationSource():
    XRays = 0
    Neutrons = 1
    
class PowderPatternSimulationConditions():
    def __init__(self):
        self.title="Default Powder Pattern"
        self.lamb = 1.54056
        self.u_resolution = 0.0002
        self.v_resolution = -0.0002
        self.w_resolution = 0.012
        self.x_resolution = 0.0015
        self.theta_min = 1.00
        self.theta_step = 0.05
        self.theta_max = 135 # or 120 ? #int(2.0*asind(stlmax*1.54056))
        self.job = PowderPatternSimulationSource.Neutrons
        self.lorentzian_size = 1900.0
        self.bkg = 50.0
        # Gs ?
        
class PowderPatternSimulator():       
    def __init__(self):
        self.y = None
        self.x = None

    def compute(self, file, simulation_conditions=None):
        # Check if file exists
        file_exists = False
        if file[-4:] == ".cif" or file[-4:] == ".cfl":
            if os.path.exists(file):
                file_exists = True
                if file[-4:] == ".cif":
                    mode = "CIF"
                else:
                    mode = "CFL"
        elif os.path.exists(file+".cfl"):
                file = file + ".cfl"
                file_exists = True
                mode = "CFL"
        elif os.path.exists(file+".cif"):
                file = file + ".cif"
                file_exists = True
                mode = "CIF"
        
        if not file_exists:
            raise IOError("No file with name: " + file)
        else:
            if mode == "CIF" and simulation_conditions is None:
                logging.warning("No PowerPatternSimulationConditions provided, use default ones with CIF file")
                simulation_conditions = PowderPatternSimulationConditions()
                
            if mode == "CFL" and (not simulation_conditions is None):
                logging.warning("Use provided PowerPatternSimulationConditions instead of the ones within the CFL file")
                mode = "CF2"
            elif mode == "CFL":
                simulation_conditions = PowderPatternSimulationConditions()
                
            # Compute powder pattern
            ret = powder_generation.crysfml_powder_compute_powder_pattern(file,
                                                                          mode,
                                                                          simulation_conditions.title,
                                                                          simulation_conditions.lamb,
                                                                          simulation_conditions.u_resolution,
                                                                          simulation_conditions.v_resolution,
                                                                          simulation_conditions.w_resolution,
                                                                          simulation_conditions.x_resolution,
                                                                          simulation_conditions.theta_min,
                                                                          simulation_conditions.theta_step,
                                                                          simulation_conditions.theta_max,
                                                                          simulation_conditions.job,
                                                                          simulation_conditions.lorentzian_size,
                                                                          simulation_conditions.bkg)
            # Get back data 
            self.y = ret["Powder pattern"]#.copy()
            self.x = np.arange(ret["xmin"], ret["xmax"]+ret["xstep"]/2, ret["xstep"])