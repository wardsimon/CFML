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
    def __init__(self, lattpar, lattangle, address=None):
        if address is not None:
            self.__address = address
        else :
            self.__address = crysfml_api.crystal_metrics_set_crystal_cell(lattpar,lattangle)["address"]

    #def __del__(self):
        ## TODO

    def printDescription(self):
        crysfml_api.crystal_metrics_write_crystal_cell(self.__address)

    def get_lattpar(self):
        return crysfml_api.crystal_metrics_get_cell(self.__address)["cell"]
    
    def get_lattangle(self):
        return crysfml_api.crystal_metrics_get_ang(self.__address)["ang"]

    def get_lattpar_refcode(self):
        return crysfml_api.crystal_metrics_get_lcell(self.__address)["lcell"]
    
    def get_lattangle_refcode(self):
        return crysfml_api.crystal_metrics_get_lang(self.__address)["lang"]
    
    def get_lattpar_std_dev(self):
        return crysfml_api.crystal_metrics_get_cell_std(self.__address)["cell_std"]
    
    def get_lattangle_std_dev(self):
        return crysfml_api.crystal_metrics_get_ang_std(self.__address)["ang_std"]
    
    def get_reciprocal_lattpar(self):
        return crysfml_api.crystal_metrics_get_rcell(self.__address)["rcell"]
    
    def get_reciprocal_lattangle(self):
        return crysfml_api.crystal_metrics_get_rang(self.__address)["rang"]
    
    def get_direct_metric_tensor(self):
        return crysfml_api.crystal_metrics_get_GD(self.__address)["GD"]
    
    def get_reciprocal_metric_tensor(self):
        return crysfml_api.crystal_metrics_get_GR(self.__address)["GR"]
    
    def get_crystal_to_orth_matrix(self):
        return crysfml_api.crystal_metrics_get_Cr_Orth_Cel(self.__address)["Cr_Orth_Cel"]
    
    def get_orth_to_crystal_matrix(self):
        return crysfml_api.crystal_metrics_get_Orth_Cr_Cel(self.__address)["Orth_Cr_Cel"]
    
    def get_BL_matrix(self):
        return crysfml_api.crystal_metrics_get_BL_M(self.__address)["BL_M"]
    
    def get_inv_BL_matrix(self):
        return crysfml_api.crystal_metrics_get_BL_Minv(self.__address)["BL_Minv"]
    
    def get_direct_cell_vol(self):
        return crysfml_api.crystal_metrics_get_cellvol(self.__address)["cellvol"]
    
    def get_reciprocal_cell_vol(self):
        return crysfml_api.crystal_metrics_get_rcellvol(self.__address)["rcellvol"]
    
    def get_direct_cell_vol_std_dev(self):
        return crysfml_api.crystal_metrics_get_stdvol(self.__address)["stdvol"]
    
    def get_cartesian_frame(self):
        return crysfml_api.crystal_metrics_get_CartType(self.__address)["CartType"]

    lattpar = property(get_lattpar)                             # Direct cell parameters
    lattangle = property(get_lattangle)
    lattpar_refcode = property(get_lattpar_refcode)                 # Code number for refinement in optimization procedures
    lattangle_refcode = property(get_lattangle_refcode)
    lattpar_std_dev = property(get_lattpar_std_dev)                 # Standard deviations of the cell parameters
    lattangle_std_dev = property(get_lattangle_std_dev)
    reciprocal_lattpar = property(get_reciprocal_lattpar)           # Reciprocal cell parameters
    reciprocal_lattangle = property(get_reciprocal_lattangle)
    direct_metric_tensor = property(get_direct_metric_tensor)       # Direct and reciprocal Metric Tensors
    reciprocal_metric_tensor = property(get_reciprocal_metric_tensor)
    crystal_to_orth_matrix = property(get_crystal_to_orth_matrix)   # P-Matrix transforming direct Crytal cell to Orthonormal basis
    orth_to_crystal_matrix = property(get_orth_to_crystal_matrix)  
    BL_matrix = property(get_BL_matrix)                             # Busing-Levy B-matrix (transforms hkl to a
    inv_BL_matrix = property(get_inv_BL_matrix)                     #    Cartesian system with x//a*, y in (a*,b*) and z//c
    direct_cell_vol = property(get_direct_cell_vol)                 # Direct and Reciprocal Cell volumes
    reciprocal_cell_vol = property(get_reciprocal_cell_vol)
    direct_cell_vol_std_dev = property(get_direct_cell_vol_std_dev) # Standard deviation of the cell volume
    cartesian_frame = property(get_cartesian_frame)                 # Cartesian Frame type: if CartType='A' the Cartesian Frame has x // a.
    
