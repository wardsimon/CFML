# **************************************************************************
#
# CrysFML API Crystal
#
# @file      Src/API_Crystal_Metrics.py
# @brief     Crystal properties based on CrysFML
#
# @homepage  https://code.ill.fr/scientific-software/crysfml
# @license   GNU LGPL (see LICENSE)
# @copyright Institut Laue Langevin 2020-now
# @authors   Scientific Computing Group at ILL (see AUTHORS)
#
# **************************************************************************

import CFML_api.crysfml_api
import CFML_api.FortranBindedClass

class Cell(CFML_api.FortranBindedClass):
    """ A class used to describe the crystal cell type(Crystal_Cell_type) in CFML.

    ... 
    Attributes
    ----------
    lattpar : ndarray(dtype='float32', ndim=1)
        Array containing the lattice parameters (a,b,c)
    lattangle : ndarray(dtype='float32', ndim=1)
        Array containing the lattice angles (alpha, beta, gamma)

    Methods
    -------
    print_description
        Prints the lattice cell description
    
    """
    
    def __init__(self, lattpar=None, lattangle=None):
        if lattpar is not None and lattangle is not None:
            self._set_fortran_address(CFML_api.crysfml_api.crystal_metrics_set_crystal_cell(lattpar,lattangle)["address"])
    
    def __del__(self):
        CFML_api.crysfml_api.crystal_metrics_del_crystal_cell(self.get_fortran_address())

    def print_description(self):
        CFML_api.crysfml_api.crystal_metrics_write_crystal_cell(self.get_fortran_address())

    def set_lattpar(self, lattpar):
        lattangle = CFML_api.crysfml_api.crystal_metrics_get_ang(self.get_fortran_address())["ang"]
        CFML_api.crysfml_api.crystal_metrics_del_crystal_cell(self.get_fortran_address())
        self._set_fortran_address(CFML_api.crysfml_api.crystal_metrics_set_crystal_cell(lattpar,lattangle)["address"])

    def set_lattangle(self, lattangle):
        lattpar = CFML_api.crysfml_api.crystal_metrics_get_cell(self.get_fortran_address())["cell"]
        CFML_api.crysfml_api.crystal_metrics_del_crystal_cell(self.get_fortran_address())
        self._set_fortran_address(CFML_api.crysfml_api.crystal_metrics_set_crystal_cell(lattpar,lattangle)["address"])

    def get_lattpar(self):
        return CFML_api.crysfml_api.crystal_metrics_get_cell(self.get_fortran_address())["cell"]

    def get_lattangle(self):
        return CFML_api.crysfml_api.crystal_metrics_get_ang(self.get_fortran_address())["ang"]
        
    def get_lattpar_refcode(self):
        return CFML_api.crysfml_api.crystal_metrics_get_lcell(self.get_fortran_address())["lcell"]

    def get_lattangle_refcode(self):
        return CFML_api.crysfml_api.crystal_metrics_get_lang(self.get_fortran_address())["lang"]
    
    def get_lattpar_std_dev(self):
        return CFML_api.crysfml_api.crystal_metrics_get_cell_std(self.get_fortran_address())["cell_std"]
    
    def get_lattangle_std_dev(self):
        return CFML_api.crysfml_api.crystal_metrics_get_ang_std(self.get_fortran_address())["ang_std"]
    
    def get_reciprocal_lattpar(self):
        return CFML_api.crysfml_api.crystal_metrics_get_rcell(self.get_fortran_address())["rcell"]
    
    def get_reciprocal_lattangle(self):
        return CFML_api.crysfml_api.crystal_metrics_get_rang(self.get_fortran_address())["rang"]
    
    def get_direct_metric_tensor(self):
        return CFML_api.crysfml_api.crystal_metrics_get_GD(self.get_fortran_address())["GD"]
    
    def get_reciprocal_metric_tensor(self):
        return CFML_api.crysfml_api.crystal_metrics_get_GR(self.get_fortran_address())["GR"]
    
    def get_crystal_to_orth_matrix(self):
        return CFML_api.crysfml_api.crystal_metrics_get_Cr_Orth_Cel(self.get_fortran_address())["Cr_Orth_Cel"]
    
    def get_orth_to_crystal_matrix(self):
        return CFML_api.crysfml_api.crystal_metrics_get_Orth_Cr_Cel(self.get_fortran_address())["Orth_Cr_Cel"]
    
    def get_BL_matrix(self):
        return CFML_api.crysfml_api.crystal_metrics_get_BL_M(self.get_fortran_address())["BL_M"]
    
    def get_inv_BL_matrix(self):
        return CFML_api.crysfml_api.crystal_metrics_get_BL_Minv(self.get_fortran_address())["BL_Minv"]
    
    def get_direct_cell_vol(self):
        return CFML_api.crysfml_api.crystal_metrics_get_cellvol(self.get_fortran_address())["cellvol"]
    
    def get_reciprocal_cell_vol(self):
        return CFML_api.crysfml_api.crystal_metrics_get_rcellvol(self.get_fortran_address())["rcellvol"]
    
    def get_direct_cell_vol_std_dev(self):
        return CFML_api.crysfml_api.crystal_metrics_get_stdvol(self.get_fortran_address())["stdvol"]
    
    def get_cartesian_frame(self):
        return CFML_api.crysfml_api.crystal_metrics_get_CartType(self.get_fortran_address())["CartType"]

    lattpar = property(get_lattpar, set_lattpar)                    # Direct cell parameters
    lattangle = property(get_lattangle, set_lattangle)
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
    
