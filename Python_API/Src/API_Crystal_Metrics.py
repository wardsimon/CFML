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
        return crysfml_api.crystallographic_symmetry_get_cell(self.__address)["cell"]
    
    def get_lattangle(self):
        return crysfml_api.crystallographic_symmetry_get_angl(self.__address)["angl"]
    
###
    # Type, public :: Crystal_Cell_Type
    #    real(kind=cp),dimension(3)   :: cell, ang
    #    integer,      dimension(3)   :: lcell, lang
    #    real(kind=cp),dimension(3)   :: cell_std, ang_std
    #    real(kind=cp),dimension(3)   :: rcell, rang
    #    real(kind=cp),dimension(3,3) :: GD,GR
    #    real(kind=cp),dimension(3,3) :: Cr_Orth_cel
    #    real(kind=cp),dimension(3,3) :: Orth_Cr_cel
    #    real(kind=cp),dimension(3,3) :: BL_M
    #    real(kind=cp),dimension(3,3) :: BL_Minv
    #    real(kind=cp)                :: CellVol
    #    real(kind=cp)                :: StdVol
    #    real(kind=cp)                :: RCellVol
    #    character (len=2)            :: CartType
    # End Type Crystal_Cell_Type      


    # def printDescription(self):
    #     crysfml_api.crystallographic_symmetry_write_spacegroup(self.__address)
        
    # def getLatticeTranslation(self):
    #     return crysfml_api.crystallographic_symmetry_get_latt_trans(self.__address)["Array"]

    # lattice_translation = property(getLatticeTranslation)
###
