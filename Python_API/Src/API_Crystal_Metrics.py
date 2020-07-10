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
    def __init__(self, lattpar=None, lattangle=None, address=None):
        if address is not None:
            self.__address = address
        elif lattpar is not None:
            self.__address = crysfml_api.crystal_metrics_set_crystal_cell(lattpar,lattangle)["address"]


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
      
###
