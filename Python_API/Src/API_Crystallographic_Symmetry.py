# **************************************************************************
#
# CrysFML API Crystallographic Symmetry
#
# @file      Src/API_Crystallographic_Symmetry.py
# @brief     Symmetry groups properties based on CrysFML
#
# @homepage  https://code.ill.fr/scientific-software/crysfml
# @license   GNU LGPL (see LICENSE)
# @copyright Institut Laue Langevin 2020-now
# @authors   Scientific Computing Group at ILL (see AUTHORS)
#
# **************************************************************************

import CFML_api.crysfml_api
import CFML_api.FortranBindedClass

class CentricType():
    CENTRIC = 0
    ACENTRIC = 1
    CENTRIC_EDGE_AT_ORIGIN = 2
    
class SymmetryOperators(CFML_api.FortranBindedClass):
    pass

class WyckoffType(CFML_api.FortranBindedClass):
    pass

class SpaceGroup(CFML_api.FortranBindedClass):
    """ A class used to describe a Crystallographic Space Group type(Space_group_type) in CFML

    ...
    Attributes
    ----------
    group_id : int
       Space group number listed in the International Tables for Crystallography

    Methods
    -------
    print_description
        Prints the lattice cell description

    """
    def __init__(self, group_id=None):
        """
        Parameters
        ----------
        group_id : int
            Space group number
        """
        
        if group_id is not None:
            self._set_fortran_address(CFML_api.crysfml_api.crystallographic_symmetry_set_spacegroup(group_id)["address"])
    
    def __del__(self):
        CFML_api.crysfml_api.crystallographic_symmetry_del_spacegroup(self.get_fortran_address())
    
    def print_description(self):
        """ Prints the lattice cell description """
        
        CFML_api.crysfml_api.crystallographic_symmetry_write_spacegroup(self.get_fortran_address())

    @property
    def lattice_translation(self):
        """
        Lattice translation matrix
        """
        return CFML_api.crysfml_api.crystallographic_symmetry_get_latt_trans(self.get_fortran_address())["latt_trans"]

    @property
    def is_hexa(self):
        """
        Is space group hexagonal
        """
        return CFML_api.crysfml_api.crystallographic_symmetry_get_hexa(self.get_fortran_address())["hexa"]

    @property
    def space_group_number(self):
        """
        Number of the space group
        """
        return CFML_api.crysfml_api.crystallographic_symmetry_get_numspg(self.get_fortran_address())["numspg"]

    @property
    def space_group_symbol(self):
        """
        Hermann-Mauguin Symbol (str)
        """
        return CFML_api.crysfml_api.crystallographic_symmetry_get_spg_symb(self.get_fortran_address())["spg_symb"]
    
    @property
    def hall_symbol(self):
        """
        Hall symbol
        """
        return crysfml_api.crystallographic_symmetry_get_hall(self.__address)["hall"]
    
    @property
    def generalized_hall_symbol(self):
        """
        Generalised Hall symbol
        """
        return crysfml_api.crystallographic_symmetry_get_ghall(self.__address)["ghall"]
    
    @property
    def crystal_system(self):
        """
        Crystal System
        """
        return crysfml_api.crystallographic_symmetry_get_CrystalSys(self.__address)["CrystalSys"]
    
    @property
    def laue_class(self):
        """
        Laue class
        """
        return crysfml_api.crystallographic_symmetry_get_Laue(self.__address)["Laue"]
    
    @property
    def point_group(self):
        """
        Point group
        """
        return crysfml_api.crystallographic_symmetry_get_PG(self.__address)["PG"]
    
    @property
    def extra_information(self):
        """
        Extra Information
        """
        return crysfml_api.crystallographic_symmetry_get_Info(self.__address)["Info"]
    
    @property
    def space_group_setting(self):
        """
        Information about the SG setting (IT,KO,ML,ZA,Table,Standard,UnConventional)
        """
        return crysfml_api.crystallographic_symmetry_get_SG_setting(self.__address)["SG_setting"]
    
    @property
    def space_group_lattice_type(self):
        """
        Lattice type
        """
        return crysfml_api.crystallographic_symmetry_get_SPG_lat(self.__address)["SPG_lat"]
    
    @property
    def space_group_lattice_type_symbol(self):
        """
        Lattice type Symbol
        """
        return crysfml_api.crystallographic_symmetry_get_SPG_latsy(self.__address)["SPG_latsy"]
    
    @property
    def number_of_lattice_points(self):
        """
        Number of lattice points in a cell
        """
        return crysfml_api.crystallographic_symmetry_get_NumLat(self.__address)["NumLat"]
    
    @property
    def bravais(self):
        """
        String with Bravais symbol + translations
        """
        return crysfml_api.crystallographic_symmetry_get_bravais(self.__address)["bravais"]
    
    @property
    def centre(self):
        """
        Alphanumeric information about the center of symmetry
        """
        return crysfml_api.crystallographic_symmetry_get_centre(self.__address)["centre"]
    
    @property
    def centric(self):
        """
        Centric or Acentric [ =0 Centric(-1 no at origin),=1 Acentric,=2 Centric(-1 at origin)]
        """
        return crysfml_api.crystallographic_symmetry_get_centred(self.__address)["centred"]
    
    @property
    def inversion_centre(self):
        """
        Fractional coordinates of the inversion centre
        """
        return crysfml_api.crystallographic_symmetry_get_centre_coord(self.__address)["centre_coord"]
    
    @property
    def number_of_symops(self):
        """
        Number of reduced set of Symmetry Operators
        """
        return crysfml_api.crystallographic_symmetry_get_numops(self.__address)["numops"]
    
    @property
    def multiplicity(self):
        """
        Multiplicity of the general position
        """
        return crysfml_api.crystallographic_symmetry_get_multip(self.__address)["multip"]
    
    @property
    def operators_minimum_number(self):
        """
        Minimum number of operators to generate the Group
        """
        return crysfml_api.crystallographic_symmetry_get_num_gen(self.__address)["num_gen"]
    
    @property
    def symmetry_operators(self):
        """
        Symmetry operators (192)
        """
        return crysfml_api.crystallographic_symmetry_get_SymOp(self.__address)["SymOp"]
    
    @property
    def symmetry_operators_as_text(self):
        """
        String form of symmetry operators
        """
        return crysfml_api.crystallographic_symmetry_get_SymopSymb(self.__address)["SymopSymb"]
    
    @property
    def wyckoff_info(self):
        """
        Wyckoff Information
        """
        return crysfml_api.crystallographic_symmetry_get_wyckoff(self.__address)["wyckoff"]
    
    @property
    def asymetric_unit(self):
        """
        Asymmetric unit in real space
        """
        return crysfml_api.crystallographic_symmetry_get_R_asym_unit(self.__address)["R_asym_unit"]
