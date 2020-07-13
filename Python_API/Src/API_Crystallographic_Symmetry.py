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

import CFML_api.crysfml_api as crysfml_api

class CentricType():
    CENTRIC = 0
    ACENTRIC = 1
    CENTRIC_EDGE_AT_ORIGIN = 2
    
class SymmetryOperators():
    def __init__(self, address):
        self.__address = None
        if address is not None:
            self.__address = address

class WyckoffType():
    def __init__(self, address):
        self.__address = None
        if address is not None:
            self.__address = address

class SpaceGroup():
    """
       integer                                       :: NumSpg=0         ! Number of the Space Group
       character(len=20)                             :: SPG_Symb=" "     ! Hermann-Mauguin Symbol
       character(len=16)                             :: Hall=" "         ! Hall symbol
       character(len=90)                             :: gHall=" "        ! Generalised Hall symbol
       character(len=12)                             :: CrystalSys=" "   ! Crystal system
       character(len= 5)                             :: Laue=" "         ! Laue Class
       character(len= 5)                             :: PG=" "           ! Point group
       character(len= 5)                             :: Info=" "         ! Extra information
       character(len=90)                             :: SG_setting=" "   ! Information about the SG setting (IT,KO,ML,ZA,Table,Standard,UnConventional)
       logical                                       :: Hexa=.false.     !
       character(len= 1)                             :: SPG_lat=" "      ! Lattice type
       character(len= 2)                             :: SPG_latsy=" "    ! Lattice type Symbol
       integer                                       :: NumLat=0         ! Number of lattice points in a cell
       real(kind=cp), allocatable,dimension(:,:)     :: Latt_trans       ! Lattice translations (3,12)
       character(len=51)                             :: Bravais=" "      ! String with Bravais symbol + translations
       character(len=80)                             :: Centre=" "       ! Alphanumeric information about the center of symmetry
       integer                                       :: Centred=0        ! Centric or Acentric [ =0 Centric(-1 no at origin),=1 Acentric,=2 Centric(-1 at origin)]
       real(kind=cp), dimension(3)                   :: Centre_coord=0.0 ! Fractional coordinates of the inversion centre
       integer                                       :: NumOps=0         ! Number of reduced set of S.O.
       integer                                       :: Multip=0         ! Multiplicity of the general position
       integer                                       :: Num_gen          ! Minimum number of operators to generate the Group
       type(Sym_Oper_Type), allocatable,dimension(:) :: SymOp            ! Symmetry operators (192)
       character(len=50),   allocatable,dimension(:) :: SymopSymb        ! Strings form of symmetry operators
       type(Wyckoff_Type)                            :: Wyckoff          ! Wyckoff Information
       real(kind=cp),dimension(3,2)                  :: R_Asym_Unit=0.0  ! Asymmetric unit in real space
    """
    def __init__(self, group_id=None, address=None):
        if address is not None:
            self.__address = address
        elif group_id is not None:
            self.__address = crysfml_api.crystallographic_symmetry_set_spacegroup(group_id)["address"]
    
    def print_description(self):
        crysfml_api.crystallographic_symmetry_write_spacegroup(self.__address)
    
    def get_lattice_translation(self):
        return crysfml_api.crystallographic_symmetry_get_latt_trans(self.__address)["latt_trans"]
    
    def is_hexa(self):
        return crysfml_api.crystallographic_symmetry_get_hexa(self.__address)["hexa"]
    
    def get_number_of_space_group(self):
        return crysfml_api.crystallographic_symmetry_get_numspg(self.__address)["numspg"]
    
    def get_space_group_symbol(self):
        return crysfml_api.crystallographic_symmetry_get_spg_symb(self.__address)["spg_symb"]
    
    def get_hall_symbol(self):
        return crysfml_api.crystallographic_symmetry_get_hall(self.__address)
    
    def get_generalized_hall_symbolp(self):
        return crysfml_api.crystallographic_symmetry_get_ghall(self.__address)
    
    def get_crystal_system(self):
        return crysfml_api.crystallographic_symmetry_get_crystalsys(self.__address)
    
    def get_laue_class(self):
        return crysfml_api.crystallographic_symmetry_get_laue(self.__address)
    
    def get_point_group(self):
        return crysfml_api.crystallographic_symmetry_get_pg(self.__address)
    
    def get_extra_information(self):
        return crysfml_api.crystallographic_symmetry_get_info(self.__address)

    def get_space_group_setting(self):
        return crysfml_api.crystallographic_symmetry_get_sg_setting(self.__address)
    
    def get_space_group_lattice_type(self):
        return crysfml_api.crystallographic_symmetry_get_spg_lat(self.__address)
    
    def get_space_group_lattice_type_symbol(self):
        return crysfml_api.crystallographic_symmetry_get_spg_latsy(self.__address)
    
    def get_number_of_lattice_points(self):
        return crysfml_api.crystallographic_symmetry_get_numlat(self.__address)
    
    def get_bravais(self):
        return crysfml_api.crystallographic_symmetry_get_bravais(self.__address)
    
    def get_centre(self):
        return crysfml_api.crystallographic_symmetry_get_centre(self.__address)
    
    def get_centric(self):
        return CentricType(crysfml_api.crystallographic_symmetry_get_centric(self.__address))

    def get_inversion_centre(self):
        return crysfml_api.crystallographic_symmetry_get_centre_coord(self.__address)
    
    def get_number_of_reduced_set(self):
        return crysfml_api.crystallographic_symmetry_get_numops(self.__address)
    
    def get_multiplicity(self):
        return crysfml_api.crystallographic_symmetry_get_multip(self.__address)
    
    def get_operators_minimum_number(self):
        return crysfml_api.crystallographic_symmetry_get_num_gen(self.__address)
    
    def get_symmetry_operators(self):
        return SymmetryOperators(crysfml_api.crystallographic_symmetry_get_symop(self.__address))
    
    def get_symmetry_operators_as_text(self):
        return crysfml_api.crystallographic_symmetry_get_sympopsymb(self.__address)
    
    def get_wyckoff_type(self):
        return WyckoffType(crysfml_api.crystallographic_symmetry_get_wycoff(self.__address))
        
    def get_asymetric_unit(self):
        return crysfml_api.crystallographic_symmetry_get_r_asym_unit(self.__address)
        
        
    number_of_space_group = property(get_number_of_space_group)
    space_group_symbol = property(get_space_group_symbol)
    hall_symbol = property(get_hall_symbol)
    generalized_hall_symbolp = property(get_generalized_hall_symbolp)
    crystal_system = property(get_crystal_system)
    laue_class = property(get_laue_class)
    point_group = property(get_point_group)
    extra_information = property(get_extra_information)
    space_group_setting = property(get_space_group_setting)
    space_group_lattice_type = property(get_space_group_lattice_type)
    space_group_lattice_type_symbol = property(get_space_group_lattice_type_symbol)
    number_of_lattice_points = property(get_number_of_lattice_points)
    lattice_translation = property(get_lattice_translation)
    bravais = property(get_bravais)
    centre = property(get_centre)
    centric = property(get_centric)
    inversion_centre = property(get_inversion_centre)
    number_of_reduced_set = property(get_number_of_reduced_set)
    multiplicity = property(get_multiplicity)
    operators_minimum_number = property(get_operators_minimum_number)
    symmetry_operators = property(get_symmetry_operators)
    symmetry_operators_as_text = property(get_symmetry_operators_as_text)
    wyckoff_type = property(get_wyckoff_type)
    asymetric_unit = property(get_asymetric_unit)
