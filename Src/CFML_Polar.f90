!!-------------------------------------------------------
!!---- Crystallographic Fortran Modules Library (CrysFML)
!!-------------------------------------------------------
!!---- The CrysFML project is distributed under LGPL. In agreement with the
!!---- Intergovernmental Convention of the ILL, this software cannot be used
!!---- in military applications.
!!----
!!---- Copyright (C) 1999-2012  Institut Laue-Langevin (ILL), Grenoble, FRANCE
!!----                          Universidad de La Laguna (ULL), Tenerife, SPAIN
!!----                          Laboratoire Leon Brillouin(LLB), Saclay, FRANCE
!!----
!!---- Authors: Juan Rodriguez-Carvajal (ILL)
!!----          Javier Gonzalez-Platas  (ULL)
!!----
!!---- Contributors: Laurent Chapon     (ILL)
!!----               Marc Janoschek     (Los Alamos National Laboratory, USA)
!!----               Oksana Zaharko     (Paul Scherrer Institute, Switzerland)
!!----               Tierry Roisnel     (CDIFX,Rennes France)
!!----               Eric Pellegrini    (ILL)
!!----
!!---- This library is free software; you can redistribute it and/or
!!---- modify it under the terms of the GNU Lesser General Public
!!---- License as published by the Free Software Foundation; either
!!---- version 3.0 of the License, or (at your option) any later version.
!!----
!!---- This library is distributed in the hope that it will be useful,
!!---- but WITHOUT ANY WARRANTY; without even the implied warranty of
!!---- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!!---- Lesser General Public License for more details.
!!----
!!---- You should have received a copy of the GNU Lesser General Public
!!---- License along with this library; if not, see <http://www.gnu.org/licenses/>.
!!----
!!---- This particular module has been developed by:
!!---- Juan Rodriguez-Carvajal & Marc Janoschek & Oksana Zaharko (OZ)
!!----
!!---- MODULE: CFML_Polarimetry
!!----
!!----   INFO: Subroutines and Functions to calculate the polarisation tensor
!!----   as it will be measured. It uses matrices defined in CFML_Crystal_Metrics in
!!----   order to calculate the polar tensor with respect to the coordinate
!!----   frame defined in the Blume equations (Phys. Rev. Vol. 130 p.1670-1676,
!!----   1963, see also the definitions below in magn_Inter_Vec_PF). As input
!!----   the nuclear structure factor, the magnetic interaction vector with
!!----   respect to the crystal frame and the matrices defined in CFML_Crystal_Metrics
!!----   for the crystal frame are needed.
!!----
!!---- HISTORY
!!----    Updates:June     - 2012: Corrections, Blume Equations in Crystallographic Frame (JRC)
!!----            January  - 2012: General revision of variables introduced by OZ (JRC)
!!----            November - 2011: New subroutines and calculations with domains (OZ)
!!----            April    - 2008:
!!----            December - 2006: Added function Write_Polar_line for more convenient
!!----                             output of matrices of many reflections in one file
!!----            April    - 2005: Created by MJ and revised by JRC
!!----
!!---- DEPENDENCIES
!!--++    Use CFML_GlobalDeps,                 only: cp, tpi
!!--++    Use CFML_Crystal_Metrics,            only: Set_Crystal_Cell, Crystal_Cell_type, Cart_Vector
!!--++    Use CFML_Math_3D,                    only: Cross_Product
!!--++    use CFML_Crystal_Metrics,            only: Set_Crystal_Cell, Crystal_Cell_type, Cart_Vector
!!--++    use CFML_Magnetic_Structure_Factors, only: MagHD_Type
!!--++    use CFML_Magnetic_Symmetry,          only: Magnetic_domain_type
!!----
!!---- VARIABLES
!!----    POLAR_INFO_TYPE
!!----    POLAR_MATRIX_TYPE
!!----    POLAR_CALC_TYPE
!!----    POLAR_CALC_LIST_TYPE
!!----    Polar_CalcMulti_List_type
!!----    POLAR_OBS_TYPE
!!----    POLAR_OBS_LIST_TYPE
!!----    Polar_ObsMulti_List_type
!!----    Polar_Calc_sVs_type
!!----    Polar_Calc_sVs_List_type
!!----    Polar_CalcMulti_sVs_List_type

!!----
!!---- PROCEDURES
!!----    Functions:
!!--++       IM_NM_Y                   [Private]
!!--++       IM_NM_Z                   [Private]
!!--++       MAG_Y                     [Private]
!!--++       MAG_Z                     [Private]
!!--++       MAGN_INTER_VEC_PF         [Private]
!!--++       MM                        [Private]
!!--++       NUC_CONTR                 [Private]
!!--++       REAL_NM_Y                 [Private]
!!--++       REAL_NM_Z                 [Private]
!!--++       TCHIRAL                   [Private]
!!----
!!----    Subroutines:
!!----       CALC_POLAR_DOM  (OZ adapted from SET_POLAR_INFO by including domain information, Feb 2009)
!!----       SET_POLAR_INFO  (Useful only for theoretical monodomains)
!!----       WRITE_POLAR_INFO
!!----       WRITE_POLAR_LINE
!!----       Calc_Polar_Dom_Efficiency
!!----       Calc_Polar_CrSec
!!----       Calc_Polar
!!----       Get_Pol_Tensor_Pc
!!
 Module CFML_Polarimetry
    !---- Used External Modules ----!
    Use CFML_GlobalDeps,                 only: cp, tpi
    Use CFML_Crystal_Metrics,            only: Set_Crystal_Cell, Crystal_Cell_type, Cart_Vector
    Use CFML_Math_3D,                    only: Cross_Product,Tensor_Product,Mat_Cross,Invert_A
    Use CFML_Math_General,               only: atan2d,sind,cosd
    use CFML_Magnetic_Structure_Factors, only: MagHD_Type
    use CFML_Magnetic_Symmetry,          only: Magnetic_domain_type
    use CFML_Geometry_SXTAL,             only: Phi_mat,Chi_mat, Psi_mat,Get_Angs_NB

    !---- Variables ----!
    implicit none

    private

    !---- List of public overloaded operators ----!

    !---- List of public functions ----!

    !---- List of public subroutines ----!
    public  :: Calc_Polar_Dom, Set_Polar_Info, Write_Polar_Info, Write_Polar_line, &
               Calc_Polar_Dom_Efficiency, Calc_Polar_CrSec, Calc_Polar, Get_Pol_Tensor_Pc

    !---- List of private Operators ----!

    !---- List of private functions ----!
    private :: Magn_Inter_Vec_PF, Nuc_Contr, Mag_Y, Mag_Z, Real_Nm_Y, Real_Nm_Z, &
               Im_Nm_Y, Im_Nm_Z, Tchiral, Mm

    !---- List of private subroutines ----!


    !---- Definitions ----!

    !!----
    !!---- TYPE :: Polar_calc_type
    !!--..
    !!---- Type, public :: Polar_calc_type
    !!----     real(kind=cp), dimension (3)        :: H     ! Scattering vector in hkl
    !!----     real(kind=cp), dimension (3)        :: SPV   ! Second vector in Scattering plane apart of scattering vector to define plane
    !!----     Type(Crystal_Cell_Type)             :: Cell  ! Unit Cell of Crystal
    !!----     real(kind=cp)                       :: P     ! magnitude of initial polarisation vector
    !!----     complex(kind=cp), dimension (3,2,24):: MiV   ! magnetic interaction vector
    !!----     complex(kind=cp)                    :: NSF   ! nuclear structure factor
    !!----     real(kind=cp)                       :: NC    ! nuclear scattering contribution
    !!----     real(kind=cp), dimension (2,24)     :: MY    ! magnetic contribution along y
    !!----     real(kind=cp), dimension (2,24)     :: MZ    ! magnetic contribution along z
    !!----     real(kind=cp), dimension (2,24)     :: RY    ! real part of nuclear magnetic interference term along y
    !!----     real(kind=cp), dimension (2,24)     :: RZ    ! real part of nuclear magnetic interference term along z
    !!----     real(kind=cp), dimension (2,24)     :: IY    ! imaginary part of nuclear magnetic interference term along y
    !!----     real(kind=cp), dimension (2,24)     :: IZ    ! imaginary part of nuclear magnetic interference term along y
    !!----     real(kind=cp), dimension (2,24)     :: TC    ! chiral contribution
    !!----     real(kind=cp), dimension (2,24)     :: MM    ! magnetic-magnetic interference term
    !!----     real(kind=cp), dimension (3,2,24)   :: CS    ! the three different elastic cross-sections depending on the direction of the initial polar vector
    !!----     real(kind=cp), dimension (3,3)      :: Pij   ! the polarisation tensor
    !!---- End Type Polar_calc_type
    !!----
    !!---- Update: February 2009 (OZ)
    !!

    !!
    Type, public :: Polar_calc_type
        real(kind=cp), dimension (3)        :: H=0.0
        real(kind=cp), dimension (3)        :: SPV=0.0
        type(Crystal_Cell_Type)             :: Cell
        real(kind=cp)                       :: P=0.0
        complex(kind=cp), dimension (3,2,24):: MiV=0.0
        complex(kind=cp)                    :: NSF=0.0
        real(kind=cp)                       :: NC=0.0
        real(kind=cp), dimension (2,24)     :: MY=0.0
        real(kind=cp), dimension (2,24)     :: MZ=0.0
        real(kind=cp), dimension (2,24)     :: RY=0.0
        real(kind=cp), dimension (2,24)     :: RZ=0.0
        real(kind=cp), dimension (2,24)     :: IY=0.0
        real(kind=cp), dimension (2,24)     :: IZ=0.0
        real(kind=cp), dimension (2,24)     :: TC=0.0
        real(kind=cp), dimension (2,24)     :: MM=0.0
        real(kind=cp), dimension (3,2,24)   :: CS=0.0
        real(kind=cp), dimension (3,3)      :: Pij=0.0
    End Type Polar_calc_type

    !!----
    !!---- TYPE :: Polar_Calc_List_type
    !!--..
    !!----     integer                                         :: NRef  ! Number of Reflections
    !!----     type(Polar_calc_type),allocatable, dimension(:) :: Polari ! Observed Polarisation tensor for the Reflection List
    !!---- End Type Polar_Calc_List_type
    !!----
    !!---- Update: Februar 2009 OZ
    !!
    Type, public :: Polar_Calc_List_type
       integer                                         :: NRef  ! Number of Reflections
       type(Polar_calc_type),allocatable, dimension(:) :: Polari ! Observed Polarisation tensor for the Reflection List
    End Type Polar_Calc_List_type

    !!
    !!---- TYPE :: Polar_CalcMulti_List_type
    !!
    !!---- Created: February 2012 OZ
    !!
    Type, public :: Polar_CalcMulti_List_type
       integer                                              :: Nset  ! Number of Datasets
       type(Polar_Calc_List_type),allocatable, dimension(:) :: Polarilist ! Calculated Polarisation tensors for NRef
    End Type Polar_CalcMulti_List_type

    !!----
    !!---- TYPE :: POLAR_INFO_TYPE
    !!--..
    !!---- Type, public :: Polar_Info_type
    !!----     real(kind=cp), dimension (3)    :: H     ! Scattering vector in hkl
    !!----     real(kind=cp), dimension (3)    :: SPV   ! Second vector in Scattering plane apart of scattering vector to define plane
    !!----     type(crystal_cell_type)         :: Cell  ! Unit Cell of Crystal
    !!----     real(kind=cp)                   :: P     ! magnitude of initial polarisation vector
    !!----     complex(kind=cp), dimension (3) :: MiV   ! magnetic interaction vector
    !!----     complex(kind=cp)                :: NSF   ! nuclear structure factor
    !!----     real(kind=cp)                   :: NC    ! nuclear scattering contribution
    !!----     real(kind=cp)                   :: MY    ! magnetic contribution along y
    !!----     real(kind=cp)                   :: MZ    ! magnetic contribution along z
    !!----     real(kind=cp)                   :: RY    ! real part of nuclear magnetic interference term along y
    !!----     real(kind=cp)                   :: RZ    ! real part of nuclear magnetic interference term along z
    !!----     real(kind=cp)                   :: IY    ! imaginary part of nuclear magnetic interference term along y
    !!----     real(kind=cp)                   :: IZ    ! imaginary part of nuclear magnetic interference term along y
    !!----     real(kind=cp)                   :: TC    ! chiral contribution
    !!----     real(kind=cp)                   :: MM    ! magnetic-magnetic interference term
    !!----     real(kind=cp), dimension (3)    :: CS    ! the three different elastic cross-sections depending on the direction of the initial polar vector
    !!----     real(kind=cp), dimension (3,3)  :: Pij   ! the polarisation tensor
    !!---- End Type Polar_Info_type
    !!----
    !!---- Update: April 2008
    !!
    Type, public :: Polar_Info_type
       real(kind=cp), dimension (3)     :: H
       real(kind=cp), dimension (3)     :: SPV
       type(crystal_cell_type)          :: Cell
       real(kind=cp)                    :: P
       complex(kind=cp), dimension (3)  :: MiV
       complex(kind=cp)                 :: NSF
       real(kind=cp)                    :: NC
       real(kind=cp)                    :: MY
       real(kind=cp)                    :: MZ
       real(kind=cp)                    :: RY
       real(kind=cp)                    :: RZ
       real(kind=cp)                    :: IY
       real(kind=cp)                    :: IZ
       real(kind=cp)                    :: TC
       real(kind=cp)                    :: MM
       real(kind=cp), dimension (3)     :: CS
       real(kind=cp), dimension (3,3)   :: PIJ
    End Type Polar_Info_type

    !!----
    !!---- TYPE :: Polar_obs_type
    !!--..
    !!----     real(kind=cp), dimension(3)  :: H      ! H +/- k
    !!----     real(kind=cp), dimension(3)  :: SPV
    !!----     real(kind=cp)                :: Pin
    !!----     real(kind=cp), dimension(3,3):: oPij   ! the observed polarisation tensor
    !!----     real(kind=cp), dimension(3,3):: soPij  ! the Sigma of polarisation tensor
    !!----     real(kind=cp), dimension(3,3):: woPij  ! the weight 1/Sigma**2
    !!---- End Type Polar_obs_type
    !!----
    !!---- Updated: February 2012 OZ
    !!
    Type, public :: Polar_obs_type
       Real(Kind=Cp), Dimension(3)     :: H
       Real(Kind=Cp), Dimension(3)     :: SPV
       Real(Kind=Cp)                   :: P
       Real(Kind=Cp), Dimension(3,3)   :: oPij
       Real(Kind=Cp), Dimension(3,3)   :: soPij
       Real(Kind=Cp), Dimension(3,3)   :: woPij
    End Type Polar_obs_type

    !!----
    !!---- TYPE :: Polar_Obs_List_type
    !!--..
    !!----     integer                                         :: NRef   ! Number of Reflections
    !!----     type(Polar_obs_type),allocatable, dimension(:)  :: Polaro ! Observed Polarisation tensor for the Reflection List
    !!---- End Type Polar_Obs_List_type
    !!----
    !!---- Update: Februar 2009 OZ
    !!
    Type, public :: Polar_Obs_List_type
       integer                                         :: NRef
       type(Polar_obs_type),allocatable, dimension(:)  :: Polaro
    End Type Polar_Obs_List_type

    !!
    !!---- TYPE :: Polar_ObsMulti_List_type
    !!
    !!---- Created: November 2011 OZ
    !!
    Type, public :: Polar_ObsMulti_List_type
       integer                                              :: Nset       ! Number of Datasets
       type(Polar_Obs_List_type),allocatable, dimension(:)  :: Polarolist ! Observed Polarisation tensor for the Reflection List
    End Type Polar_ObsMulti_List_type

    !!
    !!---- TYPE :: Polar_Calc_sVs_type
    !!
    !!---- Created: November 2011 OZ
    !!
    Type, public :: Polar_Calc_sVs_type
        real(kind=cp), dimension (3)      :: H     ! Scattering vector in hkl
        real(kind=cp), dimension (3)      :: SPV   ! Second vector in Scattering plane
        Type(Crystal_Cell_Type)           :: Cell  ! Unit Cell of Crystal
        real(kind=cp)                     :: P     ! Polarisation
        real(kind=cp), dimension (3,3)    :: Pij   ! Calculated Polarisation tensor
    End Type Polar_Calc_sVs_type

    !!----
    !!---- TYPE :: Polar_Calc_sVs_List_type
    !!
    !!---- Created: November 2011 OZ
    !!
    Type, public :: Polar_Calc_sVs_List_type
       integer                                             :: NRef      ! Number of Reflections
       type(Polar_Calc_sVs_type),allocatable, dimension(:) :: PolarisVs ! Calculated Polarisation tensor for the Reflection List
    End Type Polar_Calc_sVs_List_type

    !!
    !!---- TYPE :: Polar_CalcMulti_sVs_List_type
    !!
    !!---- Created: November 2011 OZ
    !!
    Type, public :: Polar_CalcMulti_sVs_List_type
       integer                                                  :: Nset          ! Number of Datasets
       type(Polar_Calc_sVs_List_type),allocatable, dimension(:) :: PolarisVslist ! Calculated Polarisation tensors for NRef
    End Type Polar_CalcMulti_sVs_List_type

 Contains

    !-------------------!
    !---- Functions ----!
    !-------------------!

    !!--++
    !!--++ Function Im_Nm_Y(Nsf, MiV_Pf) Result(I_Nm_Y,B_Q)
    !!--++    complex(kind=cp),               intent(in) :: NSF      !  In  -> Nuclear Structure Factor
    !!--++    complex(kind=cp), dimension(3), intent(in) :: MiV_PF   !  In  -> Magnetic Interaction Vector in polarisation frame
    !!--++    character(len=*),optional,      intent(in) :: B_Q      !  In  -> The calculation is done with original Blume equation Q_B= -Q_cryst
    !!--++    real(kind=cp)                              :: I_NM_Y   !  Out -> Imaginary part of nuclear-magnetic interference contribution along Y
    !!--++
    !!--++    (Private)
    !!--++    Calculates the imaginary part of the nuclear-magnetic interference contribution along Y
    !!--++    to scattering in the polarisation coordinate frame according to the Blume equations
    !!--++
    !!--++ Updated: April - 2005, June-2012 (JRC)
    !!
    Function Im_Nm_Y(Nsf, MiV_Pf,B_Q) Result(I_Nm_Y)
       !---- Argument ----!
       complex(kind=cp),               intent( in)  :: NSF
       complex(kind=cp), dimension(3), intent( in)  :: MiV_PF
       character(len=*),optional,      intent(in)   :: B_Q
       real(kind=cp)                                :: I_Nm_Y
       !---- Local variables ----!
       real(kind=cp) :: s

       s=-1.0
       if(present(B_Q)) s=1.0
       I_NM_Y = s*Aimag(NSF * Conjg(MiV_PF(2)) - Conjg(NSF) * MiV_PF(2))

       return
    End Function  Im_Nm_Y

    !!--++
    !!--++ Real Function Im_Nm_Z(Nsf, MiV_Pf,B_Q) Result(I_Nm_Z)
    !!--++    complex(kind=cp),              intent(in) :: NSF     !  In -> Nuclear Structure Factor
    !!--++    complex(kind=cp), dimension(3),intent(in) :: MiV_PF  !  In -> Magnetic Interaction Vector in polarisation frame
    !!--++    character(len=*),optional,     intent(in) :: B_Q     !  In  -> The calculation is done with original Blume equation Q_B= -Q_cryst
    !!--++    real(kind=cp)                             :: I_NM_Z  !  Out-> Imaginary part of nuclear magnetic interference contribution along Z
    !!--++
    !!--++    (Private)
    !!--++    Calculates the imaginary part of the nuclear magnetic interference contribution along Z
    !!--++    to scattering in the polarisation coordinate frame according to the Blume equations
    !!--++
    !!--++ Updated: April - 2005, June-2012 (JRC)
    !!
    Function Im_Nm_Z(Nsf, MiV_Pf,B_Q) Result(I_Nm_Z)
       !---- Argument ----!
       complex(kind=cp),              intent(in) :: NSF
       complex(kind=cp), dimension(3),intent(in) :: MiV_PF
       character(len=*),optional,     intent(in) :: B_Q
       real(kind=cp)                             :: I_NM_Z
       !---- Local variables ----!
       real(kind=cp) :: s

       s=-1.0
       if(present(B_Q)) s=1.0

       I_NM_Z = s*Aimag(NSF * Conjg(MiV_PF(3)) - Conjg(NSF) * MiV_PF(3))

       return
    End Function  Im_Nm_Z

    !!--++
    !!--++ Function Mag_Y(MiV_Pf) Result(My)
    !!--++    complex(kind=cp), dimension(3), intent( in) :: MiV_PF !  In  -> Magnetic Interaction Vector in polarisation frame
    !!--++    real(kind=cp)                               :: MY     !  Out -> Magnetic contribution along Y
    !!--++
    !!--++    (Private)
    !!--++    Calculates the magnetic contribution along Y to scattering in the polarisation
    !!--++    coordinate frame according to the Blume equations
    !!--++
    !!--++ Update: April - 2005
    !!
    Function Mag_Y(MiV_Pf) Result(My)
       !---- Argument ----!
       complex(kind=cp), dimension(3), intent( in)  :: MiV_PF
       real(kind=cp)                                :: MY

       !---- Local variables ----!

       MY = MiV_PF(2) * Conjg(MiV_PF(2))

       return
    End Function  Mag_Y

    !!--++
    !!--++ Function Mag_Z(MiV_Pf) Result(Mz)
    !!--++    complex(kind=cp), dimension(3), intent( in):: MiV_PF  !  In  -> Magnetic Interaction Vector in polarisation frame
    !!--++    real(kind=cp)                     :: MZ      !  Out -> Magnetic contribution along Z
    !!--++
    !!--++    (Private)
    !!--++    Calculates the magnetic contribution along Z to scattering in the polarisation
    !!--++    coordinate frame according to the Blume equations
    !!--++
    !!--++ Update: April - 2005
    !!
    Function Mag_Z(MiV_Pf) Result(Mz)
       !---- Argument ----!
       complex(kind=cp), dimension(3), intent( in) :: MiV_PF
       real(kind=cp)                               :: MZ

       !---- Local variables ----!

       MZ = MiV_PF(3) * Conjg(MiV_PF(3))

       return
    End Function  Mag_Z

    !!--++
    !!--++ Function Magn_Inter_Vec_Pf(MiV,H,Spv, Cell) Result(MiV_Pf)
    !!--++    complex(kind=cp),dimension(3), intent( in) :: MiV            !  In -> Magnetic Interaction Vector in Crystal Cartesian Coordinates
    !!--++    real(kind=cp),   dimension(3), intent( in) :: H              !  In -> Scattering Vector in hkl
    !!--++    real(kind=cp),   dimension(3), intent( in) :: SPV            !  In -> Second Scattering plane vector in hkl
    !!--++    Type (Crystal_Cell_Type),      intent(in)  :: Cell           !  In -> Cell variable which holds transformation matrices
    !!--++    complex(kind=cp), dimension(3)             :: MiV_PF         !  Out -> Magnetic Interaction Vector in polarisation coordinate frame
    !!--++
    !!--++    (Private)
    !!--++    Calculates the magnetic interaction vector in the polarisation coordinate frame according to the Blume equations
    !!--++    and therefore depending on the scattering vector.
    !!--++
    !!--++    Polarisation coordinate frame according to Blume
    !!--++    X  || scattering vector Q     (where Q is the scattering Vector in cartesian real space coordinates, it will be calculated from H and matrices in Cell)
    !!--++    Y _|_ scattering vector Q in scattering plane
    !!--++    Z _|_ scattering vector Q out of scattering plane (ATTENTION: This choice is not non-ambiguous, there are always two possible choices
    !!--++                                                       for a right handed coordinate frame which will fullfil this condition!!!)
    !!--++
    !!--++                           Y
    !!--++                          /|\
    !!--++                           |
    !!--++                   Q       | Z
    !!--++            ____________ _\o_____\ X
    !!--++            \             /      /
    !!--++             \           /
    !!--++              \         /
    !!--++               \       /
    !!--++                \     /  K_f
    !!--++             K_i \   /
    !!--++                 _\|/_
    !!--++
    !!--++    Therefore the right handed coordinate frame will be explicitly chosen like this:
    !!--++    X := Q/|Q|               where Q is the scattering Vector in cartesian real space coordinates
    !!--++    Z := (Q x SV)/|(Q x SV)| where SV is a second vector in the scattering plane in cartesian real space coordinates
    !!--++    Y := (Z x X)
    !!--++
    !!--++    ATTENTION: Be aware that the choice of SV with respect to Q will decide which of the two possible right handed coordinates fullfilling
    !!--++               the conditions above will be used!!!
    !!--++
    !!--++
    !!--++ Updated: April - 2005, June-2012 (JRC)
    !!
    Function Magn_Inter_Vec_Pf(MiVC,H,Spv, Cell) Result(MiV_Pf)
       !---- Argument ----!
       complex(kind=cp), dimension(3), intent(in) :: MiVC   !Must be provided in Crystallographic
       real(kind=cp),    dimension(3), intent(in) :: H      !Cartesian frame
       real(kind=cp),    dimension(3), intent(in) :: SPV
       Type (Crystal_Cell_Type),       intent(in) :: Cell
       complex(kind=cp), dimension(3)             :: MiV_PF

       !---- Local variables ----!
       real(kind=cp), dimension (3)            :: QSV,Q,SV,X,Y,Z
       real(kind=cp), dimension (3,3)          :: M

       Q = Cart_Vector("R",H,Cell)             !Here Q is the crystallographic Q=-Q(Blume)
       SV = Cart_Vector("R",SPV,Cell)

       X = Q / sqrt(dot_product(Q,Q))          !Crystallographic Blume-frame
       QSV= Cross_Product(Q,SV)                !components w.r.t. Cartesian
       Z = QSV / sqrt(dot_product(QSV,QSV))    !default (set in Cell) frame
       Y = Cross_Product(Z,X)


       M(1,:) = 0.0    ! X-Component of Magnetic Interaction Vector is always equal to ZERO because
                       ! of MRI = Q x (M(Q) x Q); where M(Q) is the Fourier Transform of Magnetic Density of the sample
       M(2,:) = Y
       M(3,:) = Z


       MiV_PF = Matmul(M, MiVC) !conversion of MiVC to Blume frame

       return
    End Function  Magn_Inter_Vec_Pf

    !!--++
    !!--++ Function Mm(Nsf, MiV_Pf) Result(Mmc)
    !!--++    complex(kind=cp), dimension(3), intent( in) :: MiV_PF !  In  -> Magnetic Interaction Vector in polarisation frame
    !!--++    real(kind=cp)                               :: MMC    !  Out -> magnetic-magnetic interference term
    !!--++
    !!--++    (Private)
    !!--++    Calculates the magnetic-magnetic interference contribution to scattering in the
    !!--++    polarisation coordinate frame according to the Blume equations
    !!--++
    !!--++ Update April - 2005
    !!
    Function Mm(MiV_PF) Result(Mmc)
       !---- Argument ----!
       complex(kind=cp), dimension(3), intent( in) :: MiV_PF
       real(kind=cp)                               :: MMC

       !---- Local variables ----!

       MMC = Real(MiV_PF(2) * Conjg(MiV_PF(3)) + Conjg(MiV_PF(2)) * MiV_PF(3))

       return
    End Function  Mm

    !!--++
    !!--++ Function Nuc_Contr(Nsf, MiV_Pf) Result(Nsc)
    !!--++    complex(kind=cp),               intent( in):: NSF     !  In  -> Nuclear Structure Factor
    !!--++    complex(kind=cp), dimension(3), intent( in):: MiV_PF  !  In  -> Magnetic Interaction Vector in polarisation frame
    !!--++    real(kind=cp)                              :: NSC     !  Out -> nuclear scattering contribution
    !!--++
    !!--++    (Private)
    !!--++    Calculates the nuclear contribution to scattering in the polarisation coordinate frame according to the Blume equations
    !!--++
    !!--++ Update: April - 2005
    !!
    Function Nuc_Contr(Nsf) Result(Nsc)
       !---- Argument ----!
       complex(kind=cp), intent( in)   :: NSF
       real(kind=cp)                   :: NSC

       !---- Local variables ----!

       NSC = NSF * Conjg(NSF)

       return
    End Function  Nuc_Contr

    !!--++
    !!--++ Function Real_Nm_Y(Nsf, MiV_Pf) Result(R_Nm_Y)
    !!--++    complex(kind=cp),               intent(in)  :: NSF     !  In  -> Nuclear Structure Factor
    !!--++    complex(kind=cp), dimension(3), intent(in)  :: MiV_PF  !  In  -> Magnetic Interaction Vector in polarisation frame
    !!--++    real(kind=cp)                               :: R_NM_Y  !  Out -> real part of nuclear magnetic interference contribution along Y
    !!--++
    !!--++    (Private)
    !!--++    Calculates the real part of the nuclear-magnetic interference contribution along Y to scattering in the polarisation coordinate frame according to the Blume equations
    !!--++
    !!--++ Update: April - 2005
    !!
    Function Real_Nm_Y(Nsf, MiV_Pf) Result(R_Nm_Y)
       !---- Argument ----!
       complex(kind=cp),               intent(in)  :: NSF
       complex(kind=cp), dimension(3), intent(in)  :: MiV_PF
       real(kind=cp)                               :: R_Nm_Y

       !---- Local variables ----!

       R_Nm_Y = Real(NSF * Conjg(MiV_PF(2)) + Conjg(NSF) * MiV_PF(2))

       return
    End Function  Real_Nm_Y

    !!--++
    !!--++ Function Real_Nm_Z(Nsf, MiV_Pf) Result(R_Nm_Z)
    !!--++    complex(kind=cp), intent( in)               :: NSF            !  In  -> Nuclear Structure Factor
    !!--++    complex(kind=cp), dimension(3), intent( in) :: MiV_PF         !  In  -> Magnetic Interaction Vector in polarisation frame
    !!--++    real(kind=cp)                               :: R_NM_Z         !  Out -> nuclear real part of magnetic interference contribution along Z
    !!--++
    !!--++    (Private)
    !!--++    Calculates the real part of the nuclear-magnetic interference contribution along Z
    !!--++    to scattering in the polarisation coordinate frame according to the Blume equations
    !!--++
    !!--++ Update: April - 2005
    !!
    Function Real_Nm_Z(Nsf, MiV_Pf) Result(R_Nm_Z)
       !---- Argument ----!
       complex(kind=cp),               intent( in) :: NSF
       complex(kind=cp), dimension(3), intent( in) :: MiV_PF
       real(kind=cp)                               :: R_Nm_Z

       !---- Local variables ----!

       R_Nm_Z = Real(NSF * Conjg(MiV_PF(3)) + Conjg(NSF) * MiV_PF(3))

       return
    End Function  Real_Nm_Z

    !!--++
    !!--++ Function Tchiral(MiV_Pf,B_Q) Result(Tc)
    !!--++    complex(kind=cp), dimension(3), intent( in):: MiV_PF !  In  -> Magnetic Interaction Vector in polarisation frame
    !!--++    character(len=*),optional,      intent(in) :: B_Q    !  In  -> The calculation is done with original Blume equation Q_B= -Q_cryst
    !!--++    real(kind=cp)                              :: TC     !  Out -> chiral contribution
    !!--++
    !!--++    (Private)
    !!--++    Calculates the chiral contribution to scattering in the polarisation coordinate frame
    !!--++    according to the Blume equations
    !!--++
    !!--++ Update: April - 2005
    !!
    Function Tchiral(MiV_Pf,B_Q) Result(Tc)
       !---- Argument ----!
       complex(kind=cp), dimension(3), intent(in) :: MiV_PF
       character(len=*),optional,      intent(in) :: B_Q
       real(kind=cp)                              :: TC
       !---- Local variables ----!
       real(kind=cp) :: s

       s=-1.0
       if(present(B_Q)) s=1.0

       TC = -s * Aimag(MiV_PF(2) * Conjg(MiV_PF(3)) - Conjg(MiV_PF(2)) * MiV_PF(3))

       return
    End Function  Tchiral

    !---------------------!
    !---- Subroutines ----!
    !---------------------!

    !!----
    !!---- Subroutine Calc_Polar(frame,wave,Cell, UB, Pin, NSF, Mag_dom, Mh, Pf,ok,mess,B_Q)
    !!----    real(kind=cp),                intent(in) :: wave    ! In -> Cell variable
    !!----    character(len=3),             intent(in) :: frame   ! In -> Frame used for polarisation
    !!----    type (Crystal_Cell_Type),     intent(in) :: Cell    ! In -> Busing-Levy UB-matrix
    !!----    Real(kind=cp), dimension(3,3),intent(in) :: UB      ! In -> Nuclear Structure Factor
    !!----    Real(kind=cp), dimension(3),  intent(in) :: Pin     ! In -> Incident polarisation (Cartesian or spherical)
    !!----    complex(kind=cp),             intent(in) :: NSF     ! In -> Nuclear Structure Factor
    !!----    type(Magnetic_Domain_type),   intent(in) :: Mag_Dom ! In -> Magnetic domains information
    !!----    type(MagHD_Type),             intent(in) :: Mh      ! In -> Contains Magnetic structure factor, MiV, domain info, ...
    !!----    Real(kind=cp), dimension(3),  intent(out):: Pf      ! Out ->Final polarisation in frame given by "frame"
    !!----    logical,                      intent(out):: ok      ! Out -> .true. if everything gone ok
    !!----    Character (len=*),            intent(out):: mess    ! Out ->Error message
    !!----    character(len=*), optional,   intent(in) :: B_Q     ! Original Blume equations are used Q=Q_BM
    !!----
    !!----
    !!----    We use here the crystallographic convention for the scattering vector: Q=kf-ki
    !!----    So that the equations written below have the sign of the terms containing explicitly
    !!----    the imaginary unit opposite to the terms of the original equations.
    !!----    This soubroutine calculates the final polarization in two reference frames: BL and BM
    !!----    BL-> Busing-Levy and BM-> Blume-Maleyev. The incident polarisation can be provided
    !!----    in Cartesian components or giving the module and two angles of rotation (Theta,Chi)
    !!----    corresponding to nutator and precesion coils in MAD, by default they are in Cartesian
    !!----    when frame="BL" or frame="BM". The Theta rotation corresponds to a clockwise
    !!----    rotation around the "yBL-axis" and the Chi rotation to a clockwise rotation
    !!----    around the "xBL-axis). The product of the two matrices R_chi.R_theta applied to (0,0,pol)
    !!----    provides the incident polarisation in the BL system.
    !!----    For providing "angular coordinates" frame="BLS", frame="BMS" or
    !!----    frame="MAD" (BL for Pi and scattered beam for Pf)
    !!----
    !!----    Calling N the nuclear structure factor (complex scalar) and M the magnetic
    !!----    interaction vector (complex vector), the B-M equations are:
    !!----
    !!----    I = N.N* + M.M* + (N.M* + N*.M)Pi - i(M* x M)Pi
    !!----    Pf.I = N.N* Pi - M.M* Pi + (Pi.M*) M + (Pi.M) M* -i(N*.M - N.M*)xPi +
    !!----           + N.M* + N*.M + i(M* x M)
    !!----
    !!----    Defining the Chiral vector as T=i(M* x M) and the nuclear-magnetic
    !!----    interference vector as W= 2 N.M* = Wr + i Wi, calling Inuc=N.N* and
    !!----    Imag=M.M*, the equations are written in the following form:
    !!----
    !!----    I = Inuc + Imag + Wr.Pi - T.Pi
    !!----    Pf.I = (Inuc-Imag) Pi + (Pi.M*) M + (Pi.M) M* - Wi x Pi + Wr + T
    !!----
    !!----    The Busing-Levy UB-matrix is provided as an input argument.
    !!----    Alternative to the use of partial functions
    !!----
    !!---- Created: June - 2012 JRC
    !!
    Subroutine Calc_Polar(frame,wave,Cell,UB, Pin, NSF, Mag_dom, Mh, Pf,ok,mess,B_Q)
       !---- Arguments ----!
       character(len=3),             intent(in)    :: frame
       real(kind=cp),                intent(in)    :: wave
       type (Crystal_Cell_Type),     intent(in)    :: Cell
       Real(kind=cp), dimension(3,3),intent(in)    :: UB
       Real(kind=cp), dimension(3),  intent(in)    :: Pin
       complex(kind=cp),             intent(in)    :: NSF
       type(Magnetic_Domain_type),   intent(in)    :: Mag_Dom
       type(MagHD_Type),             intent(in out):: Mh
       Real(kind=cp), dimension(3),  intent(   out):: Pf
       logical,                      intent(   out):: ok
       Character(len=*),             intent(   out):: mess
       Character(len=*), optional,   intent(in)    :: B_Q
!
!       !---- Local variables ----!
       real(kind=cp)                  :: s                 ! sign for T and Wi terms
       real(kind=cp)                  :: I_inv,Inuc,Imag   ! Inverse of elastic cross section
       real(kind=cp)                  :: gamma,omega,nu    ! Normal beam angles
       real(kind=cp)                  :: pol,ptheta,pchi   ! Module of the incident pol. and angles of nutator + rotation
       complex(kind=cp), dimension (3):: MiV, W            ! MiV for one domain and in polarisation frame,Nuclear-Magnetic interference
       integer                        :: nd,ich,nch, ierr
       Real(kind=cp), dimension(3)    :: z1,z4, Pic,T,Wr,Wi
       Real(kind=cp), dimension(3,3)  :: Um,Rot,Rot_omega,ubinv,r_pth,r_chi,BL2BM
       Real(kind=cp)                  :: suma

       ok=.true.
       nch=1
       if(Mag_Dom%chir) nch=2
       suma=0.0
       do nd=1,Mag_Dom%nd
         do ich=1,nch
           suma=suma+dot_product(Mh%MiVC(:,ich,nd),Mh%MiVC(:,ich,nd))
         end do
       end do
       if(abs(NSF) < 0.0001 .and. suma < 0.0001) then
         mess="Error: the provided reflection is Nuclear and Magnetically forbidden!"
         ok=.false.
         return
       end if
       s=1.0_cp
       if(present(B_Q)) s=-1.0_cp
       z1=Matmul(UB,Mh%h) !Cartesian coordinates of reflection hkl when all motors are at zero
       ubinv=invert_A(UB)
       call Get_Angs_NB(wave,z1,gamma,omega,nu,ierr)  !Getting normal beam angles from wavelength and z1
       if(ierr /= 0) then
         mess="Error calculating the normal beam angles"
         ok=.false.
         return
       else
         if(abs(nu) > 1.0) then !The UB matrix does not correspond to horizontal
           ok=.false.           !plane scattering geometry => gamma=2theta
           mess="The nu-angle is incompatible with horizontal-plane scattering geometry"
           return
         end if
       end if
       call Phi_mat(omega,Rot_omega)  !Rotation of the omega motor to put reflection in diffraction position
                                      !This is counter-clockwise (c s 0|-s c 0|0 0 1)
       z4=matmul(Rot_omega,z1)        !Crystallographic scattering vector in Lab-system when
                                      !the reflection is in diffraction position.
       Rot=Matmul(UB,Matmul(Cell%GD,Cell%Orth_Cr_cel)) !Conversion to BL (in diffraction position) from Crystal Cartesian

       Rot=Matmul(Rot_omega,Rot) !Rotation matrix putting the MiVC in the BL frame

       if(frame(1:2) == "BM" .or. frame(1:3) == "MAD")  then  !Crystallographic Blume-Maleyev frame
         z4 = z4 / sqrt(dot_product(z4,z4)) !Unit vector along Q in BL system
         Um(1,:) = 0.0
         Um(2,:) = Cross_Product((/0.0_cp,0.0_cp,1.0_cp/),z4)
         Um(3,:) = (/0.0_cp,0.0_cp,1.0_cp/)
         Rot=Matmul(Um,Rot)
       else if(frame(1:2) /= "BL" .and. frame(1:3) /= "MAD") then
         ok=.false.
         mess="Undefined reference frame"
         return
       end if
       !Put now the incident polarisation in Cartesian coordinates if angular
       !data have been provided.
       if(frame == "BLS" .or. frame == "BMS" .or. frame == "MAD") then
         pol=Pin(1); ptheta=Pin(2); pchi=Pin(3)
         call chi_mat(ptheta,R_pth) !This is clockwise (c 0 s|0 1 0|-s 0 c)
         call psi_mat(pchi,R_chi)   !This is counter-clockwise (1 0 0|0 c -s|0 s c)
         R_Chi=transpose(R_Chi)     !Put chi clockwise
         Pic=matmul(R_chi,matmul(R_pth,(/0.0_cp,0.0_cp,pol/))) !BL-system incident polarisation
         if(frame == "BMS" .or. frame == "MAD") then     !Transform to BM system
           call Phi_mat(gamma*0.5+90.0,BL2BM) !Active
           Pic=matmul(transpose(BL2BM),Pic)      !Incident polarisation in BM system
         end if
       else
         Pic=Pin
       end if
       ! Loop over domains
       Pf=0.0
       Inuc=real(Conjg(NSF)*NSF)
       do nd=1,Mag_Dom%nd
         do ich=1,nch
           MiV=Matmul(Rot,Mh%MiVC(:,ich,nd))  !Convert the MiVC to the frame BM or BL
           Imag=dot_Product(MiV,MiV)          ! dot_product => MiV*.MiV
           T=-s*aimag(Cross_Product(Conjg(MiV),MiV)) !Chiral Vector
           W=2.0*NSF*Conjg(MiV)       !Nuclear-Magnetic Interaction vector
           Wr=real(W); Wi=s*aimag(W)  !Real and Imaginary parts
           I_inv=1.0/(Inuc+Imag+dot_product(Wr,Pic)-dot_product(T,Pic))
           Pf=Pf+ I_inv * Mag_Dom%Pop(ich,nd)*( (Inuc-Imag)*Pic +                  &
              dot_product(Pic,MiV)*Conjg(MiV) + dot_product(Pic,Conjg(MiV))* MiV + &
              T + Wr - Cross_Product(Wi,Pic) )
         end do !loop over chiral domains
       end do !loop over S-domains

       !At the end of the calculation Pf is in Cartesian components with respect to
       !BM or BL system depending if the user wants to get the angles of the scattered
       !beam nutator and precesion we have to pass first to the Cartesian coordinates
       !in the scattered beam system. Otherwise the calculated Pf is returned as is.
       if(frame == "BLS") then   !put polarisation in scattered beam system
         call Phi_mat(gamma,Rot)
         Pf = matmul(transpose(Rot) ,Pf)
       else if(frame == "BMS" .or. frame == "MAD" ) then
         call Phi_mat(gamma,Rot)
         Pf = matmul(transpose(Rot),matmul(BL2BM, Pf) )
       else
         return
       end if
       !Calculate now the angles pchi and ptheta to bring back the polarisation
       !to the z-direction
       pol=sqrt(dot_product(Pf,Pf)) !magnitude of scattered polarisation
       if(abs(pf(3)) < 1.0e-5) then
         pchi=90.0
       else
         pchi=atan2d(pf(2),pf(3))
       end if
       ptheta=atan2d(pf(1),pf(2)*sind(pchi)+pf(3)*cosd(pchi))
       pf=(/pol,ptheta,pchi/)

       return
    End Subroutine Calc_Polar

    !!---- Subroutine Get_Pol_Tensor_Pc(frame,wave,Cell,UB,Pin, NSF, Mag_dom, Mh, Pol_tens, Pc,B_Q)
    !!----    character(len=3),             intent(in)    :: frame   ! In -> Frame used for polarisation
    !!----    real(kind=cp),                intent(in)    :: wave
    !!----    type (Crystal_Cell_Type),     intent(in)    :: Cell
    !!----    Real(kind=cp), dimension(3,3),intent(in)    :: UB
    !!----    Real(kind=cp), dimension(3),  intent(in)    :: Pin
    !!----    complex(kind=cp),             intent(in)    :: NSF
    !!----    type(Magnetic_Domain_type),   intent(in)    :: Mag_Dom
    !!----    type(MagHD_Type),             intent(in)    :: Mh
    !!----    real(kind=cp), dimension(3,3),intent(out)   :: Pol_tens
    !!----    real(kind=cp), dimension(3),  intent(out)   :: Pc
    !!----    character(len=*), optional,   intent(in)    :: B_Q  !Original Blume equations are used Q=Q_BM
    !!----
    !!----    This subroutine provides the tensor [P] and vector Pc
    !!----    in the crystallographic Blume-Maleyev frame. The input polarisation
    !!----    may be provided either in BM Cartesian coordinates or using the module
    !!----    and two angles of rotation (Theta,Chi)corresponding to nutator
    !!----    and precesion coils in MAD
    !!----    The tensorial form (for whatever reference frame) can be written as:
    !!----    Pf = [P] Pi + Pc => Pc= (Wr + T)/I
    !!----      [P] =  ([D] + [S] + [A])/I
    !!----      [D] =  (Inuc- Imag*) E
    !!----      [S] =  [M o M* + M* o M]  -> o indicates tensorial product
    !!----      [A] =  [+i(N.M* - N*.M)]cross = [-Wi]cross
    !!----
    Subroutine Get_Pol_Tensor_Pc(frame,wave,Cell,UB,Pin, NSF, Mag_dom, Mh, Pol_tens, Pc,ok,mess,B_Q)
       character(len=3),             intent(in)    :: frame
       real(kind=cp),                intent(in)    :: wave
       type (Crystal_Cell_Type),     intent(in)    :: Cell
       Real(kind=cp), dimension(3,3),intent(in)    :: UB
       Real(kind=cp), dimension(3),  intent(in)    :: Pin
       complex(kind=cp),             intent(in)    :: NSF
       type(Magnetic_Domain_type),   intent(in)    :: Mag_Dom
       type(MagHD_Type),             intent(in)    :: Mh
       real(kind=cp), dimension(3,3),intent(out)   :: Pol_tens
       real(kind=cp), dimension(3),  intent(out)   :: Pc
       logical,                      intent(   out):: ok
       Character(len=*),             intent(   out):: mess
       Character(len=*), optional,   intent(in)    :: B_Q  !Original Blume equations are used Q=Q_BM
       !--- Local variables
       real(kind=cp)                  :: s                 ! sign for T and Wi terms
       real(kind=cp)                  :: I_inv,Inuc,Imag   ! inverse of elastic cross section
       real(kind=cp)                  :: gamma,omega,nu    ! Normal beam angles
       real(kind=cp)                  :: pol,ptheta,pchi   ! Incident polarisation
       complex(kind=cp), dimension (3):: MiV, W     !MiV for one domain and in polarisation frame
       integer                        :: nd,ich,nch, ierr
       Real(kind=cp), dimension(3)    :: z1,z4, Pic,T,Wr,Wi
       Real(kind=cp), dimension(3,3)  :: Um,Rot,Rot_omega,DD,AA,SS,r_pth,r_chi,BL2BM
       Real(kind=cp), dimension(3,3), parameter  :: Identity= reshape ( (/1.0,0.0,0.0, &
                                                                          0.0,1.0,0.0, &
                                                                          0.0,0.0,1.0/),(/3,3/) )
       ok=.true.
       s=1.0_cp
       if(present(B_Q)) s=-1.0_cp
       z1=Matmul(UB,Mh%h)
       call Get_Angs_NB(wave,z1,gamma,omega,nu,ierr)
       if(ierr /= 0) then
         mess="Error calculating the normal beam angles"
         ok=.false.
         return
       else
         if(abs(nu) > 1.0) then !The UB matrix does not correspond to horizontal
           ok=.false.           !plane scattering geometry
           mess="The nu-angle is incompatible with horizontal-plane scattering geometry"
         end if
       end if
       call Phi_mat(omega,Rot_omega)  !Rotation of the omega motor to put reflection in diffraction position
       z4=matmul(Rot_omega,z1) !Crystallographic scattering vector in Lab-system when
                               !the reflection is in diffraction position.
       Rot=Matmul(UB,Matmul(Cell%GD,Cell%Orth_Cr_cel)) !Conversion to BL (in diffraction position) from Crystal Cartesian
       Rot=Matmul(Rot_omega,Rot) !Rotation matrix putting the MiVC in the BL frame
       z4 = z4 / sqrt(dot_product(z4,z4))
       Um(1,:) = 0.0_cp
       Um(2,:) = Cross_Product((/0.0_cp,0.0_cp,1.0_cp/),z4)
       Um(3,:) = (/0.0_cp,0.0_cp,1.0_cp/)
       Rot=Matmul(Um,Rot)
       !Put now the incident polarisation in Cartesian coordinates w.r.t. BM frame
       !if angular data have been provided.
       if(frame == "BLS" .or. frame == "BMS" .or. frame == "MAD") then
         pol=Pin(1); ptheta=Pin(2); pchi=Pin(3)
         call chi_mat(ptheta,R_pth) !This is clockwise (c 0 s|0 1 0|-s 0 c)
         call psi_mat(pchi,R_chi)   !This is counter-clockwise (1 0 0|0 c -s|0 s c)
         R_Chi=transpose(R_Chi)     !Put chi clockwise
         Pic=matmul(R_chi,matmul(R_pth,(/0.0_cp,0.0_cp,pol/))) !BL-system incident polarisation
         call Phi_mat(gamma*0.5_cp+90.0_cp,BL2BM) !Active
         Pic=matmul(transpose(BL2BM),Pic)   !Incident polarisation in BM system
       else
         Pic=Pin
       end if
       nch=1
       if(Mag_Dom%chir) nch=2
       ! Loop over domains
       Pol_tens = 0.0_cp
       Pc       = 0.0_cp
       Inuc=real(Conjg(NSF)*NSF)
       do nd=1,Mag_Dom%nd
         do ich=1,nch
          MiV=Matmul(Rot,Mh%MiVC(:,ich,nd)) !Convert the MiVC to the frame BM
          Imag=dot_Product(MiV,MiV)
          T=-s*aimag(Cross_Product(Conjg(MiV),MiV)) !Chiral Vector
          W=2.0_cp*NSF*Conjg(MiV)     !Nuclear-Magnetic Interaction vector
          Wr=real(W); Wi=s*aimag(W)  !Real and Imaginary parts
          I_inv=1.0_cp/(Inuc+Imag+dot_product(Wr,Pic)-dot_product(T,Pic))
          DD=(Inuc-Imag)*Identity
          SS=real(Tensor_Product(MiV,Conjg(MiV))+Tensor_Product(Conjg(MiV),MiV))
          AA=-Mat_Cross(Wi)
          Pol_tens=Pol_tens+ I_inv*Mag_Dom%Pop(ich,nd)*( DD + SS + AA )
          Pc=Pc+ I_inv*Mag_Dom%Pop(ich,nd)*(T+Wr)
         end do !loop over chiral domains
       end do !loop over S-domains
       return
    End Subroutine Get_Pol_Tensor_Pc

    !!----
    !!---- Subroutine Calc_Polar_Dom(Cell, H, SPV, Pin, NSF, Mag_dom, Mh, Polari,ok,mess,B_Q)
    !!----    Type (Crystal_Cell_Type),    intent(in)    :: Cell  !  In -> Cell variable
    !!----    real(kind=cp), dimension (3),intent(in)    :: H     !  In -> Scattering vector in hkl
    !!----    real(kind=cp), dimension(3), intent(in)    :: SPV   !  In -> Second Scattering plane vector in hkl
    !!----    real(kind=cp), intent( in)                 :: Pin   !  In -> magnitude of initial polarisation
    !!----    complex(kind=cp), intent( in)              :: NSF   !  In -> Nuclear Scattering Factor
    !!----    Type(Magnetic_Domain_type),  intent(in)    :: Mag_Dom
    !!----    Type(MagHD_Type),            intent(in out):: Mh
    !!----    Type (Polar_calc_type),      intent( out)  :: Polari !  Out ->type with all information about polarisation in
    !!----                                                                 one point hkl
    !!----    Logical,                     intent(out)   :: ok
    !!----    character(len=*),            intent(out)   :: mess
    !!----    character(len=*), optional,  intent(in)    :: B_Q  !Original Blume equations are used Q=Q_BM
    !!----
    !!----    Calculates Polarization matrix for domain case
    !!----
    !!---- Created: March - 2009 OZ, Updated: June-2012 (JRC)
    !!
    Subroutine Calc_Polar_Dom(Cell, H, SPV, Pin, NSF, Mag_dom, Mh, Polari,ok,mess,B_Q)
       !---- Arguments ----!
       type (Crystal_Cell_Type),    intent(in)       :: Cell
       real(kind=cp), dimension (3),intent(in)       :: H
       real(kind=cp), dimension(3), intent(in)       :: SPV
       real(kind=cp),               intent(in)       :: Pin
       complex(kind=cp),            intent(in)       :: NSF
       type(Magnetic_Domain_type),  intent(in)       :: Mag_Dom
       type(MagHD_Type),            intent(in out)   :: Mh
       type(Polar_calc_type),       intent(out)      :: Polari
       logical,                     intent(out)      :: ok
       character(len=*),            intent(out)      :: mess
       character(len=*), optional,  intent(in)       :: B_Q

       !---- Local variables ----!
       real(kind=cp), dimension (3)   :: sigma        ! elastic cross for different incident polarisation directions
       real(kind=cp)                  :: nc,my,mz,rnmy,rnmz,inmy,inmz,tc,mmc,a,suma !the different contribution to cross-section
       complex(kind=cp), dimension (3):: MiV, MiV_PF       !MiV for one domain and in polarisation frame
       integer                        :: nd,ich,nch


       A = tpi**3/Cell%CellVol
       ok=.true.
       !First store given info in Polari
       Polari%H = H
       Polari%SPV = SPV
       Polari%Cell = Cell
       Polari%P = Pin
       Polari%NSF = NSF

       Polari%Pij(:,:) = 0.0_cp

       nch=1
       if(Mag_Dom%chir) nch=2
       ! Loop over domains
       suma=0.0
       do nd=1,Mag_Dom%nd
         do ich=1,nch
           suma=suma+dot_product(Mh%MiVC(:,ich,nd),Mh%MiVC(:,ich,nd))      !Use the Cartesian components before calling Magn_Inter_Vec_PF
         end do
       end do
       if(abs(NSF) < 0.0001 .and. suma < 0.0001) then
         mess="Error: the provided reflection is Nuclear and Magnetically forbidden!"
         ok=.false.
         return
       end if
       !Calculate the rest and also store it in Polari
       ! Loop over domains
       do nd=1,Mag_Dom%nd
         do ich=1,nch
           MiV=Mh%MiVC(:,ich,nd)      !Use the Cartesian components before calling Magn_Inter_Vec_PF
           !magnetic interaction in polarisation frame
           MiV_PF = Magn_Inter_Vec_PF(MiV,H,SPV, Cell)
           Polari%MiV(:,ich,nd) = MiV_PF

           !the different contributions to the scattering cross-section
           nc = nuc_contr(NSF)
           Polari%NC = A * nc

           my = mag_y(MiV_PF)
           Polari%MY(ich,nd) = A * my

           mz = mag_z(MiV_PF)
           Polari%MZ(ich,nd) = A * mz

           rnmy = real_nm_y(NSF, MiV_PF)
           Polari%RY(ich,nd) = A * rnmy

           rnmz = real_nm_z(NSF, MiV_PF)
           Polari%RZ(ich,nd) = A * rnmz


           if(present(B_Q)) then
             inmy = im_nm_y(NSF, MiV_PF,"B_Q")
             inmz = im_nm_z(NSF, MiV_PF,"B_Q")
             tc   = tchiral(MiV_PF,"B_Q")
           else
             inmy = im_nm_y(NSF, MiV_PF)
             inmz = im_nm_z(NSF, MiV_PF)
             tc   = tchiral(MiV_PF)
           end if

           Polari%IY(ich,nd) = A * inmy
           Polari%IZ(ich,nd) = A * inmz
           Polari%TC(ich,nd) = tc

           mmc = mm(MiV_PF)
           Polari%MM(ich,nd) = A * mmc

           !scattering cross-section for the different initial polarisation vectors
           sigma = (/ nc + my + mz - Pin * tc, nc + my + mz + Pin * rnmy, nc + my + mz + Pin * rnmz /)
           Polari%CS(:,ich,nd) = sigma

           !summing of the polarisation matrix
           Polari%Pij(1,1) = Polari%Pij(1,1) + Mag_Dom%Pop(ich,nd)*((nc - my - mz)* Pin + tc)/sigma(1)
           Polari%Pij(1,2) = Polari%Pij(1,2) + Mag_Dom%Pop(ich,nd)*(inmz * Pin + tc)/sigma(2)
           Polari%Pij(1,3) = Polari%Pij(1,3) + Mag_Dom%Pop(ich,nd)*(-inmy * Pin + tc)/sigma(3)
           Polari%Pij(2,1) = Polari%Pij(2,1) + Mag_Dom%Pop(ich,nd)*(-inmz * Pin + rnmy)/sigma(1)
           Polari%Pij(2,2) = Polari%Pij(2,2) + Mag_Dom%Pop(ich,nd)*((nc + my - mz) * Pin + rnmy)/sigma(2)
           Polari%Pij(2,3) = Polari%Pij(2,3) + Mag_Dom%Pop(ich,nd)*(mmc * Pin + rnmy)/sigma(3)
           Polari%Pij(3,1) = Polari%Pij(3,1) + Mag_Dom%Pop(ich,nd)*(inmy * Pin + rnmz)/sigma(1)
           Polari%Pij(3,2) = Polari%Pij(3,2) + Mag_Dom%Pop(ich,nd)*(mmc * Pin + rnmz)/sigma(2)
           Polari%Pij(3,3) = Polari%Pij(3,3) + Mag_Dom%Pop(ich,nd)*((nc - my + mz) * Pin + rnmz)/sigma(3)

         end do !loop over chiral domains
       end do !loop over S-domains

       return
    End Subroutine Calc_Polar_Dom

    !!----
    !!---- Subroutine Set_Polar_Info(Cell, H, Spv, Pin, Nsf, MiV, Polari,B_Q)
    !!----    Type (Crystal_Cell_Type),       intent(in)  :: Cell     !  In -> Cell variable
    !!----    real(kind=cp), DIMENSION (3),   intent(in)  :: H        !  In -> Scattering vector in hkl
    !!----    real(kind=cp), dimension(3),    intent(in)  :: SPV      !  In -> Second Scattering plane vector in hkl
    !!----    real(kind=cp),                  intent(in)  :: Pin      !  In -> magnitude of initial polarisation
    !!----    complex(kind=cp),               intent(in)  :: NSF      !  In -> Nuclear Scattering Factor
    !!----    complex(kind=cp), dimension(3), intent(in)  :: MiV      !  In -> Magnetic interaction vector
    !!----    Type (Polar_Info_type),         intent(out) :: Polari   !  Out ->type with all information about polarisation in
    !!----                                                                     one point hkl
    !!----    character(len=*), optional,     intent(in)  :: B_Q  !Original Blume equations are used Q=Q_BM
    !!----
    !!----    Initializes the polarisation info type
    !!----
    !!---- Updated: April - 2008, June-2012 (JRC)
    !!
    Subroutine Set_Polar_Info(Cell, H, Spv, Pin, Nsf, MiV, Polari,B_Q)
       !---- Arguments ----!
       Type (Crystal_Cell_Type),       intent(in)  :: Cell
       real(kind=cp), DIMENSION (3),   intent(in)  :: H
       real(kind=cp), dimension(3),    intent(in)  :: SPV
       real(kind=cp),                  intent(in)  :: Pin
       complex(kind=cp),               intent(in)  :: NSF
       complex(kind=cp), dimension(3), intent(in)  :: MiV
       Type (Polar_Info_type),         intent(out) :: Polari
       character(len=*), optional,     intent(in)  :: B_Q

       !---- Local variables ----!
       real(kind=cp), DIMENSION (3)     :: sigma        ! elastic cross for different inicdent polarisation directions
       real(kind=cp)                    :: nc, my, mz, rnmy, rnmz, inmy, inmz, tc, mmc, A !the different contribution to cross-section
       complex(kind=cp), DIMENSION (3)  :: MiV_PF       !MiV in polarisation frame


       A = tpi**3/Cell%CellVol

       !First store given info in Polari
       Polari%H = H
       Polari%SPV = SPV
       Polari%Cell = Cell
       Polari%P = Pin
       Polari%NSF = NSF
       if(abs(dot_product(MiV,MiV)) < 0.0001 .and. abs(NSF) < 0.0001) then
         Polari%Pij=0.0
         return
       end if
       !Calculate the rest and also store it in Polari

       !magnetic interaction in polarisation frame
       MiV_PF = Magn_Inter_Vec_PF(MiV,H,SPV, Cell) !It is assumed that MiV is provided in Cartesian components
       Polari%MiV = MiV_PF
       !the different contributions to the scattering cross-section
       nc = nuc_contr(NSF)
       Polari%NC = A * nc

       my = mag_y(MiV_PF)
       Polari%MY = A * my

       mz = mag_z(MiV_PF)
       Polari%MZ = A * mz

       rnmy = real_nm_y(NSF, MiV_PF)
       Polari%RY = A * rnmy

       rnmz = real_nm_z(NSF, MiV_PF)
       Polari%RZ = A * rnmz

       if(present(B_Q)) then
         inmy = im_nm_y(NSF, MiV_PF,"B_Q")
         inmz = im_nm_z(NSF, MiV_PF,"B_Q")
         tc = tchiral(MiV_PF,"B_Q")
       else
         inmy = im_nm_y(NSF, MiV_PF)
         inmz = im_nm_z(NSF, MiV_PF)
         tc = tchiral(MiV_PF)
       end if

       Polari%IY = A * inmy
       Polari%IZ = A * inmz
       Polari%TC = tc

       mmc = mm(MiV_PF)
       Polari%MM = A * mmc

       !scattering cross-section for the different initial polarisation vectors
       sigma = (/ nc + my + mz - Pin * tc, nc + my + mz + Pin * rnmy, nc + my + mz + Pin * rnmz /)
       Polari%CS = sigma
       write(*,*) sigma
       if(dot_product(sigma,sigma) > 0.001) then
         !the polar matrix
         Polari%Pij(1,1) = ((nc - my - mz)* Pin + tc)/sigma(1)
         Polari%Pij(1,2) = (inmz * Pin + tc)/sigma(2)
         Polari%Pij(1,3) = (-inmy * Pin + tc)/sigma(3)
         Polari%Pij(2,1) = (-inmz * Pin + rnmy)/sigma(1)
         Polari%Pij(2,2) = ((nc + my - mz) * Pin + rnmy)/sigma(2)
         Polari%Pij(2,3) = (mmc * Pin + rnmy)/sigma(3)
         Polari%Pij(3,1) = (inmy * Pin + rnmz)/sigma(1)
         Polari%Pij(3,2) = (mmc * Pin + rnmz)/sigma(2)
         Polari%Pij(3,3) = ((nc - my + mz) * Pin + rnmz)/sigma(3)
       else
         Polari%Pij=0.0
       end if
       return
    End Subroutine Set_Polar_Info

    !!----
    !!----  Subroutine Write_Polar_Info(Polari, Mag_Dom, Lun, info)
    !!----    !---- Arguments ----!
    !!----    type(Polar_calc_type),     intent(in)  :: Polari     !  in -> Type with all information about polarisation in one point hkl
    !!----    type(Magnetic_Domain_type),intent(in)  :: Mag_Dom    !  in -> Magnetic domains
    !!----    integer,         optional, intent(in)  :: Lun        !  In -> Unit to write
    !!----    character(len=*),optional, intent(in)  :: info       !  in -> if info "P" also print information about coordinate frame                                                          !        if info "C" also print information about crystal
    !!----                                                         !        if info "B" also print information about both
    !!----
    !!----    Outputs the polarisation info type in nice form
    !!----
    !!---- Created: April - 2005
    !!---- Modified March - 2009 OZ for multidomain case
    !!

    Subroutine Write_Polar_Info(Polari, Mag_Dom, Lun, info)
       !---- Arguments ----!
       type(Polar_calc_type),     intent(in)  :: Polari
       type(Magnetic_Domain_type),intent(in)  :: Mag_Dom
       integer,         optional, intent(in)  :: Lun
       character(len=*),optional, intent(in)  :: info
       !---- Local variables ----!
       integer            :: iunit,nd,nch,ich,lfmt_out
       character(32)   :: tit
       character(132)  :: fmt_out

       iunit=6
       if (present(lun)) iunit=lun

       Write(unit=iunit,fmt="(/,a)")            "        Polar information:"
       Write(unit=iunit,fmt="(a,/)")            "        -------------------"
       Write(unit=iunit,fmt="(a,/)")            " => Initial parameters:"
       Write(unit=iunit,fmt="(3(a,f12.3), a)")  "    Scattering vector  (Qh,Qh,Ql) =  (", Polari%h(1) ,", ", Polari%h(2) , &
                                                ", ", Polari%h(3), ")"
       Write(unit=iunit,fmt="(3(a,f12.3), a)")  "    Add. in-plane vector      SPV =  (", Polari%SPV(1) ,", ", Polari%SPV(2) ,&
                                                ", ", Polari%SPV(3), ")"
       Write(unit=iunit,fmt="(a,f12.3)")        "              Polarisation degree =   ", Polari%P
       Write(unit=iunit,fmt="(a,/)")            " "
       IF (present(info)) THEN
         IF (info == "C" .OR. info == "c" .OR. info == "B" .OR. info == "b") THEN
           Write(unit=iunit,fmt="(a,/)")        " => Crystal information:"
           Write(unit=iunit,fmt="(3(a,f12.4))") "      a = ", Polari%Cell%cell(1),"      b = ", Polari%Cell%cell(2), "      c = ",&
                                                Polari%Cell%cell(3)
           Write(unit=iunit,fmt="(3(a,f12.3))") "  alpha = ", Polari%Cell%ang(1) ,"   beta = ", Polari%Cell%ang(2) , "  gamma = ",&
                                                Polari%Cell%ang(3)
           Write(unit=iunit,fmt="(a,f12.4)")    "                     Direct Cell Volume = ",Polari%Cell%CellVol
           Write(unit=iunit,fmt="(a,/)")        ""
        End if
         IF (info == "P" .OR. info == "p" .OR. info == "B" .OR. info == "b") THEN
           Write(unit=iunit,fmt="(a,/)")  " => Polarisation coordinate frame:"
           Write(unit=iunit,fmt="(a,/)")  " "
           Write(unit=iunit,fmt="(a,/)")  "    Polarisation coordinate frame according to Blume"
           Write(unit=iunit,fmt="(a,a,/)")"    X  || scattering vector Q    ",&
                                          " (where Q is the scattering Vector in cartesian real space coordinates)"
           Write(unit=iunit,fmt="(a,/)")  "    Y _|_ scattering vector Q in scattering plane"
           Write(unit=iunit,fmt="(a,/)")  "    Z _|_ scattering vector Q out of scattering plane"
           Write(unit=iunit,fmt="(a,/)")  "    (ATTENTION: This choice is not non-ambiguous, there are always two possible choices"
           Write(unit=iunit,fmt="(a,/)")  "    for a right handed coordinate frame which will fullfils this condition!!!)"
           Write(unit=iunit,fmt="(a,/)")  " "
           Write(unit=iunit,fmt="(a,/)")  "    Therefore the right handed coordinate frame is explicitly chosen like this:"
           Write(unit=iunit,fmt="(a,a,/)")"    X := Q/|Q|                 where Q is the scattering Vector in ", &
                                          "cartesian real space coordinates"
           Write(unit=iunit,fmt="(a,a,/)")"    Z := (Q x SVP)/|(Q x SVP)| where SVP is a second vector in the scattering ",&
                                          "plane in cartesian real space coordinates"
           Write(unit=iunit,fmt="(a,/)")  "    Y := (Z x X)"
           Write(unit=iunit,fmt="(a,/)")  " "
           Write(unit=iunit,fmt="(a)")    "                            Y"
           Write(unit=iunit,fmt="(a)")    "                          /|\"
           Write(unit=iunit,fmt="(a)")    "                           |"
           Write(unit=iunit,fmt="(a)")    "                   Q       | Z "
           Write(unit=iunit,fmt="(a)")    "            ____________ _\o_____\ X"
           Write(unit=iunit,fmt="(a)")    "            \             /      /"
           Write(unit=iunit,fmt="(a)")    "             \           / "
           Write(unit=iunit,fmt="(a)")    "              \         /"
           Write(unit=iunit,fmt="(a)")    "               \       /"
           Write(unit=iunit,fmt="(a)")    "                \     /  K_f"
           Write(unit=iunit,fmt="(a)")    "             K_i \   /"
           Write(unit=iunit,fmt="(a)")    "                 _\|/_"
           Write(unit=iunit,fmt="(a,/)")  ""
         END IF
       END IF
       Write(unit=iunit,fmt="(a,/)")         " => Interaction potentials:"
       Write(unit=iunit,fmt="(2(a,f7.3))")   "       NSF = ", real(Polari%NSF)," + i ", AIMAG(Polari%NSF)

       nch=1
       if(Mag_Dom%chir) nch=2
       do nd=1,Mag_Dom%nd
         do ich=1,nch
           Write(unit=iunit,fmt="(2(a,i2),2(a,3f7.3),a)") "   Domain # =",nd," Chiral Dom. =",ich," MiV = (",&
                         real(Polari%MiV(:,ich,nd)), ") + i(",AIMAG(Polari%MiV(:,ich,nd)),")"
         end do
       end do
       Write(unit=iunit,fmt="(a,/)")         " "
       Write(unit=iunit,fmt="(a,/)")         " => Different contributions to the cross-section:"
       Write(unit=iunit,fmt="(a,f12.3)")     "         Nuclear Contribution = ", Polari%NC
       Write(unit=iunit,fmt="(a,/)")         "             Different S-domains (Pop=100%), after : chiral counterparts)"

       nch=1
       if(Mag_Dom%chir) nch=2

! format for variable number of domains

       fmt_out(1:3)='(a,'
       write(fmt_out(4:5),'(i2.2)') Mag_Dom%nd

       tit="             Magnetic along y = "
       if(nch==1) then
         fmt_out(6:17)='(f12.3))'
         lfmt_out=17
         write(unit=iunit,fmt=fmt_out(1:lfmt_out))   tit, Polari%MY(1,1:Mag_Dom%nd)
       else if(nch==2) then
         fmt_out(6:16)='(f12.3)'
         fmt_out(17:19)=',a,'
         write(fmt_out(20:21),'(i2.2)') Mag_Dom%nd
         fmt_out(22:29)='(f12.3))'
         lfmt_out=29
         write(unit=iunit,fmt=fmt_out(1:lfmt_out))   tit, Polari%MY(1,1:Mag_Dom%nd),' :',Polari%MY(2,1:Mag_Dom%nd)
       end if

       tit="             Magnetic along z = "
       if(nch==1) then
         fmt_out(6:17)='(f12.3))'
         lfmt_out=17
         write(unit=iunit,fmt=fmt_out(1:lfmt_out))   tit, Polari%MZ(1,1:Mag_Dom%nd)
       else if(nch==2) then
         fmt_out(6:16)='(f12.3)'
         fmt_out(17:19)=',a,'
         write(fmt_out(20:21),'(i2.2)') Mag_Dom%nd
         fmt_out(22:29)='(f12.3))'
         lfmt_out=29
         write(unit=iunit,fmt=fmt_out(1:lfmt_out))   tit, Polari%MZ(1,1:Mag_Dom%nd),' :',Polari%MZ(2,1:Mag_Dom%nd)
       end if

       tit="  Real nuclear magnetic al. y = "
       if(nch==1) then
         fmt_out(6:17)='(f12.3))'
         lfmt_out=17
         write(unit=iunit,fmt=fmt_out(1:lfmt_out))   tit, Polari%RY(1,1:Mag_Dom%nd)
       else if(nch==2) then
         fmt_out(6:16)='(f12.3)'
         fmt_out(17:19)=',a,'
         write(fmt_out(20:21),'(i2.2)') Mag_Dom%nd
         fmt_out(22:29)='(f12.3))'
         lfmt_out=29
         write(unit=iunit,fmt=fmt_out(1:lfmt_out))   tit, Polari%RY(1,1:Mag_Dom%nd),' :',Polari%RY(2,1:Mag_Dom%nd)
       end if

       tit="  Real nuclear magnetic al. z = "
       if(nch==1) then
         fmt_out(6:17)='(f12.3))'
         lfmt_out=17
         write(unit=iunit,fmt=fmt_out(1:lfmt_out))   tit, Polari%RZ(1,1:Mag_Dom%nd)
       else if(nch==2) then
         fmt_out(6:16)='(f12.3)'
         fmt_out(17:19)=',a,'
         write(fmt_out(20:21),'(i2.2)') Mag_Dom%nd
         fmt_out(22:29)='(f12.3))'
         lfmt_out=29
         write(unit=iunit,fmt=fmt_out(1:lfmt_out))   tit, Polari%RZ(1,1:Mag_Dom%nd),' :',Polari%RZ(2,1:Mag_Dom%nd)
       end if

       tit="  Imag nuclear magnetic al. y = "
       if(nch==1) then
         fmt_out(6:17)='(f12.3))'
         lfmt_out=17
         write(unit=iunit,fmt=fmt_out(1:lfmt_out))   tit, Polari%IY(1,1:Mag_Dom%nd)
       else if(nch==2) then
         fmt_out(6:16)='(f12.3)'
         fmt_out(17:19)=',a,'
         write(fmt_out(20:21),'(i2.2)') Mag_Dom%nd
         fmt_out(22:29)='(f12.3))'
         lfmt_out=29
         write(unit=iunit,fmt=fmt_out(1:lfmt_out))   tit, Polari%IY(1,1:Mag_Dom%nd),' :',Polari%IY(2,1:Mag_Dom%nd)
       end if

       tit="  Imag nuclear magnetic al. z = "
       if(nch==1) then
         fmt_out(6:17)='(f12.3))'
         lfmt_out=17
         write(unit=iunit,fmt=fmt_out(1:lfmt_out))   tit, Polari%IZ(1,1:Mag_Dom%nd)
       else if(nch==2) then
         fmt_out(6:16)='(f12.3)'
         fmt_out(17:19)=',a,'
         write(fmt_out(20:21),'(i2.2)') Mag_Dom%nd
         fmt_out(22:29)='(f12.3))'
         lfmt_out=29
         write(unit=iunit,fmt=fmt_out(1:lfmt_out))   tit, Polari%IZ(1,1:Mag_Dom%nd),' :',Polari%IZ(2,1:Mag_Dom%nd)
       end if

       tit="          Chiral contribution = "
       if(nch==1) then
         fmt_out(6:17)='(f12.3))'
         lfmt_out=17
         write(unit=iunit,fmt=fmt_out(1:lfmt_out))   tit, Polari%TC(1,1:Mag_Dom%nd)
       else if(nch==2) then
         fmt_out(6:16)='(f12.3)'
         fmt_out(17:19)=',a,'
         write(fmt_out(20:21),'(i2.2)') Mag_Dom%nd
         fmt_out(22:29)='(f12.3))'
         lfmt_out=29
         write(unit=iunit,fmt=fmt_out(1:lfmt_out))   tit, Polari%TC(1,1:Mag_Dom%nd),' :',Polari%TC(2,1:Mag_Dom%nd)
       end if

       tit="            Magnetic Magnetic = "
       if(nch==1) then
         fmt_out(6:17)='(f12.3))'
         lfmt_out=17
         write(unit=iunit,fmt=fmt_out(1:lfmt_out))   tit, Polari%MM(1,1:Mag_Dom%nd)
       else if(nch==2) then
         fmt_out(6:16)='(f12.3)'
         fmt_out(17:19)=',a,'
         write(fmt_out(20:21),'(i2.2)') Mag_Dom%nd
         fmt_out(22:29)='(f12.3))'
         lfmt_out=29
         write(unit=iunit,fmt=fmt_out(1:lfmt_out))   tit, Polari%MM(1,1:Mag_Dom%nd),' :',Polari%MM(2,1:Mag_Dom%nd)
       end if

       Write(unit=iunit,fmt="(a,/)")   " "
       Write(unit=iunit,fmt="(a,/)")   " => Cross-section for initial polar vector:"
       Write(unit=iunit,fmt="(a,/)")   "             Different S-domains (Pop=100%), after : chiral counterparts)"


       tit="                      along x = "
       if(nch==1) then
         fmt_out(6:17)='(f12.3))'
         lfmt_out=17
         write(unit=iunit,fmt=fmt_out(1:lfmt_out))   tit, Polari%CS(1,1,1:Mag_Dom%nd)
       else if(nch==2) then
         fmt_out(6:16)='(f12.3)'
         fmt_out(17:19)=',a,'
         write(fmt_out(20:21),'(i2.2)') Mag_Dom%nd
         fmt_out(22:29)='(f12.3))'
         lfmt_out=29
         write(unit=iunit,fmt=fmt_out(1:lfmt_out))   tit, Polari%CS(1,1,1:Mag_Dom%nd),' :',Polari%CS(1,2,1:Mag_Dom%nd)
       end if
       tit="                      along y = "
       if(nch==1) then
         fmt_out(6:17)='(f12.3))'
         lfmt_out=17
         write(unit=iunit,fmt=fmt_out(1:lfmt_out))   tit, Polari%CS(2,1,1:Mag_Dom%nd)
       else if(nch==2) then
         fmt_out(6:16)='(f12.3)'
         fmt_out(17:19)=',a,'
         write(fmt_out(20:21),'(i2.2)') Mag_Dom%nd
         fmt_out(22:29)='(f12.3))'
         lfmt_out=29
         write(unit=iunit,fmt=fmt_out(1:lfmt_out))   tit, Polari%CS(2,1,1:Mag_Dom%nd),' :',Polari%CS(2,2,1:Mag_Dom%nd)
       end if
       tit="                      along z = "
       if(nch==1) then
         fmt_out(6:17)='(f12.3))'
         lfmt_out=17
         write(unit=iunit,fmt=fmt_out(1:lfmt_out))   tit, Polari%CS(3,1,1:Mag_Dom%nd)
       else if(nch==2) then
         fmt_out(6:16)='(f12.3)'
         fmt_out(17:19)=',a,'
         write(fmt_out(20:21),'(i2.2)') Mag_Dom%nd
         fmt_out(22:29)='(f12.3))'
         lfmt_out=29
         write(unit=iunit,fmt=fmt_out(1:lfmt_out))   tit, Polari%CS(3,1,1:Mag_Dom%nd),' :',Polari%CS(3,2,1:Mag_Dom%nd)
       end if



       Write(unit=iunit,fmt="(a,/)")           " "
       Write(unit=iunit,fmt="(a,/)")           " => Polarisation tensor as it will be measured:"
       Write(unit=iunit,fmt="(a,/)")           " "
       Write(unit=iunit,fmt="(3(a,f12.4), a)") "         /", Polari%Pij(1,1),"  ", Polari%Pij(1,2) , "  ", Polari%Pij(1,3), "  \"
       Write(unit=iunit,fmt="(3(a,f12.4), a)") "  PT  = | ", Polari%Pij(2,1),"  ", Polari%Pij(2,2) , "  ", Polari%Pij(2,3), "   |"
       Write(unit=iunit,fmt="(3(a,f12.4), a)") "         \", Polari%Pij(3,1),"  ", Polari%Pij(3,2) , "  ", Polari%Pij(3,3), "  /"

       return
    End Subroutine Write_Polar_Info

    !!----
    !!---- Subroutine Write_Polar_Line(Polari, Lun)
    !!----    Type (Polar_calc_type), intent( in)     :: Polrari !  in ->type with all information about polarization in one point hkl
    !!----    integer, optional,      intent(in)      :: lun     !  In -> Unit to write
    !!----
    !!----    Outputs the polarization info type in line form, so you can write it to a file
    !!----
    !!---- Update: May - 2005
    !!
    Subroutine Write_Polar_Line(Polari, Lun)
       !---- Arguments ----!
       Type (Polar_calc_type), intent( in)     :: Polari !
       integer, optional,      intent(in)      :: Lun

       !---- Local variables ----!
       integer            :: iunit

       iunit=6
       if (present(lun)) iunit=lun

       Write(unit=iunit,fmt="(/,a)")        "     H         K         L         NSF^2     NSF_r     NSF_i"
       Write(unit=iunit,fmt="(6(f10.6))")   Polari%H(1), Polari%H(2), Polari%H(3), Polari%NC, real(Polari%NSF), AIMAG(Polari%NSF)
       Write(unit=iunit,fmt="(/,a)")        "       Pix       Piy       Piz       Pfx       Pfy       Pfz"
       Write(unit=iunit,fmt="(6(f10.3))")   Polari%P, 0.0, 0.0, Polari%Pij(1,1), Polari%Pij(2,1), Polari%Pij(3,1)
       Write(unit=iunit,fmt="(6(f10.3))")   0.0,Polari%P, 0.0, Polari%Pij(1,2), Polari%Pij(2,2), Polari%Pij(3,2)
       Write(unit=iunit,fmt="(6(f10.3))")   0.0, 0.0, Polari%P, Polari%Pij(1,3), Polari%Pij(2,3), Polari%Pij(3,3)
       Write(unit=iunit,fmt="(a,/)")        "-------------------------------------------------------------------------"

       return
    End Subroutine Write_Polar_line

    !!----
    !!----    Subroutine Calc_Polar_Dom_Efficiency(Cell,H,SPV,Pin,NSF,Mag_dom,Mh,Polari)
    !!----      type (Crystal_Cell_Type),    intent(in)     :: Cell
    !!----      REAL(kind=cp), DIMENSION (3),intent(in)     :: H
    !!----      Real(kind=cp), dimension(3), intent(in)     :: SPV
    !!----      Real(kind=cp),               intent(in)     :: Pin
    !!----      complex(kind=cp),            intent(in)     :: NSF
    !!----      type(Magnetic_Domain_type),  intent(in)     :: Mag_Dom
    !!----      type(MagHD_Type),            intent(in out) :: Mh
    !!----      type(Polar_calc_type),       intent(out)    :: Polari
    !!----
    !!----   Calculates Polarization matrix from polarised cross-sections, accounts for partial polarization
    !!----
    !!----   Created: November - 2011 OZ
    !!

    Subroutine Calc_Polar_Dom_Efficiency(Cell,H,SPV,Pin,NSF,Mag_dom,Mh,Polari)

      Type (Crystal_Cell_Type),    intent(in)    :: Cell
      Real(kind=cp), dimension (3),intent(in)    :: H
      Real(kind=cp), dimension(3), intent(in)    :: SPV
      Real(kind=cp),               intent(in)    :: Pin
      complex(kind=cp),            intent(in)    :: NSF
      Type(Magnetic_Domain_type),  intent(in)    :: Mag_Dom
      Type(MagHD_Type),            intent(in out):: Mh
      Type(Polar_calc_type),       intent(out)   :: Polari

    !!---- Local variables ----!
      integer                           :: nd,ich,nch,i,j
      complex(kind=cp), dimension(2,2)  :: ScatAmp
      complex(kind=cp), dimension(2,2)  :: Spin_Px,Spin_Py,Spin_Pz
      Real(kind=cp)                     :: coef, Pinm, Pf !, A
      complex(kind=cp), dimension(3)    :: MiV, MiV_PF  !MiV for one domain and in polarisation frame
      complex(kind=cp), dimension(3,3,4):: sVs          ! 1dim in x,y,z 2dim out x,y,z 3dim sign 1++ 2+- 3-+ 4--
      Real(kind=cp), dimension(3,3,4)   :: CrSec
      Real(kind=cp), dimension(3,3)     :: Ipp,Ipm,Imp,Imm
      real(kind=cp), parameter          :: eps=0.00001_cp

       ! Shortcut for TASP as two benders have same efficiency
       Pf=Pin
       ! Neutron spin states
       coef=1.0_cp/sqrt(2.0_cp)

       Spin_Px(1,:)=(/coef*(1.0_cp, 0.0_cp),coef*( 1.0_cp, 0.0_cp)/) !up
       Spin_Px(2,:)=(/coef*(1.0_cp, 0.0_cp),coef*(-1.0_cp, 0.0_cp)/) !down

       Spin_Py(1,:)=(/coef*(1.0_cp, 0.0_cp),coef*( 0.0_cp, 1.0_cp)/) !up
       Spin_Py(2,:)=(/coef*(1.0_cp, 0.0_cp),coef*( 0.0_cp,-1.0_cp)/) !down

       Spin_Pz(1,:)=(/(1.0_cp, 0.0_cp),( 0.0_cp, 0.0_cp)/) !up
       Spin_Pz(2,:)=(/(0.0_cp, 0.0_cp),(-1.0_cp, 0.0_cp)/) !down

       !A = tpi**3/Cell%CellVol   !not used here
       !First store given info in Polari
       Polari%H = H
       Polari%SPV = SPV
       Polari%Cell = Cell
       Polari%P = Pin
       Polari%NSF = NSF
       Polari%Pij(:,:) = 0.0_cp
       Ipp = 0.0_cp
       Ipm = 0.0_cp
       Imp = 0.0_cp
       Imm = 0.0_cp

       nch=1
       if(Mag_Dom%chir) nch=2
       do nd=1,Mag_Dom%nd
         do ich=1,nch
           MiV=Mh%MiVC(:,ich,nd) !use Cartesian components
           ! MiV=Mh%MiV(:,ich,nd) !this does not use Cartesian components as was in the original sulbroutine
           ! Magnetic Interaction Vector in polarisation frame
           MiV_PF = Magn_Inter_Vec_PF(MiV,H,SPV, Cell)
           ! Partial Scattering Amplitudes
           ! U++ U-+
           ! U+- U--
           ScatAmp(1,1)=NSF+MiV_PF(3)
           ScatAmp(1,2)=MiV_PF(1)-(0.0_cp,1.0_cp)*MiV_PF(2)
           ScatAmp(2,1)=MiV_PF(1)+(0.0_cp,1.0_cp)*MiV_PF(2)
           ScatAmp(2,2)=NSF-MiV_PF(3)
           ! Matrix elements between two spin states <spin|Potential|spin>
           do i=1,2
             sVs(1,1,i)=dot_product(Spin_Px(i,:),matmul(ScatAmp,Spin_Px(1,:))) ! sVs xx,-xx
             sVs(2,1,i)=dot_product(Spin_Py(i,:),matmul(ScatAmp,Spin_Px(1,:))) ! sVs yx,-yx
             sVs(3,1,i)=dot_product(Spin_Pz(i,:),matmul(ScatAmp,Spin_Px(1,:))) ! sVs zx,xzx

             sVs(1,2,i)=dot_product(Spin_Px(i,:),matmul(ScatAmp,Spin_Py(1,:))) ! sVs xy,-xy
             sVs(2,2,i)=dot_product(Spin_Py(i,:),matmul(ScatAmp,Spin_Py(1,:))) ! sVs yy,-yy
             sVs(3,2,i)=dot_product(Spin_Pz(i,:),matmul(ScatAmp,Spin_Py(1,:))) ! sVs zy,-zy

             sVs(1,3,i)=dot_product(Spin_Px(i,:),matmul(ScatAmp,Spin_Pz(1,:))) ! sVs xz,-xz
             sVs(2,3,i)=dot_product(Spin_Py(i,:),matmul(ScatAmp,Spin_Pz(1,:))) ! sVs yz,-yz
             sVs(3,3,i)=dot_product(Spin_Pz(i,:),matmul(ScatAmp,Spin_Pz(1,:))) ! sVs zz,-zz
           end do

           do i=1,2
             sVs(1,1,i+2)=dot_product(Spin_Px(i,:),matmul(ScatAmp,Spin_Px(2,:))) ! sVs x-x,-x-x
             sVs(2,1,i+2)=dot_product(Spin_Py(i,:),matmul(ScatAmp,Spin_Px(2,:))) ! sVs y-x,-y-x
             sVs(3,1,i+2)=dot_product(Spin_Pz(i,:),matmul(ScatAmp,Spin_Px(2,:))) ! sVs z-x,-z-x

             sVs(1,2,i+2)=dot_product(Spin_Px(i,:),matmul(ScatAmp,Spin_Py(2,:))) ! sVs x-y,-x-y
             sVs(2,2,i+2)=dot_product(Spin_Py(i,:),matmul(ScatAmp,Spin_Py(2,:))) ! sVs y-y,-y-y
             sVs(3,2,i+2)=dot_product(Spin_Pz(i,:),matmul(ScatAmp,Spin_Py(2,:))) ! sVs z-y,-z-y

             sVs(1,3,i+2)=dot_product(Spin_Px(i,:),matmul(ScatAmp,Spin_Pz(2,:))) ! sVs x-z,-x-z
             sVs(2,3,i+2)=dot_product(Spin_Py(i,:),matmul(ScatAmp,Spin_Pz(2,:))) ! sVs y-z,-y-z
             sVs(3,3,i+2)=dot_product(Spin_Pz(i,:),matmul(ScatAmp,Spin_Pz(2,:))) ! sVs z-z,-z-z
           end do

           ! Scattering Cross-Sections and Polarization Matrices
           ! 1st direction is scattered, 2nd - incoming
           CrSec=CONJG(sVs)*sVs

           if(Pin > 0.0_cp) then
             do i=1,3
               do j=1,3
                Ipp(i,j) = Ipp(i,j) + Mag_Dom%Pop(ich,nd)*                      &
                         ( CrSec(i,j,1)*Pin*Pf + CrSec(i,j,2)*Pin*(1.0_cp-Pf) + &
                           CrSec(i,j,3)*(1.0_cp-Pin)*Pf + CrSec(i,j,4)*(1.0_cp-Pin)*(1.0_cp-Pf) )
                Imp(i,j) = Imp(i,j) + Mag_Dom%Pop(ich,nd)*                      &
                         ( CrSec(i,j,2)*Pin*Pf + CrSec(i,j,1)*Pin*(1.0_cp-Pf) + &
                           CrSec(i,j,4)*(1.0_cp-Pin)*Pf + CrSec(i,j,3)*(1.0_cp-Pin)*(1.0_cp-Pf) )
               end do
             end do
           end if

           if(Pin < 0.0_cp) then
             Pinm=-Pin
             Pf=Pinm
             do i=1,3
               do j=1,3
                Imm(i,j) = Imm(i,j) + Mag_Dom%Pop(ich,nd)* &
                         ( CrSec(i,j,4)*Pinm*Pf + CrSec(i,j,3)*Pinm*(1.0_cp-Pf) + &
                           CrSec(i,j,2)*(1.0_cp-Pinm)*Pf + CrSec(i,j,1)*(1.0_cp-Pinm)*(1.0_cp-Pf) )
                Ipm(i,j) = Ipm(i,j) + Mag_Dom%Pop(ich,nd)* &
                         ( CrSec(i,j,3)*Pinm*Pf + CrSec(i,j,4)*Pinm*(1.0_cp-Pf) + &
                           CrSec(i,j,1)*(1.0_cp-Pinm)*Pf + CrSec(i,j,2)*(1.0_cp-Pinm)*(1.0_cp-Pf) )
               end do
             end do
           end if

         end do !loop over chiral domains
       end do !loop over S-domains

       if(Pin > 0.0_cp) then
         do i=1,3
           do j=1,3
             if(Ipp(i,j)+Imp(i,j) <= eps) then
               Polari%Pij(i,j) = 0.0_cp
             else
               Polari%Pij(i,j)=(Ipp(i,j)-Imp(i,j))/(Ipp(i,j)+Imp(i,j))
             end if
           end do
         end do
       end if

       if(Pin < 0.0_cp) then
         do i=1,3
           do j=1,3
             if(Imm(i,j)+Ipm(i,j) <= eps) then
              Polari%Pij(i,j) = 0.0_cp
             else
              Polari%Pij(i,j)=-(Imm(i,j)-Ipm(i,j))/(Imm(i,j)+Ipm(i,j))
             end if
           end do
         end do
       end if

       return
    End Subroutine Calc_Polar_Dom_Efficiency

    !!----
    !!---- Subroutine Calc_Polar_CrSec(Cell,H,SPV,Pin,NSF,Mag_dom,Mh,Ipp,Ipm,Imp,Imm)
    !!----
    !!----  Type (Crystal_Cell_Type),     intent(in)    :: Cell
    !!----  Real(kind=cp), dimension (3), intent(in)    :: H
    !!----  Real(kind=cp), dimension(3),  intent(in)    :: SPV
    !!----  Real(kind=cp),                intent(in)    :: Pin
    !!----  complex(kind=cp),             intent(in)    :: NSF
    !!----  Type(Magnetic_Domain_type),   intent(in)    :: Mag_Dom
    !!----  Type(MagHD_Type),             intent(in out):: Mh
    !!----  Real(kind=cp), dimension(3,3),intent(out)   :: Ipp,Ipm,Imp,Imm
    !!----
    !!---- Calculates polarised cross-sections, accounts for partial polarization
    !!---- is useful for MultiSourceData refinement
    !!----

    Subroutine Calc_Polar_CrSec(Cell,H,SPV,Pin,NSF,Mag_dom,Mh,Ipp,Ipm,Imp,Imm)

      Type (Crystal_Cell_Type),     intent(in)     :: Cell
      Real(Kind=Cp), dimension (3), intent(in)     :: H
      Real(kind=cp), dimension(3),  intent(in)     :: SPV
      Real(kind=cp),                intent(in)     :: Pin
      complex(kind=cp),             intent(in)     :: NSF
      Type(Magnetic_Domain_type),   intent(in)     :: Mag_Dom
      Type(MagHD_Type),             intent(in out) :: Mh
      Real(kind=cp), dimension(3,3),intent(out)    :: Ipp,Ipm,Imp,Imm

    !!---- Local variables ----!
      integer                           :: nd,ich,nch,i,j
      complex(kind=cp), dimension(2,2)  :: ScatAmp
      complex(kind=cp), dimension(2,2)  :: Spin_Px,Spin_Py,Spin_Pz
      Real(kind=cp)                     :: coef, Pinm, Pf !, A
      complex(kind=cp), dimension(3)    :: MiV, MiV_PF   !MiV for one domain and in polarisation frame
      complex(kind=cp), dimension(3,3,4):: sVs ! 1dim in x,y,z 2dim out x,y,z 3dim sign 1++ 2+- 3-+ 4--
      Real(kind=cp), dimension(3,3,4)   :: CrSec

      ! shortcut for TASP as two benders have same efficiency
      Pf=Pin
      ! Neutron spin states
      coef=1.0_cp/sqrt(2.0_cp)

      Spin_Px(1,:)=(/coef*(1.0_cp, 0.0_cp),coef*( 1.0_cp, 0.0_cp)/) !up
      Spin_Px(2,:)=(/coef*(1.0_cp, 0.0_cp),coef*(-1.0_cp, 0.0_cp)/) !down

      Spin_Py(1,:)=(/coef*(1.0_cp, 0.0_cp),coef*( 0.0_cp, 1.0_cp)/) !up
      Spin_Py(2,:)=(/coef*(1.0_cp, 0.0_cp),coef*( 0.0_cp,-1.0_cp)/) !down

      Spin_Pz(1,:)=(/(1.0_cp, 0.0_cp),( 0.0_cp, 0.0_cp)/) !up
      Spin_Pz(2,:)=(/(0.0_cp, 0.0_cp),(-1.0_cp, 0.0_cp)/) !down

      !A = tpi**3/Cell%CellVol   !not used here

      Ipp = 0.0_cp
      Ipm = 0.0_cp
      Imp = 0.0_cp
      Imm = 0.0_cp

      nch=1
      if(Mag_Dom%chir) nch=2
      do nd=1,Mag_Dom%nd
        do ich=1,nch
          MiV=Mh%MiVC(:,ich,nd) !MiV must be provided in Cartesian Crystallographic Frame
          ! Magnetic Interaction Vector in polarisation frame
          MiV_PF = Magn_Inter_Vec_PF(MiV,H,SPV, Cell)
          ! Partial Scattering Amplitudes
          ! U++ U-+
          ! U+- U--
          ScatAmp(1,1)=NSF+MiV_PF(3)
          ScatAmp(1,2)=MiV_PF(1)-(0.0_cp,1.0_cp)*MiV_PF(2)
          ScatAmp(2,1)=MiV_PF(1)+(0.0_cp,1.0_cp)*MiV_PF(2)
          ScatAmp(2,2)=NSF-MiV_PF(3)
          ! Matrix elements between two spin states <spin|Potential|spin>
          do i=1,2
            sVs(1,1,i)=dot_product(Spin_Px(i,:),matmul(ScatAmp,Spin_Px(1,:))) ! sVs xx,-xx
            sVs(2,1,i)=dot_product(Spin_Py(i,:),matmul(ScatAmp,Spin_Px(1,:))) ! sVs yx,-yx
            sVs(3,1,i)=dot_product(Spin_Pz(i,:),matmul(ScatAmp,Spin_Px(1,:))) ! sVs zx,xzx

            sVs(1,2,i)=dot_product(Spin_Px(i,:),matmul(ScatAmp,Spin_Py(1,:))) ! sVs xy,-xy
            sVs(2,2,i)=dot_product(Spin_Py(i,:),matmul(ScatAmp,Spin_Py(1,:))) ! sVs yy,-yy
            sVs(3,2,i)=dot_product(Spin_Pz(i,:),matmul(ScatAmp,Spin_Py(1,:))) ! sVs zy,-zy

            sVs(1,3,i)=dot_product(Spin_Px(i,:),matmul(ScatAmp,Spin_Pz(1,:))) ! sVs xz,-xz
            sVs(2,3,i)=dot_product(Spin_Py(i,:),matmul(ScatAmp,Spin_Pz(1,:))) ! sVs yz,-yz
            sVs(3,3,i)=dot_product(Spin_Pz(i,:),matmul(ScatAmp,Spin_Pz(1,:))) ! sVs zz,-zz
          end do

          do i=1,2
            sVs(1,1,i+2)=dot_product(Spin_Px(i,:),matmul(ScatAmp,Spin_Px(2,:))) ! sVs x-x,-x-x
            sVs(2,1,i+2)=dot_product(Spin_Py(i,:),matmul(ScatAmp,Spin_Px(2,:))) ! sVs y-x,-y-x
            sVs(3,1,i+2)=dot_product(Spin_Pz(i,:),matmul(ScatAmp,Spin_Px(2,:))) ! sVs z-x,-z-x

            sVs(1,2,i+2)=dot_product(Spin_Px(i,:),matmul(ScatAmp,Spin_Py(2,:))) ! sVs x-y,-x-y
            sVs(2,2,i+2)=dot_product(Spin_Py(i,:),matmul(ScatAmp,Spin_Py(2,:))) ! sVs y-y,-y-y
            sVs(3,2,i+2)=dot_product(Spin_Pz(i,:),matmul(ScatAmp,Spin_Py(2,:))) ! sVs z-y,-z-y

            sVs(1,3,i+2)=dot_product(Spin_Px(i,:),matmul(ScatAmp,Spin_Pz(2,:))) ! sVs x-z,-x-z
            sVs(2,3,i+2)=dot_product(Spin_Py(i,:),matmul(ScatAmp,Spin_Pz(2,:))) ! sVs y-z,-y-z
            sVs(3,3,i+2)=dot_product(Spin_Pz(i,:),matmul(ScatAmp,Spin_Pz(2,:))) ! sVs z-z,-z-z
          end do

         ! Scattering Cross-Sections and Polarization Matrices
         ! 1st direction is scattered, 2nd - incoming
         CrSec=CONJG(sVs)*sVs

         if(Pin > 0.0_cp) then
           do i=1,3
             do j=1,3
              Ipp(i,j) = Ipp(i,j) + Mag_Dom%Pop(ich,nd)* &
                       ( CrSec(i,j,1)*Pin*Pf + CrSec(i,j,2)*Pin*(1.0_cp-Pf) + &
                         CrSec(i,j,3)*(1.0_cp-Pin)*Pf + CrSec(i,j,4)*(1.0_cp-Pin)*(1.0_cp-Pf) )
              Imp(i,j) = Imp(i,j) + Mag_Dom%Pop(ich,nd)* &
                       ( CrSec(i,j,2)*Pin*Pf + CrSec(i,j,1)*Pin*(1.0_cp-Pf) + &
                         CrSec(i,j,4)*(1.0_cp-Pin)*Pf + CrSec(i,j,3)*(1.0_cp-Pin)*(1.0_cp-Pf) )
             end do
           end do
         end if

         if(Pin < 0.0_cp) then
           Pinm=-Pin
           Pf=Pinm
           do i=1,3
             do j=1,3
              Imm(i,j) = Imm(i,j) + Mag_Dom%Pop(ich,nd)*                        &
                       ( CrSec(i,j,4)*Pinm*Pf + CrSec(i,j,3)*Pinm*(1.0_cp-Pf) + &
                         CrSec(i,j,2)*(1.0_cp-Pinm)*Pf + CrSec(i,j,1)*(1.0_cp-Pinm)*(1.0_cp-Pf) )
              Ipm(i,j) = Ipm(i,j) + Mag_Dom%Pop(ich,nd)*                        &
                       ( CrSec(i,j,3)*Pinm*Pf + CrSec(i,j,4)*Pinm*(1.0_cp-Pf) + &
                         CrSec(i,j,1)*(1.0_cp-Pinm)*Pf + CrSec(i,j,2)*(1.0_cp-Pinm)*(1.0_cp-Pf) )
             end do
           end do
         end if

        end do !loop over chiral domains
       end do !loop over S-domains

      return
    End Subroutine Calc_Polar_CrSec

 End Module CFML_Polarimetry
