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
!!----
!!---- MODULE: CFML_Magnetic_Structure_Factors
!!----   INFO: Main module for Structure Factors Calculations
!!----
!!---- HISTORY
!!----    Update: 07/03/2011
!!----
!!----
!!---- DEPENDENCIES
!!--++     Use CFML_GlobalDeps,                  only: sp, tpi
!!--++     Use CFML_Math_General,                only: atan2d, sort
!!--++     Use CFML_String_Utilities,            only: L_Case,U_Case, Get_LogUnit
!!--++     Use CFML_Scattering_Chemical_Tables,  only: Set_Magnetic_Form, Remove_Magnetic_Form, num_mag_form, &
!!--++                                                 Magnetic_Form
!!--++     Use CFML_Crystal_Metrics,             only: Crystal_Cell_type
!!--++     Use CFML_Crystallographic_Symmetry,   only: Space_Group_Type, Set_spaceGroup
!!--++     Use CFML_Magnetic_Symmetry,           only: ApplyMSO, MagSymm_k_type, Magnetic_Group_Type, Magnetic_Domain_type
!!--++     Use CFML_Reflections_Utilities,       only: HKL_R, HKL_Gen, Get_MaxNumRef, Reflect_Type, Reflection_List_Type, hkl_s
!!--++     Use CFML_Atom_TypeDef,                only: Matom_list_type
!!--++     Use CFML_Propagation_Vectors,         only: K_Equiv_Minus_K
!!----
!!---- VARIABLES
!!--..    Types
!!--++    HR_TYPE                      [Private]
!!----    MAGH_TYPE
!!----    MAGH_LIST_TYPE
!!----    MAGHD_TYPE
!!----    MAGHD_LIST_TYPE
!!--++    AJH                          [Private]
!!--++    BJH                          [Private]
!!----    ERR_MSFAC
!!----    ERR_MSFAC_MESS
!!--++    HR                           [Private]
!!--++    HT                           [Private]
!!--++    MFI                          [Private]
!!--++    MFR                          [Private]
!!--++    MSF_INITIALIZED              [Private]
!!----    PN
!!--++    TH                           [Private]
!!----
!!---- PROCEDURES
!!----    Functions:
!!--++       MFJ                         [Private]
!!----
!!----    Subroutines:
!!----       CALC_MAG_INTERACTION_VECTOR
!!----       CALC_MAGNETIC_STRF_MIV
!!----       CALC_MAGNETIC_STRF_MIV_DOM
!!----       CALC_MAGNETIC_STRF_TENSOR
!!--++       CALC_TABLE_MAB              [Private]
!!--++       CALC_TABLE_TH               [Private]
!!--++       CREATE_TABLE_HR_HT          [Private]
!!--++       CREATE_TABLE_MFR            [Private]
!!----       GEN_SATELLITES
!!----       INIT_MAG_STRUCTURE_FACTORS
!!----       MAG_STRUCTURE_FACTORS
!!----       MODIFY_MSF
!!--++       SET_FIXED_TABLES            [Private]
!!--++       SUM_MAB                     [Private]
!!----       WRITE_MAG_STRUCTURE_FACTORS
!!----
!!
 Module CFML_Magnetic_Structure_Factors
    !---- Use Modules ----!
    Use CFML_GlobalDeps,                  only: cp, sp, dp, tpi
    Use CFML_Math_General,                only: atan2d, sort
    Use CFML_String_Utilities,            only: L_Case,U_Case, Get_LogUnit
    Use CFML_Scattering_Chemical_Tables,  only: Set_Magnetic_Form, Remove_Magnetic_Form, num_mag_form, &
                                                Magnetic_Form
    Use CFML_Crystal_Metrics,             only: Crystal_Cell_type
    Use CFML_Crystallographic_Symmetry,   only: Space_Group_Type, Set_spaceGroup
    Use CFML_Magnetic_Symmetry,           only: ApplyMSO, MagSymm_k_type, Magnetic_Group_Type, Write_Magnetic_Structure, &
                                                Magnetic_Domain_type
    Use CFML_Reflections_Utilities,       only: HKL_R, HKL_Gen, Get_MaxNumRef, Reflect_Type, Reflection_List_Type, hkl_s
    Use CFML_Atom_TypeDef,                only: Matom_list_type, Allocate_mAtom_list
    Use CFML_Propagation_Vectors,         only: K_Equiv_Minus_K

    !---- Variables ----!
    implicit none

    private

    !---- List of public functions ----!

    !---- List of public subroutines ----!
    public :: Calc_Mag_Interaction_Vector, Gen_satellites, Init_Mag_Structure_Factors, &
              Mag_Structure_Factors, Modify_MSF, Write_Mag_Structure_Factors,          &
              Calc_Magnetic_StrF_MiV, Calc_Magnetic_StrF_MiV_dom, Init_Err_MSfac,      &
              Calc_Magnetic_Strf_Tensor

    !---- List of private functions ----!
    private :: mFj

    !---- List of private subroutines ----!
    private :: Calc_Table_MAB, Create_Table_mFR,  &
               Create_Table_HR_HT, Set_Fixed_Tables, Calc_Table_TH, Sum_MAB

    !---- Definitions ----!

    !!--++
    !!--++ TYPE :: HR_TYPE
    !!--..
    !!--++    Type, Private :: HR_Type
    !!--++       real(kind=cp), dimension(3) :: H
    !!--++    End Type HR_Type
    !!--++
    !!--++    (Private)
    !!--++    Define the scatering vector vector  H+k
    !!--++
    !!--++ Update: April - 2005
    !!
    Type, Private :: HR_Type
       real(kind=cp), dimension(3) :: H
    End Type HR_Type

    !!----
    !!---- TYPE :: MAGH_TYPE
    !!--..
    !!----    Type, Public  :: MagH_Type
    !!----       logical                         :: keqv_minus  !True if k equivalent to -k
    !!----       integer                         :: mult        !Multiplicity of the reflection (useful for powder calculations)
    !!----       integer                         :: num_k       !number of the propagation vector vk
    !!----       real(kind=cp)                   :: signp       !+1 for -vk   and -1 for +vk
    !!----       real(kind=cp)                   :: s           !sinTheta/Lambda
    !!----       real(kind=cp)                   :: sqMiV       !Square of the Magnetic Interaction vector
    !!----       real(kind=cp),    dimension(3)  :: H           ! H +/- k
    !!----       complex(kind=cp), dimension(3)  :: MsF         !magnetic structure factor
    !!----       complex(kind=cp), dimension(3,3):: TMsF        !tensorial magnetic structure factor
    !!----       complex(kind=cp), dimension(3)  :: MiV         !magnetic interaction vector
    !!----       complex(kind=cp), dimension(3)  :: MiVC        !magnetic interaction vector in Cartesian components
    !!----    End Type  MagH_Type
    !!----
    !!----    Define the scatering vector vector  H+k and the sign -1 for H+k and +1 for H-k.
    !!----    Includes the magnetic interaction vector MiV = Mper = M
    !!----
    !!---- Updated: April-2005, June - 2012, June -2014
    !!
    Type, Public  :: MagH_Type
       logical                          :: keqv_minus  !True if k equivalent to -k
       integer                          :: mult        !Multiplicity of the reflection (useful for powder calculations)
       integer                          :: num_k       !number of the propagation vector vk
       real(kind=cp)                    :: signp       !+1 for -vk   and -1 for +vk
       real(kind=cp)                    :: s           !sinTheta/Lambda
       real(kind=cp)                    :: sqMiV       !Square of the Magnetic Interaction vector
       real(kind=cp),    dimension(3)   :: H           ! H +/- k
       complex(kind=cp), dimension(3)   :: MsF         !Magnetic structure factor w.r.t. unitary Crystal Frame
       complex(kind=cp), dimension(3,3) :: TMsF        !tensorial magnetic structure factor
       complex(kind=cp), dimension(3)   :: MiV         !Magnetic interaction vector w.r.t. unitary Crystal Frame
       complex(kind=cp), dimension(3)   :: MiVC        !Magnetic interaction vector in Cartesian components w.r.t. Crystal Frame
    End Type  MagH_Type

    !!----
    !!----  TYPE :: MAGH_LIST_TYPE
    !!--..
    !!----     Type, Public  :: MagH_List_Type
    !!----        integer                                   :: Nref
    !!----        Type(MagH_Type),allocatable, dimension(:) :: Mh
    !!----     End Type MagH_List_Type
    !!----
    !!----     Define a list of magnetic reflections containing the
    !!----     scatering vector, the magnetic structure factor and
    !!----     the magnetic interaction vector.
    !!----
    !!---- Update: April - 2005
    !!
    Type, Public  :: MagH_List_Type
       integer                                   :: Nref
       Type(MagH_Type),allocatable, dimension(:) :: Mh
    End Type MagH_List_Type

    !!----
    !!----  TYPE :: MAGHD_TYPE
    !!--..
    !!----     Type, Public  :: MagHD_Type
    !!----        logical                            :: keqv_minus  !True if k equivalent to -k
    !!----        integer                            :: num_k       !number of the propagation vector vk
    !!----        real(kind=cp)                      :: signp       !+1 for -vk   and -1 for +vk
    !!----        real(kind=cp)                      :: s           !sinTheta/Lambda
    !!----        real(kind=cp)                      :: sqAMiV      !Square of the Average Magnetic Interaction vector
    !!----        real(kind=cp)                      :: sqMiV       !Average of the Square of Magnetic Interaction vectors
    !!----        real(kind=cp),   dimension(3)      :: H           ! H +/- k
    !!----        complex(kind=cp),dimension(3,2,24) :: MsF         !Magnetic structure factors of each domain (second dimension for chirality domains)
    !!----        complex(kind=cp),dimension(3,2,24) :: MiV         !Magnetic interaction vector of each domain
    !!----        complex(kind=cp),dimension(3,2,24) :: MiVC        !Magnetic interaction vector of each domain w.r.t. to Cartesian Crystal Frame
    !!----        complex(kind=cp),dimension(3)      :: AMiV        !Average Magnetic interaction vector = 1/nd Sum{ pop(i) Miv(:,i)} in Cartesian Frame
    !!----     End Type  MagHD_Type
    !!----
    !!----    Define the scatering vector vector  H+k and the sign -1 for H+k and +1 for H-k.
    !!----    Includes the average magnetic interaction vector AMiV(:) = 1/nd Sum[i]{ pop(i) MiVC(:,i)}
    !!----    This type should be used whenever magnetic domains are present (single crystal work)
    !!----
    !!---- Updated: November - 2006, June 2012 (JRC)
    !!
    Type, Public  :: MagHD_Type
       logical                            :: keqv_minus
       integer                            :: num_k     !number of the propagation vector vk
       real(kind=cp)                      :: signp     !+1 for -vk   and -1 for +vk
       real(kind=cp)                      :: s         !sinTheta/Lambda
       real(kind=cp)                      :: sqAMiV    !Square of the Average Magnetic Interaction vector
       real(kind=cp)                      :: sqMiV     !Average of the Square of Magnetic Interaction vectors
       real(kind=cp),   dimension(3)      :: H         ! H +/- k
       complex(kind=cp),dimension(3,2,24) :: MsF       !Magnetic structure factors of each domain (second dimension for chirality domains)
       complex(kind=cp),dimension(3,2,24) :: MiV       !Magnetic interaction vector of each domain
       complex(kind=cp),dimension(3,2,24) :: MiVC      !Magnetic interaction vector of each domain in Cartesian Crystal space
       complex(kind=cp),dimension(3)      :: AMiV      !Average Magnetic interaction vector = 1/nd Sum{ pop(i) Miv(:,i)}
    End Type  MagHD_Type

    !!----
    !!----  MAGHD_LIST_TYPE
    !!----    Type, Public  :: MagHD_List_Type
    !!----       integer                                    :: Nref
    !!----       Type(MagHD_Type),allocatable, dimension(:) :: Mh
    !!----    End Type MagHD_List_Type
    !!----
    !!----    Define a list of magnetic reflections containing the
    !!----    scatering vector, the magnetic structure factor and
    !!----    the magnetic interaction vector for each of the domains.
    !!----
    !!---- Update: February - 2009 Oksana Zaharko
    !!
    Type, Public  :: MagHD_List_Type
       integer                                    :: Nref
       Type(MagHD_Type),allocatable, dimension(:) :: Mh
    End Type MagHD_List_Type

    !!--++
    !!--++ AJH
    !!--++     real(kind=cp), dimension(:,:,:), allocatable, private :: Ajh
    !!--++
    !!--++     (Private)
    !!--++     Array for Aj(h). The dimensions are Ajh(3,Natoms,Nref)
    !!--++     where the magnetic structure factor vector is
    !!--++           M(h)=p Sum_j [Fj(h){Aj(h)+i Bj(h)}]
    !!--++
    !!--++ Update: April - 2005
    !!
    real(kind=cp), dimension(:,:,:), allocatable, private :: AJH

    !!--++
    !!--++ BJH
    !!--++     real(kind=cp), dimension(:,:,:), allocatable, private :: Bjh
    !!--++
    !!--++     (Private)
    !!--++     Array for Bj(h). The dimensions are Bjh(3,Natoms,Nref)
    !!--++     where the magnetic structure factor vector is
    !!--++           M(h)=pSum_j[Fj(h){Aj(h)+i Bj(h)}]
    !!--++
    !!--++ Update: April - 2005
    !!
    real(kind=cp), dimension(:,:,:), allocatable, private :: BJH

    !!----
    !!---- ERR_MSFAC_MESS
    !!----    character(len=150), public :: err_msfac_mess
    !!----
    !!----    String containing information about the last error in Magnetic_CFML_Structure_Factors
    !!----
    !!---- Update: April - 2005
    !!
    character(len=150), public :: err_msfac_mess

    !!----
    !!---- ERR_MSFAC
    !!----    logical, public :: err_msfac
    !!----
    !!----    Logical Variable indicating an error in Magnetic_CFML_Structure_Factors
    !!----
    !!---- Update: April - 2005
    !!
    logical, public :: err_msfac

    !!--++
    !!--++ HR
    !!--++     Type(HR_Type), dimension(:,:), allocatable, private :: Hr
    !!--++
    !!--++     (Private)
    !!--++     Array for HR Calculations. The dimension are HR(Nsymop,NRef)
    !!--++     Transformed H by the rotational part of the crystallographic
    !!--++     symmetry operators:  HR = matmul (H, Rs)
    !!--++
    !!--++ Update: April - 2005
    !!
    type(HR_Type), dimension(:,:), allocatable, private :: HR

    !!--++
    !!--++ HT
    !!--++     real(kind=cp), dimension(:,:), allocatable, private :: Ht
    !!--++
    !!--++     (Private)
    !!--++     Array for HT Calculations. The dimension are HT(Nsymop,Nref)
    !!--++     Scalar products of H.t, real sccatering vector by the translational
    !!--++     parts of the crystallographic symmetry operators: HT= dot_product(H,Ts)
    !!--++
    !!--++ Update: april - 2005
    !!
    real(kind=cp), dimension(:,:), allocatable, private :: HT

    !!--++
    !!--++ MFI
    !!--++     real(kind=cp), dimension(:), allocatable, private :: AFPP
    !!--++
    !!--++     (Private)
    !!--++     Array for imaginary part of magnetic form factor.
    !!--++     The dimension is: mFI(Natoms,NRef)
    !!--++
    !!--++ Update: April - 2005
    !!
    real(kind=cp), dimension(:,:), allocatable, private :: mFI

    !!--++
    !!--++ MFR
    !!--++     real(kind=cp), dimension(:,:), allocatable, private :: mFR
    !!--++
    !!--++     (Private)
    !!--++     Array for Magnetic form Factors. The dimensions are mFR(Natoms,NRef)
    !!--++
    !!--++ Update: April - 2005
    !!
    real(kind=cp), dimension(:,:), allocatable, private :: mFR

    !!--++
    !!--++ MSF_INITIALIZED
    !!--++    logical, private :: MSF_Initialized
    !!--++
    !!--++    (Private)
    !!--++    Logical Variable indicating if the module has been initialized.
    !!--++
    !!--++ Update: april - 2005
    !!
    logical, private :: MSF_Initialized=.false.

    !!----
    !!---- PN
    !!----    real(kind=dp), parameter, public :: pn=0.2695420113693928312
    !!----
    !!----    pn=Constant=  1/2 * gamma (µN) * r0 (units of 10^-12 cm)
    !!----    gamma: magnetic moment of neutrons in nuclear magnetons = 1.91304272(45) µN
    !!----    r0   : Classical radius of the electron = 0.28179403267(27) × 10-12 cm
    !!----
    !!---- Update: May - 2015
    !!
    real(kind=dp), parameter, public :: pn=0.2695420113693928312

    !!--++
    !!--++ TH
    !!--++    real(kind=cp), dimension(:,:), allocatable, private :: Th
    !!--++
    !!--++    (Private)
    !!--++    Array for TH Calculations. The dimension are TH(Natoms,Nref)
    !!--++
    !!--++ Update: April - 2005
    !!
    real(kind=cp), dimension(:,:), allocatable, private :: TH

 Contains

    !-------------------!
    !---- Functions ----!
    !-------------------!

    !!--++
    !!--++ Function Mfj(S,Coeff)
    !!--++    real(kind=cp),             intent(in) :: s !sin Theta/Lambda
    !!--++    real(kind=cp),dimension(7),intent(in) :: coeff
    !!--++
    !!--++    (Private)
    !!--++    Magnetic form factor calculation according to:
    !!--++    Fj(s)=Sum_i{1,6,2}[Coeff(i)*exp(-Coeff(i+1)*s*s)] + Coeff(7)
    !!--++
    !!--++ Update: April - 2005
    !!
    Function Mfj(S,Coeff) Result(Res)
       !---- Arguments ----!
       real(kind=cp),             intent(in) :: s  !sin Theta/Lambda
       real(kind=cp),dimension(7),intent(in) :: coeff
       real(kind=cp)                         :: res

       !---- Local variables ----!
       integer :: i

       res=coeff(7)
       do i=1,6,2
          res=res + coeff(i)*exp(-coeff(i+1)*s*s)
       end do

       return
    End Function Mfj

    !---------------------!
    !---- Subroutines ----!
    !---------------------!

    !!----
    !!---- Subroutine Calc_Mag_Interaction_Vector(Reflex,Cell)
    !!----    type(MagH_List_Type),     intent(in out) :: Reflex
    !!----    type(Crystal_Cell_type),  intent(in)     :: Cell
    !!----
    !!----    Calculate the Magnetic Interaction vector from Magnetic
    !!----    Structure factors, reflections and cell parameters.
    !!----    The components are given with respect to the crystallographic
    !!----    unitary direct cell system: {e1,e2,e3} and with respect to
    !!----    the Cartesian frame defined in Cell.
    !!----
    !!---- Updated: April - 2005, June-2012 (JRC,"mode" removed)
    !!
    Subroutine Calc_Mag_Interaction_Vector(Reflex,Cell)
       !---- Argument ----!
       type(MagH_List_Type),     intent(in out) :: Reflex
       type(Crystal_Cell_type),  intent(in)     :: Cell

       !---- Local variables ----!
       integer                       :: j
       real(kind=cp)                 :: s
       real(kind=cp),dimension(3)    :: ed,er
       complex(kind=cp),dimension(3) :: M, MiV,MiVC

       !---- Calculation of the Magnetic Interaction vector (in unitary and Cartesian Crystal Frames)----!
       do j=1,reflex%nref
          s  = 2.0*reflex%Mh(j)%s  !1/d=r*, M = M// + Mp   => Mp = M - M// = M - (M.e)e
          er = reflex%Mh(j)%h/s    !unitary vector referred to the reciprocal basis
          ed = matmul(cell%GR,er)  !  "        "       "             direct    "
          M  = reflex%Mh(j)%MsF / Cell%cell    !Magnetic structure factor in basis {a,b,c}
          MiV = M - dot_product(er,M) * ed     !Magnetic interaction vector in basis {a,b,c}
          reflex%Mh(j)%MiV =  MiV * Cell%cell  !Magnetic interaction vector in basis {e1,e2,e3}
          MiVC  = matmul(Cell%Cr_Orth_cel,MiV) !Magnetic interaction vector in Cartesian components
          reflex%Mh(j)%MiVC =  MiVC
          reflex%Mh(j)%sqMiV= dot_product(MiVC, MiVC)
       end do

       return
    End Subroutine Calc_Mag_Interaction_Vector

    !!----
    !!---- Subroutine Calc_Magnetic_Strf_Miv(Cell,Mgp,Atm,Mh)
    !!----    type(Crystal_Cell_type),  intent(in)     :: Cell
    !!----    type(MagSymm_k_Type),     intent(in)     :: MGp
    !!----    type(Matom_list_type),    intent(in)     :: Atm
    !!----    type(MagH_Type),          intent(in out) :: Mh
    !!----
    !!----    Calculate the Magnetic Interaction vector from Magnetic
    !!----    Structure factors, reflections and cell parameters.
    !!----    Whatever kind of settings for symmetry operators is allowed.
    !!----    The components are given with respect to the crystallographic
    !!----    unitary direct cell system: {e1,e2,e3} and with respect to the
    !!----    Cartesian frame defined in Cell.
    !!----
    !!---- Updated: April - 2005, June 2012, November 2014 (JRC)
    !!
    Subroutine Calc_Magnetic_Strf_Miv(Cell,Mgp,Atm,Mh)
       !---- Arguments ----!
       type(Crystal_Cell_type),  intent(in)     :: Cell
       type(MagSymm_k_Type),     intent(in)     :: MGp
       type(Matom_list_type),    intent(in)     :: Atm
       type(MagH_Type),          intent(in out) :: Mh

       !---- Local Variables ----!
       integer                            :: i,j,k,nvk,m, n
       real(kind=cp)                      :: arg,anis,onh,ph,s,b,ht,mFF,tho, isig, x
       real(kind=cp),    dimension(3)     :: h,ed,er
       real(kind=cp),    dimension(6)     :: beta
       real(kind=cp),    dimension(3,3)   :: Mcos,Msin
       real(kind=cp),    dimension(3)     :: ar,ai,br,bi,Ajh,Bjh,aa,bb
       complex(kind=cp)                   :: ci
       complex(kind=cp), dimension(3)     :: Mc, MiV, Sk, GMh

       s=Mh%s
       if (.not. Mh%keqv_minus ) then
          onh=0.5
       else
          onh=1.0
       end if
       nvk= Mh%num_k
       aa=0.0; bb=0.0
       isig=Mh%signp

       if (MGp%nirreps == 0) then

          do i=1,Atm%natoms
             m= Atm%Atom(i)%imat(nvk)
             if(m == 0) cycle  !Calculate only with contributing atoms

             !---- Isotropic Debye-Waller factor * occupation * p=0.5*re*gamma * Magnetic form-factors mFF
             b=atm%atom(i)%biso
             j=atm%atom(i)%ind(2)  !pointer to the magnetic form factor coefficients
             mFF=mfj(s,Magnetic_Form(j)%SctM)
             tho= pn*atm%atom(i)%occ*mFF*exp(-b*s*s)

             Mcos=0.0 ; Msin=0.0

             do k=1,MGp%NumOps
                h=Hkl_R(Mh%h,MGp%symop(k))
                ht=dot_product(Mh%h,MGp%SymOp(k)%Tr)
                ph= isig * (Atm%atom(i)%Mphas(m) + MGp%MSymOp(k,m)%Phas)
                arg=tpi*(dot_product(h,Atm%atom(i)%x)+ht + ph)
                anis=1.0
                if (Atm%atom(i)%thtype == "aniso") then
                   beta=Atm%atom(i)%u(1:6)
                   anis=     h(1)*h(1)*beta(1)+     h(2)*h(2)*beta(2)+    h(3)*h(3)*beta(3) &
                        +2.0*h(1)*h(2)*beta(4)+ 2.0*h(1)*h(3)*beta(5)+2.0*h(2)*h(3)*beta(6)
                   anis=exp(-anis)
                end if
                Mcos(:,:)=Mcos(:,:)+cos(arg)*anis*MGp%MSymOp(k,m)%Rot(:,:)
                Msin(:,:)=Msin(:,:)+sin(arg)*anis*MGp%MSymOp(k,m)%Rot(:,:)
             end do ! symmetry

             ar =       onh*matmul(Mcos,Atm%atom(i)%SkR(:,nvk)/Cell%cell)*Cell%cell  !The introduction of Cell%cell
             ai = -isig*onh*matmul(Mcos,Atm%atom(i)%SkI(:,nvk)/Cell%cell)*Cell%cell  !is for handling the possibility
             br =       onh*matmul(Msin,Atm%atom(i)%SkR(:,nvk)/Cell%cell)*Cell%cell  !of using non conventional settings for
             bi = -isig*onh*matmul(Msin,Atm%atom(i)%SkI(:,nvk)/Cell%cell)*Cell%cell  !symmetry operators
             Ajh(:) = ar(:) - bi(:)
             Bjh(:) = br(:) + ai(:)

             aa(:)= aa(:) + tho*ajh(:)
             bb(:)= bb(:) + tho*bjh(:)

          end do ! Atoms
          Mh%MsF(:)=cmplx(aa(:),bb(:)) * MGp%Num_Lat * MGp%Centred

       else  !Now magnetic structure described in terms of basis functions (No magnetic rotation matrices are provided)

          Mh%MsF(:)=cmplx(0.0,0.0)
          do i=1,Atm%natoms
             m= Atm%Atom(i)%imat(nvk)
             if(m == 0) cycle  !Calculate only with contributing atoms
             !---- Isotropic Debye-Waller factor * occupation * p=0.5*re*gamma * Magnetic form-factors mFF
             b=atm%atom(i)%biso
             j=atm%atom(i)%ind(2)  !pointer to the magnetic form factor coefficients
             mFF=mfj(s,Magnetic_Form(j)%SctM)
             tho= pn*atm%atom(i)%occ*mFF*exp(-b*s*s)

             GMh(:)=cmplx(0.0,0.0)

             do k=1,MGp%NumOps
                h=Hkl_R(Mh%h,MGp%symop(k))
                ht=dot_product(Mh%h,MGp%SymOp(k)%Tr)
                ph= isig*Atm%atom(i)%Mphas(m)
                arg=tpi*(dot_product(h,Atm%atom(i)%x)+ht + ph)
                anis=1.0
                if (Atm%atom(i)%thtype == "aniso") then
                   beta=Atm%atom(i)%u(1:6)
                   anis=     h(1)*h(1)*beta(1)+     h(2)*h(2)*beta(2)+    h(3)*h(3)*beta(3) &
                        +2.0*h(1)*h(2)*beta(4)+ 2.0*h(1)*h(3)*beta(5)+2.0*h(2)*h(3)*beta(6)
                   anis=exp(-anis)
                end if
                Sk(:) = cmplx(0.0,0.0)
                do n=1,abs(MGp%nbas(m)) !cannot be greater than 12 at present
                   x=real(MGp%icomp(n,m))
                   ci=cmplx(1.0-x,-isig*x)
                   Sk(:)=Sk(:)+ Atm%atom(i)%cbas(n,nvk)*ci* cmplx(  Real(MGp%basf(:,n,k,m)), -isig*aimag(MGp%basf(:,n,k,m))  )
                end do
                GMh(:)=GMh(:) + anis*Sk(:)*CMPLX(COS(arg),SIN(arg))  !Atomic contribution to geometric Magnetic Structure Factor
             end do ! symmetry
             Mh%MsF(:)=Mh%MsF(:) + tho*onh*GMh(:)
          end do ! Atoms
       end if
       Mh%MsF(:)=Mh%MsF(:) * MGp%Num_Lat * MGp%Centred
       !---- Calculation of the Magnetic Interaction vector ----!
       s  = 2.0*Mh%s            !1/d=r*, M = M// + Mp   => Mp = M - M// = M - (M.e)e
       er = Mh%h/s              !unitary vector referred to the reciprocal basis
       ed = matmul(cell%GR,er)  !  "        "       "             direct    "
       Mc  = Mh%MsF / Cell%cell                !Magnetic structure factor in basis {a,b,c}
       MiV = Mc - dot_product(er,Mc) * ed      !Magnetic interaction vector in basis {a,b,c}
       Mh%MiV  =  MiV * Cell%cell              !Magnetic interaction vector in basis {e1,e2,e3}
       Mh%MiVC = matmul(Cell%Cr_Orth_cel,MiV)  !Magnetic interaction vector in Cartesian components
       Mh%sqMiV= dot_product(Mh%MiVC, Mh%MiVC)

       return
    End Subroutine Calc_Magnetic_StrF_MiV

    !!----
    !!---- Subroutine Calc_Magnetic_Strf_Miv_Dom(Cell,Mgp,Atm,Mag_Dom,Mh)
    !!----    type(Crystal_Cell_type),   intent(in)     :: Cell
    !!----    type(MagSymm_k_Type),      intent(in)     :: MGp
    !!----    type(Matom_list_type),     intent(in)     :: Atm
    !!----    type(Magnetic_Domain_type),intent(in)     :: Mag_Dom
    !!----    type(MagHD_Type),          intent(in out) :: Mh
    !!----
    !!----    Calculate the Magnetic Interaction vector from Magnetic
    !!----    Structure factors, reflections and cell parameters.
    !!----    Whatever kind of settings for symmetry operators is allowed.
    !!----    The components are given with respect to the crystallographic
    !!----    unitary direct cell system: {e1,e2,e3} and with respect to
    !!----    Cartesian frame defined in Cell.
    !!----    In this subroutine the presence of magnetic domains is
    !!----    taken into account
    !!----
    !!---- Updated: September - 2010, July-2012, November 2014 (JRC)
    !!
    Subroutine Calc_Magnetic_Strf_Miv_Dom(Cell,Mgp,Atm,Mag_Dom,Mh)
       !---- Arguments ----!
       type(Crystal_Cell_type),   intent(in)      :: Cell
       type(MagSymm_k_Type),      intent(in)      :: MGp
       type(Matom_list_type),     intent(in)      :: Atm
       type(Magnetic_Domain_type),intent(in)      :: Mag_Dom
       type(MagHD_Type),          intent(in out)  :: Mh

       !---- Local Variables ----!
       integer                            :: i,j,k,nvk,m, n, nd, ich, nch
       real(kind=cp)                      :: arg,anis,onh,ph,s,b,ht,mFF,tho, isig, x
       real(kind=cp),    dimension(3)     :: h,ed,er,h_dom,xpos
       real(kind=cp),    dimension(6)     :: beta
       real(kind=cp),    dimension(3,3)   :: Mcos,Msin
       real(kind=cp),    dimension(3)     :: ar,ai,br,bi,Ajh,Bjh,aa,bb,Skr,Ski

       complex(kind=cp)                   :: ci
       complex(kind=cp), dimension(3)     :: Mc, MiV, Sk, GMh
       real(kind=cp),dimension(2), parameter :: ch=(/1.0,-1.0/)

       s=Mh%s
       if (.not. Mh%keqv_minus ) then
          onh=0.5
       else
          onh=1.0
       end if
       nvk= Mh%num_k

       isig=Mh%signp
       nch=1
       if (Mag_Dom%chir) nch=2

       if (MGp%nirreps == 0) then

          do nd=1,Mag_Dom%nd
             Mh%MsF(:,:,nd)=cmplx(0.0,0.0)
             do ich=1,nch
                aa=0.0; bb=0.0
                if(Mag_Dom%twin) then
                  h_dom=matmul(Mh%h,real(Mag_Dom%Dmat(:,:,nd)))
                else
                  h_dom=Mh%h
                end if
                do i=1,Atm%natoms
                   m= Atm%Atom(i)%imat(nvk)
                   if (m == 0) cycle  !Calculate only with contributing atoms
                   if(Mag_Dom%twin) then
                     Skr=Atm%atom(i)%SkR(:,nvk)
                     Ski=Atm%atom(i)%SkI(:,nvk)
                   else
                     Skr= matmul(Mag_Dom%Dmat(:,:,nd),Atm%atom(i)%SkR(:,nvk))
                     Ski= matmul(Mag_Dom%Dmat(:,:,nd),ch(ich)*Atm%atom(i)%SkI(:,nvk))
                   end if
                   !---- Isotropic Debye-Waller factor * occupation * p=0.5*re*gamma * Magnetic form-factors mFF
                   b=atm%atom(i)%biso
                   j=atm%atom(i)%ind(2)  !pointer to the magnetic form factor coefficients
                   mFF=mfj(s,Magnetic_Form(j)%SctM)
                   tho= pn*atm%atom(i)%occ*mFF*exp(-b*s*s)

                   Mcos=0.0
                   Msin=0.0
                   xpos=Atm%atom(i)%x
                   if(Mag_Dom%trans) xpos=matmul(Mag_Dom%Dmat(:,:,nd),xpos) + Mag_Dom%Dt(:,nd)

                   do k=1,MGp%NumOps
                      h=Hkl_R(h_dom,MGp%symop(k))
                      ht=dot_product(h_dom,MGp%SymOp(k)%Tr)
                      ph= isig * (Atm%atom(i)%Mphas(m) + MGp%MSymOp(k,m)%Phas)
                      arg=tpi*(dot_product(h,xpos)+ht + ph)
                      anis=1.0
                      if (Atm%atom(i)%thtype == "aniso") then
                         beta=Atm%atom(i)%u(1:6)
                         anis=     h(1)*h(1)*beta(1)+     h(2)*h(2)*beta(2)+    h(3)*h(3)*beta(3) &
                              +2.0*h(1)*h(2)*beta(4)+ 2.0*h(1)*h(3)*beta(5)+2.0*h(2)*h(3)*beta(6)
                         anis=exp(-anis)
                      end if
                      Mcos(:,:)=Mcos(:,:)+cos(arg)*anis*MGp%MSymOp(k,m)%Rot(:,:)
                      Msin(:,:)=Msin(:,:)+sin(arg)*anis*MGp%MSymOp(k,m)%Rot(:,:)
                   end do ! symmetry

                   ar =       onh*matmul(Mcos,SkR(:)/Cell%cell)*Cell%cell  !The introduction of Cell%cell
                   ai = -isig*onh*matmul(Mcos,SkI(:)/Cell%cell)*Cell%cell  !is for handling the possibility
                   br =       onh*matmul(Msin,SkR(:)/Cell%cell)*Cell%cell  !of using non conventional settings for
                   bi = -isig*onh*matmul(Msin,SkI(:)/Cell%cell)*Cell%cell  !symmetry operators
                   Ajh(:) = ar(:) - bi(:)
                   Bjh(:) = br(:) + ai(:)

                   aa(:)= aa(:) + tho*ajh(:)
                   bb(:)= bb(:) + tho*bjh(:)

                end do ! Atoms
                Mh%MsF(:,ich,nd)=cmplx(aa(:),bb(:)) * MGp%Num_Lat * MGp%Centred
             end do ! Chirality Domains
          end do !Domains

       else  !Now magnetic structure described in terms of basis functions (No magnetic rotation matrices are provided)

          do nd=1,Mag_dom%nd
             Mh%MsF(:,:,nd)=cmplx(0.0,0.0)
             do ich=1,nch
                aa=0.0; bb=0.0
                if(Mag_Dom%twin) then
                  h_dom=matmul(Mh%h,real(Mag_Dom%Dmat(:,:,nd)))
                else
                  h_dom=Mh%h
                end if

                do i=1,Atm%natoms
                   m= Atm%Atom(i)%imat(nvk)
                   if (m == 0) cycle  !Calculate only with contributing atoms
                   xpos=Atm%atom(i)%x
                   if(Mag_Dom%trans) xpos=matmul(Mag_Dom%Dmat(:,:,nd),xpos) + Mag_Dom%Dt(:,nd)

                   !---- Isotropic Debye-Waller factor * occupation * p=0.5*re*gamma * Magnetic form-factors mFF
                   b=atm%atom(i)%biso
                   j=atm%atom(i)%ind(2)  !pointer to the magnetic form factor coefficients
                   mFF=mfj(s,Magnetic_Form(j)%SctM)
                   tho= pn*atm%atom(i)%occ*mFF*exp(-b*s*s)

                   GMh(:)=cmplx(0.0,0.0)

                   do k=1,MGp%NumOps
                      h=Hkl_R(h_dom,MGp%symop(k))
                      ht=dot_product(h_dom,MGp%SymOp(k)%Tr)
                      ph= isig*Atm%atom(i)%Mphas(m)
                      arg=tpi*(dot_product(h,xpos)+ht + ph)
                      anis=1.0
                      if (Atm%atom(i)%thtype == "aniso") then
                         beta=Atm%atom(i)%u(1:6)
                         anis=     h(1)*h(1)*beta(1)+     h(2)*h(2)*beta(2)+    h(3)*h(3)*beta(3) &
                              +2.0*h(1)*h(2)*beta(4)+ 2.0*h(1)*h(3)*beta(5)+2.0*h(2)*h(3)*beta(6)
                         anis=exp(-anis)
                      end if
                      Sk(:) = cmplx(0.0,0.0)
                      do n=1,abs(MGp%nbas(m)) !cannot be greater than 12 at present
                         x=real(MGp%icomp(n,m))
                         ci=cmplx(1.0-x,-isig*x)
                         Sk(:)=Sk(:)+ Atm%atom(i)%cbas(n,nvk)*ci* cmplx(Real(MGp%basf(:,n,k,m)), -isig*aimag(MGp%basf(:,n,k,m)))
                      end do
                      Sk(:) = cmplx (real(Sk),ch(ich)*aimag(Sk))
                      if(.not. Mag_Dom%twin) then
                        Sk(:)=matmul(real(Mag_Dom%Dmat(:,:,nd)),Sk(:))
                      end if
                      GMh(:)=GMh(:) + anis*Sk(:)*CMPLX(COS(arg),SIN(arg))  !Atomic contribution to geometric Magnetic Structure Factor
                   end do ! symmetry
                   Mh%MsF(:,ich,nd)=Mh%MsF(:,ich,nd) + tho*onh*GMh(:)
                end do ! Atoms
                Mh%MsF(:,ich,nd)=Mh%MsF(:,ich,nd) * MGp%Num_Lat * MGp%Centred
             end do ! Chirality Domains
          end do ! Domains
       end if

       !---- Calculation of the Magnetic Interaction vectors ----!
       s  = 2.0*Mh%s            !1/d=r*  M = M// + Mp   => Mp = M - M// = M - (M.e)e
       er = Mh%h/s              !unitary vector referred to the reciprocal basis
       ed = matmul(cell%GR,er)  !  "        "       "             direct    "
       Mh%AMiV(:) = cmplx(0.0,0.0)
       Mh%sqMiV = 0.0
       do nd=1,Mag_dom%nd
          do ich=1,nch
             Mc  = Mh%MsF(:,ich,nd) / Cell%cell    !Magnetic structure factor in basis {a,b,c}
             MiV = Mc - dot_product(er,Mc) * ed    !Magnetic interaction vector in basis {a,b,c}
             Mh%MiVC(:,ich,nd)  = matmul(Cell%Cr_Orth_cel,MiV)  !Magnetic interaction vector in Cartesian components
             Mh%MiV(:,ich,nd) =  MiV * Cell%cell  !Magnetic interaction vector in basis {e1,e2,e3}
             Mh%AMiV(:)=Mh%AMiV(:)+ Mh%MiVC(:,ich,nd) * Mag_Dom%Pop(ich,nd)
             Mh%sqMiV = Mh%sqMiV + dot_product(Mh%MiVC(:,ich,nd), Mh%MiVC(:,ich,nd))* Mag_Dom%Pop(ich,nd)
          end do ! Chirality Domains
       end do ! Domains
       Mh%sqAMiV= dot_product(Mh%AMiV, Mh%AMiV)

       return
    End Subroutine Calc_Magnetic_Strf_Miv_Dom

    !!----
    !!---- Subroutine Calc_Magnetic_Strf_Tensor(SpG,Atm,Mh)
    !!----    type(Space_Group_Type),   intent(in)     :: SpG
    !!----    type(Matom_list_type),    intent(in)     :: Atm
    !!----    type(MagH_Type),          intent(in out) :: Mh
    !!----
    !!----    Calculate the Tensorial Magnetic Structure factor of the
    !!----    reflection provided in Mh. Only reasonable settings for symmetry
    !!----    operators are allowed to get correct values in this subroutine.
    !!----    The components are given with respect to the crystallographic
    !!----    unitary direct cell system: {e1,e2,e3} and with respect to the
    !!----    Cartesian frame defined in Cell.
    !!----
    !!---- Created: June - 2014 (JRC)
    !!
    Subroutine Calc_Magnetic_Strf_Tensor(SpG,Atm,Mh)
       !---- Arguments ----!
       type(Space_Group_Type),   intent(in)     :: SpG
       type(Matom_list_type),    intent(in)     :: Atm
       type(MagH_Type),          intent(in out) :: Mh

       !---- Local Variables ----!
       integer                            :: i,j,k
       real(kind=cp)                      :: arg,anis,s,b,ht,mFF,tho
       real(kind=cp),    dimension(3)     :: h
       real(kind=cp),    dimension(6)     :: beta
       real(kind=cp),    dimension(3,3)   :: Mcos,Msin,chi,chit

       s=Mh%s
       Mh%TMsF=0.0

       do i=1,Atm%natoms
          !---- Isotropic Debye-Waller factor * occupation * p=0.5*re*gamma * Magnetic form-factors mFF
          b=atm%atom(i)%biso
          j=atm%atom(i)%ind(2)  !pointer to the magnetic form factor coefficients
          mFF=mfj(s,Magnetic_Form(j)%SctM)
          tho= pn*atm%atom(i)%occ*mFF*exp(-b*s*s)
          chi=reshape((/atm%atom(i)%chi(1),atm%atom(i)%chi(4), atm%atom(i)%chi(5), &
                        atm%atom(i)%chi(4),atm%atom(i)%chi(2), atm%atom(i)%chi(6), &
                        atm%atom(i)%chi(6),atm%atom(i)%chi(6), atm%atom(i)%chi(3) /),(/3,3/))
          Mcos=0.0 ; Msin=0.0
          do k=1,SpG%Numops
             h=Hkl_R(Mh%h,SpG%symop(k))
             ht=dot_product(Mh%h,SpG%SymOp(k)%Tr)
             arg=tpi*(dot_product(h,Atm%atom(i)%x)+ht)
             anis=1.0
             if (Atm%atom(i)%thtype == "aniso") then
                beta=Atm%atom(i)%u(1:6)
                anis=     h(1)*h(1)*beta(1)+     h(2)*h(2)*beta(2)+    h(3)*h(3)*beta(3) &
                     +2.0*h(1)*h(2)*beta(4)+ 2.0*h(1)*h(3)*beta(5)+2.0*h(2)*h(3)*beta(6)
                anis=exp(-anis)
             end if
             chit=matmul(SpG%SymOp(k)%Rot,chi)
             chit=matmul(chit,transpose(SpG%SymOp(k)%Rot))
             Mcos(:,:)=Mcos(:,:)+cos(arg)*anis*chit
             if(SpG%Centred /= 2) Msin(:,:)=Msin(:,:)+sin(arg)*anis*chit
             !write(*,"(a,10f10.4)") "  arg, Mcos_part: ",arg,cos(arg)*anis*chit
          end do ! symmetry
          Mh%TMsF=Mh%TMsF+tho*cmplx(Mcos,Msin)
       end do ! Atoms
       Mh%TMsF=Mh%TMsF*SpG%Centred*SpG%NumLat
       return
    End Subroutine Calc_Magnetic_StrF_Tensor

    !!--++
    !!--++ Subroutine Calc_Table_MAB(Cell,Mlist,Atm,Mgp)
    !!--++    type(Crystal_Cell_type),intent(in) :: Cell
    !!--++    type(MagH_List_Type),   intent(in) :: MList
    !!--++    type(Matom_list_type),  intent(in) :: Atm
    !!--++    type(MagSymm_k_type),   intent(in) :: MGp
    !!--++
    !!--++    (Private)
    !!--++    Calculate Table with Aj(h) and Bj(h) values
    !!--++
    !!--++ Update: April - 2005, November 2014 (JRC)
    !!
    Subroutine Calc_Table_Mab(Cell,Mlist,Atm,Mgp)
       !---- Arguments ----!
       type(Crystal_Cell_type),intent(in) :: Cell
       type(MagH_List_Type),   intent(in) :: MList
       type(Matom_list_type),  intent(in) :: Atm
       type(MagSymm_k_type),   intent(in) :: MGp

       !---- Local Variables ----!
       integer                       :: i,j,k,n,nvk,m, NRef
       real(kind=cp)                 :: arg,anis,onh,ph, x, isig
       real(kind=cp),dimension(3)    :: h
       complex(kind=cp),dimension(3) :: Sk,GMh
       complex(kind=cp)              :: ci
       real(kind=cp),dimension(6)    :: beta
       real(kind=cp),dimension(3,3)  :: Mcos,Msin
       real(kind=cp), dimension(3)   :: ar,ai,br,bi

       Ajh=0.0
       Bjh=0.0
       Nref=MList%Nref
       if (MGp%Centred == 2) then
          do j=1,Nref
             if (.not. Mlist%Mh(j)%keqv_minus ) then
                onh=0.5
             else
                onh=1.0
             end if
             nvk= Mlist%Mh(j)%num_k
             isig=Mlist%Mh(j)%signp

             do i=1,Atm%natoms
                m= Atm%Atom(i)%imat(nvk)
                if(m == 0) cycle
                Mcos=0.0
                do k=1,MgP%NumOps
                   h=hr(k,j)%h
                   ph= isig * (Atm%atom(i)%Mphas(m) + MGp%MSymOp(k,m)%Phas)
                   arg=tpi*(dot_product(h,Atm%atom(i)%x)+ht(k,j) + ph)
                   anis=1.0
                   if(Atm%atom(i)%thtype == "aniso") then
                      beta=Atm%atom(i)%u(1:6)
                      anis=     h(1)*h(1)*beta(1)+     h(2)*h(2)*beta(2)+    h(3)*h(3)*beta(3) &
                           +2.0*h(1)*h(2)*beta(4)+ 2.0*h(1)*h(3)*beta(5)+2.0*h(2)*h(3)*beta(6)
                      anis=exp(-anis)
                   end if
                   Mcos(:,:)=Mcos(:,:)+cos(arg)*anis*MGp%MSymOp(k,m)%Rot(:,:)
                end do ! symmetry
                Ajh(:,i,j) = onh*matmul(Mcos,Atm%atom(i)%SkR(:,nvk)/Cell%cell)*REAL(MGp%mcentred,kind=cp)*Cell%cell
             end do ! Atoms
          end do ! Reflections

       else

          if (MGp%nirreps == 0) then
             do j=1,Nref
                if (.not. Mlist%Mh(j)%keqv_minus ) then
                   onh=0.5
                else
                   onh=1.0
                end if
                isig=Mlist%Mh(j)%signp
                nvk= Mlist%Mh(j)%num_k

                do i=1,Atm%natoms
                   m= Atm%Atom(i)%imat(nvk)
                   if (m == 0) cycle

                   Mcos=0.0 ; Msin=0.0
                   do k=1,MGp%NumOps
                      h=hr(k,j)%h
                      ph = isig * (Atm%atom(i)%Mphas(m) + MGp%MSymOp(k,m)%Phas)
                      arg=tpi*(dot_product(h,Atm%atom(i)%x)+ht(k,j) + ph)
                      anis=1.0
                      if (Atm%atom(i)%thtype == "aniso") then
                         beta=Atm%atom(i)%u(1:6)
                         anis=     h(1)*h(1)*beta(1)+     h(2)*h(2)*beta(2)+    h(3)*h(3)*beta(3) &
                              +2.0*h(1)*h(2)*beta(4)+ 2.0*h(1)*h(3)*beta(5)+2.0*h(2)*h(3)*beta(6)
                         anis=exp(-anis)
                      end if
                      Mcos(:,:)=Mcos(:,:)+cos(arg)*anis*MGp%MSymOp(k,m)%Rot(:,:)
                      Msin(:,:)=Msin(:,:)+sin(arg)*anis*MGp%MSymOp(k,m)%Rot(:,:)
                   end do ! symmetry
                   ar =       onh*matmul(Mcos,Atm%atom(i)%SkR(:,nvk)/Cell%cell)*Cell%cell
                   ai = -isig*onh*matmul(Mcos,Atm%atom(i)%SkI(:,nvk)/Cell%cell)*Cell%cell
                   br =       onh*matmul(Msin,Atm%atom(i)%SkR(:,nvk)/Cell%cell)*Cell%cell
                   bi = -isig*onh*matmul(Msin,Atm%atom(i)%SkI(:,nvk)/Cell%cell)*Cell%cell

                   Ajh(:,i,j) = ar(:) - bi(:)
                   Bjh(:,i,j) = br(:) + ai(:)
                end do ! Atoms
             end do ! Reflections

          else   !now magnetic structure described in terms of basis functions

             do j=1,Nref
                if (.not. Mlist%Mh(j)%keqv_minus ) then
                   onh=0.5
                else
                   onh=1.0
                end if
                nvk= Mlist%Mh(j)%num_k
                isig=Mlist%Mh(j)%signp

                do i=1,Atm%natoms
                   m= Atm%Atom(i)%imat(nvk)
                   if(m == 0) cycle
                   GMh=cmplx(0.0,0.0)
                   do k=1,MGp%NumOps
                      h=hr(k,j)%h
                      ph= isig* Atm%atom(i)%Mphas(m)
                      arg=tpi*(dot_product(h,Atm%atom(i)%x)+ht(k,j) + ph)
                      anis=1.0
                      if (Atm%atom(i)%thtype == "aniso") then
                         beta=Atm%atom(i)%u(1:6)
                         anis=     h(1)*h(1)*beta(1)+     h(2)*h(2)*beta(2)+    h(3)*h(3)*beta(3) &
                              +2.0*h(1)*h(2)*beta(4)+ 2.0*h(1)*h(3)*beta(5)+2.0*h(2)*h(3)*beta(6)
                         anis=exp(-anis)
                      end if
                      Sk(:) = cmplx(0.0,0.0)
                      do n=1,abs(MGp%nbas(m)) !cannot be greater than 12 at present
                         x=real(MGp%icomp(n,m),kind=cp)
                         ci=cmplx(1.0-x,-isig*x)
                         Sk(:)=Sk(:)+ Atm%atom(i)%cbas(n,nvk)*ci*cmplx(  Real(MGp%basf(:,n,k,m)), -isig*aimag(MGp%basf(:,n,k,m))  )
                      end do
                      GMh(:)=GMh(:) + onh*anis*Sk(:)*CMPLX(COS(arg),SIN(arg))  !Atomic contribution to geometric Magnetic Structure Factor
                   end do ! symmetry

                   Ajh(:,i,j) = real(GMh)
                   Bjh(:,i,j) = aimag(GMh)

                end do ! Atoms
             end do ! Reflections
          end if
       end if

       return
    End Subroutine Calc_Table_Mab

    !!--++
    !!--++ Subroutine Calc_Table_Th(Reflex,Atm)
    !!--++    type(MagH_List_Type),   intent(in) :: Reflex
    !!--++    type(Matom_list_type),  intent(in) :: Atm
    !!--++
    !!--++    (Private)
    !!--++    Calculate the Table of Isotropic Thermal contribution and occupation
    !!--..    TH(Natoms,Nref) multiplied by pn= 0.2695420113693928312
    !!--++
    !!--++ Update: April - 2005
    !!
    Subroutine Calc_Table_Th(Reflex,Atm)
       !---- Argument ----!
       type(MagH_List_Type),   intent(in) :: Reflex
       type(Matom_list_type),  intent(in) :: Atm

       !---- Local variables ----!
       integer          :: i,j
       real(kind=cp)    :: b,s

       !---- Isotropic Debye-Waller factor * occupation * p=0.5*re*gamma ----!
       do j=1,reflex%nref
          s=reflex%Mh(j)%s
          do i=1,atm%natoms
             b=atm%atom(i)%biso
             th(i,j)= pn*atm%atom(i)%occ*exp(-b*s*s)
          end do
       end do

       return
    End Subroutine Calc_Table_Th

    !!--++
    !!--++ Subroutine Create_Table_Hr_Ht(Reflex,Grp)
    !!--++    type(MagH_List_Type),   intent(in) :: Reflex
    !!--++    type(MagSymm_k_type),   intent(in) :: Grp
    !!--++
    !!--++    (Private)
    !!--++    Calculate a Table with HR and HT values
    !!--..       Hr(Grp%Numops,Reflex%Nref)
    !!--..       HT(Grp%Numops,Reflex%Nref)
    !!--++
    !!--++ Update: April - 2005
    !!
    Subroutine Create_Table_Hr_Ht(Reflex,Grp)
       !---- Arguments ----!
       type(MagH_List_Type),   intent(in) :: Reflex
       type(MagSymm_k_type),   intent(in) :: Grp

       !---- Local Variables ----!
       integer :: i,j

       do j=1,reflex%nref
          do i=1,grp%NumOps
             hr(i,j)%h=Hkl_R(reflex%Mh(j)%h,Grp%symop(i))
             ht(i,j)=dot_product(real(reflex%Mh(j)%h,kind=cp),Grp%SymOp(i)%Tr)
          end do
       end do

       return
    End Subroutine Create_Table_HR_HT

    !!--++
    !!--++ Subroutine Create_Table_mFR(Reflex,Atm,Lun)
    !!--++    type(MagH_List_Type),   intent(in) :: Reflex
    !!--++    type(Matom_list_type),  intent(in) :: Atm
    !!--++    integer, optional,      intent(in) :: lun
    !!--++
    !!--++    (Private)
    !!--++    Calculate a Table of Magnetic form Factors
    !!--..     mFR(Natoms,Nref)
    !!--++
    !!--++ Update: April - 2005
    !!
    Subroutine Create_Table_mFR(Reflex,Atm,lun)
       !---- Arguments ----!
       type(MagH_List_Type),   intent(in) :: Reflex
       type(Matom_list_type),  intent(in) :: Atm
       integer, optional,      intent(in) :: lun

       !---- Local Variables ----!
       character(len=4)               :: symbcar
       integer                        :: i,j, k,n,L
       integer, dimension(atm%natoms) :: ix,jx,ia

       !---- Init ----!
       err_msfac=.false.

       !---- Load form factor values for Magnetic Scattering ----!
       call Set_Magnetic_Form()

       !---- Find Species in Magnetic_Form ----!
       ix=0
       jx=0
       n=0
       do i=1,atm%natoms
          symbcar=u_case(atm%atom(i)%SfacSymb)
          do j=1,num_mag_form
             if (symbcar /= Magnetic_Form(j)%Symb) cycle
             ix(i)=j
             if(any(jx == j) ) exit
             n=n+1
             jx(n)=j
             ia(n)=i
             exit
          end do
       end do

       if (present(lun)) then
          write(unit=lun,fmt="(/,a)") "  INFORMATION FROM TABULATED MAGNETIC FORM FACTORS"
          write(unit=lun,fmt="(a,/)") "  ================================================"
       End if

       if (any(ix==0)) then
          err_msfac=.true.
          err_msfac_mess="The Species "//symbcar//" was not found"
          return
       else
          !---- Fill mFR Table ----!
          do j=1,reflex%nref
             do i=1,atm%natoms
                mFR(i,j)=mfj(reflex%Mh(j)%s,Magnetic_Form(ix(i))%SctM)
             end do
          end do
       end if

       !---- Printing Information ----!
       if (present(lun)) then
          write(unit=lun,fmt="(/,a,/)")    "   MAGNETIC FORM FACTOR COEFFICIENTS: {A(i),B(i),I=1,3},C  "
          write(unit=lun,fmt="(a,i3)")     "   Number of chemically different species: ",n
          write(unit=lun,fmt="(/,a)") &
               "     Atom     a1       b1       a2       b2       a3       b3       c      "
          do k=1,n
             j = jx(k)
             i = ia(k)
             write(unit=lun,fmt="(a,7F9.5)")    &
                           "     "//atm%atom(i)%SfacSymb, (Magnetic_Form(j)%SctM(L), L=1,7)
          end do
          write(unit=lun,fmt="(/,/)")
       end if

       call Remove_Magnetic_Form()

       return
    End Subroutine Create_Table_mFR

    !!----
    !!---- Subroutine Gen_Satellites(Cell,Grp,Smax,H,Ord,Powder,hkl)
    !!----    type(Crystal_Cell_type),               intent(in)     :: cell
    !!----    type(MagSymm_k_Type),                  intent(in)     :: Grp
    !!----    real(kind=cp),                         intent(in)     :: smax
    !!----    type(MagH_List_Type),                  intent(in out) :: H
    !!----    logical, optional,                     intent(in)     :: ord
    !!----    logical, optional,                     intent(in)     :: powder
    !!----    type (Reflection_List_Type), optional, intent(in)     :: hkl
    !!----
    !!----    Generates half reciprocal sphere of integer reflections and
    !!----    add satellites according to the information given in Grp.
    !!----    If Ord is given the reflections are reordered
    !!----    by increasing sinTheta/Lambda.
    !!----    If Powder is given, the unique reflections for a powder pattern
    !!----    are generated. The extinctions are obtained from a calculation of
    !!----    a random magnetic structure respecting the symmetry provided in Grp.
    !!----    If hkl is provided, the call to HKL_GEN is avoided.
    !!----    The subroutine constructs partially the object H.
    !!----
    !!---- Created:    April - 2005
    !!---- Updated: December - 2011, November 2014 (JRC)
    !!
    Subroutine Gen_Satellites(Cell,Grp,Smax,H,Ord,Powder,hkl)
       !---- Arguments ----!
       type(Crystal_Cell_type),               intent(in)     :: cell
       type(MagSymm_k_Type),                  intent(in)     :: Grp
       real(kind=cp),                         intent(in)     :: smax
       type(MagH_List_Type),                  intent(in out) :: H
       logical, optional,                     intent(in)     :: ord
       logical, optional,                     intent(in)     :: powder
       type (Reflection_List_Type), optional, intent(in)     :: hkl

       !---- Local variables ----!
       integer                                       :: i,j,k, numref,num_ref,n, ng, lu, addk,nmat,ik
       real(kind=cp)                                 :: smgen, maxr, epv,s1,s2,sqmiv1,sqmiv2, epm
       type(Space_Group_Type)                        :: G
       type(Reflect_Type), dimension(:), allocatable :: reflex
       type(MagH_List_Type)                          :: Hloc
       character(len=20)                             :: Symb
       real(kind=cp),dimension(3)                    :: kv,hr,ks
       logical                                       :: keq
       integer,            dimension(:), allocatable :: ind
       logical,            dimension(:), allocatable :: treated
       type(Matom_list_type)                         :: Am   !useful only the the powder case
       character(len=*),parameter, dimension(8)      :: label=(/"Ho","Er","Gd","Dy","Mn","Fe","Co","Ni"/)
       character(len=*),parameter, dimension(8)      :: sfac= (/"JHO3","JER3","MGD3","JDY3","MMN3","MFE3","MCO2","MNI2"/)

       epv=0.0001
       epm=0.01
       Symb=Grp%latt//" -1"
       Call set_spacegroup(symb,G)

       ! Determine the higher reciprocal cell parameter to add it to the given smax
       maxr=maxval(Cell%rcell)
       smgen=smax+maxr  !generate reflections up to smgen
       if(present(hkl)) then
         numref=hkl%nref
       else
         numref= Get_MaxNumRef(smgen, Cell%CellVol)
       end if
       if(allocated(reflex)) deallocate(reflex)
       allocate(reflex(numref))
       if(present(hkl)) then
         num_ref=numref
         do i=1,num_ref
          reflex(i)%h=hkl%Ref(i)%h
          reflex(i)%s=hkl%Ref(i)%s
          reflex(i)%Mult=hkl%Ref(i)%Mult
         end do
       else
         call hkl_gen(Cell,G,.true.,0.0_cp,smgen,num_ref,reflex)
       end if

       !calculate the total number of satellites
       n=Grp%nkv
       ng=0
       do i=1,n
          kv=Grp%kvec(:,i)
          ng=ng+1
          if( .not. K_Equiv_Minus_K(kv,Grp%latt) ) ng=ng+1
       end do

       !Adding the 000 satellites is an additional 000 reflection, allocate also
       !the reflections of the form (-h,k,0) + kvect with h<=1 that have to be explicitly generated
       !when kvect is equivalent to -kvect
       addk = 2*(nint(smgen/Cell%rcell(2)) + 1)
       Hloc%nref= (num_ref+1+addk) * ng

       if(allocated(Hloc%mh)) deallocate(Hloc%mh)
       allocate(Hloc%mh(Hloc%nref))

       ng=0
       !Generate the 000 satellites, except for k=(000)
       do i=1,n
          kv=Grp%kvec(:,i)
          ng=ng+1
          keq = K_Equiv_Minus_K(kv,Grp%latt)
          hloc%Mh(ng)%keqv_minus= keq
          hloc%Mh(ng)%num_k= i
          hloc%Mh(ng)%signp=-1.0
          ks=kv
          hloc%Mh(ng)%s = hkl_s(ks,Cell)
          if (hloc%Mh(ng)%s < 0.0000001) then
             ng=ng-1
             cycle
          end if
          hloc%Mh(ng)%h=ks
       end do

       !rest of reflections
       outer: do j=1,num_ref
          hr=real(reflex(j)%h,kind=cp)
          do i=1,n
             kv=Grp%kvec(:,i)
             ng=ng+1
             if (ng > Hloc%nref) exit outer
             ks=hr+kv
             hloc%Mh(ng)%s = hkl_s(ks,Cell)
             if ( hloc%Mh(ng)%s > smax) then     !avoid reflection with s>smax
                ng=ng-1
                cycle
             end if
             keq = K_Equiv_Minus_K(kv,Grp%latt)
             hloc%Mh(ng)%keqv_minus= keq
             hloc%Mh(ng)%num_k= i
             hloc%Mh(ng)%signp=-1.0
             hloc%Mh(ng)%h=ks
             if ( .not. keq) then
                ng=ng+1
                if (ng > Hloc%nref) exit outer
                ks=hr-kv
                hloc%Mh(ng)%s = hkl_s(ks,Cell)
                if ( hloc%Mh(ng)%s > smax) then  !avoid reflection with s>smax
                   ng=ng-1
                   cycle
                end if
                hloc%Mh(ng)%keqv_minus= keq
                hloc%Mh(ng)%num_k= i
                hloc%Mh(ng)%signp= 1.0
                hloc%Mh(ng)%h=ks
             else if ( abs(hr(3)) < epv .and. abs(hr(1)) <= 1.0+epv) then   !complete with the reflections (-h,k,0)+kvec
                ng=ng+1
                if (ng > Hloc%nref) exit outer
                ks=(/-hr(1),-hr(2),0.0_cp/)+kv
                hloc%Mh(ng)%s = hkl_s(ks,Cell)
                if( hloc%Mh(ng)%s > smax) then     !avoid reflection with s>smax
                  ng=ng-1
                  cycle
                end if
                keq = K_Equiv_Minus_K(kv,Grp%latt)
                hloc%Mh(ng)%keqv_minus= keq
                hloc%Mh(ng)%num_k= i
                hloc%Mh(ng)%signp=-1.0
                hloc%Mh(ng)%h=ks
             end if
          end do
       end do outer

       Hloc%nref=ng  !update to the real number of effectively generated reflections

       if (present(powder)) then   !assumed ord=.true. even if not present

          !Generates a fictive magnetic structure with the provided symmetry to test the equivalent
          !reflections, get the real multiplicity and eliminate the systematic absences.
          if (allocated(ind)) deallocate(ind)
          allocate(ind(ng))
          if (allocated(treated)) deallocate(treated)
          allocate(treated(ng))

          treated=.false.
          call sort(hloc%Mh(:)%s,ng,ind)

          ! Re-order the local reflections
          call Get_LogUnit(lu)
          open(unit=lu,status="scratch",form="unformatted",action="readwrite")
          do i=1,ng
             j=ind(i)
             write(unit=lu) hloc%Mh(j)%keqv_minus,hloc%Mh(j)%num_k,hloc%Mh(j)%signp,hloc%Mh(j)%s,hloc%Mh(j)%h
          end do
          rewind(unit=lu)

          do j=1,ng
             read(unit=lu)  hloc%Mh(j)%keqv_minus,hloc%Mh(j)%num_k,hloc%Mh(j)%signp,hloc%Mh(j)%s,hloc%Mh(j)%h
             hloc%Mh(j)%mult=2
          end do
          close(unit=lu)

          !Calculate the magnetic structure factors for an arbitrary structure of the same symmetry as that
          !provided. First construct a magnetic atom list with Ho+3 form factors
          nmat = Grp%nmsym
          Call Allocate_mAtom_list(nmat,Am)
          call RANDOM_SEED()
          do i=1,nmat
             Am%atom(i)%lab=label(i)       !Label
             Am%atom(i)%SfacSymb=sfac(i)   !Formfactor label
             call random_number(ks)
             Am%atom(i)%x= ks              !Fract. coord.
             Am%atom(i)%Biso=0.3           !Is. Temp. Fact.
             Am%atom(i)%occ=1.0            !occupation
             call random_number(maxr)
             Am%atom(i)%nvk= max(1,nint(maxr*n))
             do ik=1,n
                call random_number(maxr)
                Am%atom(i)%imat(ik)= max(1,nint(maxr*n))
                call random_number(ks)
                Am%atom(i)%Skr(:,ik)= ks(:)*8.0
                Am%atom(i)%Ski(:,ik)= 0.0
                call random_number(maxr)
                Am%atom(i)%mphas(ik)= maxr
             end do
          end do

          call Mag_Structure_Factors(Cell,Am,Grp,hloc)
          call Calc_Mag_Interaction_vector(hloc,cell)
          maxr=maxval(hloc%Mh(:)%sqMiV)
          epm=maxr*0.00001

          !Lines for debugging
          call Get_LogUnit(lu)
          open(unit=lu,file="powder_test.sfa",status="replace",action="write")
          call Write_Magnetic_Structure(lu,Grp,Am)
          call Write_Mag_Structure_Factors(lu,hloc,Grp)
          !End Lines for debugging
          !Start analysis
          ind(:) = 0

          do i=1,ng
             if(treated(i)) cycle
             sqmiv1= hloc%Mh(i)%sqMiV
             ind(i)=2
             s1= hloc%Mh(i)%s
             treated(i) = .true.
             if(sqmiv1 < epm ) ind(i)=0
             do j=i+1,ng
                s2= hloc%Mh(j)%s
                if( abs(s1-s2) > epv) exit
                sqmiv2= hloc%Mh(j)%sqMiV
                if( abs(sqmiv1-sqmiv2) > epm) exit
                !Passing here give an equivalent reflection
                ind(j) = 0
                if(sqmiv2 > epm) ind(i)=ind(i)+2
                treated(j)=.true.
             end do
          end do

          ! Determine the number of independent reflections
          k=0
          do i=1,ng
             if(ind(i) /= 0) k=k+1
          end do

          !Allocate the strictly needed magnetic reflections
          H%nref=k
          if(allocated(H%mh)) deallocate(H%mh)
          allocate(H%mh(k))
          k=0
          do j=1,ng
             if (ind(j) /= 0) then
                k=k+1
                h%Mh(k)%keqv_minus = hloc%Mh(j)%keqv_minus
                h%Mh(k)%num_k      = hloc%Mh(j)%num_k
                h%Mh(k)%signp      = hloc%Mh(j)%signp
                h%Mh(k)%s          = hloc%Mh(j)%s
                h%Mh(k)%h          = hloc%Mh(j)%h
                h%Mh(k)%mult       = ind(j)
                h%Mh(k)%sqMiV      = hloc%Mh(j)%sqMiV
                h%Mh(k)%MiV        = cmplx(0.0,0.0)
                h%Mh(k)%MsF        = cmplx(0.0,0.0)
             end if
          end do
          call Write_Mag_Structure_Factors(lu,h,Grp)

       else  !  present(powder)=.false.

          H%nref=ng
          if(allocated(H%mh)) deallocate(H%mh)
          allocate(H%mh(ng))

          if (present(ord)) then
             if (allocated(ind)) deallocate(ind)
             allocate(ind(ng))

             if (ord) then !Reordering reflections by increasing sinTheta/Lambda
                call sort(hloc%Mh(:)%s,ng,ind)
             else
                do j=1,ng
                   ind(j)=j
                end do
             end if

             do j=1,ng
                i=ind(j)
                h%Mh(j)%keqv_minus = hloc%Mh(i)%keqv_minus
                h%Mh(j)%num_k      = hloc%Mh(i)%num_k
                h%Mh(j)%signp      = hloc%Mh(i)%signp
                h%Mh(j)%s          = hloc%Mh(i)%s
                h%Mh(j)%h          = hloc%Mh(i)%h
                h%Mh(j)%mult       = 2
                h%Mh(j)%sqMiV      = 0.0
                h%Mh(j)%MiV        = cmplx(0.0,0.0)
                h%Mh(j)%MsF        = cmplx(0.0,0.0)
             end do

          else   ! present(ord)=.false.
             h%nref=ng-1
             do j=1,h%nref
                h%Mh(j)%keqv_minus = hloc%Mh(j)%keqv_minus
                h%Mh(j)%num_k      = hloc%Mh(j)%num_k
                h%Mh(j)%signp      = hloc%Mh(j)%signp
                h%Mh(j)%s          = hloc%Mh(j)%s
                h%Mh(j)%h          = hloc%Mh(j)%h
                h%Mh(j)%mult       = 2
                h%Mh(j)%sqMiV      = 0.0
                h%Mh(j)%MiV        = cmplx(0.0,0.0)
                h%Mh(j)%MsF        = cmplx(0.0,0.0)
             end do
          end if

       end if

       return
    End Subroutine Gen_Satellites

    !!----
    !!---- Subroutine Init_Err_MSfac()
    !!----
    !!----    Initialize the errors flags in this Module
    !!----
    !!---- Update: April - 2008
    !!
    Subroutine Init_Err_MSfac()

       err_msfac=.false.
       err_msfac_mess=" "

       return
    End Subroutine Init_Err_MSfac

    !!----
    !!---- Subroutine Init_Mag_Structure_Factors(Reflex,Atm,Grp,Lun)
    !!----    type(MagH_List_Type),           intent(in) :: Reflex
    !!----    type(Matom_list_type),          intent(in) :: Atm
    !!----    type(MagSymm_k_Type),           intent(in) :: Grp
    !!----    integer,              optional, intent(in) :: lun
    !!----
    !!----    Allocates and initializes arrays for Magnetic Structure Factors calculations.
    !!----    A calculation of fixed tables is also performed.
    !!----
    !!---- Update: April - 2005
    !!
    Subroutine Init_Mag_Structure_Factors(Reflex,Atm,Grp,lun)
       !---Arguments ---!
       type(MagH_List_Type),           intent(in) :: Reflex
       type(Matom_list_type),          intent(in) :: Atm
       type(MagSymm_k_Type),           intent(in) :: Grp
       integer,              optional, intent(in) :: lun

       !--- Local variables ---!
       integer :: Natm, Multr
       integer :: ierr

       call init_err_msfac()

       Natm = Atm%natoms
       Multr= Grp%Numops

       !---- Magnetic Scattering factor tables ----!
       if (allocated(mFR)) deallocate(mFR)
       allocate(mFR(Natm,Reflex%Nref),stat=ierr)
       if (ierr /=0) then
          err_msfac=.true.
          err_msfac_mess=" Unable to allocate mFR"
          return
       end if
       mFR=0.0

       !---- HR Table ----!
       if (allocated(HR)) deallocate(HR)
       allocate(HR(Multr,Reflex%Nref),stat=ierr)
       if (ierr /=0) then
          err_msfac=.true.
          err_msfac_mess=" Unable to allocate HR"
          return
       end if
       HR=HR_Type(0)

       !---- HT Table ----!
       if (allocated(HT)) deallocate(HT)
       allocate(HT(Multr,Reflex%Nref),stat=ierr)
       if (ierr /=0) then
          err_msfac=.true.
          err_msfac_mess=" Unable to allocate HTR"
          return
       end if
       HT=0.0

       if (allocated(TH)) deallocate(TH)
       allocate(TH(Natm,Reflex%Nref),stat=ierr)
       if (ierr /=0) then
          err_msfac=.true.
          err_msfac_mess=" Unable to allocate HTR"
          return
       end if
       TH=0.0

       if (allocated(Ajh)) deallocate(Ajh)
       allocate(Ajh(3,Natm,Reflex%Nref), stat=ierr)
       if (ierr /=0) then
          err_msfac=.true.
          err_msfac_mess=" Unable to allocate Aj(h)"
          return
       end if
       Ajh=0.0

       if (allocated(Bjh)) deallocate(Bjh)
       allocate(Bjh(3,Natm,Reflex%Nref), stat=ierr)
       if (ierr /=0) then
          err_msfac=.true.
          err_msfac_mess=" Unable to allocate Bj(h)"
          return
       end if
       Bjh=0.0

       if (present(lun)) then
          call Set_Fixed_Tables(Reflex,Atm,Grp,lun=lun)
       else
          call Set_Fixed_Tables(Reflex,Atm,Grp)
       end if

       if (.not. err_msfac) MSF_Initialized=.true.

       return
    End Subroutine Init_Mag_Structure_Factors

    !!----
    !!---- Subroutine Mag_Structure_Factors(Cell,Atm,Grp,Reflex)
    !!----    !---- Arguments ----!
    !!----    type(Crystal_Cell_type),  intent(in)     :: Cell
    !!----    type(Matom_list_type),    intent(in)     :: Atm
    !!----    type(MagSymm_k_Type),     intent(in)     :: Grp
    !!----    type(MagH_List_Type),     intent(in out) :: Reflex
    !!----
    !!----    Calculate the Magnetic Structure Factors from a list of magnetic Atoms
    !!----    and a set of reflections. A call to Init_Mag_Structure_Factors
    !!----    is a pre-requisite for using this subroutine. In any case
    !!----    the subroutine calls Init_Mag_Structure_Factors if SF_initialized=.false.
    !!----    The argument "Cell" has been added in order to consider whatever kind of
    !!----    settings for symmetry operators.
    !!----
    !!---- Update: April - 2005, November 2014 (JRC)
    !!
    Subroutine Mag_Structure_Factors(Cell,Atm,Grp,Reflex)
       !---- Arguments ----!
       type(Crystal_Cell_type),  intent(in)     :: Cell
       type(Matom_list_type),    intent(in)     :: Atm
       type(MagSymm_k_Type),     intent(in)     :: Grp
       type(MagH_List_Type),     intent(in out) :: Reflex


       call init_err_msfac()
       if(.not. MSF_Initialized) call Init_Mag_Structure_Factors(Reflex,Atm,Grp)

       !---- Table TH ----!
       Call Calc_Table_TH(Reflex,Atm)

       !---- Table AB ----!
       call Calc_Table_MAB(Cell,Reflex,Atm,Grp)

       !---- Final Calculation ----!
       call Sum_MAB(Reflex,Atm%Natoms,Grp%Centred)

       return
    End Subroutine Mag_Structure_Factors

    !!----
    !!---- Subroutine Modify_MSF(Reflex,Atm,Grp,List,Nlist)
    !!----    !---- Arguments ----!
    !!----    type(MagH_List_Type),         intent(in) :: Reflex
    !!----    type(Matom_list_type),        intent(in) :: Atm
    !!----    type(MagSymm_k_Type),         intent(in) :: Grp
    !!----    integer,dimension(:),         intent(in) :: List
    !!----    integer,                      intent(in) :: NList
    !!----
    !!----    Recalculation of Magnetic Structure Factors because a
    !!----    list of Atoms parameters were modified. The "List" variable
    !!----    contains the numbers in the list of the atoms to be changed.
    !!----
    !!---- Update: April - 2005
    !!
    Subroutine Modify_MSF(Reflex,Atm,Grp,List,Nlist)
       !---- Arguments ----!
       type(MagH_List_Type),         intent(in out) :: Reflex
       type(Matom_list_type),        intent(in)     :: Atm
       type(MagSymm_k_Type),         intent(in)     :: Grp
       integer,dimension(:),         intent(in)     :: List
       integer,                      intent(in)     :: NList

       !---- Local variables ----!
       integer                       :: i,j,k,ii,nvk,m, n
       real(kind=cp)                 :: arg,onh,anis,ph, isig, x
       real(kind=cp),dimension(3,3)  :: Mcos,Msin
       complex(kind=cp),dimension(3) :: Sk,GMh
       complex(kind=cp)              :: ci
       real(kind=cp),dimension(3)    :: ar,ai,br,bi,h
       real(kind=cp),dimension(6)    :: beta

       if (Grp%Centred == 2) then

          do j=1,Reflex%Nref
            if (.not. Reflex%Mh(j)%keqv_minus ) then
               onh=0.5
            else
               onh=1.0
            end if
            nvk= Reflex%Mh(j)%num_k
            isig=Reflex%Mh(j)%signp
             do ii=1,Nlist
                i=list(ii)
                m= Atm%Atom(i)%imat(nvk)
                arg=0.0
                Mcos=0.0
                do k=1,grp%NumOps
                   h=hr(k,j)%h
                   ph= isig*(Atm%atom(i)%Mphas(m) + Grp%MSymOp(k,m)%Phas)
                   arg=tpi*(dot_product(h,Atm%atom(i)%x)+ht(k,j) + ph)
                   anis=1.0
                   if (Atm%atom(i)%thtype == "aniso") then
                      beta=Atm%atom(i)%u(1:6)
                      anis=     h(1)*h(1)*beta(1)+     h(2)*h(2)*beta(2)+    h(3)*h(3)*beta(3) &
                           +2.0*h(1)*h(2)*beta(4)+ 2.0*h(1)*h(3)*beta(5)+2.0*h(2)*h(3)*beta(6)
                      anis=exp(-anis)
                   end if
                   Mcos(:,:)=Mcos(:,:)+cos(arg)*anis*Grp%MSymOp(k,m)%Rot(:,:)
                end do ! symmetry
                ar =  onh*matmul(Mcos,Atm%atom(i)%SkR(:,nvk))
                Ajh(:,i,j)= ar(:)
             end do ! NList
          end do ! Reflections

       else

          if (Grp%nirreps == 0) then
             do j=1,Reflex%Nref
                if (.not. Reflex%Mh(j)%keqv_minus ) then
                   onh=0.5
                else
                   onh=1.0
                end if
                nvk= Reflex%Mh(j)%num_k
                isig=Reflex%Mh(j)%signp
                do ii=1,Nlist
                   i=list(ii)
                   m= Atm%Atom(i)%imat(nvk)
                   arg=0.0
                   Mcos=0.0 ; Msin=0.0
                   do k=1,grp%NumOps
                      h=hr(k,j)%h
                      ph=isig*(Atm%atom(i)%Mphas(m) + Grp%MSymOp(k,m)%Phas)
                      arg=tpi*(dot_product(h,Atm%atom(i)%x)+ht(k,j) + ph)
                      anis=1.0
                      if (Atm%atom(i)%thtype == "aniso") then
                         beta=Atm%atom(i)%u(1:6)
                         anis=     h(1)*h(1)*beta(1)+     h(2)*h(2)*beta(2)+    h(3)*h(3)*beta(3) &
                              +2.0*h(1)*h(2)*beta(4)+ 2.0*h(1)*h(3)*beta(5)+2.0*h(2)*h(3)*beta(6)
                         anis=exp(-anis)
                      end if
                      Mcos(:,:)=Mcos(:,:)+cos(arg)*anis*Grp%MSymOp(k,m)%Rot(:,:)
                      Msin(:,:)=Msin(:,:)+sin(arg)*anis*Grp%MSymOp(k,m)%Rot(:,:)
                   end do ! symmetry
                   ar =       onh*matmul(Mcos,Atm%atom(i)%SkR(:,nvk))
                   ai = -isig*onh*matmul(Mcos,Atm%atom(i)%SkI(:,nvk))
                   br =       onh*matmul(Msin,Atm%atom(i)%SkR(:,nvk))
                   bi = -isig*onh*matmul(Msin,Atm%atom(i)%SkI(:,nvk))
                   Ajh(:,i,j)= ar(:) - bi(:)
                   Bjh(:,i,j)= br(:) + ai(:)
                end do ! NList
             end do ! Reflections

          else

             do j=1,Reflex%Nref
                if (.not. Reflex%Mh(j)%keqv_minus ) then
                   onh=0.5
                else
                   onh=1.0
                end if
                nvk= Reflex%Mh(j)%num_k
                isig=Reflex%Mh(j)%signp

                do ii=1,Nlist
                   i=List(ii)
                   m= Atm%Atom(i)%imat(nvk)
                   if(m == 0) cycle
                   GMh=cmplx(0.0,0.0)
                   do k=1,Grp%NumOps
                      h=hr(k,j)%h
                      ph= isig* Atm%atom(i)%Mphas(m)
                      arg=tpi*(dot_product(h,Atm%atom(i)%x)+ht(k,j) + ph)
                      anis=1.0
                      if (Atm%atom(i)%thtype == "aniso") then
                         beta=Atm%atom(i)%u(1:6)
                         anis=     h(1)*h(1)*beta(1)+     h(2)*h(2)*beta(2)+    h(3)*h(3)*beta(3) &
                              +2.0*h(1)*h(2)*beta(4)+ 2.0*h(1)*h(3)*beta(5)+2.0*h(2)*h(3)*beta(6)
                         anis=exp(-anis)
                      end if
                      Sk(:) = cmplx(0.0,0.0)
                      do n=1,abs(Grp%nbas(m)) !cannot be greater than 12 at present
                         x=real(Grp%icomp(n,m),kind=cp)
                         ci=cmplx(1.0-x,-isig*x)
                         Sk=Sk+Atm%atom(i)%cbas(n,nvk)*ci*cmplx(Real(Grp%basf(:,n,k,m),kind=cp),-isig*aimag(Grp%basf(:,n,k,m)))
                      end do
                      GMh(:)=GMh(:) + onh*anis*Sk(:)*CMPLX(COS(arg),SIN(arg))  !Atomic contribution to geometric Magnetic Structure Factor
                   end do ! symmetry
                   Ajh(:,i,j) = real(GMh)
                   Bjh(:,i,j) = aimag(GMh)
                end do ! Nlist
             end do ! Reflections
          end if
       end if

       !---- Recalculation of MSF ----!
       call Sum_MAB(Reflex,Atm%Natoms,Grp%Centred)

       return
    End Subroutine Modify_MSF

    !!--++
    !!--++ Subroutine Set_Fixed_Tables(Reflex,Atm,Grp,lun)
    !!--++    type(MagH_List_Type),         intent(in) :: Reflex
    !!--++    type(Matom_list_type),        intent(in) :: Atm
    !!--++    type(MagSymm_k_Type),         intent(in) :: Grp
    !!--++    integer, optional,            intent(in) :: lun
    !!--++
    !!--++    (Private)
    !!--++    Calculates arrays that are fixed during all further
    !!--++    calculations
    !!--++
    !!--++ Update: April - 2005
    !!
    Subroutine Set_Fixed_Tables(Reflex,Atm,Grp,lun)
       !---- Arguments ----!
       type(MagH_List_Type),         intent(in) :: Reflex
       type(Matom_list_type),        intent(in) :: Atm
       type(MagSymm_k_Type),         intent(in) :: Grp
       integer, optional,            intent(in) :: lun

       !---- Local variables ----!

       !---- Table HR - HT ----!
       call Create_Table_HR_HT(Reflex,Grp)

       !---- Table mFR ----!
       if (present(lun)) then
          call Create_Table_mFR(Reflex,Atm,lun=lun)
       else
          call Create_Table_mFR(Reflex,Atm)
       end if

       !---- Modify the scattering factor tables to include the
       !---- multipliers factors concerning centre of symmetry and
       !---- centred translations
       if (Grp%mCentred == 2) mFR=2.0*mFR
       if (Grp%Num_Lat  > 1)  mFR=Grp%Num_Lat*mFR

       return
    End Subroutine Set_Fixed_Tables

    !!--++
    !!--++ Subroutine Sum_MAB(Reflex,Natm,icent)
    !!--++    type(MagH_List_Type),   intent(in out) :: Reflex
    !!--++    integer,                intent(in)     :: Natm
    !!--++    integer,                intent(in)     :: icent
    !!--++
    !!--++    (Private)
    !!--++    Calculate the Final Sum for Structure Factors calculations
    !!--++
    !!--++ Update: April - 2005
    !!
    Subroutine Sum_MAB(Reflex,Natm,icent)
       !---- Arguments ----!
       type(MagH_List_Type), intent(in out)  :: Reflex
       integer,              intent(in)      :: Natm
       integer,              intent(in)      :: icent

       !---- Local Variables ----!
       integer                      :: i,j
       real(kind=cp), dimension(3)  :: aa,bb


       !---- Fj(h)*Aj(h) ----!
       if (icent == 2) then    !Calculation for centrosymmetric structures
          do j=1,reflex%nref
             aa=0.0
             do i=1,Natm
                aa(:)= aa(:) + mFR(i,j)*th(i,j)*ajh(:,i,j)
             end do
             Reflex%Mh(j)%MsF(:)=cmplx(aa(:))
          end do

       else       !Calculation for non-centrosymmetric structures
          !---- Final Sum ----!
          do j=1,reflex%Nref
             aa=0.0
             bb=0.0
             do i=1,Natm
                aa(:)= aa(:) + mFR(i,j)*th(i,j)*ajh(:,i,j)
                bb(:)= bb(:) + mFR(i,j)*th(i,j)*bjh(:,i,j)
             end do
             Reflex%Mh(j)%MsF(:)=cmplx(aa(:),bb(:))
          end do
       end if

       return
    End Subroutine Sum_MAB

    !!----
    !!---- Subroutine Write_Structure_Factors(Lun,Reflex,Grp)
    !!----    integer,               intent(in) :: Lun
    !!----    type(MagH_List_Type),  intent(in) :: Reflex
    !!----    type(MagSymm_k_Type),  intent(in) :: Grp
    !!----
    !!----    Writes in logical unit=lun the list of structure factors
    !!----    contained in the array hkl
    !!----
    !!---- Update: April - 2005
    !!
    Subroutine Write_Mag_Structure_Factors(Lun,Reflex,Grp)
       !---- Argument ----!
       integer,               intent(in) :: lun
       type(MagH_List_Type),  intent(in) :: Reflex
       type(MagSymm_k_Type),  intent(in) :: Grp

       !---- Local Variables ----!
       integer                     :: i,nv, sig, mul
       real(kind=cp)               :: sqMiV, dspc
       integer,dimension(3)        :: h
       real(kind=cp), dimension(3) :: vk,hr

       write(unit=lun,fmt="(/,a)")   "    LIST OF REFLECTIONS AND MAGNETIC STRUCTURE FACTORS"
       write(unit=lun,fmt="(a,/)")   "    =================================================="
       write(unit=lun,fmt="(a,i6)") " => Total number of reflections  : ",reflex%Nref
       write(unit=lun,fmt="(a,i6)") " => Number of propagation vectors: ",Grp%nkv
       do i=1,Grp%nkv
          write(unit=lun,fmt="(a,i2,a,3f8.4,a)") " => Propagation vectors #",i," = (",Grp%kvec(:,i)," )"
       end do

       write(unit=lun,fmt="(/,a,/)") &
        "    Hr      Kr      Lr       H   K   L   nvk   Mult   dspc      |MiV|^2     Mrx      Mry      Mrz      "// &
        "Mix      Miy      Miz     MiVrx    MiVry    MiVrz    MiVix    MiViy    MiViz"
       do i=1,reflex%Nref
          hr=reflex%Mh(i)%h
          sig=reflex%Mh(i)%signp
          mul=reflex%Mh(i)%mult
          nv =reflex%Mh(i)%Num_k
          vk=Grp%kvec(:,nv)
          h=nint(hr+sig*vk)
          sqMiV= reflex%Mh(i)%sqMiV
          dspc=0.5/reflex%Mh(i)%S
          write(unit=lun,fmt="(3f8.3,tr2,3i4,i5,i6,f9.4,f13.5,12f9.4)") hr,h, -sig*nv, mul, dspc,sqMiV, &
               real(reflex%Mh(i)%MsF),aimag(reflex%Mh(i)%MsF), real(reflex%Mh(i)%MiV),aimag(reflex%Mh(i)%MiV)
       end do

       return
    End Subroutine Write_Mag_Structure_Factors

 End Module CFML_Magnetic_Structure_Factors
