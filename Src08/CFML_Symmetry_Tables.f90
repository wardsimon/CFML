!!-------------------------------------------------------
!!---- Crystallographic Fortran Modules Library (CrysFML)
!!-------------------------------------------------------
!!---- The CrysFML project is distributed under LGPL. In agreement with the
!!---- Intergovernmental Convention of the ILL, this software cannot be used
!!---- in military applications.
!!----
!!---- Copyright (C) 1999-2019  Institut Laue-Langevin (ILL), Grenoble, FRANCE
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
!!---- MODULE: CFML_Symmetry_Tables
!!----   INFO: Tabulated information on Crystallographic Symmetry
!!----
!!---- HISTORY
!!----    Update: 04/03/2011
!!----
!!
 Module CFML_Symmetry_Tables
    !---- Use modules ----!
    Use CFML_GlobalDeps
    Use CFML_Strings,    only: u_case, pack_string, string_count

    !---- Variables ----!
    implicit none

    private

    !---- List of public subroutines ----!
    public :: get_compact_hm_symbol, Get_HM_Symbol_from_Compact_HM, get_it_generators, &
              get_spaceG_Symbols 
    public :: set_shubnikov_info,set_spgr_info, set_system_equiv, set_wyckoff_info
    public :: remove_shubnikov_info, remove_spgr_info, remove_system_equiv, remove_wyckoff_info

    !---- List of private subroutines ----!
    private :: set_IT_gen
    
    !---- Types Definitions ----!
    
    !!----
    !!---- TYPE :: SHUB_SPGR_INFO_TYPE
    !!----
    !!
    Type, public :: Shub_Spgr_Info_Type
       character(len=8) :: ID_BNS =" "   !  ID Number of BNS     
       character(len=15) :: BNS    =" "   !  BNS symbol
       character(len=12) :: ID_OG  =" "   !  ID number of OG
       character(len=15) :: OG     =" "   !  OG symbol
       character(len=25) :: STD    =" "   !  Proposed standard symbol
       character(len=25) :: MHall  =" "   !  Magnetic Hall symbol
   End Type Shub_Spgr_Info_Type
    
    !!----
    !!---- TYPE :: SPGR_INFO_TYPE
    !!--..
    !!
    Type, public :: Spgr_Info_Type
       integer                 :: N         =0      ! Number of the Spacegroup  
       character (len=12)      :: HM        =" "    ! Hermann-Mauguin           
       character (len=16)      :: Hall      =" "    ! Hall                      
       integer                 :: Laue      =0      ! Laue Group                
       integer                 :: Pg        =0      ! Point group               
       integer, dimension(6)   :: Asu       =0      ! Asymmetric unit * 24      
       character (len= 5)      :: Inf_Extra =" "    ! Extra information         
    End Type Spgr_Info_Type
    
    !!----
    !!---- TYPE :: TABLE_EQUIV_TYPE
    !!--..
    !!
    Type, public :: Table_Equiv_Type
       character(len= 6)      :: SC=" "             ! Schoenflies
       character(len=17)      :: ML=" "             ! Miller & Love
       character(len=18)      :: KO=" "             ! Kovalev
       character(len=32)      :: BC=" "             ! Bradley & Cracknell
       character(len=18)      :: ZA=" "             ! Zak
    End Type Table_Equiv_Type
    
    !!----
    !!---- TYPE :: WYCK_INFO_TYPE
    !!--..
    !!
    Type, public :: Wyck_Info_Type
       character (len=12)               :: HM      =" "
       integer                          :: Norbit  =0
       character (len=15),dimension(26) :: Corbit  =" "
    End Type Wyck_Info_Type
    
    !---- Private Variables ----!
    character(len=120), private, dimension(230) :: it_spg_gen=" "          ! Vector containing Generators of the Space group

    !---- Public variables ----!
    Type(Shub_Spgr_Info_Type), allocatable, dimension(:), public :: Shubnikov_Info ! General Info about Shubnikov Magnetic grousp
    Type(Spgr_Info_Type),      allocatable, dimension(:), public :: Spgr_Info      ! General Info about Space Groups
    Type(Table_Equiv_Type),    allocatable, dimension(:), public :: System_Equiv
    Type(Wyck_Info_Type),      allocatable, dimension(:), public :: Wyckoff_Info
    
    !---- Parameters ----!
    integer,       parameter, public  :: CUBIC = 554            ! Cubic parameter index
    integer,       parameter, public  :: HEXAG = 527            ! Index parameter for hexagonal Groups
    integer,       parameter, public  :: MONOC =  15            ! Index parameter for Monoclinic Groups
    integer,       parameter, public  :: ORTHOR  = 163          ! Index parameter for Orthorhombic Groups
    integer,       parameter, public  :: TETRA = 410            ! Index parameter for Tetragonal Groups
    integer,       parameter, public  :: TRIGO = 495            ! Index parameter for Trigonal Groups
    integer,       parameter, public  :: NUM_SPGR_INFO = 612    ! Total dimension of SPGR_INFO
    integer,       parameter, public  :: NUM_SHUBNIKOV = 1651   ! Total dimension of SHUBNIKOV_INFO


    !> Lattice Traslations
    real(kind=cp), dimension(3,2), parameter, public :: LTR_A =reshape ( (/0.0,0.0,0.0, 0.0,0.5,0.5/), (/3,2/) )
    real(kind=cp), dimension(3,2), parameter, public :: LTR_B =reshape ( (/0.0,0.0,0.0, 0.5,0.0,0.5/), (/3,2/) )
    real(kind=cp), dimension(3,2), parameter, public :: LTR_C =reshape ( (/0.0,0.0,0.0, 0.5,0.5,0.0/), (/3,2/) )
    real(kind=cp), dimension(3,2), parameter, public :: LTR_I =reshape ( (/0.0,0.0,0.0, 0.5,0.5,0.5/), (/3,2/) )
    real(kind=cp), dimension(3,4), parameter, public :: &
                   LTR_F =reshape( (/0.0,0.0,0.0, 0.0,0.5,0.5, 0.5,0.0,0.5, 0.5,0.5,0.0 /),(/3,4/) )
    real(kind=cp), dimension(3,3), parameter, public :: &
                   LTR_R =reshape( (/0.0,0.0,0.0, 2.0/3.0,1.0/3.0,1.0/3.0,  1.0/3.0,2.0/3.0,2.0/3.0/),(/3,3/) )

    
    !> Bradley & Cracknell Notation
    character(len=*), dimension(24), parameter, public  :: BC_D6h =[            &  
       "  E  "," C+_3"," C-_3"," C_2 "," C-_6"," C+_6","C'_23","C'_21","C'_22", &
       "C`_23","C`_21","C`_22","  I  "," S-_6"," S+_6"," s_h "," S+_3"," S-_3", &
       " s_v3"," s_v1"," s_v2"," s_d3"," s_d1"," s_d2"]

    character(len=*), dimension(48), parameter, public :: BC_Oh =[              &
       "  E  "," C_2z"," C_2y"," C_2x","C+_31","C+_34","C+_33","C+_32","C-_31", &
       "C-_33","C-_32","C-_34"," C_2a"," C_2b","C-_4z","C+_4z","C-_4x"," C_2d", &
       " C_2f","C+_4x","C+_4y"," C_2c","C-_4y"," C_2e","  I  "," s_z "," s_y ", &
       " s_x ","S-_61","S-_64","S-_63","S-_62","S+_61","S+_63","S+_62","S+_64", &
       " s_da"," s_db","S+_4z","S-_4z","S+_4x"," s_dd"," s_df","S-_4x","S-_4y", &
       " s_dc","S+_4y"," s_de"]

    !> Magnetic Array
    character(len=*), dimension(72), parameter, public :: DEPMAT = [        &    
       "( Dx, Dy, Dz)      ","(-Dx,-Dy, Dz)      ","(-Dx, Dy,-Dz)      ",   &
       "( Dx,-Dy,-Dz)      ","( Dz, Dx, Dy)      ","( Dz,-Dx,-Dy)      ",   &
       "(-Dz,-Dx, Dy)      ","(-Dz, Dx,-Dy)      ","( Dy, Dz, Dx)      ",   &
       "(-Dy, Dz,-Dx)      ","( Dy,-Dz,-Dx)      ","(-Dy,-Dz, Dx)      ",   &
       "( Dy, Dx,-Dz)      ","(-Dy,-Dx,-Dz)      ","( Dy,-Dx, Dz)      ",   &
       "(-Dy, Dx, Dz)      ","( Dx, Dz,-Dy)      ","(-Dx, Dz, Dy)      ",   &
       "(-Dx,-Dz,-Dy)      ","( Dx,-Dz, Dy)      ","( Dz, Dy,-Dx)      ",   &
       "( Dz,-Dy, Dx)      ","(-Dz, Dy, Dx)      ","(-Dz,-Dy,-Dx)      ",   &
       "(-Dx,-Dy,-Dz)      ","( Dx, Dy,-Dz)      ","( Dx,-Dy, Dz)      ",   &
       "(-Dx, Dy, Dz)      ","(-Dz,-Dx,-Dy)      ","(-Dz, Dx, Dy)      ",   &
       "( Dz, Dx,-Dy)      ","( Dz,-Dx, Dy)      ","(-Dy,-Dz,-Dx)      ",   &
       "( Dy,-Dz, Dx)      ","(-Dy, Dz, Dx)      ","( Dy, Dz,-Dx)      ",   &
       "(-Dy,-Dx, Dz)      ","( Dy, Dx, Dz)      ","(-Dy, Dx,-Dz)      ",   &
       "( Dy,-Dx,-Dz)      ","(-Dx,-Dz, Dy)      ","( Dx,-Dz,-Dy)      ",   &
       "( Dx, Dz, Dy)      ","(-Dx, Dz,-Dy)      ","(-Dz,-Dy, Dx)      ",   &
       "(-Dz, Dy,-Dx)      ","( Dz,-Dy,-Dx)      ","( Dz, Dy, Dx)      ",   &
       "( Dx   ,    Dy, Dz)","(   -Dy, Dx-Dy, Dz)","(-Dx+Dy,-Dx   , Dz)",   &
       "(-Dx   ,   -Dy, Dz)","(    Dy,-Dx+Dy, Dz)","( Dx-Dy, Dx   , Dz)",   &
       "(    Dy, Dx   ,-Dz)","( Dx-Dy,   -Dy,-Dz)","(-Dx   ,-Dx+Dy,-Dz)",   &
       "(   -Dy,-Dx   ,-Dz)","(-Dx+Dy,    Dy,-Dz)","( Dx   , Dx-Dy,-Dz)",   &
       "(-Dx   ,   -Dy,-Dz)","(    Dy,-Dx+Dy,-Dz)","( Dx-Dy, Dx   ,-Dz)",   &
       "( Dx   ,    Dy,-Dz)","(   -Dy, Dx-Dy,-Dz)","(-Dx+Dy,-Dx   ,-Dz)",   &
       "(   -Dy,-Dx   , Dz)","(-Dx+Dy,    Dy, Dz)","( Dx   , Dx-Dy, Dz)",   &
       "(    Dy, Dx   , Dz)","( Dx-Dy,   -Dy, Dz)","(-Dx   ,-Dx+Dy, Dz)"]
       
    character(len=*), dimension(72), parameter, public :: MAGMAT = [        &
       "( Mx, My, Mz)      ","(-Mx,-My, Mz)      ","(-Mx, My,-Mz)      ",   &
       "( Mx,-My,-Mz)      ","( Mz, Mx, My)      ","( Mz,-Mx,-My)      ",   &
       "(-Mz,-Mx, My)      ","(-Mz, Mx,-My)      ","( My, Mz, Mx)      ",   &
       "(-My, Mz,-Mx)      ","( My,-Mz,-Mx)      ","(-My,-Mz, Mx)      ",   &
       "( My, Mx,-Mz)      ","(-My,-Mx,-Mz)      ","( My,-Mx, Mz)      ",   &
       "(-My, Mx, Mz)      ","( Mx, Mz,-My)      ","(-Mx, Mz, My)      ",   &
       "(-Mx,-Mz,-My)      ","( Mx,-Mz, My)      ","( Mz, My,-Mx)      ",   &
       "( Mz,-My, Mx)      ","(-Mz, My, Mx)      ","(-Mz,-My,-Mx)      ",   &
       "(-Mx,-My,-Mz)      ","( Mx, My,-Mz)      ","( Mx,-My, Mz)      ",   &
       "(-Mx, My, Mz)      ","(-Mz,-Mx,-My)      ","(-Mz, Mx, My)      ",   &
       "( Mz, Mx,-My)      ","( Mz,-Mx, My)      ","(-My,-Mz,-Mx)      ",   &
       "( My,-Mz, Mx)      ","(-My, Mz, Mx)      ","( My, Mz,-Mx)      ",   &
       "(-My,-Mx, Mz)      ","( My, Mx, Mz)      ","(-My, Mx,-Mz)      ",   &
       "( My,-Mx,-Mz)      ","(-Mx,-Mz, My)      ","( Mx,-Mz,-My)      ",   &
       "( Mx, Mz, My)      ","(-Mx, Mz,-My)      ","(-Mz,-My, Mx)      ",   &
       "(-Mz, My,-Mx)      ","( Mz,-My,-Mx)      ","( Mz, My, Mx)      ",   &
       "( Mx   ,    My, Mz)","(   -My, Mx-My, Mz)","(-Mx+My,-Mx   , Mz)",   &
       "(-Mx   ,   -My, Mz)","(    My,-Mx+My, Mz)","( Mx-My, Mx   , Mz)",   &
       "(    My, Mx   ,-Mz)","( Mx-My,   -My,-Mz)","(-Mx   ,-Mx+My,-Mz)",   &
       "(   -My,-Mx   ,-Mz)","(-Mx+My,    My,-Mz)","( Mx   , Mx-My,-Mz)",   &
       "(-Mx   ,   -My,-Mz)","(    My,-Mx+My,-Mz)","( Mx-My, Mx   ,-Mz)",   &
       "( Mx   ,    My,-Mz)","(   -My, Mx-My,-Mz)","(-Mx+My,-Mx   ,-Mz)",   &
       "(   -My,-Mx   , Mz)","(-Mx+My,    My, Mz)","( Mx   , Mx-My, Mz)",   &
       "(    My, Mx   , Mz)","( Mx-My,   -My, Mz)","(-Mx   ,-Mx+My, Mz)"]   


    !> International Symbols For Point Group Elements Of 6/mmm (D6h)
    character(len=*), dimension(24), parameter, public :: IntSymD6h =[       &
       "  1           "," 3+ ( 0, 0, z)"," 3- ( 0, 0, z)","  2 ( 0, 0, z)",  &
       " 6- ( 0, 0, z)"," 6+ ( 0, 0, z)","  2 ( x, x, 0)","  2 ( x, 0, 0)",  &
       "  2 ( 0, y, 0)","  2 ( x,-x, 0)","  2 ( x,2x, 0)","  2 (2x, x, 0)",  &
       " -1           ","-3+ ( 0, 0, z)","-3- ( 0, 0, z)","  m ( x, y, 0)",  &
       "-6- ( 0, 0, z)","-6+ ( 0, 0, z)","  m ( x,-x, z)","  m ( x,2x, z)",  &
       "  m (2x, x, z)","  m ( x, x, z)","  m ( x, 0, z)","  m ( 0, y, z)"]

    !> International Symbols For Point Group Elements Of M3M (Oh)
    character(len=*), dimension(48), parameter, public :: IntSymOh = [       &
       "  1           ","  2 ( 0, 0, z)","  2 ( 0, y, 0)","  2 ( x, 0, 0)",  &
       " 3+ ( x, x, x)"," 3+ (-x, x,-x)"," 3+ ( x,-x,-x)"," 3+ (-x,-x, x)",  &
       " 3- ( x, x, x)"," 3- ( x,-x,-x)"," 3- (-x,-x, x)"," 3- (-x, x,-x)",  &
       "  2 ( x, x, 0)","  2 ( x,-x, 0)"," 4- ( 0, 0, z)"," 4+ ( 0, 0, z)",  &
       " 4- ( x, 0, 0)","  2 ( 0, y, y)","  2 ( 0, y,-y)"," 4+ ( x, 0, 0)",  &
       " 4+ ( 0, y, 0)","  2 ( x, 0, x)"," 4- ( 0, y, 0)","  2 (-x, 0, x)",  &
       " -1           ","  m ( x, y, 0)","  m ( x, 0, z)","  m ( 0, y, z)",  &
       "-3+ ( x, x, x)","-3+ (-x, x,-x)","-3+ ( x,-x,-x)","-3+ (-x,-x, x)",  &
       "-3- ( x, x, x)","-3- ( x,-x,-x)","-3- (-x,-x, x)","-3- (-x, x,-x)",  &
       "  m ( x,-x, z)","  m ( x, x, z)","-4- ( 0, 0, z)","-4+ ( 0, 0, z)",  &
       "-4- ( x, 0, 0)","  m ( x, y,-y)","  m ( x, y, y)","-4+ ( x, 0, 0)",  &
       "-4+ ( 0, y, 0)","  m (-x, y, x)","-4- ( 0, y, 0)","  m ( x, y, x)"]

    !> Kovalev Notation
    character(len=*), dimension(24), parameter, public :: KOV_D6H=[        &
       " h1"," h3"," h5"," h4"," h6"," h2","h11"," h9"," h7"," h8","h12",  &
       "h10","h13","h15","h17","h16","h18","h14","h23",                    &
       "h21","h19","h20","h24","h22"]

    character(len=*), dimension(48), parameter, public :: KOV_OH=[                &
       " h1"," h4"," h3"," h2"," h9","h10","h12","h11"," h5"," h7"," h6"," h8",   &
       "h16","h13","h15","h14","h20","h18","h17","h19","h24","h23",               &
       "h22","h21","h25","h28","h27","h26","h33","h34","h36","h35",               &
       "h29","h31","h30","h32","h40","h37","h39","h38","h44","h42",               &
       "h41","h43","h48","h47","h46","h45"]

    !> Lattice Traslations
    character(len=*), dimension(8) , parameter, public  :: LATT =[     &
       "  P: { 000 }                                       ",          &
       "  A: { 000;  0  1/2 1/2 }+                         ",          &
       "  B: { 000; 1/2  0  1/2 }+                         ",          &
       "  C: { 000; 1/2 1/2  0  }+                         ",          &
       "  I: { 000; 1/2 1/2 1/2 }+                         ",          &
       "  R: { 000; 2/3 1/3 1/3; 1/3 2/3 2/3   }+          ",          &
       "  F: { 000;  0  1/2 1/2; 1/2  0  1/2; 1/2 1/2  0 }+",          &
       "  Z: { 000;  Unconventional Z-centering vectors  }+"]

    !> Laue symbols
    character(len=*), dimension(16), parameter, public :: LAUE_CLASS=[  &
       "-1   ","2/m  ","mmm  ","4/m  ","4/mmm","-3 R ","-3m R","-3   ", &
       "-3m1 ","-31m ","6/m  ","6/mmm","m-3  ","m-3m ","m3   ","m3m  "]

    !> Litvin points Operators
    character(len=*), dimension(48), parameter, public :: LITVIN_POINT_OP_LABEL=[               &
       "1       ","2x      ","2y      ","2z      ","3xyz-1  ","3xy-z   ","3-xyz   ","3x-yz   ", &
       "3xyz    ","3x-yz-1 ","3xy-z-1 ","3-xyz-1 ","2-xy    ","4z      ","4z-1    ","2xy     ", &
       "2-yz    ","2yz     ","4x      ","4x-1    ","2-xz    ","4y-1    ","2xz     ","4y      ", &
       "-1      ","mx      ","my      ","mz      ","-3xyz-1 ","-3xy-z  ","-3-xyz  ","-3x-yz  ", &
       "-3xyz   ","-3x-yz-1","-3xy-z-1","-3-xyz-1","m-xy    ","-4z     ","-4z-1   ","mxy     ", &
       "m-yz    ","myz     ","-4x     ","-4x-1   ","m-xz    ","-4y-1   ","mxz     ","-4y     "]

    character(len=*), dimension(48), parameter, public :: LITVIN_POINT_OP=[  &
       "x,y,z   ", "x,-y,-z ", "-x,y,-z ", "-x,-y,z ", "y,z,x   ",           &
       "y,-z,-x ", "-y,z,-x ", "-y,-z,x ", "z,x,y   ", "z,-x,-y ",           &
       "-z,x,-y ", "-z,-x,y ", "-y,-x,-z", "-y,x,z  ", "y,-x,z  ",           &
       "y,x,-z  ", "-x,-z,-y", "-x,z,y  ", "x,-z,y  ", "x,z,-y  ",           &
       "-z,-y,-x", "-z,y,x  ", "z,-y,x  ", "z,y,-x  ", "-x,-y,-z",           &
       "-x,y,z  ", "x,-y,z  ", "x,y,-z  ", "-y,-z,-x", "-y,z,x  ",           &
       "y,-z,x  ", "y,z,-x  ", "-z,-x,-y", "-z,x,y  ", "z,-x,y  ",           &
       "z,x,-y  ", "y,x,z   ", "y,-x,-z ", "-y,x,-z ", "-y,-x,z ",           &
       "x,z,y   ", "x,-z,-y ", "-x,z,-y ", "-x,-z,y ", "z,y,x   ",           &
       "z,-y,-x ", "-z,y,-x ", "-z,-y,x "]

    character(len=*), dimension(24), parameter, public :: LITVIN_POINT_OP_HEX_LABEL=[  &
       "1    ","6z   ","3z   ","2z   ","3z-1 ","6z-1 ","2x   ","21   ",                &
       "2xy  ","22   ","2y   ","23   ","-1   ","-6z  ","-3z  ","mz   ",                &
       "-3z-1","-6z-1","mx   ","m1   ","mxy  ","m2   ","my   ","m3   "]

    character(len=*), dimension(24), parameter, public :: LITVIN_POINT_OP_HEX=[       &
       "x,y,z     ","x-y,x,z   ","-y,x-y,z  ","-x,-y,z   ","-x+y,-x,z ","y,-x+y,z  ", &
       "x-y,-y,-z ","x,x-y,-z  ","y,x,-z    ","-x+y,y,-z ","-x,-x+y,-z","-y,-x,-z  ", &
       "-x,-y,-z  ","-x+y,-x,-z","y,-x+y,-z ","x,y,-z    ","x-y,x,-z  ","-y,x-y,-z ", &
       "-x+y,y,z  ","-x,-x+y,z ","-y,-x,z   ","x-y,-y,z  ","x,x-y,z   ","y,x,z     "]

    !> Miller & Love Notation
    character(len=*), dimension(24), parameter, public :: ML_D6H=[                &
       " 1"," 3"," 5"," 4"," 6"," 2"," 9"," 7","11","12","10"," 8","13","15","17",&
       "16","18","14","21","19","23","24","22","20"]

    character(len=*), dimension(48), parameter, public :: ML_OH=[                 &
       " 1"," 4"," 3"," 2"," 9","10","12","11"," 5"," 7"," 6"," 8","16","13","15",&
       "14","20","18","17","19","24","23","22","21","25","28","27","26","33","34",&
       "36","35","29","31","30","32","40","37","39","38","44","42","41","43","48",&
       "47","46","45"]

    !> Matrix Types For Rotational Operators In Conventional Basis: 1->24 Oh, 25->36 D6h
    Integer,  dimension(36,3,3), parameter, public :: MOD6 = reshape ( [     &
       1,-1,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,                   &
      -1, 1, 0, 0, 0, 0, 1, 0,-1,-1, 0, 1, 0, 1,-1, 0,-1, 1,                   &
       0, 0, 0, 0, 1,-1,-1, 1, 0, 0, 0, 0, 1,-1,-1, 1, 0, 0,                   &
       0, 0, 0, 0, 0, 0, 0, 1,-1, 0,-1, 1, 1, 0,-1,-1, 0, 1,                   &
       0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 1, 0, 0, 0, 0, 0, 0,                   &
       0, 0,-1, 1, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                   &
       0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 0, 0,                   &
       0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 1,-1, 1,-1, 0,-1, 1, 0,                   &
       1,-1, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                   &
       0, 0, 1,-1, 1,-1, 1,-1, 0,-1, 1, 0, 0,-1, 1, 0, 1,-1,                   &
       0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1,                   &
      -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                   &
       0, 0, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                   &
       0, 0, 1, 1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                   &
       0, 0, 0, 0, 0, 0, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0, 1, 1,                   &
      -1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                   &
       1, 1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1, 1, 1, 0, 0,                   &
       0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1], [36,3,3] )

    !> Point Group Symbols
    character(len=*), dimension(42), parameter, public :: POINT_GROUP=[   &
       "1    ","-1   ","2    ","m    ","2/m  ","222  ","mm2  ","m2m  ",   &
       "2mm  ","mmm  ","4    ","-4   ","4/m  ","422  ","4mm  ","-42m ",   &
       "-4m2 ","4/mmm","3    ","-3   ","32   ","3m   ","-3m  ","312  ",   &
       "31m  ","-31m ","6    ","-6   ","6/m  ","622  ","6mm  ","-62m ",   &
       "-6m2 ","6/mmm","23   ","m-3  ","432  ","-43m ","m-3m ","m3   ",   &
       "m3m  ","-3m1 "]
       
    !> System type
    character(len=*), dimension(7) , parameter, public:: SYS_CRY =[    &
       "Triclinic   ","Monoclinic  ","Orthorhombic","Tetragonal  ",    &
       "Trigonal    ","Hexagonal   ","Cubic       "]   

    !> X_D6H
    character(len=*), dimension(24), parameter, public   :: X_D6H = [        &
       "( x  ,   y, z)","(  -y, x-y, z)","(-x+y,-x  , z)","(-x  ,  -y, z)",  &
       "(   y,-x+y, z)","( x-y, x  , z)","(   y, x  ,-z)","( x-y,  -y,-z)",  &
       "(-x  ,-x+y,-z)","(  -y,-x  ,-z)","(-x+y,   y,-z)","( x  , x-y,-z)",  &
       "(-x  ,  -y,-z)","(   y,-x+y,-z)","( x-y, x  ,-z)","( x  ,   y,-z)",  &
       "(  -y, x-y,-z)","(-x+y,-x  ,-z)","(  -y,-x  , z)","(-x+y,   y, z)",  &
       "( x  , x-y, z)","(   y, x  , z)","( x-y,  -y, z)","(-x  ,-x+y, z)"]

    character(len=*), dimension(48), parameter, public  :: X_OH = [                   &
       "( x, y, z)","(-x,-y, z)","(-x, y,-z)","( x,-y,-z)","( z, x, y)","( z,-x,-y)", &
       "(-z,-x, y)","(-z, x,-y)","( y, z, x)","(-y, z,-x)","( y,-z,-x)","(-y,-z, x)", &
       "( y, x,-z)","(-y,-x,-z)","( y,-x, z)","(-y, x, z)","( x, z,-y)","(-x, z, y)", &
       "(-x,-z,-y)","( x,-z, y)","( z, y,-x)","( z,-y, x)","(-z, y, x)","(-z,-y,-x)", &
       "(-x,-y,-z)","( x, y,-z)","( x,-y, z)","(-x, y, z)","(-z,-x,-y)","(-z, x, y)", &
       "( z, x,-y)","( z,-x, y)","(-y,-z,-x)","( y,-z, x)","(-y, z, x)","( y, z,-x)", &
       "(-y,-x, z)","( y, x, z)","(-y, x,-z)","( y,-x,-z)","(-x,-z, y)","( x,-z,-y)", &
       "( x, z, y)","(-x, z,-y)","(-z,-y, x)","(-z, y,-x)","( z,-y,-x)","( z, y, x)"]

    !> Zak Notation
    character(len=*), dimension(24), parameter, public :: ZAK_D6H =[            &
       "   E   "," C(z)_3","C(2z)_3","  C_2  ","C(5z)_6"," C(z)_6","  U(xy)",   &
       "  U(x) ","  U(y) ","  U(3) ","  U(2) ","  U(1) ","   I   ","S(5z)_6",   &
       " S(z)_6","  s(z) "," S(z)_3","S(2z)_3"," s(xy) ","  s(x) ","  s(y) ",   &
       "  s(3) ","  s(2) ","  s(1) "]

    character(len=*), dimension(48), parameter, public :: ZAK_OH =[             &
       "     E     ","    U(z)   ","    U(y)   ","    U(x)   ","  C(xyz)_3 ",   &
       " C(-xy-z)_3"," C(x-y-z)_3"," C(-x-yz)_3"," C(2xyz)_3 ","C(2x-y-z)_3",   &
       " C(2x-yz)_3","C(-2xy-z)_3","    U(xy)  ","   U(-xy)  ","   C(3z)_4 ",   &
       "   C(z)_4  ","   C(3x)_4 ","    U(yz)  ","   U(y-z)  ","   C(x)_4  ",   &
       "   C(y)_4  ","    U(xz)  ","   C(3y)_4 ","   U(x-z)  ","      I    ",   &
       "    s(z)   ","    s(y)   ","    s(x)   "," S(5xyz)_6 ","S(-5xy-z)_6",   &
       "S(5x-y-z)_6","S(-5x-yz)_6","  S(xyz)_6 "," S(x-y-z)_6"," S(-x-yz)_6",   &
       " S(-xy-z)_6","    s(xy)  ","   s(-xy)  ","   S(z)_4  ","  S(3z)_4  ",   &
       "   S(x)_4  ","    s(yz)  ","   s(y-z)  ","  S(3x)_4  ","  S(3y)_4  ",   &
       "    s(xz)  ","   S(y)_4  ","   s(x-z)  "]

 
    Interface
       Module Function Get_Compact_HM_Symbol(HM) Result(C_HM)
          !---- Arguments ----!
          character(len=*), intent(in) :: HM
          character(len=:), allocatable :: C_HM
       End Function Get_Compact_HM_Symbol
       
       Module Function Get_HM_Symbol_from_Compact_HM(C_HM) Result(HM)
          !---- Arguments ----!
          character(len=*), intent(in) :: C_HM
          character(len=:), allocatable :: HM
       End Function Get_HM_Symbol_from_Compact_HM
       
       Module Function Get_IT_Generators(Spg) Result(StrGen)
          !---- Arguments ----!
          character(len=*), intent(in)  :: spg      
          character(len=:), allocatable :: Strgen  
       End Function Get_IT_Generators
       
       Module Subroutine Get_SpaceG_Symbols(Str, HM, Hall, IT, C_HM)
          !---- Arguments ----!
          character(len=*),           intent(in)  :: Str   
          character(len=*), optional, intent(out) :: HM    
          character(len=*), optional, intent(out) :: Hall  
          integer,          optional, intent(out) :: IT    
          character(len=*), optional, Intent(out) :: C_HM
       End Subroutine Get_SpaceG_Symbols
       
       Module Subroutine Set_It_Gen()
          !---- Arguments ----!
       End Subroutine Set_It_Gen
       
       Module Subroutine Set_Shubnikov_Info()
          !---- Arguments ----!
       End Subroutine Set_Shubnikov_Info
       
       Module Subroutine Set_Spgr_Info()
          !---- Arguments ----!
       End Subroutine Set_Spgr_Info
       
       Module Subroutine Set_System_Equiv()
          !---- Arguments ----!
       End Subroutine Set_System_Equiv  
       
       Module Subroutine Set_Wyckoff_Info()
          !---- Arguments ----!
       End Subroutine Set_Wyckoff_Info
       
       Module Subroutine Remove_Shubnikov_Info()
          !---- Arguments ----!
       End Subroutine Remove_Shubnikov_Info 
       
       Module Subroutine Remove_Spgr_Info()
          !---- Arguments ----!
       End Subroutine Remove_Spgr_Info 
       
       Module Subroutine Remove_System_Equiv()
          !---- Arguments ----!
       End Subroutine Remove_System_Equiv
       
       Module Subroutine Remove_Wyckoff_Info()
          !---- Arguments ----!
       End Subroutine Remove_Wyckoff_Info
           
    End Interface
 
 End Module CFML_Symmetry_Tables
