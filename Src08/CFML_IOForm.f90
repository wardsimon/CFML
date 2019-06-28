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
!!----               Ross Angel         (University of Pavia)
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
!!---- MODULE: CFML_IO_Formats
!!----   INFO: Creation/Conversion for several formats
!!----
!!
Module CFML_IOForm

    !---- Use modules ----!
    Use CFML_GlobalDeps,                only: CP, PI, EPS, TAB, Err_CFML, Clear_Error
    Use CFML_Strings,                   only: l_case, u_case, get_num, cut_string, get_words, &
                                              get_numstd, Read_Key_Value, Read_Key_ValueSTD
    Use CFML_Atoms,                     only: Atm_Type, Atm_Std_Type, Matm_std_type, Atm_Ref_Type, &
                                              AtList_Type, Allocate_Atom_List
    Use CFML_Metrics,                   only: Cell_Type, Cell_G_Type, Set_Crystal_Cell, U_equiv
    Use CFML_gSpaceGroups,              only: SpG_Type

    !---- Variables ----!
    implicit none

    private


    !---- Public Functions ----!
    
    !---- Public subroutines ----!

    !---- Definitions ----!


    !!----
    !!---- TYPE :: FLIST_TYPE
    !!--..
    !!
    Type, public :: Fil_Type
       character(len=:),               allocatable :: Fname      ! Name of file
       integer,                                    :: Iunit =0   ! Logical unit
       integer                                     :: nlines=0   ! Number of lines
       character(len=:), dimension(:), allocatable :: line       ! String contains
    End Type Fil_Type


    !---- Overloaded Zone ----!
    
    
    !---- Interface zone ----!
    Interface
       Module Subroutine Read_CFL_Atom(lines,n_ini, n_end, At_List)
          !---- Arguments ----!
          character(len=*), dimension(:), intent(in)     :: lines
          integer,                        intent(in out) :: n_ini
          integer,                        intent(in)     :: n_end
          Type (AtList_Type),             intent(out)    :: At_List
       End Subroutine Read_CFL_Atom
       
       Module Subroutine Read_CFL_Cell(lines, n_ini, n_end, Cell)
          !---- Arguments ----!
          character(len=*), dimension(:),  intent(in)     :: lines   
          integer,                         intent(in out) :: n_ini   
          integer,                         intent(in)     :: n_end   
          class(Cell_Type),                intent(out)    :: Cell    
       End Subroutine Read_CFL_Cell
       
       Module Subroutine Read_Cif_Atom(lines,n_ini,n_end, At_List)
          !---- Arguments ----!
          character(len=*), dimension(:),   intent(in)      :: lines
          integer,                          intent(in out)  :: n_ini
          integer,                          intent(in)      :: n_end
          type (AtList_type),               intent(out)     :: At_List
       End Subroutine Read_Cif_Atom
       
       Module Subroutine Read_Cif_Cell(lines, n_ini, n_end, Cell)
          !---- Arguments ----!
          character(len=*), dimension(:),  intent(in)     :: lines   
          integer,                         intent(in out) :: n_ini   
          integer,                         intent(in)     :: n_end   
          class(Cell_Type),                intent(out)    :: Cell    
       End Subroutine Read_Cif_Cell
       
       Module Subroutine Read_Cif_ChemName(lines,N_ini,N_End,ChemName)
          !---- Arguments ----!
          character(len=*),  dimension(:), intent(in) :: lines
          integer,           intent(in out)           :: n_ini
          integer,           intent(in)               :: n_end
          character(len=*),  intent(out)              :: ChemName
       End Subroutine Read_Cif_ChemName
       
       Module Subroutine Read_Cif_Z(lines, n_ini, n_end, Z)
          !---- Arguments ----!
          character(len=*), dimension(:),  intent(in)     :: lines   
          integer,                         intent(in out) :: n_ini   
          integer,                         intent(in)     :: n_end   
          real(kind=cp),                   intent(out)    :: Z       
       End Subroutine Read_Cif_Z
       
       Module Subroutine Read_Cif_Wave(lines, n_ini, n_end, Wave)
          !---- Arguments ----!
          character(len=*), dimension(:),  intent(in)     :: lines   
          integer,                         intent(in out) :: n_ini   
          integer,                         intent(in)     :: n_end   
          real(kind=cp),                   intent(out)    :: Wave    
       End Subroutine Read_Cif_Wave
       
       Module Subroutine Read_Cif_Cont(lines,N_Ini,N_End,N_Elem_Type,Elem_Type,N_Elem)
          !---- Arguments ----!
          character(len=*), dimension(:),      intent(in)      :: lines
          integer,                             intent(in out)  :: n_ini
          integer,                             intent(in)      :: n_end
          integer,                             intent(out)     :: n_elem_type
          character(len=*), dimension(:),      intent(out)     :: elem_type
          real(kind=cp), dimension(:),optional,intent(out)     :: n_elem
       End Subroutine Read_Cif_Cont
       
       Module Subroutine Read_Cif_Pressure(lines,N_ini,N_End, P, SigP)
          !---- Arguments ----!
          character(len=*),  dimension(:), intent(in) :: lines
          integer,           intent(in out)           :: n_ini
          integer,           intent(in)               :: n_end
          real(kind=cp),     intent(out)              :: p
          real(kind=cp),     intent(out)              :: sigp
       End Subroutine Read_Cif_Pressure
       
       Module Subroutine Read_Cif_Title(lines,N_Ini,N_End,Title)
          !---- Arguments ----!
          character(len=*),  dimension(:), intent(in) :: lines
          integer,           intent(in out)           :: n_ini
          integer,           intent(in)               :: n_end
          character(len=*),  intent(out)              :: title
       End Subroutine Read_Cif_Title
       
       Module Subroutine Read_Cif_Temp(lines,N_Ini,N_End,T,SigT)
          !---- Arguments ----!
          character(len=*),  dimension(:), intent(in) :: lines
          integer,           intent(in out)           :: n_ini
          integer,           intent(in)               :: n_end
          real(kind=cp),     intent(out)              :: T
          real(kind=cp),     intent(out)              :: sigT
       End Subroutine Read_Cif_Temp
       
       Module Subroutine Read_Cif_Hall(lines, N_Ini, N_End, Hall)
          !---- Arguments ----!
          character(len=*), dimension(:), intent(in) :: lines
          integer,          intent(in out)           :: n_ini
          integer,          intent(in)               :: n_end
          character(len=*), intent(out)              :: Hall
       End Subroutine Read_Cif_Hall
       
       Module Subroutine Read_Cif_HM(lines, N_Ini, N_End, Spgr_Hm)
          !---- Arguments ----!
          character(len=*),  dimension(:), intent(in) :: lines
          integer,           intent(in out)           :: n_ini
          integer,           intent(in)               :: n_end
          character(len=*),  intent(out)              :: spgr_hm
       End Subroutine Read_Cif_HM
       
       Module Subroutine Read_Cif_Symm(lines,N_Ini,N_End, N_Oper, Oper_Symm)
          !---- Arguments ----!
          character(len=*), dimension(:), intent(in)     :: lines
          integer,                        intent(in out) :: n_ini
          integer,                        intent(in)     :: n_end
          integer,                        intent(out)    :: n_oper
          character(len=*), dimension(:), intent(out)    :: oper_symm
       End Subroutine Read_Cif_Symm
       
       Module Subroutine Read_Shx_Atom(lines, n_ini, n_end, n_fvar, fvar, elem_type, cell, At_List)
          !---- Arguments ----!
          character(len=*), dimension(:), intent(in)      :: lines
          integer,                        intent(in out)  :: n_ini
          integer,                        intent(in)      :: n_end
          integer,                        intent(in)      :: n_fvar
          real(kind=cp), dimension(:),    intent(in)      :: fvar
          character(len=*), dimension(:), intent(in)      :: elem_type
          class(Cell_G_Type),             intent(in)      :: Cell
          type (AtList_type),             intent(out)     :: At_List
       End Subroutine Read_Shx_Atom
       
       Module Subroutine Read_Shx_Cell(lines, n_ini, n_end, Cell)
          !---- Arguments ----!
          character(len=*), dimension(:),  intent(in)     :: lines   
          integer,                         intent(in out) :: n_ini   
          integer,                         intent(in)     :: n_end   
          class(Cell_Type),                intent(out)    :: Cell    
       End Subroutine Read_Shx_Cell
       
       Module Subroutine Read_Shx_Wave(lines, n_ini, n_end, Wave)
          !---- Arguments ----!
          character(len=*), dimension(:),  intent(in)     :: lines   
          integer,                         intent(in out) :: n_ini   
          integer,                         intent(in)     :: n_end   
          real(kind=cp),                   intent(out)    :: Wave    
       End Subroutine Read_Shx_Wave
       
       Module Subroutine Read_Shx_Z(lines, n_ini, n_end, Z)
          !---- Arguments ----!
          character(len=*), dimension(:),  intent(in)     :: lines   
          integer,                         intent(in out) :: n_ini   
          integer,                         intent(in)     :: n_end   
          real(kind=cp),                   intent(out)    :: Z       
       End Subroutine Read_Shx_Z
       
       Module Subroutine Read_Shx_Fvar(Lines,n_ini,n_end, n_fvar, fvar)
          !---- Arguments ----!
          character(len=*), dimension(:), intent(in)    :: lines
          integer,                        intent(in out):: n_ini
          integer,                        intent(in)    :: n_end
          integer,                        intent(out)   :: n_fvar
          real(kind=cp), dimension(:),    intent(out)   :: fvar
       End Subroutine Read_Shx_Fvar
       
       Module Subroutine Read_Shx_Titl(lines,n_ini,n_end,Title)
          !---- Arguments ----!
          character(len=*),dimension(:), intent(in)     :: lines
          integer,                       intent(in out) :: n_ini
          integer,                       intent(in)     :: n_end
          character(len=*),              intent(out)    :: title
       End Subroutine Read_Shx_Titl
       
       Module Subroutine Read_Shx_Latt(lines,n_ini,n_end,latt)
          !---- Arguments ----!
          character(len=*), dimension(:), intent(in) :: lines
          integer,           intent(in out)          :: n_ini
          integer,           intent(in)              :: n_end
          integer,           intent(out)             :: latt
       End Subroutine Read_Shx_Latt
       
       Module Subroutine Read_Shx_Cont(lines,n_ini,n_end, n_elem_type, elem_type, n_elem)
          !---- Arguments ----!
          character(len=*), dimension(:),           intent(in)      :: lines
          integer,                                  intent(in out)  :: n_ini
          integer,                                  intent(in)      :: n_end
          integer,                                  intent(out)     :: n_elem_type
          character(len=*), dimension(:),           intent(out)     :: elem_type
          real(kind=cp),    dimension(:), optional, intent(out)     :: n_elem
       End Subroutine Read_Shx_Cont
       
       Module Subroutine Read_Shx_Symm(lines,n_ini,n_end,n_oper,oper_symm)
          !---- Arguments ----!
          character(len=*), dimension(:), intent(in) :: lines
          integer,          intent(in out)           :: n_ini
          integer,          intent(in)               :: n_end
          integer,          intent(out)              :: n_oper
          character(len=*), dimension(:),intent(out) :: oper_symm
       End Subroutine Read_Shx_Symm 
       
       Module Subroutine Write_Shx_Template(filename,code,title,lambda,z,cell,spg,At_List)
          !---- Arguments ----!
          character(len=*),        intent(in) :: filename
          integer,                 intent(in) :: code        
          character(len=*),        intent(in) :: title
          real(kind=cp),           intent(in) :: lambda
          integer,                 intent(in) :: z
          class(cell_Type),        intent(in) :: cell
          class(SpG_Type),         intent(in) :: SpG
          type(atlist_type),       intent(in) :: at_List
       End Subroutine Write_Shx_Template   
          
    End Interface   



End Module CFML_IOForm

