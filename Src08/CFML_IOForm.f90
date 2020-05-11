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
   Use CFML_GlobalDeps,        only: SP, CP, PI, EPS, TAB, Err_CFML, Clear_Error
   Use CFML_Maths,             only: Get_Eps_Math
   Use CFML_Rational
   Use CFML_Strings,           only: l_case, u_case, get_num, cut_string, get_words, &
                                     get_numstd, Read_Key_Value, Read_Key_ValueSTD,  &
                                     string_numstd, Number_Lines, Reading_Lines,     &
                                     FindFMT, Init_FindFMT, String_Array_Type,       &
                                     File_type, Reading_File, Get_Transf, Get_Extension

   Use CFML_Atoms,             only: Atm_Type, Atm_Std_Type, Matm_std_type, Atm_Ref_Type, &
                                     AtList_Type, Allocate_Atom_List, Init_Atom_Type

   Use CFML_Metrics,           only: Cell_Type, Cell_G_Type, Set_Crystal_Cell, U_equiv, &
                                     get_U_from_Betas, get_Betas_from_U, get_Betas_from_B

   Use CFML_gSpaceGroups,      only: SpG_Type, SuperSpaceGroup_Type, Kvect_Info_Type,   &
                                     Change_Setting_SpaceG, Set_SpaceGroup, Get_Multip_Pos,&
                                     Get_Orbit, Get_Moment_Ctr, Get_TFourier_Ctr

   Use CFML_DiffPatt,          only: DiffPat_Type, DiffPat_E_Type

   !---- Variables ----!
   implicit none

   private


   !---- Public subroutines ----!

   public :: Read_Xtal_Structure, &
             Write_Cif_Template, Write_SHX_Template

   !--------------------!
   !---- PARAMETERS ----!
   !--------------------!
   character(len=*), parameter :: DIGCAR="0123456789+-"  ! Digit character and signs
   integer,          parameter :: MAX_PHASES=30          ! Number of Maximum Phases
   real(kind=cp),    parameter :: EPSV=0.0001_cp         ! Small real value to be used for decisions


   !---------------!
   !---- TYPES ----!
   !---------------!

   !!----
   !!---- TYPE :: INTERVAL_TYPE
   !!--..
   !!
   Type, public :: Interval_Type
      real(kind=cp) :: Mina=0.0_cp  !low limit
      real(kind=cp) :: Maxb=0.0_cp  !high limit
   End Type Interval_Type

   !!----
   !!---- TYPE :: JOB_INFO_TYPE
   !!--..
   !!
   Type, public :: Job_Info_type
      character(len=120)                            :: Title       =" "       ! Title
      integer                                       :: Num_Phases  =0         ! Number of phases
      integer                                       :: Num_Patterns=0         ! Number of patterns
      integer                                       :: Num_cmd     =0         ! Number of command lines
      character(len=16),  dimension(:), allocatable :: Patt_typ               ! Type of Pattern
      character(len=128), dimension(:), allocatable :: Phas_nam               ! Name of phases
      character(len=128), dimension(:), allocatable :: cmd                    ! Command lines: text for actions
      type(interval_type),dimension(:), allocatable :: range_stl              ! Range in sinTheta/Lambda
      type(interval_type),dimension(:), allocatable :: range_q                ! Range in 4pi*sinTheta/Lambda
      type(interval_type),dimension(:), allocatable :: range_d                ! Range in d-spacing
      type(interval_type),dimension(:), allocatable :: range_2theta           ! Range in 2theta-spacing
      type(interval_type),dimension(:), allocatable :: range_Energy           ! Range in Energy
      type(interval_type),dimension(:), allocatable :: range_tof              ! Range in Time of Flight
      type(interval_type),dimension(:), allocatable :: Lambda                 ! Lambda
      real(kind=cp)      ,dimension(:), allocatable :: ratio                  ! ratio lambda2/lambda1
      real(kind=cp)      ,dimension(:), allocatable :: dtt1,dtt2              ! d-to-TOF coefficients
   End Type Job_Info_type

   !---- Overloaded Zone ----!
   !Interface Readn_Set_Xtal_Structure
   !   Module Procedure Readn_Set_Xtal_Structure_Split  ! For Cell, Spg, A types
   !   !Module Procedure Readn_Set_Xtal_Structure_Molcr ! For Molecular Crystal Type
   !End Interface

   !---- Interface zone ----!
   Interface
      Module Subroutine Get_Job_Info(cfl,Job_info, i_ini,i_end)
         type(File_Type),      intent(in)  :: cfl
         type(job_info_type),  intent(out) :: Job_info
         integer,              intent(in)  :: i_ini, i_end
      End Subroutine Get_Job_Info

      Module Subroutine Read_Atom(Str, Atm)
         character(len=*), intent(in)   :: Str
         class (Atm_Type), intent(out)  :: Atm
      End Subroutine Read_Atom

      Module Subroutine Read_Cell(Str,Celda,Std,Cell,CFrame)
         character(len=*),                      intent(in)  :: Str
         real(kind=cp), dimension(6),           intent(out) :: Celda
         real(kind=cp), dimension(6), optional, intent(out) :: Std
         class(Cell_Type),            optional, intent(out) :: Cell
         character(len=*),            optional, intent(in)  :: CFrame
      End Subroutine Read_Cell

      Module Subroutine Read_Modulation_Amplitudes(Str, Atm, Ulabel, Nt)
         character(len=*),    intent(in )     :: str
         class(MAtm_Std_Type),intent(in out)  :: Atm
         character(len=*),    intent(in)      :: ulabel
         integer,             intent(in)      :: nt
      End Subroutine Read_Modulation_Amplitudes

      Module Subroutine Read_Moment(Str, Atm)
         character(len=*), intent(in )     :: Str
         Class (Atm_Type), intent(in out)  :: Atm
      End Subroutine Read_Moment

      Module Subroutine Read_RngSintL(Str, v1,v2)
         character(len=*), intent(in)  :: str
         real(kind=cp),    intent(out) :: v1,v2
      End Subroutine Read_RngSintL

      Module Subroutine Read_SpaceGroup(Str,Spg)
         character(len=*), intent(in)     :: Str
         class(SpG_Type),  intent(out)    :: SpG
      End Subroutine Read_SpaceGroup

      Module Subroutine Read_Transf(str, trans, orig)
         character(len=*),                intent(in)     :: str
         real(kind=cp),dimension(3,3),    intent(out)    :: trans
         real(kind=cp),dimension(3  ),    intent(out)    :: orig
      End Subroutine Read_Transf

      Module Subroutine Read_UTherms(Str, Atm)
         character(len=*),  intent(in )     :: Str
         Class (Atm_Type),  intent(in out)  :: Atm
      End Subroutine Read_UTherms

      Module Subroutine Read_Wavelength(Str,v1,v2,ratio)
         character(len=*), intent(in)     :: Str
         real(kind=cp),    intent(out)    :: v1,v2
         real(kind=cp),    intent(out)    :: ratio
      End Subroutine Read_Wavelength

      Module Subroutine Read_CFL_Atoms(cfl, AtmList, Type_Atm, d, i_ini, i_end)
         type(File_Type),      intent(in)     :: cfl
         Type(AtList_Type),    intent(out)    :: AtmList
         character(len=*),     intent(in)     :: Type_Atm
         integer,              intent(in)     :: d
         integer, optional,    intent(in)     :: i_ini, i_end
      End Subroutine Read_CFL_Atoms

      Module Subroutine Read_CFL_Cell(cfl, Cell, CFrame, i_ini, i_end )
         type(File_Type),            intent(in)     :: cfl
         class(Cell_Type),           intent(out)    :: Cell
         character(len=*), optional, intent( in)    :: CFrame
         integer, optional,          intent(in)     :: i_ini, i_end
      End Subroutine Read_CFL_Cell

      Module Subroutine Read_CFL_KVectors(cfl, Kvec, i_ini, i_end)
         type(File_Type),         intent(in)     :: cfl
         type(kvect_info_Type),   intent(out)    :: Kvec
         integer, optional,          intent(in)     :: i_ini, i_end
      End Subroutine Read_CFL_KVectors

      Module Subroutine Read_CFL_SpG(cfl, SpG, xyz_type, i_ini, i_end)
         Type(File_Type),                 intent(in)     :: cfl
         class(SpG_Type),                 intent(out)    :: SpG
         character(len=*), optional,      intent(in)     :: xyz_type
         integer, optional,               intent(in)     :: i_ini, i_end
      End Subroutine Read_CFL_SpG

      Module Subroutine Write_CFL_Atoms(AtmList, Lun, Cell)
         Type(AtList_Type),            intent(in) :: AtmList
         integer,            optional, intent(in) :: Lun
         class(Cell_G_Type), optional, intent(in) :: Cell
      End Subroutine Write_CFL_Atoms

       Module Subroutine Write_CFL_File(Lun,Cell, SpG, Atm, Title)
          integer,                       intent(in)    :: lun
          class(Cell_G_Type),            intent(in)    :: Cell
          class(SpG_Type)  ,             intent(in)    :: SpG
          Type(AtList_Type), optional,   intent(in)    :: Atm
          character(len=*),  optional,   intent(in)    :: Title
      End Subroutine Write_CFL_File

      Module Subroutine Read_CIF_Atom(lines,n_ini,n_end, AtmList)
         character(len=*), dimension(:),   intent(in)      :: lines
         integer,                          intent(in out)  :: n_ini
         integer,                          intent(in)      :: n_end
         type (AtList_type),               intent(out)     :: AtmList
      End Subroutine Read_CIF_Atom

      Module Subroutine Read_CIF_Cell(lines, n_ini, n_end, Cell)
         character(len=*), dimension(:),  intent(in)     :: lines
         integer,                         intent(in out) :: n_ini
         integer,                         intent(in)     :: n_end
         class(Cell_Type),                intent(out)    :: Cell
      End Subroutine Read_CIF_Cell

      Module Subroutine Read_CIF_ChemName(lines,N_ini,N_End,ChemName)
         character(len=*),  dimension(:), intent(in) :: lines
         integer,           intent(in out)           :: n_ini
         integer,           intent(in)               :: n_end
         character(len=*),  intent(out)              :: ChemName
      End Subroutine Read_Cif_ChemName

      Module Subroutine Read_CIF_Z(lines, n_ini, n_end, Z)
         character(len=*), dimension(:),  intent(in)     :: lines
         integer,                         intent(in out) :: n_ini
         integer,                         intent(in)     :: n_end
         real(kind=cp),                   intent(out)    :: Z
      End Subroutine Read_CIF_Z

      Module Subroutine Read_CIF_Wave(lines, n_ini, n_end, Wave)
         character(len=*), dimension(:),  intent(in)     :: lines
         integer,                         intent(in out) :: n_ini
         integer,                         intent(in)     :: n_end
         real(kind=cp),                   intent(out)    :: Wave
      End Subroutine Read_CIF_Wave

      Module Subroutine Read_CIF_Cont(lines,N_Ini,N_End,N_Elem_Type,Elem_Type,N_Elem)
         character(len=*), dimension(:),      intent(in)      :: lines
         integer,                             intent(in out)  :: n_ini
         integer,                             intent(in)      :: n_end
         integer,                             intent(out)     :: n_elem_type
         character(len=*), dimension(:),      intent(out)     :: elem_type
         real(kind=cp), dimension(:),optional,intent(out)     :: n_elem
      End Subroutine Read_CIF_Cont

      Module Subroutine Read_CIF_Pressure(lines,N_ini,N_End, P, SigP)
         character(len=*),  dimension(:), intent(in) :: lines
         integer,           intent(in out)           :: n_ini
         integer,           intent(in)               :: n_end
         real(kind=cp),     intent(out)              :: p
         real(kind=cp),     intent(out)              :: sigp
      End Subroutine Read_CIF_Pressure

      Module Subroutine Read_CIF_Title(lines,N_Ini,N_End,Title)
         character(len=*),  dimension(:), intent(in) :: lines
         integer,           intent(in out)           :: n_ini
         integer,           intent(in)               :: n_end
         character(len=*),  intent(out)              :: title
      End Subroutine Read_CIF_Title

      Module Subroutine Read_CIF_Temp(lines,N_Ini,N_End,T,SigT)
         character(len=*),  dimension(:), intent(in) :: lines
         integer,           intent(in out)           :: n_ini
         integer,           intent(in)               :: n_end
         real(kind=cp),     intent(out)              :: T
         real(kind=cp),     intent(out)              :: sigT
      End Subroutine Read_CIF_Temp

      Module Subroutine Read_CIF_Hall(lines, N_Ini, N_End, Hall)
         character(len=*), dimension(:), intent(in) :: lines
         integer,          intent(in out)           :: n_ini
         integer,          intent(in)               :: n_end
         character(len=*), intent(out)              :: Hall
      End Subroutine Read_CIF_Hall

      Module Subroutine Read_CIF_HM(lines, N_Ini, N_End, Spgr_Hm)
         character(len=*),  dimension(:), intent(in) :: lines
         integer,           intent(in out)           :: n_ini
         integer,           intent(in)               :: n_end
         character(len=*),  intent(out)              :: spgr_hm
      End Subroutine Read_CIF_HM

      Module Subroutine Read_CIF_Symm(lines,N_Ini,N_End, N_Oper, Oper_Symm)
         character(len=*), dimension(:), intent(in)     :: lines
         integer,                        intent(in out) :: n_ini
         integer,                        intent(in)     :: n_end
         integer,                        intent(out)    :: n_oper
         character(len=*), dimension(:), intent(out)    :: oper_symm
      End Subroutine Read_CIF_Symm

      Module Subroutine Write_CIF_Powder_Profile(filename,Pat,r_facts)
         character(len=*),                      intent(in) :: filename
         class(DiffPat_Type),                   intent(in) :: Pat
         real(kind=cp), dimension(4), optional, intent(in) :: r_facts
      End Subroutine Write_CIF_Powder_Profile

      Module Subroutine Read_SHX_Atom(shx, n_fvar, fvar, elem_type, cell, AtmList)
         type(File_Type),                intent(in)      :: shx
         integer,                        intent(in)      :: n_fvar
         real(kind=cp), dimension(:),    intent(in)      :: fvar
         character(len=*), dimension(:), intent(in)      :: elem_type
         class(Cell_G_Type),             intent(in)      :: Cell
         type (AtList_type),             intent(out)     :: AtmList
      End Subroutine Read_SHX_Atom

      Module Subroutine Read_SHX_Cell(shx, Cell)
         type(File_Type),                 intent(in)     :: shx
         class(Cell_Type),                intent(out)    :: Cell
      End Subroutine Read_SHX_Cell

      Module Subroutine Read_SHX_Wave(shx, Wave)
         type(File_Type),  intent(in)     :: shx
         real(kind=cp),    intent(out)    :: Wave
      End Subroutine Read_SHX_Wave

      Module Subroutine Read_SHX_Z(shx, Z)
         type(File_Type),  intent(in)     :: shx
         real(kind=cp),    intent(out)    :: Z
      End Subroutine Read_SHX_Z

      Module Subroutine Read_SHX_Fvar(shx, n_fvar, fvar)
         type(File_Type),                intent(in)    :: shx
         integer,                        intent(out)   :: n_fvar
         real(kind=cp), dimension(:),    intent(out)   :: fvar
      End Subroutine Read_SHX_Fvar

      Module Subroutine Read_SHX_Titl(shx,Title)
         type(File_Type),  intent(in)     :: shx
         character(len=*), intent(out)    :: title
      End Subroutine Read_SHX_Titl

      Module Subroutine Read_SHX_Latt(shx,latt)
         type(File_Type),  intent(in)    :: shx
         integer,          intent(out)   :: latt
      End Subroutine Read_SHX_Latt

      Module Subroutine Read_SHX_Cont(shx, n_elem_type, elem_type, n_elem)
         type(File_Type),                           intent(in)     :: shx
         integer,                                  intent(out)     :: n_elem_type
         character(len=*), dimension(:),           intent(out)     :: elem_type
         real(kind=cp),    dimension(:), optional, intent(out)     :: n_elem
      End Subroutine Read_SHX_Cont

      Module Subroutine Read_SHX_Symm(shx,n_oper,oper_symm)
         type(File_Type),  intent(in)               :: shx
         integer,          intent(out)              :: n_oper
         character(len=*), dimension(:),intent(out) :: oper_symm
      End Subroutine Read_SHX_Symm

      Module Subroutine Write_SHX_Template(filename,code,title,lambda,z,cell,spg,AtmList)
         !---- Arguments ----!
         character(len=*),        intent(in) :: filename
         integer,                 intent(in) :: code
         character(len=*),        intent(in) :: title
         real(kind=cp),           intent(in) :: lambda
         integer,                 intent(in) :: z
         class(cell_Type),        intent(in) :: cell
         class(SpG_Type),         intent(in) :: SpG
         type(atlist_type),       intent(in) :: atmList
      End Subroutine Write_SHX_Template

      Module Subroutine Read_XTal_CFL(cfl, Cell, SpG, AtmList, Nphase, CFrame, Job_Info)
         type(File_Type),               intent(in)  :: cfl
         class(Cell_Type),              intent(out) :: Cell
         class(SpG_Type),               intent(out) :: SpG
         Type(AtList_Type),             intent(out) :: Atmlist
         Integer,             optional, intent(in)  :: Nphase
         character(len=*),    optional, intent(in)  :: CFrame
         Type(Job_Info_type), optional, intent(out) :: Job_Info
      End Subroutine Read_XTal_CFL

      Module Subroutine Read_XTal_SHX(shx, Cell, SpG, Atm)
         type(File_Type),                 intent(in)  :: shx
         class (Cell_G_Type),             intent(out) :: Cell
         class (SpG_Type),                intent(out) :: SpG
         type (AtList_type),              intent(out) :: Atm
      End Subroutine Read_XTal_SHX

      Module Subroutine Write_CIF_Template(filename, Cell, SpG, Atmlist, Type_data, Code)
         character(len=*),           intent(in) :: filename     ! Filename
         class(Cell_G_Type),         intent(in) :: Cell         ! Cell parameters
         class(SpG_Type),            intent(in) :: SpG          ! Space group information
         Type (AtList_Type),         intent(in) :: AtmList      ! Atoms
         integer,                    intent(in) :: Type_data    ! 0,2:Single crystal diffraction; 1:Powder
         character(len=*),           intent(in) :: Code         ! Code or name of the structure
      End Subroutine Write_CIF_Template

    End Interface

    Contains

    !!----
    !!---- READ_XTAL_STRUCTURE
    !!----
    !!----
    !!----
    !!----
    !!---- 09/05/2020
    !!
    Subroutine Read_Xtal_Structure(filenam, Cell, Spg, Atm, IPhase, FType)
       !---- Arguments ----!
       character(len=*),          intent( in)     :: filenam    ! Name of the file
       class (Cell_G_Type),       intent(out)     :: Cell       ! Cell object
       class (SpG_Type),          intent(out)     :: SpG        ! Space Group object
       type (Atlist_type),        intent(out)     :: Atm        ! Atom List object
       integer,         optional, intent(in)      :: IPhase     ! Number of phase
       type(File_Type), optional, intent(out)     :: FType      ! File type

       !---- Local Variables ----!
       character(len=6):: Ext
       type(File_Type) :: F

       !> Init
       call clear_error()

       !> Load filename
       F=Reading_File(trim(filenam))
       if (err_CFML%Ierr /= 0) return

       if (F%nlines ==0) then
          Err_CFML%Ierr=1
          Err_CFML%Msg="Zero lines in the file "//trim(filenam)
          return
       end if

       !> Extension
       Ext=get_extension(trim(filenam))

       select case (trim(u_case(ext)))
          case ('CFL')
             call Read_XTal_CFL(f, Cell, SpG, Atm)
          case ('CIF')
          case ('INS','RES')
             call Read_XTal_SHX(f, Cell, SpG, Atm)
          case ('PCR')
          case default
             Err_CFML%Ierr=1
             Err_CFML%Msg="The extension file is unknown. Pease, check it! Ext= "//trim(filenam)
             return
       end select

       !> End
       if (present(FType)) FType=F

    End Subroutine Read_Xtal_Structure

    !!--++
    !!--++ Subroutine Readn_Set_Xtal_Structure_Split(filenam,Cell,SpG,A,Type_Atm,Mode,Iphase,Job_Type,File_List,CFrame)
    !!--++    character(len=*),              intent( in)     :: filenam    ! In -> Name of the file
    !!--++    Type (Crystal_Cell_Type),      intent(out)     :: Cell       ! Out -> Cell object
    !!--++    Type (Space_Group_Type),       intent(out)     :: SpG        ! Out -> Space Group object
    !!--++    Type (atom_list_type),         intent(out)     :: A          ! Out -> Atom_List object
    !!--++    Character(len=*),              intent( in)     :: Type_Atm   ! In ->  Type of atoms (Atm,Atm_Std, Matm_std,Atm_ref,Matm_Ref
    !!--++    Character(len=*),    optional, intent( in)     :: Mode       ! In -> if Mode="CIF" filenam
    !!--++                                                                         is of CIF type format
    !!--++    Integer,             optional, intent( in)     :: Iphase     ! Number of the phase.
    !!--++    Type(Job_Info_type), optional, intent(out)     :: Job_Info   ! Diffraction conditions
    !!--++    Type(file_list_type),optional, intent(in out)  :: file_list  ! Complete file to be used by
    !!--++                                                                   the calling program or other procedures
    !!--++    Character(len=*),    optional, intent( in)     :: CFrame     !Cartesian Frame
    !!--++
    !!--++    Overloaded
    !!--++    Subroutine to read and input file and construct the crystal structure
    !!--++    in terms of the objects Cell, SpG and A. The optional argument Iphase is an integer
    !!--++    telling to the program to read the phase number Iphase in the case of the presence
    !!--++    of more than one phase. If absent only the first phase is read.
    !!--++
    !!--++ Update: April - 2005, Febraury 2020
    !!
    !Subroutine Readn_Set_Xtal_Structure_Split(Filenam, Cell, SpG, A,Type_Atm,Mode,Iphase,Job_Info,file_list,CFrame)
    !
    !   select case(modec)
    !
    !       !case("cif")
    !       !   if (present(iphase)) then
    !       !      if(present(CFrame)) then
    !       !        call readn_set_xtal_cif(file_dat,nlines,Cell,Spg, A,CFrame,NPhase=IPhase)
    !       !      else
    !       !        call readn_set_xtal_cif(file_dat,nlines,Cell,Spg, A,NPhase=IPhase)
    !       !      end if
    !       !   else
    !       !      if(present(CFrame)) then
    !       !        call readn_set_xtal_cif(file_dat,nlines,Cell,Spg,A,CFrame)
    !       !      else
    !       !        call readn_set_xtal_cif(file_dat,nlines,Cell,Spg,A)
    !       !      end if
    !       !   end if
    !       !
    !       !case("pcr")
    !       !   if (present(iphase)) then
    !       !      if(present(CFrame)) then
    !       !        call readn_set_xtal_pcr(file_dat,nlines,Cell,Spg, A,CFrame,NPhase=IPhase)
    !       !      else
    !       !        call readn_set_xtal_pcr(file_dat,nlines,Cell,Spg, A,NPhase=IPhase)
    !       !      end if
    !       !   else
    !       !      if(present(CFrame)) then
    !       !        call readn_set_xtal_pcr(file_dat,nlines,Cell,Spg,A,CFrame)
    !       !      else
    !       !        call readn_set_xtal_pcr(file_dat,nlines,Cell,Spg,A)
    !       !      end if
    !       !   end if
    !       !
    !       case default
    !          !---- CFL Format ----!
    !          if (present(Job_Info)) then
    !             if (present(iphase)) then
    !                if(present(CFrame)) then
    !                  call readn_set_xtal_cfl(file_dat,nlines,Cell,Spg,A,Type_Atm,CFrame,NPhase=IPhase,Job_Info=Job_Info)
    !                else
    !                  call readn_set_xtal_cfl(file_dat,nlines,Cell,Spg,A,Type_Atm,NPhase=IPhase,Job_Info=Job_Info)
    !                end if
    !             else
    !                if(present(CFrame)) then
    !                  call readn_set_xtal_cfl(file_dat,nlines,Cell,Spg,A,Type_Atm,CFrame,Job_Info=Job_Info)
    !                else
    !                  call readn_set_xtal_cfl(file_dat,nlines,Cell,Spg,A,Type_Atm,Job_Info=Job_Info)
    !                end if
    !             end if
    !          else
    !             if (present(iphase)) then
    !                if(present(CFrame)) then
    !                  call readn_set_xtal_cfl(file_dat,nlines,Cell,Spg,A,Type_Atm,CFrame,NPhase=IPhase)
    !                else
    !                  call readn_set_xtal_cfl(file_dat,nlines,Cell,Spg,A,Type_Atm,NPhase=IPhase)
    !                end if
    !             else
    !                if(present(CFrame)) then
    !                  call readn_set_xtal_cfl(file_dat,nlines,Cell,Spg,A,Type_Atm,CFrame)
    !                else
    !                  call readn_set_xtal_cfl(file_dat,nlines,Cell,Spg,A,Type_Atm)
    !                end if
    !             end if
    !          end if
    !
    !   end select
    !
    !End Subroutine Readn_Set_Xtal_Structure_Split

End Module CFML_IOForm

