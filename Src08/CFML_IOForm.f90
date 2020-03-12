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
    Use CFML_GlobalDeps,        only: CP, PI, EPS, TAB, Err_CFML, Clear_Error
    Use CFML_Strings,           only: l_case, u_case, get_num, cut_string, get_words, &
                                      get_numstd, Read_Key_Value, Read_Key_ValueSTD,  &
                                      string_numstd, Number_Lines, Reading_Lines, FindFMT, &
                                      Init_FindFMT, String_Array_Type, File_type
    Use CFML_Atoms,             only: Atm_Type, Atm_Std_Type, Matm_std_type, Atm_Ref_Type, &
                                      AtList_Type, Allocate_Atom_List
    Use CFML_Metrics,           only: Cell_Type, Cell_G_Type, Set_Crystal_Cell, U_equiv, &
                                      get_U_from_Betas, get_Betas_from_U, get_Betas_from_B
    Use CFML_gSpaceGroups,      only: SpG_Type, SuperSpaceGroup_Type, Kvect_Info_Type,   &
                                      Change_Setting_SpaceG, Set_SpaceGroup, Get_Multip_Pos,&
                                      Get_Orbit, Get_Moment_Ctr, Get_TFourier_Ctr
    Use CFML_Maths,             only: Get_Eps_Math

    Use CFML_Rational
    !---- Variables ----!
    implicit none

    private


    !---- Public Functions ----!

    !---- Public subroutines ----!
    public :: Readn_Set_Xtal_Structure, Read_CFL_Cell, Read_CFL_SpG, Read_CFL_Atoms,Write_Atom_List, &
              Read_Kinfo, Check_Symmetry_Constraints
    real(kind=cp), parameter :: EPSV=0.0001_cp     ! Small real value to be used for decisions
    !---- Definitions ----!

    !!----
    !!---- TYPE :: INTERVAL_TYPE
    !!--..
    !!---- Type, public :: interval_type
    !!----    real(kind=cp) :: mina  !low limit
    !!----    real(kind=cp) :: maxb  !high limit
    !!---- End Type interval_type
    !!----
    !!---- Update: February - 2005
    !!
    Type, public :: interval_type
       real(kind=cp) :: mina  !low limit
       real(kind=cp) :: maxb  !high limit
    End Type interval_type

    !!----
    !!---- TYPE :: JOB_INFO_TYPE
    !!--..
    !!---- Type, public :: Job_Info_type
    !!----    character(len=120)                            :: Title          ! Title
    !!----    integer                                       :: Num_Phases     ! Number of phases
    !!----    integer                                       :: Num_Patterns   ! Number of patterns
    !!----    integer                                       :: Num_cmd        ! Number of command lines
    !!----    character(len=16),  dimension(:), allocatable :: Patt_typ       ! Type of Pattern
    !!----    character(len=128), dimension(:), allocatable :: Phas_nam       ! Name of phases
    !!----    character(len=128), dimension(:), allocatable :: cmd            ! Command lines: text for actions
    !!----    type(interval_type),dimension(:), allocatable :: range_stl      ! Range in sinTheta/Lambda
    !!----    type(interval_type),dimension(:), allocatable :: range_q        ! Range in 4pi*sinTheta/Lambda
    !!----    type(interval_type),dimension(:), allocatable :: range_d        ! Range in d-spacing
    !!----    type(interval_type),dimension(:), allocatable :: range_2theta   ! Range in 2theta-spacing
    !!----    type(interval_type),dimension(:), allocatable :: range_Energy   ! Range in Energy
    !!----    type(interval_type),dimension(:), allocatable :: range_tof      ! Range in Time of Flight
    !!----    type(interval_type),dimension(:), allocatable :: Lambda         ! Lambda
    !!----    real(kind=cp)      ,dimension(:), allocatable :: ratio          ! ratio lambda2/lambda1
    !!----    real(kind=cp)      ,dimension(:), allocatable :: dtt1,dtt2      ! d-to-TOF coefficients
    !!---- End Type Job_Info_type
    !!----
    !!---- Update: February - 2005
    !!
    Type, public :: Job_Info_type
       character(len=120)                            :: Title
       integer                                       :: Num_Phases
       integer                                       :: Num_Patterns
       integer                                       :: Num_cmd
       character(len=16),  dimension(:), allocatable :: Patt_typ
       character(len=128), dimension(:), allocatable :: Phas_nam
       character(len=128), dimension(:), allocatable :: cmd
       type(interval_type),dimension(:), allocatable :: range_stl
       type(interval_type),dimension(:), allocatable :: range_q
       type(interval_type),dimension(:), allocatable :: range_d
       type(interval_type),dimension(:), allocatable :: range_2theta
       type(interval_type),dimension(:), allocatable :: range_Energy
       type(interval_type),dimension(:), allocatable :: range_tof
       type(interval_type),dimension(:), allocatable :: Lambda
       real(kind=cp)      ,dimension(:), allocatable :: ratio
       real(kind=cp)      ,dimension(:), allocatable :: dtt1,dtt2
    End Type Job_Info_type

    !---- Overloaded Zone ----!
    Interface Read_CFL_SpG
      Module Procedure Read_CFL_SpG_lines
      Module Procedure Read_CFL_SpG_FileTyp
    End Interface

    Interface Read_CFL_Cell
      Module Procedure Read_CFL_Cell_lines
      Module Procedure Read_CFL_Cell_FileTyp
    End Interface

    Interface Readn_Set_Xtal_Structure
       Module Procedure Readn_Set_Xtal_Structure_Split  ! For Cell, Spg, A types
       !Module Procedure Readn_Set_Xtal_Structure_Molcr ! For Molecular Crystal Type
    End Interface

    !---- Interface zone ----!
    Interface

       Module Subroutine Check_Symmetry_Constraints(SpG,Atm)
         class(SpG_Type),   intent(in)     :: SpG
         type(AtList_Type), intent(in out) :: Atm
       End Subroutine Check_Symmetry_Constraints

       Module Subroutine Write_Atom_List(A, Iunit)
          !---- Arguments ----!
          type(AtList_Type), intent(in) :: A        ! Atom list object
          integer, optional, intent(in) :: IUnit    ! Logical unit
       End Subroutine Write_Atom_List

       Module Subroutine Get_Job_Info(file_dat,i_ini,i_end,Job_info)
          !---- Arguments ----!
          character(len=*), dimension(:), intent( in) :: file_dat
          integer,                        intent( in) :: i_ini,i_end
          type(job_info_type),            intent(out) :: Job_info
       End Subroutine Get_Job_Info

       Module Subroutine Read_CFL_Atoms(lines,n_ini, n_end, At_List,Type_Atm,d)
          !---- Arguments ----!
          character(len=*), dimension(:), intent(in)     :: lines
          integer,                        intent(in out) :: n_ini
          integer,                        intent(in)     :: n_end
          Type (AtList_Type),             intent(out)    :: At_List
          character(len=*),               intent(in)     :: Type_Atm
          integer,                        intent(in)     :: d
       End Subroutine Read_CFL_Atoms

       Module Subroutine Read_CFL_Cell_Lines(lines, n_ini, n_end, Cell, CFrame)
          !---- Arguments ----!
          character(len=*), dimension(:),  intent(in)     :: lines   ! Containing information
          integer,                         intent(in out) :: n_ini   ! Index to start
          integer,                         intent(in)     :: n_end   ! Index to Finish
          class(Cell_Type),                intent(out)    :: Cell    ! Cell object
          Character(len=*), optional,      intent( in)    :: CFrame
       End Subroutine Read_CFL_Cell_Lines

       Module Subroutine Read_CFL_Cell_FileTyp(cfl,Cell, CFrame)
          !---- Arguments ----!
          type(File_Type),                 intent(in)     :: cfl     ! File_type object
          class(Cell_Type),                intent(out)    :: Cell    ! Cell object
          Character(len=*), optional,      intent( in)    :: CFrame
       End Subroutine Read_CFL_Cell_FileTyp

       Module Subroutine Read_CFL_SpG_Lines(lines, n_ini, n_end, SpG, xyz_type)
          !---- Arguments ----!
          character(len=*), dimension(:),  intent(in)     :: lines
          integer,                         intent(in out) :: n_ini
          integer,                         intent(in)     :: n_end
          class(SpG_Type),                 intent(out)    :: SpG
          character(len=*), optional,      intent(in)     :: xyz_type
       End Subroutine Read_CFL_SpG_Lines

       Module Subroutine Read_CFL_SpG_FileTyp(cfl, SpG, xyz_type)
          !---- Arguments ----!
          type(File_Type),                 intent(in)     :: cfl     ! File_type object
          class(SpG_Type),                 intent(out)    :: SpG
          character(len=*), optional,      intent(in)     :: xyz_type
       End Subroutine Read_CFL_SpG_FileTyp

       Module Subroutine Read_kinfo(cfl, kinfo)
          !---- Arguments ----!
          type(File_Type),                 intent(in)     :: cfl     ! File_type object
          type(kvect_info_Type),           intent(out)    :: kinfo
       End Subroutine Read_kinfo

       Module Subroutine Readn_Set_XTal_CFL(file_dat,nlines,Cell,SpG,A,Type_Atm,CFrame,NPhase,Job_Info,xyz_type)
          !---- Arguments ----!
          character(len=*),dimension(:),intent(in)   :: file_dat
          integer,                      intent(in)   :: nlines
          class (Cell_Type),            intent(out)  :: Cell
          class (SPG_Type),             intent(out)  :: SpG
          Type (AtList_Type),           intent(out)  :: A
          character(len=*),             intent(in)   :: Type_Atm
          character(len=*),    optional,intent(in)   :: CFrame
          Integer,             optional,intent(in)   :: Nphase
          Type(Job_Info_type), optional,intent(out)  :: Job_Info
          character(len=*),    optional,intent(in)   :: xyz_type
       End Subroutine Readn_Set_XTal_CFL

      !Module Subroutine Read_Cif_Atom(lines,n_ini,n_end, At_List)
       !   !---- Arguments ----!
       !   character(len=*), dimension(:),   intent(in)      :: lines
       !   integer,                          intent(in out)  :: n_ini
       !   integer,                          intent(in)      :: n_end
       !   type (AtList_type),               intent(out)     :: At_List
       !End Subroutine Read_Cif_Atom
       !
       !Module Subroutine Read_Cif_Cell(lines, n_ini, n_end, Cell)
       !   !---- Arguments ----!
       !   character(len=*), dimension(:),  intent(in)     :: lines
       !   integer,                         intent(in out) :: n_ini
       !   integer,                         intent(in)     :: n_end
       !   class(Cell_Type),                intent(out)    :: Cell
       !End Subroutine Read_Cif_Cell
       !
       !Module Subroutine Read_Cif_ChemName(lines,N_ini,N_End,ChemName)
       !   !---- Arguments ----!
       !   character(len=*),  dimension(:), intent(in) :: lines
       !   integer,           intent(in out)           :: n_ini
       !   integer,           intent(in)               :: n_end
       !   character(len=*),  intent(out)              :: ChemName
       !End Subroutine Read_Cif_ChemName
       !
       !Module Subroutine Read_Cif_Z(lines, n_ini, n_end, Z)
       !   !---- Arguments ----!
       !   character(len=*), dimension(:),  intent(in)     :: lines
       !   integer,                         intent(in out) :: n_ini
       !   integer,                         intent(in)     :: n_end
       !   real(kind=cp),                   intent(out)    :: Z
       !End Subroutine Read_Cif_Z
       !
       !Module Subroutine Read_Cif_Wave(lines, n_ini, n_end, Wave)
       !   !---- Arguments ----!
       !   character(len=*), dimension(:),  intent(in)     :: lines
       !   integer,                         intent(in out) :: n_ini
       !   integer,                         intent(in)     :: n_end
       !   real(kind=cp),                   intent(out)    :: Wave
       !End Subroutine Read_Cif_Wave
       !
       !Module Subroutine Read_Cif_Cont(lines,N_Ini,N_End,N_Elem_Type,Elem_Type,N_Elem)
       !   !---- Arguments ----!
       !   character(len=*), dimension(:),      intent(in)      :: lines
       !   integer,                             intent(in out)  :: n_ini
       !   integer,                             intent(in)      :: n_end
       !   integer,                             intent(out)     :: n_elem_type
       !   character(len=*), dimension(:),      intent(out)     :: elem_type
       !   real(kind=cp), dimension(:),optional,intent(out)     :: n_elem
       !End Subroutine Read_Cif_Cont
       !
       !Module Subroutine Read_Cif_Pressure(lines,N_ini,N_End, P, SigP)
       !   !---- Arguments ----!
       !   character(len=*),  dimension(:), intent(in) :: lines
       !   integer,           intent(in out)           :: n_ini
       !   integer,           intent(in)               :: n_end
       !   real(kind=cp),     intent(out)              :: p
       !   real(kind=cp),     intent(out)              :: sigp
       !End Subroutine Read_Cif_Pressure
       !
       !Module Subroutine Read_Cif_Title(lines,N_Ini,N_End,Title)
       !   !---- Arguments ----!
       !   character(len=*),  dimension(:), intent(in) :: lines
       !   integer,           intent(in out)           :: n_ini
       !   integer,           intent(in)               :: n_end
       !   character(len=*),  intent(out)              :: title
       !End Subroutine Read_Cif_Title
       !
       !Module Subroutine Read_Cif_Temp(lines,N_Ini,N_End,T,SigT)
       !   !---- Arguments ----!
       !   character(len=*),  dimension(:), intent(in) :: lines
       !   integer,           intent(in out)           :: n_ini
       !   integer,           intent(in)               :: n_end
       !   real(kind=cp),     intent(out)              :: T
       !   real(kind=cp),     intent(out)              :: sigT
       !End Subroutine Read_Cif_Temp
       !
       !Module Subroutine Read_Cif_Hall(lines, N_Ini, N_End, Hall)
       !   !---- Arguments ----!
       !   character(len=*), dimension(:), intent(in) :: lines
       !   integer,          intent(in out)           :: n_ini
       !   integer,          intent(in)               :: n_end
       !   character(len=*), intent(out)              :: Hall
       !End Subroutine Read_Cif_Hall
       !
       !Module Subroutine Read_Cif_HM(lines, N_Ini, N_End, Spgr_Hm)
       !   !---- Arguments ----!
       !   character(len=*),  dimension(:), intent(in) :: lines
       !   integer,           intent(in out)           :: n_ini
       !   integer,           intent(in)               :: n_end
       !   character(len=*),  intent(out)              :: spgr_hm
       !End Subroutine Read_Cif_HM
       !
       !Module Subroutine Read_Cif_Symm(lines,N_Ini,N_End, N_Oper, Oper_Symm)
       !   !---- Arguments ----!
       !   character(len=*), dimension(:), intent(in)     :: lines
       !   integer,                        intent(in out) :: n_ini
       !   integer,                        intent(in)     :: n_end
       !   integer,                        intent(out)    :: n_oper
       !   character(len=*), dimension(:), intent(out)    :: oper_symm
       !End Subroutine Read_Cif_Symm
       !
       !Module Subroutine Write_Cif_Powder_Profile(filename)
       !   !---- Arguments ----!
       !   character(len=*), intent(in) :: filename
       !End Subroutine Write_Cif_Powder_Profile
       !
       !Module Subroutine Read_Shx_Atom(lines, n_ini, n_end, n_fvar, fvar, elem_type, cell, At_List)
       !   !---- Arguments ----!
       !   character(len=*), dimension(:), intent(in)      :: lines
       !   integer,                        intent(in out)  :: n_ini
       !   integer,                        intent(in)      :: n_end
       !   integer,                        intent(in)      :: n_fvar
       !   real(kind=cp), dimension(:),    intent(in)      :: fvar
       !   character(len=*), dimension(:), intent(in)      :: elem_type
       !   class(Cell_G_Type),             intent(in)      :: Cell
       !   type (AtList_type),             intent(out)     :: At_List
       !End Subroutine Read_Shx_Atom
       !
       !Module Subroutine Read_Shx_Cell(lines, n_ini, n_end, Cell)
       !   !---- Arguments ----!
       !   character(len=*), dimension(:),  intent(in)     :: lines
       !   integer,                         intent(in out) :: n_ini
       !   integer,                         intent(in)     :: n_end
       !   class(Cell_Type),                intent(out)    :: Cell
       !End Subroutine Read_Shx_Cell
       !
       !Module Subroutine Read_Shx_Wave(lines, n_ini, n_end, Wave)
       !   !---- Arguments ----!
       !   character(len=*), dimension(:),  intent(in)     :: lines
       !   integer,                         intent(in out) :: n_ini
       !   integer,                         intent(in)     :: n_end
       !   real(kind=cp),                   intent(out)    :: Wave
       !End Subroutine Read_Shx_Wave
       !
       !Module Subroutine Read_Shx_Z(lines, n_ini, n_end, Z)
       !   !---- Arguments ----!
       !   character(len=*), dimension(:),  intent(in)     :: lines
       !   integer,                         intent(in out) :: n_ini
       !   integer,                         intent(in)     :: n_end
       !   real(kind=cp),                   intent(out)    :: Z
       !End Subroutine Read_Shx_Z
       !
       !Module Subroutine Read_Shx_Fvar(Lines,n_ini,n_end, n_fvar, fvar)
       !   !---- Arguments ----!
       !   character(len=*), dimension(:), intent(in)    :: lines
       !   integer,                        intent(in out):: n_ini
       !   integer,                        intent(in)    :: n_end
       !   integer,                        intent(out)   :: n_fvar
       !   real(kind=cp), dimension(:),    intent(out)   :: fvar
       !End Subroutine Read_Shx_Fvar
       !
       !Module Subroutine Read_Shx_Titl(lines,n_ini,n_end,Title)
       !   !---- Arguments ----!
       !   character(len=*),dimension(:), intent(in)     :: lines
       !   integer,                       intent(in out) :: n_ini
       !   integer,                       intent(in)     :: n_end
       !   character(len=*),              intent(out)    :: title
       !End Subroutine Read_Shx_Titl
       !
       !Module Subroutine Read_Shx_Latt(lines,n_ini,n_end,latt)
       !   !---- Arguments ----!
       !   character(len=*), dimension(:), intent(in) :: lines
       !   integer,           intent(in out)          :: n_ini
       !   integer,           intent(in)              :: n_end
       !   integer,           intent(out)             :: latt
       !End Subroutine Read_Shx_Latt
       !
       !Module Subroutine Read_Shx_Cont(lines,n_ini,n_end, n_elem_type, elem_type, n_elem)
       !   !---- Arguments ----!
       !   character(len=*), dimension(:),           intent(in)      :: lines
       !   integer,                                  intent(in out)  :: n_ini
       !   integer,                                  intent(in)      :: n_end
       !   integer,                                  intent(out)     :: n_elem_type
       !   character(len=*), dimension(:),           intent(out)     :: elem_type
       !   real(kind=cp),    dimension(:), optional, intent(out)     :: n_elem
       !End Subroutine Read_Shx_Cont
       !
       !Module Subroutine Read_Shx_Symm(lines,n_ini,n_end,n_oper,oper_symm)
       !   !---- Arguments ----!
       !   character(len=*), dimension(:), intent(in) :: lines
       !   integer,          intent(in out)           :: n_ini
       !   integer,          intent(in)               :: n_end
       !   integer,          intent(out)              :: n_oper
       !   character(len=*), dimension(:),intent(out) :: oper_symm
       !End Subroutine Read_Shx_Symm
       !
       !Module Subroutine Write_Shx_Template(filename,code,title,lambda,z,cell,spg,At_List)
       !   !---- Arguments ----!
       !   character(len=*),        intent(in) :: filename
       !   integer,                 intent(in) :: code
       !   character(len=*),        intent(in) :: title
       !   real(kind=cp),           intent(in) :: lambda
       !   integer,                 intent(in) :: z
       !   class(cell_Type),        intent(in) :: cell
       !   class(SpG_Type),         intent(in) :: SpG
       !   type(atlist_type),       intent(in) :: at_List
       !End Subroutine Write_Shx_Template

    End Interface

    Contains

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
    Subroutine Readn_Set_Xtal_Structure_Split(filenam,Cell,SpG,A,Type_Atm,Mode,Iphase,Job_Info,file_list,CFrame)
       !---- Arguments ----!
       character(len=*),             intent( in)     :: filenam
       class (Cell_Type),            intent(out)     :: Cell
       class (SpG_Type),             intent(out)     :: SpG
       type (Atlist_type),           intent(out)     :: A
       Character(len=*),             intent( in)     :: Type_Atm
       Character(len=*),    optional,intent( in)     :: Mode
       Integer,             optional,intent( in)     :: Iphase
       Type(Job_Info_type), optional,intent(out)     :: Job_Info
       Type(File_type),     optional,intent(in out)  :: file_list
       Character(len=*),    optional,intent( in)     :: CFrame

       !---- Local variables -----!
       character(len=132), allocatable, dimension(:) :: file_dat
       character(len=3)                              :: modec
       integer                                       :: nlines,i

       call Clear_Error()

       nlines=0
       if (present(file_list)) nlines=file_list%nlines

       !---- Number of Lines in the input file ----!
       if(nlines == 0) then
           nlines=Number_Lines(trim(filenam))
           if (nlines==0) then
              Err_CFML%Ierr=1
              Err_CFML%Msg="The file "//trim(filenam)//" contains nothing"
              return
           else
              if (allocated(file_dat)) deallocate( file_dat)
              allocate(file_dat(nlines))
              call reading_Lines(trim(filenam),nlines,file_dat)
           end if
           if (present(file_list)) then
              file_list%nlines=nlines
              file_list%Fname=trim(filenam)
              if (allocated(file_list%line)) deallocate(file_list%line)
              allocate(file_list%line(nlines))
              do i=1,nlines
                file_list%line(i)%str=trim(file_dat(i))
              end do
           end if
       else
           if (allocated(file_dat)) deallocate( file_dat)
           allocate(file_dat(nlines))
           do i=1,nlines
             file_dat(i)=file_list%line(i)%str
           end do
       end if


       !---- Define the type of file: CIF, CFL, RES,... ----!
       modec=" "
       if (present(mode)) modec=l_case(mode(1:3))

       select case(modec)

           !case("cif")
           !   if (present(iphase)) then
           !      if(present(CFrame)) then
           !        call readn_set_xtal_cif(file_dat,nlines,Cell,Spg, A,CFrame,NPhase=IPhase)
           !      else
           !        call readn_set_xtal_cif(file_dat,nlines,Cell,Spg, A,NPhase=IPhase)
           !      end if
           !   else
           !      if(present(CFrame)) then
           !        call readn_set_xtal_cif(file_dat,nlines,Cell,Spg,A,CFrame)
           !      else
           !        call readn_set_xtal_cif(file_dat,nlines,Cell,Spg,A)
           !      end if
           !   end if
           !
           !case("pcr")
           !   if (present(iphase)) then
           !      if(present(CFrame)) then
           !        call readn_set_xtal_pcr(file_dat,nlines,Cell,Spg, A,CFrame,NPhase=IPhase)
           !      else
           !        call readn_set_xtal_pcr(file_dat,nlines,Cell,Spg, A,NPhase=IPhase)
           !      end if
           !   else
           !      if(present(CFrame)) then
           !        call readn_set_xtal_pcr(file_dat,nlines,Cell,Spg,A,CFrame)
           !      else
           !        call readn_set_xtal_pcr(file_dat,nlines,Cell,Spg,A)
           !      end if
           !   end if
           !
           !case("shx")
           !   if(present(CFrame)) then
           !     call readn_set_xtal_shx(file_dat,nlines,Cell,Spg,A,CFrame)
           !   else
           !     call readn_set_xtal_shx(file_dat,nlines,Cell,Spg,A)
           !   end if

           case default
              !---- CFL Format ----!
              if (present(Job_Info)) then
                 if (present(iphase)) then
                    if(present(CFrame)) then
                      call readn_set_xtal_cfl(file_dat,nlines,Cell,Spg,A,Type_Atm,CFrame,NPhase=IPhase,Job_Info=Job_Info)
                    else
                      call readn_set_xtal_cfl(file_dat,nlines,Cell,Spg,A,Type_Atm,NPhase=IPhase,Job_Info=Job_Info)
                    end if
                 else
                    if(present(CFrame)) then
                      call readn_set_xtal_cfl(file_dat,nlines,Cell,Spg,A,Type_Atm,CFrame,Job_Info=Job_Info)
                    else
                      call readn_set_xtal_cfl(file_dat,nlines,Cell,Spg,A,Type_Atm,Job_Info=Job_Info)
                    end if
                 end if
              else
                 if (present(iphase)) then
                    if(present(CFrame)) then
                      call readn_set_xtal_cfl(file_dat,nlines,Cell,Spg,A,Type_Atm,CFrame,NPhase=IPhase)
                    else
                      call readn_set_xtal_cfl(file_dat,nlines,Cell,Spg,A,Type_Atm,NPhase=IPhase)
                    end if
                 else
                    if(present(CFrame)) then
                      call readn_set_xtal_cfl(file_dat,nlines,Cell,Spg,A,Type_Atm,CFrame)
                    else
                      call readn_set_xtal_cfl(file_dat,nlines,Cell,Spg,A,Type_Atm)
                    end if
                 end if
              end if

       end select

    End Subroutine Readn_Set_Xtal_Structure_Split

End Module CFML_IOForm

