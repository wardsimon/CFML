!!-------------------------------------------------------
!!---- Crystallographic Fortran Modules Library (CrysFML)
!!-------------------------------------------------------
!!---- The CrysFML project is distributed under LGPL. In agreement with the
!!---- Intergovernmental Convention of the ILL, this software cannot be used
!!---- in military applications.
!!----
!!---- Copyright (C) 1999-2018  Institut Laue-Langevin (ILL), Grenoble, FRANCE
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
!!---- MODULE: CFML_Diffraction_Patterns
!!----   INFO: Diffraction Patterns Information
!!----
!!
 Module CFML_Diffraction_Patterns
    !---- Use Modules ----!
    Use CFML_GlobalDeps,       only : cp, ops_sep, err_cfml, clear_error
    Use CFML_Math_General,     only : spline, splint, locate, second_derivative
    use CFML_String_Utilities, only : FindFmt,  Init_FindFmt , u_case, get_word, get_num, Get_NumStd

    implicit none

    private

    !---- List of public functions ----!
    public ::  calc_fwhm_peak

    !---- List of public subroutines ----!
    public ::  Allocate_Diffraction_Pattern, Deallocate_Diffraction_Pattern, &
               Read_Pattern
               !Calc_Background, Del_NoisyPoints, &
               !Read_Background_File, Read_Pattern,      &
               !Write_Pattern_XYSig,&
               !Write_Pattern_FreeFormat, Add_Diffraction_Patterns,      &
               !Write_Pattern_INSTRM5

    !---- Definitions ----!
    
    !!----
    !!---- TYPE :: DIFFPAT_TYPE
    !!----
    Type, public :: DiffPat_Type
       character(len=180)                        :: Title  =" "        ! Indentification/ Title
       character(len=20)                         :: KindRad=" "        ! Type of Radiation
       character(len=20)                         :: ScatVar=" "        ! 2Theta, TOF, Q, s, d-spacing, SinTL/L,...
       real(kind=cp)                             :: xmin   =0.0_cp     ! Maximum and Minimum values for X and Y 
       real(kind=cp)                             :: xmax   =0.0_cp  
       real(kind=cp)                             :: ymin   =0.0_cp  
       real(kind=cp)                             :: ymax   =0.0_cp  
       integer                                   :: NPts   =0          ! Number of Points
       logical                                   :: SigVar =.true.     ! .True. for sigma values / .False. for variance
       real(kind=cp), dimension(5)               :: Wave   =0.0_cp     ! Wave1, Wave2, Dtt1, Dtt2,....
       real(kind=cp), allocatable, dimension (:) :: x 
       real(kind=cp), allocatable, dimension (:) :: y 
       real(kind=cp), allocatable, dimension (:) :: sigma 
    End Type DiffPat_Type
    
    !!----
    !!---- TYPE :: DiffPat_E_Type
    !!----
    Type, public, extends (DiffPat_Type) ::  DiffPat_E_Type
       character(len=30)           :: Instr   =" "                    ! Instrument name
       character(len=80)           :: Filename=" "                    ! Filename
       character(len=512)          :: FilePath=" "                    ! Path
       real(kind=cp)               :: Monitor =0.0_cp
       real(kind=cp)               :: Norm_Mon=0.0_cp
       real(kind=cp)               :: Col_Time=0.0_cp
       real(kind=cp)               :: Tsample =298.0_cp               ! Sample temperature
       real(kind=cp)               :: Tset    =0.0_cp                 ! Wished temperature
       logical                     :: CT_Step =.false.                ! Constant step
       logical                     :: al_x    =.false.                ! Flags for variable allocations 
       logical                     :: al_y    =.false.                ! 
       logical                     :: al_sigma=.false.                ! 
       logical                     :: al_ycalc=.false.                ! 
       logical                     :: al_bgr  =.false.                ! 
       logical                     :: al_istat=.false.                ! 
       real(kind=cp), allocatable, dimension (:) :: ycalc             ! Caculated intensity 
       real(kind=cp), allocatable, dimension (:) :: bgr               ! Background
       integer,       allocatable, dimension (:) :: istat             ! Information about point "i"
       integer,       allocatable, dimension (:) :: ND                ! Number of Detectors contributing to point "i" 
    End Type DiffPat_E_Type

    !!----
    !!---- TYPE :: DIFFPAT_G_TYPE
    !!----
    Type, public, extends (DiffPat_E_Type) :: DiffPat_G_Type
       character(len=40)                           :: Legend_X=" "     !x-axis legend, eg. "Lambda (Angstroms)"
       character(len=40)                           :: Legend_Y=" "     !y-axis legend, eg. "Intensity (arb. units)"
       logical                                     :: gy      =.false. ! Flags for Graphics 
       logical                                     :: gycalc  =.false.
       logical                                     :: gsigma  =.false.
       logical                                     :: gbgr    =.false.
    End Type DiffPat_G_Type

    
    !---- Interfaces - Overlap ----!
    Interface Read_Pattern
       Module procedure Read_Pattern_Mult
       Module procedure Read_Pattern_One
    End Interface
    
    Interface
       Module Function Calc_FWHM_Peak(Pat, Xi, Yi, Ybi, RLim) Result(v)
          !---- Arguments ----!
          class(DiffPat_Type),       intent(in) :: Pat      ! Pattern object
          real(kind=cp),             intent(in) :: Xi       ! (Xi,Yi) for point i
          real(kind=cp),             intent(in) :: Yi       ! 
          real(kind=cp),             intent(in) :: Ybi      ! Background at point i
          real(kind=cp),optional,    intent(in) :: RLim     ! Limit range in X units to search the point
          real(kind=cp)                         :: V
       End Function Calc_FWHM_Peak   
       
       Module Subroutine Calc_BackGround(Pat, Ncyc, Np, Xmin, Xmax)
          !---- Arguments ----!
          class(DiffPat_E_Type),     intent(in out) :: Pat        ! Pattern object
          integer,                   intent(in)     :: NCyc       ! Number of Cycles to apply
          integer,                   intent(in)     :: Np         ! Number of extension points at L and R.
          real(kind=cp), optional,   intent(in)     :: Xmin       ! Min, max values in X
          real(kind=cp), optional,   intent(in)     :: Xmax
       End Subroutine Calc_BackGround
       
       Module Subroutine Del_NoisyPoints(Pat, NoisyP, FileInfo)
          !---- Arguments ----!
          class(DiffPat_Type),  intent(in out) :: Pat        ! Pattern object
          integer,              intent(out)    :: NoisyP     ! Noisy points
          logical, optional,    intent(in)     :: FileInfo   ! .true. For create an information file
       End Subroutine Del_NoisyPoints
       
       Module Subroutine Read_Background_File(Bck_File, Bck_Mode, Pat)
          !---- Arguments ----!
          character(len=*),         intent(in   )    :: bck_file      ! Path+Filename of Background file
          character(len=*),         intent(in   )    :: bck_mode
          class(DiffPat_E_Type),  intent(in out)   :: Pat
       End Subroutine Read_Background_File
       
       Module Subroutine Set_Background_Poly(Pat, Bkpos, Bckx, N)
          !---- Arguments ----!
          class(DiffPat_E_Type),         intent(in out) :: Pat
          real (kind=cp),                intent(in    ) :: bkpos
          real (kind=cp), dimension(:),  intent(in    ) :: bckx
          integer,                       intent(in    ) :: n
       End Subroutine Set_Background_Poly
       
       Module Subroutine Set_Background_Inter(Pat, Bcky, Bckx, N)
          !---- Arguments ----!
          class(DiffPat_E_Type),         intent(in out) :: Pat
          real (kind=cp), dimension(:),  intent(in out) :: bcky
          real (kind=cp), dimension(:),  intent(in out) :: bckx
          integer,                       intent(in    ) :: n
       End Subroutine Set_Background_Inter
       
       Module Subroutine Write_Pattern_XYSig(Filename,Pat,excl,xmin,xmax)
          !---- Arguments ----!
          character(len=*),               intent(in) :: filename     ! Path+Filename
          class(DiffPat_Type),            intent(in) :: Pat          ! Pattern object
          logical, dimension(:),optional, intent(in) :: excl         ! Exclusion zones
          real(kind=cp),        optional, intent(in) :: xmin         ! Limits
          real(kind=cp),        optional, intent(in) :: xmax    
       End Subroutine Write_Pattern_XYSig
       
       Module Subroutine Write_Pattern_FreeFormat(Filename,Pat,excl,xmin,xmax)
          !---- Arguments ----!
          character (len=*),               intent(in)     :: Filename      ! Path+Filename
          class(DiffPat_Type),             intent(in out) :: Pat           ! Pat object
          logical, dimension(:),optional,  intent(in)     :: excl
          real,                 optional,  intent(in)     :: xmin
          real,                 optional,  intent(in)     :: xmax
       End Subroutine Write_Pattern_FreeFormat   
       
       Module Subroutine Write_Pattern_INSTRM5(Filename,Pat,excl,xmin,xmax,var)
          !---- Arguments ----!
          character(len=*),                intent(in)     :: Filename
          class(DiffPat_E_Type),           intent(in out) :: Pat
          logical, dimension(:),optional,  intent(in)     :: excl
          real,                 optional,  intent(in)     :: xmin
          real,                 optional,  intent(in)     :: xmax
          character(len=*),     optional,  intent(in)     :: var
       End Subroutine Write_Pattern_INSTRM5   
       
       Module Subroutine Read_Pattern_Panalytical_CSV(Filename,Pat)
          !---- Arguments ----!
          character (len=*),   intent(in)  :: Filename      ! Path+Filename
          class(DiffPat_type), intent(out) :: pat
       End Subroutine Read_Pattern_Panalytical_CSV 
       
       Module Subroutine Read_Pattern_Panalytical_JCP(Filename, Pat)
          !---- Arguments ----!
          character(len=*),    intent(in)  :: Filename      ! Path+Filename
          class(DiffPat_type), intent(out) :: pat 
       End Subroutine Read_Pattern_Panalytical_JCP  
       
       Module Subroutine Read_Pattern_Panalytical_UDF(Filename, Pat)
          !---- Arguments ----!
          character(len=*),    intent(in)  :: Filename      ! Path+Filename
          class(DiffPat_type), intent(out) :: pat
       End Subroutine Read_Pattern_Panalytical_UDF 
       
       Module Subroutine Read_Pattern_Panalytical_XRDML(Filename,Pat)
          !---- Arguments ----!
          character(len=*),    intent(in)  :: Filename      ! Path+Filename
          class(DiffPat_type), intent(out) :: pat 
       End Subroutine Read_Pattern_Panalytical_XRDML
       
       Module Subroutine Read_Pattern_D1A_D2B(Filename,Pat)
          !---- Arguments ----!
          character(len=*),      intent(in)  :: Filename      ! Path+Filename
          class(DiffPat_E_Type), intent(out) :: Pat
       End Subroutine Read_Pattern_D1A_D2B 
       
       Module Subroutine Read_Pattern_D1A_D2B_OLD(Filename,Pat)
          !---- Arguments ----!
          character(len=*),      intent(in)  :: Filename      ! Path+Filename
          class(DiffPat_E_Type), intent(out) :: Pat
       End Subroutine Read_Pattern_D1A_D2B_OLD  
       
       Module Subroutine Read_Pattern_D1B_D20(Filename,Pat)
          !---- Arguments ----!
          character(len=*),      intent(in)  :: Filename      ! Path+Filename
          class(DiffPat_E_Type), intent(out) :: Pat
       End Subroutine Read_Pattern_D1B_D20
       
       Module Subroutine Read_Pattern_DMC(Filename,Pat)
          !---- Arguments ----!
          character(len=*),        intent(in)  :: Filename      ! Path+Filename
          class(DiffPat_Type),     intent(out) :: Pat
       End Subroutine Read_Pattern_DMC
       
       Module Subroutine Read_Pattern_Isis_M(Filename, VPat, NPat)
          !---- Arguments ----!
          Character(len=*),                   intent(in)  :: Filename
          class(DiffPat_Type), dimension(:),  intent(out) :: VPat
          integer,                            intent(out) :: Npat          ! Number of Patterns readed
       End Subroutine Read_Pattern_Isis_M
       
       Module Subroutine Read_Pattern_G41(Filename,Pat)
          !---- Arguments ----!
          character(len=*),      intent(in)  :: Filename      ! Path+Filename
          class(DiffPat_E_Type), intent(out) :: Pat
       End Subroutine Read_Pattern_G41
       
       Module Subroutine Read_Pattern_CIF(Filename,Pat)
          !---- Arguments ----!
          character(len=*),    intent(in)  :: Filename      ! Path+Filename
          class(DiffPat_Type), intent(out) :: Pat
       End Subroutine Read_Pattern_CIF  
       
       Module Subroutine Read_Pattern_NLS(Filename,Pat)
          !---- Arguments ----!
          character(len=*),    intent(in)  :: Filename      ! Path+Filename
          class(DiffPat_Type), intent(out) :: Pat 
       End Subroutine Read_Pattern_NLS  
       
       Module Subroutine Read_Pattern_GSAS(Filename, Pat, mode)
          !---- Arguments ----!
          character(len=*),           intent(in)  :: Filename      ! Path+Filename
          class(DiffPat_Type),        intent(out) :: Pat
          character(len=*), optional, intent(in)  :: Mode 
       End Subroutine Read_Pattern_GSAS 
       
       Module Subroutine Read_Pattern_Socabim(Filename,Pat)
          !---- Arguments ----!
          character(len=*),    intent(in)  :: Filename      ! Path+Filename
          class(DiffPat_Type), intent(out) :: Pat 
       End Subroutine Read_Pattern_Socabim   
       
       Module Subroutine Read_Pattern_Free(Filename,Pat,ext)
          !---- Arguments ----!
          character(len=*),           intent(in)  :: Filename      ! Path+Filename
          class(DiffPat_Type),        intent(out) :: Pat
          character(len=*), optional, intent(in)  :: ext
       End Subroutine Read_Pattern_Free  
       
       Module Subroutine Read_Pattern_Time_Variable(Filename, Pat)
          !---- Arguments ----!
          character(len=*),    intent(in)  :: Filename      ! Path+Filename
          class(DiffPat_Type), intent(out) :: Pat 
       End Subroutine Read_Pattern_Time_Variable
       
       Module Subroutine Read_Pattern_XYSigma(Filename, Pat, PDF, Header)
          !---- Arguments ----!
          character(len=*),           intent(in)  :: Filename      ! Path+Filename
          class(DiffPat_Type),        intent(out) :: Pat
          logical,          optional, intent(in)  :: PDF
          character(len=*), optional, intent (out):: Header
       End Subroutine Read_Pattern_XYSigma   
       
       Module Subroutine Add_Patterns(Patterns, N, Active, Pat, VNorm)
          !---- Arguments ----!
          class(DiffPat_Type), dimension(:), intent(in)  :: Patterns
          integer,                           intent(in)  :: N
          logical,             dimension(:), intent(in)  :: Active
          class(DiffPat_Type),               intent(out) :: Pat
          real(kind=cp), optional,           intent(in)  :: VNorm
       End Subroutine Add_Patterns
       
    End Interface

 Contains
    !!----
    !!---- SUBROUTINE ALLOCATE_DIFFRACTION_PATTERN
    !!----
    !!----    Allocate the Part of Diffractions Patterns
    !!----
    !!---- Update: December - 2005
    !!
    Subroutine Allocate_Diffraction_Pattern(Pat,Npts)
       !---- Arguments ----!
       class(DiffPat_Type), intent (in out) :: Pat      ! Pattern object
       integer, optional,   intent (in)     :: npts     ! Number of points

       !---- Local variables ----!
       integer :: n

       !> Init
       call clear_error()
        
       N=Pat%npts
       if (present(npts)) N=Npts
       
       if (n <= 0) then
          err_CFML%state=.true.
          err_CFML%Flag=2
          err_CFML%Msg=" Attempt to allocate a Diffraction Pattern with 0-dimension "
          return
       end if
       
       !> Allocating
       Pat%npts=N
       
       !> Class (DiffPat_Type)
       if (allocated(pat%x) ) deallocate(pat%x)
       allocate(pat%x(n))
       pat%x=0.0_cp
       
       if (allocated(pat%y) ) deallocate(pat%y)
       allocate(pat%y(n))
       pat%y=0.0_cp
       
       if (allocated(pat%sigma) ) deallocate(pat%sigma)
       allocate(pat%sigma(n))
       pat%sigma=0.0_cp
       
       !> class (DiffPat_E_Type)
       select type(Pat)
          class is (DiffPat_E_Type)
             if (allocated(pat%bgr) ) deallocate(pat%bgr)
             allocate(pat%bgr(n))
             pat%bgr=0.0_cp
             
             if (allocated(pat%ycalc) ) deallocate(pat%ycalc)
             allocate(pat%ycalc(n))
             pat%ycalc=0.0_cp
             
             if (allocated(pat%istat) ) deallocate(pat%istat)
             allocate(pat%istat(n))
             pat%istat=1
             
             if (allocated(pat%nd) ) deallocate(pat%nd)
             allocate(pat%nd(n))
             pat%nd=0
            
             Pat%al_x=.true.
             Pat%al_y=.true.
             Pat%al_sigma=.true.
             Pat%al_ycalc=.true.
             Pat%al_bgr  =.true.
             Pat%al_istat=.true.
            
       end select
       
       !> class (DiffPat_G_Type)
       select type(Pat)
          type is (DiffPat_G_Type)
             Pat%gy    =.false.
             Pat%gycalc=.false.
             Pat%gsigma=.false.
             Pat%gbgr  =.false.
             
       end select

       return
    End Subroutine Allocate_Diffraction_Pattern
    
    !!----
    !!---- Subroutine Deallocate_Diffraction_Pattern
    !!----
    !!----    De-Allocate components of the object "pat", of type Diffraction_Pattern_Type
    !!----    depending on the value of the MODE string. At present the following MODE
    !!----    values are available:
    !!----      "DATA " -> x,y remain allocated                  (purge sigma,ycalc,bgr,istat)
    !!----      "DATAS" -> x,y,sigma remain allocated            (purge ycalc,bgr,istat)
    !!----      "RIETV" -> x,y,sigma,ycalc,bgr remain allocated  (purge istat)
    !!----      "GRAPH" -> x,y,sigma,istat remain allocated      (purge ycalc, bgr)
    !!----      "PRF  " -> x,y,sigma,ycalc,bgr,istat, everything remains allocated (changed w.r.t. previous version)
    !!----
    !!----
    !!----
    !!---- Update: December - 2005
    !!---- Updated: December - 2005
    !!
    Subroutine Deallocate_Diffraction_Pattern(Pat,Mode)
       !---- Arguments ----!
       class(DiffPat_Type), intent (in out) :: Pat       ! Pattern object
       character(len=*),    intent (in)     :: Mode      ! Type of the Deallocation

       !---- Local Variables ----!
       character(len=5) :: car
       
       car=adjustl(u_case(mode))
       
       select case (trim(car))
          case ("DATA")
             select type(Pat)
                type is (DiffPat_Type)
                   if (allocated(Pat%sigma)) deallocate(Pat%sigma) 
                   
                class is (DiffPat_E_Type)
                   if (allocated(Pat%sigma)) deallocate(Pat%sigma) 
                   if (allocated(Pat%ycalc)) deallocate(Pat%ycalc) 
                   if (allocated(Pat%bgr))   deallocate(Pat%bgr) 
                   if (allocated(Pat%istat)) deallocate(Pat%istat) 
                   !if (allocated(Pat%ND))    deallocate(Pat%ND) 
                   
                   Pat%al_sigma=.false.
                   Pat%al_ycalc=.false.
                   Pat%al_bgr  =.false.
                   Pat%al_istat=.false.
             end select
            
          case ("DATAS")
             select type(Pat)
                class is (DiffPat_E_Type)
                   if (allocated(Pat%ycalc)) deallocate(Pat%ycalc) 
                   if (allocated(Pat%bgr))   deallocate(Pat%bgr) 
                   if (allocated(Pat%istat)) deallocate(Pat%istat) 
                   !if (allocated(Pat%ND))    deallocate(Pat%ND)
                   
                   Pat%al_ycalc=.false.
                   Pat%al_bgr  =.false.
                   Pat%al_istat=.false.
             end select
            
          case ("RIETV")
             select type(Pat)
                class is (DiffPat_E_Type)  
                   if (allocated(Pat%istat)) deallocate(Pat%istat)
                   
                   Pat%al_istat=.false.
             end select
            
          case ("GRAPH")
             select type(Pat)
                type is (DiffPat_E_Type) 
                   if (allocated(Pat%ycalc)) deallocate(Pat%ycalc) 
                   if (allocated(Pat%bgr))   deallocate(Pat%bgr)
                   
                   Pat%al_ycalc=.false.
                   Pat%al_bgr  =.false.
                    
                type is (DiffPat_G_Type)
                   if (allocated(Pat%ycalc)) deallocate(Pat%ycalc) 
                   if (allocated(Pat%bgr))   deallocate(Pat%bgr)
                   
                   Pat%al_ycalc=.false.
                   Pat%al_bgr  =.false.
                   
                   Pat%gycalc=.false.
                   Pat%gbgr  =.false.  
             end select
            
          case ("PRF")
             ! Nothing to do
       end select 

       return
    End Subroutine Deallocate_Diffraction_Pattern
    
    
    !!--++
    !!--++ Subroutine Read_Pattern_Mult
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Read patterns from a Filename
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Read_Pattern_Mult(filename, Patts, NPats, mode)
       !---- Arguments ----!
       character(len=*),                   intent (in)      :: Filename        ! Path + Name of File containing Patterns
       class(DiffPat_Type), dimension(:),  intent (out)     :: Patts           ! Pattern objects
       integer,                            intent (in out)  :: NPats           ! In: Number of Patterns to read, Out: Number of Patterns readed
       character(len=*), optional,         intent (in)      :: Mode            ! Mode: ISIS, GSAS, XYSIGMA

       !---- Local variables ----!
       logical :: esta
       integer :: i_dat, ier,i

       !> Init
       call clear_error()

       !> Mode option
       if (present(mode)) then
          select case (u_case(mode))
              case ("XYSIGMA")
                 !   call  Read_Pattern_xysigma_m(dif_pat,npat)

              case ("ISIS")
                 call Read_Pattern_isis_m(trim(filename),Patts,NPats)
                 !dif_pat%diff_kind = "neutrons_tof"
                 !dif_pat%scat_var =  "TOF"
                 !dif_pat%xax_text =  "TOF (micro-seconds)"
                 !dif_pat%yax_text =  "Intensity (arb. units)"
                 !dif_pat%instr  = " 14  - "//mode

              case ("GSAS")
                 !   call Read_Pattern_gsas_m(dif_pat,npat)      ! GSAS file

              case default
                 Err_CFML%state=.true.
                 Err_CFML%Flag=2
                 Err_CFML%Msg="Invalid Mode for Read Multi Patterns. Please, check it!"
          end select
          
       else   
          
       end if

       return
    End Subroutine Read_Pattern_Mult

    
    !!--++
    !!--++ Subroutine Read_Pattern_One
    !!--++
    !!--++    Read one pattern from a Filename. If sig is present the content of Dif_Pat%sigma
    !!--++    is the true sigma not the variance.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Read_Pattern_One(Filename, Pat, Mode, Sig, Header)
       !---- Arguments ----!
       character(len=*),            intent (in)      :: Filename       ! Path Name of File
       class(DiffPat_Type),         intent (in out)  :: Pat            ! Pattern object
       character(len=*), optional,  intent (in)      :: mode           ! Mode
       logical,          optional,  intent (in)      :: sig            ! Sigma or Variance
       character(len=*), optional,  intent (out)     :: header         ! Header

       !---- Local Variables ----!
       character(len=6)                               :: extdat !extension of panalytical file
       character(len=4)                               :: tofn
       character(len=12)                              :: typef, modem !extension of panalytical file
       logical                                        :: esta,gr
       integer                                        :: i, i_dat,ier

       !> Init
       call clear_error()


       Typef="DEFAULT"
       TofN=" "
       
       if (present(mode)) then
          Typef=trim(u_case(mode))
          if (Typef(1:7) == "GSASTOF") then
             if (len_trim(mode) > 7) then
                tofn="TOF"//mode(8:8)
                typef="GSASTOF"
             else
                tofn="TOF"
             end if
          end if
       end if     
       
       select case (trim(typef))
          case ("CIF")
             call Read_Pattern_CIF(trim(filename),Pat)
             if (err_CFML%state) return
             
             select type(Pat)
                type is (DiffPat_E_Type)
                   pat%ct_step = .false.
                   pat%instr  = "  XY  - "//trim(typef)
                  
                type is (DiffPat_G_Type)
                   pat%ct_step = .false.
                   pat%instr  =  "  XY  - "//trim(typef)
                   pat%legend_Y ="Intensity (arb. units)"
             end select
            
          case ("D1B" , "D20")  
             pat%kindrad = "Neutrons_CW"
             pat%scatvar = "2Theta"
             select type(Pat)
                class is (DiffPat_E_Type) 
                   call Read_Pattern_D1B_D20(trim(filename), Pat)
                   if (err_CFML%state) return
                   
                   pat%ct_step = .true.
                   pat%instr=trim(Typef)
             
                   select type(Pat)
                      type is (DiffPat_G_Type)
                         pat%legend_X ="2Theta(degrees)"
                         pat%legend_Y ="Intensity (arb. units)"
                   end select
             end select      
            
          case ("NLS")
             call Read_Pattern_NLS(trim(filename), Pat)
             if (err_CFML%state) return
             
             pat%kindrad = "XRay_CW"
             pat%scatvar = "2Theta"
             
             select type(Pat)
                type is (DiffPat_E_Type)
                   pat%ct_step = .true.
                   pat%instr=trim(Typef)
                  
                type is (DiffPat_G_Type)
                   pat%ct_step = .true.
                   pat%instr=trim(Typef)
                   
                   pat%legend_X ="2Theta(degrees)"
                   pat%legend_Y ="Intensity (arb. units)"
             end select
            
          case ("G41") 
             pat%kindrad = "Neutrons_CW"
             pat%scatvar = "2Theta"
             select type (Pat)
                class is (DiffPat_E_Type)
                   call Read_Pattern_G41(trim(filename), Pat)
                   if (err_CFML%state) return
             
                   pat%ct_step = .true.
                   pat%instr=trim(Typef)
                   
                   select type(pat)
                      type  is (DiffPat_G_Type)
                         pat%legend_X ="2Theta(degrees)"
                         pat%legend_Y ="Intensity (arb. units)"
                   end select      
             end select 
            
          case ("D1A","D2B","3T2","G42")  
             pat%kindrad = "Neutrons_CW"
             pat%scatvar = "2Theta"
             select type(Pat)
                class is (DiffPat_E_Type)
                   call Read_Pattern_D1A_D2B(trim(filename), Pat)
                   if (err_CFML%state) return
             
                   pat%ct_step = .true.
                   pat%instr=trim(Typef)
                   
                   select type(pat)
                      type is (DiffPat_G_Type)
                         pat%legend_X ="2Theta(degrees)"
                         pat%legend_Y ="Intensity (arb. units)"
                   end select      
             end select
            
          case ("D1AOLD", "D2BOLD","OLDD1A", "OLDD2B") 
             pat%kindrad = "Neutrons_CW"
             pat%scatvar = "2Theta"
             select type(Pat)
                class is (DiffPat_E_Type)
                   call Read_Pattern_D1A_D2B_OLD(trim(filename), Pat)
                   if (err_CFML%state) return
             
                   pat%ct_step = .true.
                   pat%instr=trim(Typef)
                   
                   select type(pat)
                      type is (DiffPat_G_Type)
                         pat%legend_X ="2Theta(degrees)"
                         pat%legend_Y ="Intensity (arb. units)"
                   end select      
             end select
            
          case ("DMC","HRPT") 
             pat%kindrad = "Neutrons_CW"
             pat%scatvar = "2Theta"
             select type(Pat)
                class is (DiffPat_E_Type)
                   call Read_Pattern_DMC(trim(filename), Pat)
                   if (err_CFML%state) return
             
                   pat%ct_step = .true.
                   pat%instr=trim(Typef)
                   
                   select type(pat)
                      type is (DiffPat_G_Type)
                         pat%legend_X ="2Theta(degrees)"
                         pat%legend_Y ="Intensity (arb. units)"
                   end select      
             end select
             
          case ("SOCABIM") 
              call Read_Pattern_SOCABIM(trim(filename), Pat)
             if (err_CFML%state) return
             
             pat%kindrad = "XRay_CW"
             pat%scatvar = "2Theta"
             select type(Pat)
                type is (DiffPat_E_Type)
                   pat%ct_step = .true.
                   pat%instr=trim(Typef)
                  
                type is (DiffPat_G_Type)
                   pat%ct_step = .true.
                   pat%instr=trim(Typef)
                   
                   pat%legend_X ="2Theta(degrees)"
                   pat%legend_Y ="Intensity (arb. units)"
             end select
            
          case ("XYSIGMA") 
             i=index(filename,".",back=.true.)
             if (i > 0) then
                extdat=u_case(filename(i:))
             else
                extdat="---"
             end if     
             select case (trim(extdat))
                case (".GR")
                   if (present(header)) then
                      call  Read_Pattern_xysigma(filename, Pat, GR, Header)
                   else
                      call  Read_Pattern_xysigma(filename, Pat, GR)
                   end if
                     
                case ("---") 
                   if (present(header)) then
                      call  Read_Pattern_xysigma(filename, Pat, Header=Header)
                   else
                      call  Read_Pattern_xysigma(filename, Pat)
                   end if 
             end select
             if (err_CFML%state) return
             
             pat%kindrad = "Unknown"
             if (pat%x(pat%npts) > 180) then
                pat%scatvar =  "TOF"
             else
                pat%scatvar =  "2Theta"
             end if  
             select type(Pat)
                type is (DiffPat_E_Type)
                   pat%instr=trim(Typef)
                  
                type is (DiffPat_G_Type)
                   pat%instr=trim(Typef)
                   
                   pat%legend_Y ="Intensity (arb. units)"
                   if (index(pat%scatvar,"TOF") > 0) then
                      pat%legend_X ="TOF(micro-seconds)"
                   else
                      pat%legend_X ="2Theta(degrees)"
                   end if  
             end select    
             
          case ("GSAS") 
              call Read_Pattern_GSAS(trim(filename), Pat)
             if (err_CFML%state) return
             
             pat%kindrad = "Constant_Wavelength"
             pat%scatvar = "2Theta"
             select type(Pat)
                type is (DiffPat_E_Type)
                   pat%instr=trim(Typef)
                  
                type is (DiffPat_G_Type)
                   pat%instr=trim(Typef)
                   
                   pat%legend_X ="2Theta(degrees)"
                   pat%legend_Y ="Intensity (arb. units)"
             end select
             
          case ("GSASTOF") 
              call Read_Pattern_GSAS(trim(filename), Pat, TOFN)
             if (err_CFML%state) return
             
             pat%kindrad = "Neutrons_TOF"
             pat%scatvar = "TOF"
             select type(Pat)
                type is (DiffPat_E_Type)
                   pat%instr=trim(Typef)
                  
                type is (DiffPat_G_Type)
                   pat%instr=trim(Typef)
                   
                   pat%legend_X ="TOF(micro-seconds)"
                   pat%legend_Y ="Intensity (arb. units)"
             end select  
             
          case ("PANALYTICAL") 
             i=index(filename,".",back=.true.)
             extdat=u_case(filename(i:))
             select case (trim(extdat))
                case (".CSV")
                   call Read_Pattern_PANalytical_CSV(trim(filename), Pat)
                   
                case (".UDF")
                   call Read_Pattern_PANalytical_UDF(trim(filename), Pat)
                   
                case (".JCP")
                   call Read_Pattern_PANalytical_JCP(trim(filename), Pat)
                
                case(".XRDML") 
                   call Read_Pattern_PANalytical_XRDML(trim(filename), Pat)
             end select
             if (err_CFML%state) return
             
             pat%kindrad = "XRays_CW"
             pat%scatvar = "2Theta"
             select type(Pat)
                type is (DiffPat_E_Type)
                   pat%instr=trim(Typef)
                  
                type is (DiffPat_G_Type)
                   pat%instr=trim(Typef)
                   
                   pat%legend_X ="2Theta(degrees)"
                   pat%legend_Y ="Intensity (arb. units)"
             end select 
             
            
          case ("TIMEVARIABLE") 
             call Read_Pattern_Time_Variable(trim(filename), Pat)
             if (err_CFML%state) return
             
             pat%kindrad = "XRays_CW"
             pat%scatvar = "2Theta"
             select type(Pat)
                type is (DiffPat_E_Type)
                   pat%instr=trim(Typef)
                  
                type is (DiffPat_G_Type)
                   pat%instr=trim(Typef)
                   
                   pat%legend_X ="2Theta(degrees)"
                   pat%legend_Y ="Intensity (arb. units)"
             end select 
             
          case default   
             i=index(filename,".",back=.true.)
             extdat=u_case(filename(i:))
             call Read_Pattern_FREE(trim(filename), Pat, extdat)
             if (err_CFML%state) return
             
             pat%kindrad = "unknown"
             if (pat%x(pat%npts) > 180.0 ) then
                pat%scatvar =  "TOF"
             else
                pat%scatvar =  "2Theta"
             end if  
             select type (Pat)
                type is (DiffPat_E_Type)
                   pat%instr   = "Free format"
                   pat%ct_step=.true.
                   
                type is (DiffPat_G_Type)  
                   pat%instr   = "Free format"
                   pat%ct_step=.true.
                   pat%Legend_Y="Intensity (arb. units)"
                   
                   if (index(pat%scatvar,'TOF')> 0) then
                      pat%Legend_X="TOF(micro-seconds)"
                   else
                      pat%Legend_X="2Theta (degrees)"
                   end if  
             end select 
             
       end select
       
       !> CHECK PLEASE
       if (present(sig)) then
          if (.not. sig) pat%sigma=sqrt(pat%sigma)
       end if
       
       return
    End Subroutine Read_Pattern_One

 End Module CFML_Diffraction_Patterns
