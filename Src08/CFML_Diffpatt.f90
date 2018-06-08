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
    public ::  Allocate_Diffraction_Pattern, Calc_Background, Del_NoisyPoints, &
               Read_Background_File, Read_Pattern,      &
               Purge_Diffraction_Pattern, Write_Pattern_XYSig,&
               Write_Pattern_FreeFormat, Add_Diffraction_Patterns,      &
               Write_Pattern_INSTRM5

    !---- Definitions ----!
    
    !!----
    !!---- TYPE :: DIFFPAT_TYPE
    !!----
    Type, public :: DiffPat_Type
       character(len=180)                        :: Title  =" "                  ! Indentification/ Title
       character(len=20)                         :: KindRad=" "                  ! Type of Radiation
       character(len=20)                         :: ScatVar=" "                  ! 2Theta, TOF, Q, s, d-spacing, SinTL/L,...
       real(kind=cp)                             :: xmin   =0.0_cp               ! Maximum and Minimum values for X and Y 
       real(kind=cp)                             :: xmax   =0.0_cp  
       real(kind=cp)                             :: ymin   =0.0_cp  
       real(kind=cp)                             :: ymax   =0.0_cp  
       integer                                   :: NPts   =0                    ! Number of Points
       logical                                   :: SigVar =.true.               ! .True. for sigma values / .False. for variance
       real(kind=cp), dimension(5)               :: Wave   =0.0_cp               ! Wave1, Wave2, Dtt1, Dtt2,....
       real(kind=cp), allocatable, dimension (:) :: x 
       real(kind=cp), allocatable, dimension (:) :: y 
       real(kind=cp), allocatable, dimension (:) :: sigma 
    End Type DiffPat_Type
    
    !!----
    !!---- TYPE :: DIFFPAT_ILL_TYPE
    !!----
    Type, public, extends (DiffPat_Type) ::  DiffPat_ILL_Type
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
    End Type DiffPat_ILL_Type

    !!----
    !!---- TYPE :: DIFFRACTION_PATTERN_TYPE
    !!----
    Type, public, extends (DiffPat_ILL_Type) :: DiffPat_G_Type
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
       Module Pure Function Calc_FWHM_Peak(Pat, Xi, Yi, Ybi, RLim) Result(v)
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
          class(DiffPat_ILL_Type),   intent(in out) :: Pat        ! Pattern object
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
          class(DiffPat_ILL_Type),  intent(in out)   :: Pat
       End Subroutine Read_Background_File
       
       Module Subroutine Set_Background_Poly(Pat, Bkpos, Bckx, N)
          !---- Arguments ----!
          class(DiffPat_ILL_Type),       intent(in out) :: Pat
          real (kind=cp),                intent(in    ) :: bkpos
          real (kind=cp), dimension(:),  intent(in    ) :: bckx
          integer,                       intent(in    ) :: n
       End Subroutine Set_Background_Poly
       
       Module Subroutine Set_Background_Inter(Pat, Bcky, Bckx, N)
          !---- Arguments ----!
          class(DiffPat_ILL_Type),       intent(in out) :: Pat
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
          class(DiffPat_ILL_Type),         intent(in out) :: Pat
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
          character(len=*),        intent(in)  :: Filename      ! Path+Filename
          class(DiffPat_ILL_Type), intent(out) :: Pat
       End Subroutine Read_Pattern_D1A_D2B 
       
       Module Subroutine Read_Pattern_D1A_D2B_OLD(Filename,Pat)
          !---- Arguments ----!
          character(len=*),        intent(in)  :: Filename      ! Path+Filename
          class(DiffPat_ILL_Type), intent(out) :: Pat
       End Subroutine Read_Pattern_D1A_D2B_OLD  
       
       Module Subroutine Read_Pattern_D1B_D20(Filename,Pat)
          !---- Arguments ----!
          character(len=*),        intent(in)  :: Filename      ! Path+Filename
          class(DiffPat_ILL_Type), intent(out) :: Pat
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
          character(len=*),        intent(in)  :: Filename      ! Path+Filename
          class(DiffPat_ILL_Type), intent(out) :: Pat
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
       
       !> class (DiffPat_ILL_Type)
       select type(Pat)
          class is (DiffPat_ILL_Type)
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
            
             al_x=.true.; al_y=.true.; al_sigma=.true.
             al_ycalc=.true.
             al_bgr  =.true.
             al_istat=.true.
            
       end select
       
       !> class (DiffPat_G_Type)
       select type(Pat)
          type is (DiffPat_G_Type)
             gy    =.false.
             gycalc=.false.
             gsigma=.false.
             gbgr  =.false.
             
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
                   
                class is (DiffPat_ILL_Type)
                   if (allocated(Pat%sigma)) deallocate(Pat%sigma) 
                   if (allocated(Pat%ycalc)) deallocate(Pat%ycalc) 
                   if (allocated(Pat%bgr))   deallocate(Pat%bgr) 
                   if (allocated(Pat%istat)) deallocate(Pat%istat) 
                   !if (allocated(Pat%ND))    deallocate(Pat%ND) 
                   
                   al_sigma=.false.
                   al_ycalc=.false.
                   al_bgr  =.false.
                   al_istat=.false.
             end select
            
          case ("DATAS")
             select type(Pat)
                class is (DiffPat_ILL_Type)
                   if (allocated(Pat%ycalc)) deallocate(Pat%ycalc) 
                   if (allocated(Pat%bgr))   deallocate(Pat%bgr) 
                   if (allocated(Pat%istat)) deallocate(Pat%istat) 
                   !if (allocated(Pat%ND))    deallocate(Pat%ND)
                   
                   al_ycalc=.false.
                   al_bgr  =.false.
                   al_istat=.false.
             end select
            
          case ("RIETV")
             select type(Pat)
                class is (DiffPat_ILL_Type)  
                   if (allocated(Pat%istat)) deallocate(Pat%istat)
                   
                   al_istat=.false.
             end select
            
          case ("GRAPH")
             select type(Pat)
                type is (DiffPat_ILL_Type) 
                   if (allocated(Pat%ycalc)) deallocate(Pat%ycalc) 
                   if (allocated(Pat%bgr))   deallocate(Pat%bgr)
                   
                   al_ycalc=.false.
                   al_bgr  =.false.
                    
                type is (DiffPat_G_Type)
                   if (allocated(Pat%ycalc)) deallocate(Pat%ycalc) 
                   if (allocated(Pat%bgr))   deallocate(Pat%bgr)
                   
                   al_ycalc=.false.
                   al_bgr  =.false.
                   
                   gycalc=.false.
                   gbgr  =.false.  
             end select
            
          case ("PRF")
             ! Nothing to do
       end select 

       return
    End Subroutine Deallocate_Diffraction_Pattern
    
    

    !!----
    !!---- Subroutine Read_Pattern(Filename, Dif_Pat, Mode, header)
    !!--<<                   or   (Filename, Dif_Pat, NumPat, Mode, header)
    !!----    character(len=*),                              intent (in)    :: Filename
    !!----    type (diffraction_pattern_type),               intent (in out):: Dif_Pat
    !!----    character(len=*), optional,                    intent (in)    :: mode
    !!----
    !!----    character(len=*),                              intent (in)    :: Filename
    !!----    type (diffraction_pattern_type), dimension(:), intent (in out):: Dif_Pat
    !!----    integer,                                       intent (out)   :: numpat
    !!----    character(len=*), optional,                    intent (in)    :: mode
    !!----    character(len=*), optional,                    intent (out)   :: header
    !!-->>
    !!----    Read one pattern from a Filename
    !!----
    !!---- Update: February - 2005
    !!

    

    

    

    

    

    


    

    

    !!--++
    !!--++ Subroutine Read_Pattern_Mult(Filename,Dif_Pat, NumPat, Mode)
    !!--++    character(len=*),                                          intent (in)      :: filename
    !!--++    type (diffraction_pattern_type), dimension(:),             intent (in out)  :: dif_pat
    !!--++    integer,                                                   intent (out)     :: numpat
    !!--++    character(len=*), optional,                                intent (in)      :: mode
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Read one pattern from a Filename
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Read_Pattern_Mult(filename, dif_pat, numpat, mode)
       !---- Arguments ----!
       character(len=*),                                          intent (in)      :: filename
       type (diffraction_pattern_type), dimension(:),             intent (in out)  :: dif_pat
       integer,                                                   intent (in out)  :: numpat
       character(len=*), optional,                                intent (in)      :: mode

       !---- Local variables ----!
       logical :: esta
       integer :: i_dat, ier,i

       call init_err_diffpatt()

       inquire(file=filename,exist=esta)
       if ( .not. esta) then
          Err_diffpatt=.true.
          ERR_DiffPatt_Mess=" The file "//trim(filename)//" doesn't exist"
          return
       else
          call get_logunit(i_dat)
          open(unit=i_dat,file=trim(filename),status="old",action="read",position="rewind",iostat=ier)
          if (ier /= 0) then
             Err_diffpatt=.true.
             ERR_DiffPatt_Mess=" Error opening the file "//trim(filename)
             return
          end if
          i=index(filename,OPS_SEP,back=.true.)
          If( i /= 0) then
            dif_pat%filename=trim(filename(i+1:))
            dif_pat%filepath=filename(1:i)
          Else
            dif_pat%filename=trim(filename)
            dif_pat%filepath="."//OPS_SEP
          End if
       end if

       if (present(mode)) then
          select case (u_case(mode))
              case ("XYSIGMA")
                 !   call  Read_Pattern_xysigma_m(dif_pat,npat)

              case ("ISIS")
                 call Read_Pattern_isis_m(i_dat,dif_pat,numpat)
                 dif_pat%diff_kind = "neutrons_tof"
                 dif_pat%scat_var =  "TOF"
                 dif_pat%xax_text =  "TOF (micro-seconds)"
                 dif_pat%yax_text =  "Intensity (arb. units)"
                 dif_pat%instr  = " 14  - "//mode

              case ("GSAS")
                 !   call Read_Pattern_gsas_m(dif_pat,npat)      ! GSAS file

              case default
                 Err_diffpatt=.true.
                 ERR_DiffPatt_Mess="Invalid Mode"
                 return
          end select
          return
       end if
       close(unit=i_dat,iostat=ier)

       if (ier/=0) then
           Err_diffpatt=.true.
           ERR_DiffPatt_Mess=" Problems closing data file"
       end if

       return
    End Subroutine Read_Pattern_Mult

    

    !!--++
    !!--++ Subroutine Read_Pattern_One(Filename,Dif_Pat, Mode,header,sig)
    !!--++    character(len=*),                intent (in)    :: filename
    !!--++    type (diffraction_pattern_type), intent(in out) :: Dif_Pat
    !!--++    character(len=*), optional,      intent (in)    :: mode
    !!--++    logical,          optional,      intent (in)    :: sig
    !!--++
    !!--++    Read one pattern from a Filename. If sig is present the content of Dif_Pat%sigma
    !!--++    is the true sigma not the variance.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Read_Pattern_One(Filename,Dif_Pat,Mode,header,sig)
       !---- Arguments ----!
       character(len=*),                intent (in)      :: filename
       type (diffraction_pattern_type), intent (in out)  :: dif_pat
       character(len=*), optional,      intent (in)      :: mode
       character(len=*), optional,      intent (out)     :: header
       logical,          optional,      intent (in)      :: sig

       !---- Local Variables ----!
       character(len=6)                               :: extdat !extension of panalytical file
       character(len=4)                               :: tofn
       character(len=12)                              :: modem !extension of panalytical file
       logical                                        :: esta,gr
       integer                                        :: i, i_dat,ier

       call init_err_diffpatt()

       inquire(file=filename,exist=esta)
       if (.not. esta) then
          Err_diffpatt=.true.
          ERR_DiffPatt_Mess=" The file "//trim(filename)//" doesn't exist"
          return
       else
          call get_logunit(i_dat)
          open(unit=i_dat,file=trim(filename),status="old",action="read",position="rewind",iostat=ier)
          if (ier /= 0) then
             Err_diffpatt=.true.
             ERR_DiffPatt_Mess=" Error opening the file: "//trim(filename)
             return
          end if
          i=index(filename,OPS_SEP,back=.true.)
          If( i /= 0) then
            dif_pat%filename=trim(filename(i+1:))
            dif_pat%filepath=filename(1:i)
          Else
            dif_pat%filename=trim(filename)
            dif_pat%filepath="."//OPS_SEP
          End if
       end if

       if (present(mode)) then
          modem=u_case(mode)
          if(modem(1:7) == "GSASTOF") then
            if(len_trim(modem) > 7) then
               tofn="TOF"//modem(8:8)
               modem="GSASTOF"
            else
               tofn="TOF"
            end if
          end if
       else
          modem="DEFAULT"
       end if

       select case (modem)

          case ("CIF")
             call Read_Pattern_CIF(i_dat,dif_pat)
             dif_pat%yax_text =  "Intensity (arb. units)"
             dif_pat%instr  = "  XY  - "//mode
             dif_pat%ct_step = .false.

          case ("D1B" , "D20")
             call Read_Pattern_d1b_d20(i_dat,dif_pat)
             dif_pat%diff_kind = "neutrons_cw"
             dif_pat%scat_var =  "2theta"
             dif_pat%xax_text =  "2theta(degrees)"
             dif_pat%yax_text =  "Intensity (arb. units)"
             dif_pat%instr  = "  3  - "//mode
             dif_pat%ct_step = .true.

          case ("NLS")                   ! Data from N.L.S (Brookhaven) Synchrotron Radiation  ,data from synchrotron source and correct data for dead time
             call Read_Pattern_nls(i_dat,dif_pat)
             dif_pat%diff_kind = "xrays_cw"
             dif_pat%scat_var =  "2theta"
             dif_pat%xax_text =  "2theta(degrees)"
             dif_pat%yax_text =  "Intensity (arb. units)"
             dif_pat%instr  = "  4  - "//mode
             dif_pat%ct_step = .true.

          case ("G41")                   ! Data from general format of two axis instruments with fixed step in twotheta
             call Read_Pattern_g41(i_dat,dif_pat)
             dif_pat%diff_kind = "neutrons_cw"
             dif_pat%scat_var =  "2theta"
             dif_pat%xax_text =  "2theta(degrees)"
             dif_pat%yax_text =  "Intensity (arb. units)"
             dif_pat%instr  = "  5  - "//mode
             dif_pat%ct_step = .true.

          case ("D1A","D2B","3T2","G42")
             call Read_Pattern_d1a_d2b(i_dat,dif_pat)     ! Data from D1A,D2B  (Files *.sum, renamed *.dat, as prepared by D1ASUM or D2BSUM programs)
             dif_pat%diff_kind = "neutrons_cw"
             dif_pat%scat_var =  "2theta"
             dif_pat%xax_text =  "2theta(degrees)"
             dif_pat%yax_text =  "Intensity (arb. units)"
             dif_pat%instr  = "  6  - "//mode
             dif_pat%ct_step = .true.

          case ("D1AOLD", "D2BOLD","OLDD1A", "OLDD2B")
             call Read_Pattern_d1a_d2b_old(i_dat,dif_pat)
             dif_pat%diff_kind = "neutrons_cw"
             dif_pat%scat_var =  "2theta"
             dif_pat%xax_text =  "2theta(degrees)"
             dif_pat%yax_text =  "Intensity (arb. units)"
             dif_pat%instr  = "  1  - "//mode
             dif_pat%ct_step = .true.

          case ("DMC","HRPT")                   ! Data from DMC,HRPT
             call Read_Pattern_dmc(i_dat,dif_pat)
             dif_pat%diff_kind = "neutrons_cw"
             dif_pat%scat_var =  "2theta"
             dif_pat%xax_text =  "2theta(degrees)"
             dif_pat%yax_text =  "Intensity (arb. units)"
             dif_pat%instr  = "  8  - "//mode
             dif_pat%ct_step = .true.

          case ("SOCABIM")
             call  Read_Pattern_socabim(i_dat,dif_pat)
             dif_pat%diff_kind = "xrays_cw"
             dif_pat%scat_var =  "2theta"
             dif_pat%xax_text =  "2theta(degrees)"
             dif_pat%yax_text =  "Intensity (arb. units)"
             dif_pat%instr  = "  9  - "//mode
             dif_pat%ct_step = .true.

          case ("XYSIGMA")            !XYSIGMA  data file
             !Determine if the patter is of G(r) type from PDFGUI
             i=index(dif_pat%filename,".",back=.true.)
             if(present(header)) then
                 if(i /= 0 .and. dif_pat%filename(i:i+2) == ".gr") then
                   call  Read_Pattern_xysigma(i_dat, dif_pat,gr,header)
                 else
                   call  Read_Pattern_xysigma(i_dat, dif_pat,header=header)
                   dif_pat%diff_kind = "unknown"
                   dif_pat%instr  = " 10  - "//mode
                 end if
             else
                 if(i /= 0 .and. dif_pat%filename(i:i+2) == ".gr") then
                   call  Read_Pattern_xysigma(i_dat, dif_pat,gr)
                 else
                   call  Read_Pattern_xysigma(i_dat, dif_pat)
                   dif_pat%diff_kind = "unknown"
                   dif_pat%instr  = " 10  - "//mode
                 end if
             end if
             if(len_trim(dif_pat%scat_var) == 0) then
               if(dif_pat%x(dif_pat%npts) > 180.0) then
                   dif_pat%scat_var =  "TOF"
               else
                   dif_pat%scat_var =  "2theta"
               end if
             end if

          case ("GSAS")
             call Read_Pattern_gsas(i_dat,dif_pat)         ! GSAS file
             dif_pat%diff_kind = "constant_wavelength"
             dif_pat%scat_var =  "2theta"
             dif_pat%xax_text =  "2theta(degrees)"
             dif_pat%yax_text =  "Intensity (arb. units)"
             dif_pat%instr  = " 12  - "//mode

          case ("GSASTOF")
             call Read_Pattern_gsas(i_dat,dif_pat,tofn)         ! GSAS file for TOF
             dif_pat%diff_kind = "neutrons_tof"
             dif_pat%scat_var =  "TOF"
             dif_pat%xax_text =  "TOF(micro-seconds)"
             dif_pat%yax_text =  "Intensity (arb. units)"
             dif_pat%instr  = " 12  - "//mode

          case ("PANALYTICAL")
             i=index(filename,".",back=.true.)
             extdat=u_case(filename(i:))
             dif_pat%diff_kind = "xrays_cw"
             dif_pat%scat_var =  "2theta"
             dif_pat%xax_text =  "2theta(degrees)"
             dif_pat%yax_text =  "Intensity (arb. units)"
             dif_pat%instr  = " 13  - "//mode

             select case (extdat)
                case(".CSV")
                   CALL Read_Pattern_PANalytical_CSV(i_dat,dif_pat)

                case(".UDF")
                   CALL Read_Pattern_PANalytical_UDF(i_dat,dif_pat)

                case(".JCP")
                   CALL Read_Pattern_PANalytical_JCP(i_dat,dif_pat)

                case(".XRDML")
                   CALL Read_Pattern_PANalytical_XRDML(i_dat,dif_pat)
             end select

          case ("TIMEVARIABLE")
             call Read_Pattern_time_variable(i_dat,dif_pat)
             dif_pat%diff_kind = "xrays_cw"
             dif_pat%scat_var =  "2theta"
             dif_pat%xax_text =  "2theta(degrees)"
             dif_pat%yax_text =  "Intensity (arb. units)"
             dif_pat%instr  = " 11  - "//mode

          case default
             i=index(filename,".",back=.true.)
             extdat=u_case(filename(i:))
             call Read_Pattern_free(i_dat,dif_pat,extdat)
             if(Err_diffpatt) return
             dif_pat%diff_kind = "unknown"
             dif_pat%instr  = "  0  - "//"Free format"
             dif_pat%ct_step = .true.
             if(len_trim(dif_pat%yax_text) == 0) dif_pat%yax_text =  "Intensity (arb. units)"
             if(len_trim(dif_pat%xax_text) == 0)  then
                if(dif_pat%x(dif_pat%npts) > 180.0 ) then
                    dif_pat%scat_var =  "TOF"
                    dif_pat%xax_text =  "TOF(micro-seconds)"
                else
                    dif_pat%scat_var =  "2theta"
                    dif_pat%xax_text =  "2theta(degrees)"
                end if
             else
                if(len_trim(dif_pat%scat_var) == 0) dif_pat%scat_var =  "2theta"
             end if
       end select

       close(unit=i_dat,iostat=ier)
       if(present(sig)) then
          dif_pat%sigma=sqrt(dif_pat%sigma)
          dif_pat%sig_var=.false.
       end if
       if (ier/=0) then
          Err_diffpatt=.true.
          ERR_DiffPatt_Mess=" Problems closing the data file: "//trim(filename)
       end if

       return
    End Subroutine Read_Pattern_One

    


    

 End Module CFML_Diffraction_Patterns
