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
!!---- You should have received a copy of the GNU of Lesser General Public
!!---- License along with this library; if not, see <http://www.gnu.org/licenses/>.
!!----
!!----
!!---- MODULE: CFML_GlobalDeps (Windows version)
!!----   INFO: Precision for CrysFML library and Operating System information
!!----         All the global variables defined in this module are implicitly public.
!!----
!!----
!!
Module CFML_GlobalDeps
   !---- Variables ----!
   Implicit None

   public

   !--------------------!
   !---- PARAMETERS ----!
   !--------------------!

   !---- Operating System ----!
   character(len=*), parameter :: OPS_NAME = "LIN"               ! O.S. Name
   character(len=*), parameter :: OPS_SEP = "/"                  ! O.S. directory separator character
   integer,          parameter :: OPS = 2                        ! O.S. Flag -> 1:Win 2:Lin 3:Mac

   !---- Compiler ----!
   !character(len=*), parameter :: COMPILER = "IFOR"              ! Intel Compiler
   character(len=*), parameter :: COMPILER = "GFOR"              ! GFortran Compiler

   !---- Precision ----!
   integer, parameter :: DP = selected_real_kind(14,150)         ! Double precision
   integer, parameter :: SP = selected_real_kind(6,30)           ! Simple precision
   integer, parameter :: LI = selected_int_kind(16)              ! Long Integer

   integer, parameter :: CP = SP                                 ! Current precision

   !---- Trigonometric ----!
   real(kind=DP), parameter :: PI = 3.141592653589793238463_DP   ! PI value
   real(kind=DP), parameter :: TO_DEG  = 180.0_DP/PI             ! Conversion from Radians to Degrees
   real(kind=DP), parameter :: TO_RAD  = PI/180.0_DP             ! Conversion from Degrees to Radians
   real(kind=DP), parameter :: TPI = 6.283185307179586476925_DP  ! 2*PI

   !---- Numeric ----!
   real(kind=DP), parameter :: DEPS=0.00000001_DP                ! Epsilon value use for comparison of real numbers (DOUBLE)
   real(kind=CP), parameter :: EPS=0.00001_CP                    ! Epsilon value use for comparison of real numbers 
   real(kind=CP), parameter :: V_EPSI=epsilon(1.0_CP)             ! Epsilon value for Current precision
   real(kind=CP), parameter :: V_HUGE=huge(1.0_CP)                ! Huge value for current precision 
   real(kind=CP), parameter :: V_TINY=tiny(1.0_CP)                ! Tiny value for current precision

   !---- Special Characters ----!
   character(len=2), parameter   :: NEWLINE = char(10)           ! Newline character
   character(len=1), parameter   :: TAB     = char(9)            ! TAB character
   
   !---- Error Flags ----!
   logical, parameter        :: CFML_DEBUG=.false.       ! For checking test
   
   Type :: Err_Type
      integer                         :: IErr =0                  ! =0: No error, < 0: Warning, > 0: Error
      character(len=180)              :: Msg=" "                  ! Text for Message
      integer                         :: nl=0                     ! number of lines
      character(len=132),dimension(5) :: Txt=" "                  ! Extra Message information
   End Type Err_Type
   Type (Err_Type)       :: Err_CFML                             ! Error Information for CFML
  
 Contains

   !!----
   !!---- FUNCTION DIRECTORY_EXISTS
   !!----
   !!----    Generic function dependent of the compiler that return
   !!----    a logical value if a directory exists or not.
   !!----
   !!---- Update: 11/07/2015
   !!
   Function Directory_Exists(DirName) Result(info)
      !---- Argument ----!
      character(len=*), intent(in) :: DirName               ! Directory name
      logical                      :: info                  ! Return Value

      !---- Local Variables ----!
      character(len=:), allocatable :: linea
      integer                       :: nlong

      !> Init value
      info=.false.
      
      !> Check
      if (len_trim(dirname)<= 0) return
      
      linea=adjustl(dirname)
      nlong=len_trim(linea)

      if (linea(nlong:nlong) /= OPS_SEP) linea=trim(linea)//OPS_SEP

      !> Compiler
      select case (trim(compiler))
         case ('IFOR')
            !inquire(directory=trim(linea), exist=info)

         case default
            inquire(file=trim(linea)//'.' , exist=info)
      end select

      return
   End Function Directory_Exists
   
   !!----
   !!---- CLEAR_ERROR
   !!----    Reset information on Error Variables for CFML
   !!----
   !!---- 27/03/2019 
   !!
   Subroutine Clear_Error()
      
      Err_CFML%IErr=0
      Err_CFML%Msg =" "
      Err_CFML%nl=0
      Err_CFML%Txt=" "
      
      return
   End Subroutine Clear_Error 

End Module CFML_GlobalDeps
