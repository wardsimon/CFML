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
!!---- MODULE: CFML_LSQ_TypeDef
!!----   INFO: Type definitions for LSQ Routines
!!----
!!---- HISTORY
!!--..    Update: 01/07/2015
!!--..
!!---- DEPENDENCIES
!!----    CFML_GlobalDeps
!!----
!!---- VARIABLES
!!----    MAX_FREE_PAR
!!--..
!!----    LSQ_CONDITIONS_TYPE
!!----    LSQ_DATA_TYPE
!!----    LSQ_STATE_VECTOR_TYPE
!!----
!!---- PROCEDURES
!!----    Subroutines:
!!----    MODIFY_CODE_STATE_VECTOR
!!----
!!
Module CFML_LSQ_TypeDef
   !---- Use Files ----!
   Use CFML_GlobalDeps,   only: Cp, Dp

   !---- Variables ----!
   implicit None

   !---- Definitions ----!

   !!----
   !!---- MAX_FREE_PAR
   !!----    integer, parameter, public  :: Max_Free_Par
   !!----
   !!----    Maximum number of free parameters (3000)
   !!----    (it may be changed at will!)
   !!----
   !!---- Update: 01/07/2015
   !!
   integer, parameter :: Max_Free_Par=3000   !Maximum number of free parameters


   !!----
   !!---- TYPE :: LSQ_CONDITIONS_TYPE
   !!--..
   !!----  Type, public :: LSQ_Conditions_type
   !!----     logical          :: constr          ! if true box constraint of percent% are applied to parameters
   !!----     logical          :: reached         ! if true convergence was reached in the algorithm
   !!----     integer          :: corrmax         ! value of correlation in % to output
   !!----     integer          :: nfev            ! number of function evaluations (output component, useful for assessing LM algorithm)
   !!----     integer          :: njev            ! number of Jacobian evaluations                 "
   !!----     integer          :: icyc            ! number of cycles of refinement or maximum number of function evaluations in LM
   !!----                                         ! In LM procedures the default value is icyc = maxfev = 100(npvar+1)
   !!----     integer          :: npvar           ! number of effective free parameters of the model
   !!----     integer          :: iw              ! indicator for weighting scheme (if iw=1 => w=1/yc)
   !!----     integer          :: nprint          ! indicator for printing during iterations, if nprint > 0 printing each nprint iterations
   !!----     real(kind=cp)    :: tol             ! tolerance value for applying stopping criterion in LM algorithm
   !!----     real(kind=cp)    :: percent         ! %value of maximum variation of a parameter w.r.t.
   !!----                                         ! the intial value before fixing it
   !!----  End Type LSQ_Conditions_type
   !!----
   !!----  Derived type encapsulating all necessary conditions for running the LSQ algorithm
   !!----
   !!---- Update: 01/07/2015
   !!
   Type :: LSQ_Conditions_Type
      logical          :: constr=.false.  ! if true box constraint of percent% are applied to parameters
      logical          :: reached=.false. ! if true convergence was reached in the algorithm
      integer          :: corrmax=50      ! value of correlation in % to output
      integer          :: nfev=0          ! number of function evaluations (output component, useful for assessing LM algorithm)
      integer          :: njev=0          ! number of Jacobian evaluations                 "
      integer          :: icyc=0          ! number of cycles of refinement or maximum number of function evaluations in LM
                                          ! In LM procedures the default value is icyc = maxfev = 100(npvar+1)
      integer          :: npvar=0         ! number of effective free parameters of the model
      integer          :: iw=0            ! indicator for weighting scheme (if iw=1 => w=1/yc)
      integer          :: nprint=0        ! indicator for printing during iterations, if nprint > 0 printing each nprint iterations
      real(kind=cp)    :: tol=0.0         ! tolerance value for applying stopping criterion in LM algorithm
      real(kind=cp)    :: percent=0.0     ! %value of maximum variation of a parameter w.r.t.
                                          ! the intial value before fixing it
   End Type LSQ_Conditions_Type

   !!----
   !!---- TYPE :: LSQ_DATA_TYPE
   !!--..
   !!----
   !!----  Type, public :: LSQ_Data_Type
   !!----     integer                                    :: nobs  !total number of observations
   !!----     integer                                    :: iw    !Indicator for type of values contained in component sig
   !!----     real(kind=cp),    dimension(:),allocatable :: x     !Vector containing a relevant quantity for each observation (x-coordinate ...)
   !!----     real(kind=cp),    dimension(:),allocatable :: y     !Vector containing the observed values
   !!----     real(kind=cp),    dimension(:),allocatable :: sw    !Vector containing the standard deviation of observations if iw=0
   !!----                                                         !or the weight factors for least squares refinement if iw=1
   !!----     real(kind=cp),    dimension(:),allocatable :: yc    !Vector containing the calculated values
   !!----  End Type LSQ_Data_type
   !!----
   !!----  Derived type encapsulating the observed and calculated data as well as the
   !!----  weighting factors, a variable related with each observed value and the
   !!----  total number of observations. It is responsibility of the calling program to
   !!----  allocate the components before calling the Marquardt_fit procedure.
   !!----
   !!---- Update: August - 2009
   !!
   Type :: LSQ_Data_Type
      integer                                    :: nobs  !total number of observations
      integer                                    :: iw    !Indicator for type of values contained in component sw
      real(kind=cp),    dimension(:),allocatable :: x     !Vector containing a relevant quantity for each observation (x-coordinate ...)
      real(kind=cp),    dimension(:),allocatable :: y     !Vector containing the observed values
      real(kind=cp),    dimension(:),allocatable :: sw    !Vector containing the standard deviation of observations if iw=0
                                                          !or the weight factors for least squares refinement if iw=1
      real(kind=cp),    dimension(:),allocatable :: yc    !Vector containing the calculated values
   End Type LSQ_Data_type

   !!----
   !!---- TYPE :: LSQ_STATE_VECTOR_TYPE
   !!--..
   !!----
   !!----  Type, public :: LSQ_State_Vector_Type
   !!----     integer                                    :: np         !total number of model parameters <= Max_Free_Par
   !!----     logical                                    :: code_comp  !If .true. the codes are interpreted as number in the LSQ list
   !!----     integer(kind=2)                            :: code_max   !Maximum code number (used in case of code_comp=.true.)
   !!----     real(kind=cp),     dimension(Max_Free_Par) :: mul        !Vector of multipliers (used in case of code_comp=.true.
   !!----     real(kind=cp),     dimension(Max_Free_Par) :: pv         !Vector of parameters
   !!----     real(kind=cp),     dimension(Max_Free_Par) :: spv        !Vector of standard deviations
   !!----     real(kind=cp),     dimension(Max_Free_Par) :: dpv        !Vector of derivatives at a particular point
   !!----     integer(kind=2),   dimension(Max_Free_Par) :: code       !pointer for selecting variable parameters
   !!----     character(len=40), dimension(Max_Free_Par) :: nampar     !Names of parameters
   !!----  End Type LSQ_State_Vector_Type
   !!----
   !!----  Derived type encapsulating the vector state defining a set of parameter
   !!----  for calculating the model function and running the LSQ algorithm.
   !!----  Now, with the introduction of code_comp and mul, the codes may be also interpreted
   !!----  as the ordinal number in the LSQ list of independent parameters. Depending on the
   !!----  the way the user program attributes codes and constraints a call to the subroutine
   !!----  Modify_Codes_State_Vector (see below)
   !!----
   !!---- Update: January- 2014
   !!
   Type :: LSQ_State_Vector_Type
      integer                                    :: np         !total number of model parameters <= Max_Free_Par
      logical                                    :: code_comp  !If .true. the codes are interpreted as number in the LSQ list
      integer(kind=2)                            :: code_max   !Maximum code number (used in case of code_comp=.true.)
      real(kind=cp),     dimension(Max_Free_Par) :: mul        !Vector of multipliers (used in case of code_comp=.true.)
      real(kind=cp),     dimension(Max_Free_Par) :: pv         !Vector of parameters
      real(kind=cp),     dimension(Max_Free_Par) :: spv        !Vector of standard deviations
      real(kind=cp),     dimension(Max_Free_Par) :: dpv        !Vector of derivatives at a particular point
      integer(kind=2),   dimension(Max_Free_Par) :: code       !pointer for selecting variable parameters
      character(len=40), dimension(Max_Free_Par) :: nampar     !Names of parameters
   End Type LSQ_State_Vector_type

 contains
    !!----
    !!---- Subroutine Modify_Codes_State_Vector(Lcodes,Multip,Lcode_max)
    !!----   integer,       dimension(:),intent (in out) :: Lcodes
    !!----   real(kind=cp), dimension(:),intent (in out) :: Multip
    !!----   integer,                    intent (in out) :: Lcode_max
    !!----
    !!----  This subroutine should be called with the arguments corresponding to
    !!----  the components  Lcodes=LSQ_State_Vector%code, Multip=LSQ_State_Vector%mul
    !!----  and Lcode_max=maxval(LSQ_State_Vector%code) before using the LSQ methods
    !!----  in which and object of LSQ_State_Vector_type exist and LSQ_State_Vector%code_comp=.true.
    !!----  This allows the use of constraints once the codes have been attributed.
    !!----
    Subroutine Modify_Codes_State_Vector(Lcodes,Multip,Lcode_max)
       !---- Arguments ----!
       integer,       dimension(:),intent (in out) :: Lcodes
       real(kind=cp), dimension(:),intent (in out) :: Multip
       integer,                    intent (in out) :: Lcode_max

       !--- Local variables ---!
       integer                       ::  k, L , j, n_given, Lcm, n_att,ndisp, nn, ni, ndispm,maxs
       integer, dimension(Lcode_max) :: dispo
       real(kind=cp), parameter      :: e=0.001

       !> Check correlated parameters and already used codes
       ndisp=0
       n_given=0
       dispo(:)=0
       Lcm=size(Lcodes)
       !> First Pass
       Do L=1, Lcode_max
          ni=0
          do k=1, Lcm
             if (L == Lcodes(k)) ni=ni+1
          end do
          if (ni == 0) then
             ndisp=ndisp+1   !number of disponible codes
             dispo(ndisp)=L  !disponible code number
          else
             n_given=n_given+1  !number of given codes
          end if
       end do !=1,Lcode_max

       if (ndisp == 0) return  !all codes have been attributed

       !> Attributing numbers to parameters (not already attributed) with multipliers
       !> different from zero. First the attribution is taken from the vector dispo() and
       !> continued, after finishing the disponible codes, from Lcode_max+1, ...
       !> If after attributing the codes ndisp /=0, then a displacement of all parameters
       !> is done and the maximum number of parameters to be refined is diminished by
       !> ndisp
       ndispm=ndisp !number of disponible codes before attributing code numbers
       ni=0
       nn=Lcode_max
       n_att=0
       do j=1,Lcm
          if (abs(Multip(j)) > e .and. Lcodes(j) == 0 ) then
             if (abs(Multip(j)) > 1.001) then
                Multip(j)=sign(1.0_cp,Multip(j))*(abs(Multip(j))-1.0)
             end if
             if (ndisp==0) then
                nn=nn+1
                Lcodes(j)=nn
                n_att=n_att+1
             else
                ni=ni+1
                Lcodes(j)=dispo(ni)
                n_att=n_att+1
                ndisp=ndisp-1
             end if
          end if
       end do

       if (ndisp == 0) then
          if (nn > Lcode_max) Lcode_max=nn
          return  !all parameters have been attributed
       end if

       maxs=n_given+n_att    !Number of refined parameters (given + attributed)

       !> Third Pass ndisp /=0 => Displacement of codes needed to avoid holes in the matrix.
       n_att=ndispm-ndisp+1
       do L=Lcode_max, maxs+1,-1
          nn=0
          do j =1,Lcm
             if (L == Lcodes(j)) then
                Lcodes(j)=dispo(n_att)
                nn=nn+1
             end if
          end do
          if (nn == 0) cycle
          n_att=n_att+1
          if(n_att > ndispm) exit
       end do !i=Lcode_max,maxs+1,-1
       Lcode_max=maxs

       return
    End Subroutine Modify_Codes_State_Vector

End Module CFML_LSQ_TypeDef
