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
!!---- MODULE: CFML_Spherical_Harmonics
!!----   INFO: Spherical Harmonics routines
!!----
!!---- HISTORY
!!----    Update: 03/03/2011
!!----
!!---- DEPENDENCIES
!!--++    Use CFML_GlobalDeps, only: cp, dp, eps, pi, to_rad
!!----
!!---- VARIABLES
!!----    ERR_SPHER
!!----    ERR_SPHER_MESS
!!----
!!---- PROCEDURES
!!----    Functions:
!!----       CUBIC_HARM_ANG
!!----       CUBIC_HARM_UCVEC
!!--++       ENVJ                      [Private]
!!----       INT_SLATER_BESSEL
!!--++       PLGNDR                    [Private]
!!----       REAL_SPHER_HARM_ANG
!!----       REAL_SPHER_HARM_UCVEC
!!----       REAL_SPHER_HARMCHARGE_UCVEC
!!--++       START1                    [Private]
!!--++       START2                    [Private]
!!----
!!----    Subroutines:
!!----       INIT_ERR_SPHER
!!----       PIKOUT_LJ_CUBIC
!!----       SPHJN
!!----
!!
 Module CFML_Spherical_Harmonics

    !---- Use Modules ----!
    Use CFML_GlobalDeps, Only: Cp, Dp, Eps, Pi, To_Rad

    implicit none

    private

    !---- List of public functions ----!
    public :: Cubic_Harm_Ang, Cubic_Harm_Ucvec, Int_Slater_Bessel, Real_Spher_Harm_Ang, &
              Real_Spher_Harm_Ucvec, Real_Spher_HarmCharge_Ucvec

    !---- List of public subroutines ----!
    public :: Init_Err_Spher, Pikout_Lj_Cubic, Sphjn

    !---- Private Functions / Routines ----!
    private :: Envj, Plgndr, Start1, Start2

    !---- Definitions ----!

    !!----
    !!---- ERR_SPHER
    !!----    logical, public    :: err_spher
    !!----
    !!----    Logical Variable indicating an error in CFML_Spherical_Harmonics module
    !!----
    !!---- Update: February - 2005
    !!
    logical, public    :: Err_Spher

    !!----
    !!---- ERR_Spher_Mess
    !!----    character(len=150), public :: ERR_Spher_Mess
    !!----
    !!----    String containing information about the last error
    !!----
    !!---- Update: February - 2005
    !!
    character(len=150), public :: Err_Spher_Mess

 Contains

    !---- Functions ----!

    !!----
    !!---- Function Cubic_Harm_Ang(L,M,Theta,Phi) Result(Klm)
    !!----    integer,         intent(in) :: l
    !!----    integer,         intent(in) :: m
    !!----    real(kind=cp),   intent(in) :: theta
    !!----    real(kind=cp),   intent(in) :: phi
    !!----    real(kind=cp)               :: Klm
    !!----
    !!----    Calculation of the cubic harmonics given in Table 3 of reference
    !!----    M.Kara & K. Kurki-Suonio, Acta Cryt. A37, 201 (1981).
    !!----    Only up to tenth order.
    !!----
    !!---- Update: February - 2005
    !!
    Pure Function Cubic_Harm_Ang(L,M,Theta,Phi) Result(Klm)
       !---- Arguments ----!
       integer,      intent (in) :: l,m
       real(kind=cp),intent (in) :: theta,phi
       real(kind=cp)             :: klm

       select case (l)
          case(0)   ! 0 1
             klm=real_spher_harm_ang(0,0,1,theta,phi)
          case(3)   ! 3 1
             klm=real_spher_harm_ang(3,2,-1,theta,phi)
          case(4)   ! 4 1
             klm=0.5*sqrt(7.0/3.0)*real_spher_harm_ang(4,0,1,theta,phi)
             klm=klm+0.5*sqrt(5.0/3.0)*real_spher_harm_ang(4,4,1,theta,phi)
          case(6)
             if (m == 1) then   ! 6 1
                klm=0.5*sqrt(0.5)*real_spher_harm_ang(6,0,1,theta,phi)
                klm=klm-0.5*sqrt(7.0/2.0)*real_spher_harm_ang(6,4,1,theta,phi)
             else               ! 6 2
                klm=0.25*sqrt(11.0)*real_spher_harm_ang(6,2,1,theta,phi)
                klm=klm-0.25*sqrt(5.0)*real_spher_harm_ang(6,6,1,theta,phi)
             end if
          case(7)   ! 7 1
             klm=0.5*sqrt(13.0/6.0)*real_spher_harm_ang(7,2,-1,theta,phi)
             klm=klm+0.5*sqrt(11.0/6.0)*real_spher_harm_ang(7,6,-1,theta,phi)
          case(8)   ! 8 1
             klm=0.125*sqrt(33.0)*real_spher_harm_ang(8,0,1,theta,phi)
             klm=klm+0.25*sqrt(7.0/3.0)*real_spher_harm_ang(8,4,1,theta,phi)
             klm=klm+0.125*sqrt(65.0/3.0)*real_spher_harm_ang(8,8,1,theta,phi)
          case(9)
             if (m == 1) then   ! 9 1
                klm=0.25*sqrt(3.0)*real_spher_harm_ang(9,2,-1,theta,phi)
                klm=klm-0.25*sqrt(13.0)*real_spher_harm_ang(9,6,-1,theta,phi)
             else               ! 9 2
                klm=0.5*sqrt(17.0/6.0)*real_spher_harm_ang(9,4,-1,theta,phi)
                klm=klm-0.5*sqrt(7.0/6.0)*real_spher_harm_ang(9,8,-1,theta,phi)
             end if
          case(10)
             if (m == 1) then   ! 10 1
                klm=0.125*sqrt(65.0/6.0)*real_spher_harm_ang(10,0,1,theta,phi)
                klm=klm-0.25*sqrt(11.0/2.0)*real_spher_harm_ang(10,4,1,theta,phi)
                klm=klm-0.125*sqrt(187.0/6.0)*real_spher_harm_ang(10,8,1,theta,phi)
             else               ! 10 2
                klm=0.125*sqrt(247.0/6.0)*real_spher_harm_ang(10,2,1,theta,phi)
                klm=klm+(1.0/16.0)*sqrt(19.0/3.0)*real_spher_harm_ang(10,6,1,theta,phi)
                klm=klm-(1.0/16.0)*sqrt(85.0)*real_spher_harm_ang(10,10,1,theta,phi)
             end if
          case default
             klm=0.0
       end select

       return
    End Function Cubic_Harm_Ang

    !!----
    !!---- Function Cubic_Harm_Ucvec(L,M,U) Result(Klm)
    !!----    integer,                    intent(in) :: l
    !!----    integer,                    intent(in) :: m
    !!----    real(kind=cp),dimension(3), intent(in) :: u
    !!----    real(kind=cp)                          :: Klm
    !!--<<
    !!----    Calculation of the cubic harmonics given in Table 3 of reference
    !!----    M.Kara & K. Kurki-Suonio, Acta Cryt. A37, 201 (1981).
    !!----    Only up to tenth order. A control of errors is included.
    !!----    For "m3m" symmetry, calculations include up to L=20 M=2 using the
    !!----    coefficients from F.M. Mueller and M.G. Priestley, Phys Rev 148, 638 (1966)
    !!-->>
    !!----
    !!---- Update: February - 2005
    !!
    Pure Function Cubic_Harm_Ucvec(L,M,U) Result(Klm)
       !---- Arguments ----!
       integer,                    intent (in) :: l,m
       real(kind=cp),dimension(3), intent (in) :: u
       real(kind=cp)                           :: Klm

       select case (l)
          case(0)   ! 0 1
             klm=real_spher_harm_ucvec(0,0,1,u)
          case(3)   ! 3 1
             klm=real_spher_harm_ucvec(3,2,-1,u)
          case(4)   ! 4 1
             klm=0.5*sqrt(7.0/3.0)*real_spher_harm_ucvec(4,0,1,u)
             klm=klm+0.5*sqrt(5.0/3.0)*real_spher_harm_ucvec(4,4,1,u)
          case(6)
             if (m == 1) then   ! 6 1
                klm=0.5*sqrt(0.5)*real_spher_harm_ucvec(6,0,1,u)
                klm=klm-0.5*sqrt(7.0/2.0)*real_spher_harm_ucvec(6,4,1,u)
             else               ! 6 2
                klm=0.25*sqrt(11.0)*real_spher_harm_ucvec(6,2,1,u)
                klm=klm-0.25*sqrt(5.0)*real_spher_harm_ucvec(6,6,1,u)
             end if
          case(7)   ! 7 1
             klm=0.5*sqrt(13.0/6.0)*real_spher_harm_ucvec(7,2,-1,u)
             klm=klm+0.5*sqrt(11.0/6.0)*real_spher_harm_ucvec(7,6,-1,u)
          case(8)   ! 8 1
             klm=0.125*sqrt(33.0)*real_spher_harm_ucvec(8,0,1,u)
             klm=klm+0.25*sqrt(7.0/3.0)*real_spher_harm_ucvec(8,4,1,u)
             klm=klm+0.125*sqrt(65.0/3.0)*real_spher_harm_ucvec(8,8,1,u)
          case(9)
             if (m == 1) then   ! 9 1
                klm=0.25*sqrt(3.0)*real_spher_harm_ucvec(9,2,-1,u)
                klm=klm-0.25*sqrt(13.0)*real_spher_harm_ucvec(9,6,-1,u)
             else               ! 9 2
                klm=0.5*sqrt(17.0/6.0)*real_spher_harm_ucvec(9,4,-1,u)
                klm=klm-0.5*sqrt(7.0/6.0)*real_spher_harm_ucvec(9,8,-1,u)
             end if
          case(10)
             if (m == 1) then   ! 10 1
                klm=0.125*sqrt(65.0/6.0)*real_spher_harm_ucvec(10,0,1,u)
                klm=klm-0.25*sqrt(11.0/2.0)*real_spher_harm_ucvec(10,4,1,u)
                klm=klm-0.125*sqrt(187.0/6.0)*real_spher_harm_ucvec(10,8,1,u)
             else               ! 10 2
                klm=0.125*sqrt(247.0/6.0)*real_spher_harm_ucvec(10,2,1,u)
                klm=klm+(1.0/16.0)*sqrt(19.0/3.0)*real_spher_harm_ucvec(10,6,1,u)
                klm=klm-(1.0/16.0)*sqrt(85.0)*real_spher_harm_ucvec(10,10,1,u)
             end if
          case(12)
             ! from here only the cubic harmonics for m3m symmetry are calculated
             ! coefficients from F.M. Mueller and M.G. Priestley, Phys Rev 148, 638 (1966)
             if (m == 1) then   ! 12 1
                klm=    0.69550265*real_spher_harm_ucvec(12,0,1,u)
                klm=klm+0.31412555*real_spher_harm_ucvec(12,4,1,u)
                klm=klm+0.34844954*real_spher_harm_ucvec(12,8,1,u)
                klm=klm+0.54422797*real_spher_harm_ucvec(12,12,1,u)
             else               ! 12 2
                klm=    0.55897937*real_spher_harm_ucvec(12,4,1,u)
                klm=klm-0.80626751*real_spher_harm_ucvec(12,8,1,u)
                klm=klm+0.19358400*real_spher_harm_ucvec(12,12,1,u)
             end if
          case(14)
             ! 12 1
             klm=    0.44009645*real_spher_harm_ucvec(14,0,1,u)
             klm=klm-0.45768183*real_spher_harm_ucvec(14,4,1,u)
             klm=klm-0.49113230*real_spher_harm_ucvec(14,8,1,u)
             klm=klm-0.59634848*real_spher_harm_ucvec(14,12,1,u)
          case(16)
             if (m == 1) then   ! 16 1
                klm=    0.68136168*real_spher_harm_ucvec(16,0,1,u)
                klm=klm+0.27586801*real_spher_harm_ucvec(16,4,1,u)
                klm=klm+0.29048987*real_spher_harm_ucvec(16,8,1,u)
                klm=klm+0.32756975*real_spher_harm_ucvec(16,12,1,u)
                klm=klm+0.51764542*real_spher_harm_ucvec(16,16,1,u)
             else               ! 16 2
                klm=    0.63704821*real_spher_harm_ucvec(16,4,1,u)
                klm=klm-0.32999033*real_spher_harm_ucvec(16,8,1,u)
                klm=klm-0.64798073*real_spher_harm_ucvec(16,12,1,u)
                klm=klm+0.25572816*real_spher_harm_ucvec(16,16,1,u)
             end if
          case(18)
             if (m == 1) then   ! 18 1
                klm=    0.45791513*real_spher_harm_ucvec(18,0,1,u)
                klm=klm-0.38645598*real_spher_harm_ucvec(18,4,1,u)
                klm=klm-0.40209462*real_spher_harm_ucvec(18,8,1,u)
                klm=klm-0.43746593*real_spher_harm_ucvec(18,12,1,u)
                klm=klm-0.53657149*real_spher_harm_ucvec(18,16,1,u)
             else               ! 18 2
               klm=    0.14872751*real_spher_harm_ucvec(18,4,1,u)
               klm=klm-0.63774601*real_spher_harm_ucvec(18,8,1,u)
               klm=klm+0.72334167*real_spher_harm_ucvec(18,12,1,u)
               klm=klm-0.21894515*real_spher_harm_ucvec(18,16,1,u)
             end if
          case(20)
             if (m == 1) then   ! 20 1
                klm=    0.67141495*real_spher_harm_ucvec(20,0,1,u)
                klm=klm+0.24982619*real_spher_harm_ucvec(20,4,1,u)
                klm=klm+0.25782846*real_spher_harm_ucvec(20,8,1,u)
                klm=klm+0.27469333*real_spher_harm_ucvec(20,12,1,u)
                klm=klm+0.31248919*real_spher_harm_ucvec(20,16,1,u)
                klm=klm+0.49719956*real_spher_harm_ucvec(20,20,1,u)
             else               ! 20 2
                klm=    0.66299538*real_spher_harm_ucvec(20,4,1,u)
                klm=klm-0.11295259*real_spher_harm_ucvec(20,8,1,u)
                klm=klm-0.42738441*real_spher_harm_ucvec(20,12,1,u)
                klm=klm-0.52810433*real_spher_harm_ucvec(20,16,1,u)
                klm=klm+0.29347435*real_spher_harm_ucvec(20,20,1,u)
             end if
          case default
             klm=0.0
       end select

       return
    End Function Cubic_Harm_Ucvec

    !!--++
    !!--++ Function Envj(N,X) Result(Y)
    !!--++    integer, intent(in)       :: n
    !!--++    real(Kind=dp), intent(in) :: x
    !!--++    real(Kind=dp)             :: y
    !!--++
    !!--++    (PRIVATE)
    !!--++
    !!--++ Update: February - 2005
    !!
    Pure Function Envj(N,X) Result(Y)
       !---- Arguments ----!
       integer,       intent(in) :: n
       real(Kind=dp), intent(in) :: x
       real(Kind=dp)             :: y

       y=0.5_dp*log10(6.28_dp*real(n,kind=dp))-n*log10(1.36_dp*x/real(n,kind=dp))

       return
    End Function Envj

    !!----
    !!---- Function Int_Slater_Bessel(N,L,Z,S) Result(Int_Slater_Besselv)
    !!----    integer,           intent(in) :: n
    !!----    integer,           intent(in) :: l
    !!----    real(kind=cp),     intent(in) :: z
    !!----    real(kind=cp),     intent(in) :: s
    !!----    real(kind=cp)                 :: Int_Slater_Besselv
    !!--<<
    !!----    Returns the integral:
    !!----    Int[0,inf]{r**(n+2)*exp(-chi*r)*j_l(s*r)} dr
    !!----    where j_l is the spherical Bessel function of order l
    !!----    only -1 <= n  and 0 <= l <= n+1
    !!-->>
    !!----
    !!---- Update: February - 2005
    !!
    Pure Function Int_Slater_Bessel(N,L,Z,S) Result(Int_Slater_Besselv)
       !---- arguments ----!
       integer,       intent(in) :: n, l
       real(kind=cp), intent(in) :: z, s
       real(kind=cp)             :: int_slater_besselv

       !---- Local Variables ----!
       real(kind=cp)  :: isb0, tmp
       integer        :: fact, i

       Int_Slater_Besselv =0.0
       if (n >= -1 .and. l<=n+1 .and. l>=0) then
          if (l == 0) then
             fact=1
          else
             fact= 1
             do i=1,l
                fact=fact*i
             end do
          end if
          if (l == n+1) then
             Int_Slater_Besselv = (2.0_cp*s/(s*s+z*z))**l*fact/(s*s+z*z)
             return
          end if
          if (l == n) then
             Int_Slater_Besselv = real(2*l+2,kind=cp)*z* (2.0_cp*s/(s*s+z*z))**l*fact/(s*s+z*z)/(s*s+z*z)
             return
          end if

          !---- We arrive here iif l>=n+1
          !---- We calculate the case n=l-1
          isb0 = (2.0_cp*s/(s*s+z*z))**l*fact/(s*s+z*z)

          !---- We calculate the case n=l
          Int_Slater_Besselv = real(2*l+2,kind=cp)*z* (2.0_cp*s/(s*s+z*z))**l*fact/(s*s+z*z)/(s*s+z*z)

          !---- And apply the recursivity up to n
          do i=l+1,n
             tmp= Int_Slater_Besselv
             Int_Slater_Besselv = (2.0_cp*(i+1)*z*Int_Slater_Besselv - &
                                  real(i+l+1,kind=cp)*real(i-l,kind=cp)*ISB0)/(s*s+z*z)
             isb0= tmp
          end do

          return
       end if

       return
    End Function Int_Slater_Bessel

    !!--++
    !!--++ Function Plgndr(L,M,X) Result(Plmx)
    !!--++    integer,       intent(in) :: l
    !!--++    integer,       intent(in) :: m
    !!--++    real(kind=cp), intent(in) :: x
    !!--++    real(kind=cp)             :: plmx
    !!--++
    !!--++    (PRIVATE)
    !!--++    Compute the Legendre Polynomial Pml(x). Here m and l are
    !!--++    integers satisfying 0 <= m <= l, while x lies in the range
    !!--++    -1 <= x <= 1.
    !!--++
    !!--++ Update: February - 2005
    !!
    Pure Function Plgndr(l,m,x) result(plmx)
       !---- Arguments ----!
       integer,      intent (in) :: l,m
       real(kind=cp),intent (in) :: x
       real(kind=cp)             :: plmx

       !---- Local variables ----!
       integer       :: i, ll
       real(kind=cp) :: fact, pll, pmm, pmmp1, somx2

       if (m < 0 .or. m > l .or. abs(x) > 1.0) then
          plmx=0.0
          return
       end if

       pmm=1.0
       if (m > 0) then                  !Compute P(m,m)
          somx2=sqrt((1.0_cp-x)*(1.0_cp+x))
          fact=1.0_cp
          do i=1,m
             pmm=-pmm*fact*somx2
             fact=fact+2.0_cp
          end do
       end if

       if (l == m) then
          plmx=pmm
       else
          pmmp1=x*real(2*m+1,kind=cp)*pmm           !Compute P(m,m+1)
          if (l == m+1) then
             plmx=pmmp1                 !Compute P(m,L), L > m+1
          else
             do ll=m+2,l
                pll=(x*real(2*ll-1,kind=cp)*pmmp1-real(ll+m-1,kind=cp)*pmm)/real(ll-m,kind=cp)
                pmm=pmmp1
                pmmp1=pll
             end do
             plmx=pll
          end if
       end if

       return
    End Function Plgndr

    !!----
    !!---- FUNCTION Real_Spher_Harm_Ang(L,M,P,Theta,Phi) Result(Ylmp)
    !!----    integer,       intent(in) :: l         !  In -> Index l >= 0
    !!----    integer,       intent(in) :: m         !  In -> Index m <= l
    !!----    integer,       intent(in) :: p         !  In -> +1: cosinus, -1: sinus
    !!----    real(kind=cp), intent(in) :: theta     !  In -> Spherical coordinate in degree
    !!----    real(kind=cp), intent(in) :: phi       !  In -> Spherical coordinate in degree
    !!----    real(kind=cp)             :: ylmp      ! Out -> Value of ylmn(theta,phi)
    !!----
    !!----    Real Spherical Harmonics: M.Kara & K. Kurki-Suonio, Acta Cryt. A37, 201 (1981)
    !!----    Input spherical coordinates Theta & Phi in degrees
    !!----
    !!---- Update: February - 2005
    !!
    Pure Function Real_Spher_Harm_Ang(l,m,p,theta,phi) result(ylmp)
       !---- Arguments ----!
       integer,      intent (in) :: l,m,p
       real(kind=cp),intent (in) :: theta,phi
       real(kind=cp)             :: ylmp

       !---- Local Variables ----!
       real(kind=dp)             :: norm
       real(kind=cp)             :: pphi,x
       integer                   :: i

       x=cos(theta*to_rad)
       if (p > 0) then
          pphi = cos(m*phi*to_rad)
       else
          pphi = sin(m*phi*to_rad)
       end if
       ylmp=plgndr(l,m,x)*pphi
       norm=real(2*l+1,kind=dp)

       do i=l-m+1,l+m
          norm=norm/real(i,kind=dp)
       end do
       if (m == 0) then
          norm=sqrt(norm/(4.0_dp*pi))
       else
          norm=sqrt(norm/(2.0_dp*pi))
       end if
       ylmp=norm*ylmp

       return
    End Function Real_Spher_Harm_Ang

    !!----
    !!---- Function Real_Spher_Harm_Ucvec(L,M,P,U) Result(Ylmp)
    !!----    integer,                    intent(in) :: l       !  In -> Index l >= 0
    !!----    integer,                    intent(in) :: m       !  In -> Index m <= l
    !!----    integer,                    intent(in) :: p       !  In -> +1: cosinus, -1: sinus
    !!----    real(kind=cp),dimension(3), intent(in) :: u       !  In -> unit vector in cartesian coordinates
    !!----    real(kind=cp)                          :: ylmp    ! Out -> Value of ylmn(u)
    !!--<<
    !!----    Real Spherical Harmonics: M.Kara & K. Kurki-Suonio, Acta Cryt. A37, 201 (1981)
    !!----    Input U: unit vector in cartesian coordinates. If the provided vector is not
    !!----    normalized to unit, this version calculates ylmp along the unit vector defined by
    !!----    the given vector that is not modified on output
    !!-->>
    !!---- Update: February - 2005, June - 2014 (JRC, simplify the calculation when l=m=0)
    !!
    Pure Function Real_Spher_Harm_Ucvec(l,m,p,u) result(ylmp)
       !---- Arguments ----!
       integer,                    intent (in) :: l,m,p
       real(kind=cp),dimension(3), intent (in) :: u
       real(kind=cp)                           :: ylmp

       !---- Local Variables ----!
       real(kind=dp)             :: norm
       real(kind=dp),parameter   :: sqrt_pi=1.7724538509055160272981674833411 !sqrt(pi)
       real(kind=cp)             :: pphi,x,ss
       integer                   :: i
       real(kind=cp),dimension(3):: v

       if(l == m .and. m == 0) then !for trivial zero-order give the value and return
         ylmp=0.5_dp/sqrt_pi
         return
       end if
       v=u
       ss=dot_product(v,v)
       if (abs(ss-1.0_cp) > 2.0_cp*eps) then  !Test the provided unit vector
          if(ss > eps) then
            v=v/sqrt(ss)             !Normalize the input vector
          else
            ylmp=0.0_cp
            return
          end if
       end if
       x=v(3)                    !costheta
       if (abs(x-1.0_cp) < eps .or. abs(x+1.0_cp) < eps) then
          pphi=0.0
       else
          pphi=atan2(v(2),v(1))  !This is the good function to obtain the true spherical angles
       end if
       if (p > 0) then
          pphi = cos(m*pphi)     !cos(m*phi)
       else
          pphi = sin(m*pphi)     !sin(m*phi)
       end if
       ylmp=plgndr(l,m,x)*pphi
       norm=real(2*l+1,kind=dp)

       do i=l-m+1,l+m  !this is not executed when l=m=0 !!!
          norm=norm/real(i,kind=dp)
       end do
       if (m == 0) then
          norm=sqrt(norm)/(2.0_dp*sqrt_pi)
       else
          norm=sqrt(norm/2.0_dp)/sqrt_pi
       end if
       ylmp=norm*ylmp

       return
    End Function Real_Spher_Harm_Ucvec

    !!----
    !!---- Function Real_Spher_HarmCharge_Ucvec(L,M,P,U) Result(Dlmp)
    !!----    integer,                    intent(in) :: l       !  In -> Index l >= 0
    !!----    integer,                    intent(in) :: m       !  In -> Index m <= l
    !!----    integer,                    intent(in) :: p       !  In -> +1: cosinus, -1: sinus
    !!----    real(kind=cp),dimension(3), intent(in) :: u       !  In -> unit vector in cartesian coordinates
    !!----    real(kind=cp)                          :: dlmp    ! Out -> Value of dlmn(u)
    !!--<<
    !!----    This function difers from the previous one (real spherical harmonics) only in the
    !!----    normalization constants, so that: Dlmp= Clmp Ylmp
    !!----    The constants Clmp are selected so that Int[abs(Dlmp) dOmega] = 2 -d(l,0)
    !!----    They have been calculated using Gaussian/Legendre quadrature.
    !!-->>
    !!---- Update: April - 2009
    !!
    Pure Function Real_Spher_HarmCharge_Ucvec(L,M,P,U) Result(Dlmp)
       !---- Arguments ----!
       integer,                    intent (in) :: l,m,p
       real(kind=cp),dimension(3), intent (in) :: u
       real(kind=cp)                           :: Dlmp

       !---- Local Variables ----!
       integer                   :: i
       real(kind=dp),dimension(0:10,0:10,2) :: Clmp
       if(l > 10 .or. m > 10 .or. m > l) then
          Dlmp=0.0_cp
          return
       end if
       !Assign values to the normalization constants
       Clmp = 0.0
       Clmp(  0,  0,  1)=     0.28209481_dp
       Clmp(  1,  1,  2)=     0.65147001_dp
       Clmp(  1,  0,  1)=     0.65147001_dp
       Clmp(  1,  1,  1)=     0.65146995_dp
       Clmp(  2,  2,  2)=     0.68646848_dp
       Clmp(  2,  1,  2)=     0.68646836_dp
       Clmp(  2,  0,  1)=     0.65552908_dp
       Clmp(  2,  1,  1)=     0.68646836_dp
       Clmp(  2,  2,  1)=     0.68646848_dp
       Clmp(  3,  3,  2)=     0.71929127_dp
       Clmp(  3,  2,  2)=     0.69189525_dp
       Clmp(  3,  1,  2)=     0.70087820_dp
       Clmp(  3,  0,  1)=     0.65613419_dp
       Clmp(  3,  1,  1)=     0.70087814_dp
       Clmp(  3,  2,  1)=     0.69189525_dp
       Clmp(  3,  3,  1)=     0.71929145_dp
       Clmp(  4,  4,  2)=     0.74899846_dp
       Clmp(  4,  3,  2)=     0.70616245_dp
       Clmp(  4,  2,  2)=     0.69879586_dp
       Clmp(  4,  1,  2)=     0.70847464_dp
       Clmp(  4,  0,  1)=     0.65620995_dp
       Clmp(  4,  1,  1)=     0.70847464_dp
       Clmp(  4,  2,  1)=     0.69879586_dp
       Clmp(  4,  3,  1)=     0.70616263_dp
       Clmp(  4,  4,  1)=     0.74899834_dp
       Clmp(  5,  5,  2)=     0.77591366_dp
       Clmp(  5,  4,  2)=     0.72266090_dp
       Clmp(  5,  3,  2)=     0.70547515_dp
       Clmp(  5,  2,  2)=     0.70407319_dp
       Clmp(  5,  1,  2)=     0.71306634_dp
       Clmp(  5,  0,  1)=     0.65617186_dp
       Clmp(  5,  1,  1)=     0.71306634_dp
       Clmp(  5,  2,  1)=     0.70407319_dp
       Clmp(  5,  3,  1)=     0.70547533_dp
       Clmp(  5,  4,  1)=     0.72266084_dp
       Clmp(  5,  5,  1)=     0.77591383_dp
       Clmp(  6,  6,  2)=     0.80047959_dp
       Clmp(  6,  5,  2)=     0.73945135_dp
       Clmp(  6,  4,  2)=     0.71554828_dp
       Clmp(  6,  3,  2)=     0.70703226_dp
       Clmp(  6,  2,  2)=     0.70801049_dp
       Clmp(  6,  1,  2)=     0.71609598_dp
       Clmp(  6,  0,  1)=     0.65610999_dp
       Clmp(  6,  1,  1)=     0.71609592_dp
       Clmp(  6,  2,  1)=     0.70801049_dp
       Clmp(  6,  3,  1)=     0.70703244_dp
       Clmp(  6,  4,  1)=     0.71554822_dp
       Clmp(  6,  5,  1)=     0.73945147_dp
       Clmp(  6,  6,  1)=     0.80047959_dp
       Clmp(  7,  7,  2)=     0.82308108_dp
       Clmp(  7,  6,  2)=     0.75586903_dp
       Clmp(  7,  5,  2)=     0.72696197_dp
       Clmp(  7,  4,  2)=     0.71344930_dp
       Clmp(  7,  3,  2)=     0.70896560_dp
       Clmp(  7,  2,  2)=     0.71099448_dp
       Clmp(  7,  1,  2)=     0.71822095_dp
       Clmp(  7,  0,  1)=     0.65604913_dp
       Clmp(  7,  1,  1)=     0.71822089_dp
       Clmp(  7,  2,  1)=     0.71099448_dp
       Clmp(  7,  3,  1)=     0.70896578_dp
       Clmp(  7,  4,  1)=     0.71344924_dp
       Clmp(  7,  5,  1)=     0.72696209_dp
       Clmp(  7,  6,  1)=     0.75586903_dp
       Clmp(  7,  7,  1)=     0.82308108_dp
       Clmp(  8,  8,  2)=     0.84402764_dp
       Clmp(  8,  7,  2)=     0.77168250_dp
       Clmp(  8,  6,  2)=     0.73882592_dp
       Clmp(  8,  5,  2)=     0.72154719_dp
       Clmp(  8,  4,  2)=     0.71311980_dp
       Clmp(  8,  3,  2)=     0.71081281_dp
       Clmp(  8,  2,  2)=     0.71331161_dp
       Clmp(  8,  1,  2)=     0.71977973_dp
       Clmp(  8,  0,  1)=     0.65599412_dp
       Clmp(  8,  1,  1)=     0.71977967_dp
       Clmp(  8,  2,  1)=     0.71331161_dp
       Clmp(  8,  3,  1)=     0.71081293_dp
       Clmp(  8,  4,  1)=     0.71311975_dp
       Clmp(  8,  5,  1)=     0.72154731_dp
       Clmp(  8,  6,  1)=     0.73882592_dp
       Clmp(  8,  7,  1)=     0.77168250_dp
       Clmp(  8,  8,  1)=     0.84402764_dp
       Clmp(  9,  9,  2)=     0.86356676_dp
       Clmp(  9,  8,  2)=     0.78682709_dp
       Clmp(  9,  7,  2)=     0.75072247_dp
       Clmp(  9,  6,  2)=     0.73046678_dp
       Clmp(  9,  5,  2)=     0.71902418_dp
       Clmp(  9,  4,  2)=     0.71349597_dp
       Clmp(  9,  3,  2)=     0.71247119_dp
       Clmp(  9,  2,  2)=     0.71514851_dp
       Clmp(  9,  1,  2)=     0.72096473_dp
       Clmp(  9,  0,  1)=     0.65595061_dp
       Clmp(  9,  1,  1)=     0.72096467_dp
       Clmp(  9,  2,  1)=     0.71514851_dp
       Clmp(  9,  3,  1)=     0.71247137_dp
       Clmp(  9,  4,  1)=     0.71349585_dp
       Clmp(  9,  5,  1)=     0.71902430_dp
       Clmp(  9,  6,  1)=     0.73046678_dp
       Clmp(  9,  7,  1)=     0.75072247_dp
       Clmp(  9,  8,  1)=     0.78682709_dp
       Clmp(  9,  9,  1)=     0.86356682_dp
       Clmp( 10, 10,  2)=     0.88188988_dp
       Clmp( 10,  9,  2)=     0.80130792_dp
       Clmp( 10,  8,  2)=     0.76244813_dp
       Clmp( 10,  7,  2)=     0.73974365_dp
       Clmp( 10,  6,  2)=     0.72590035_dp
       Clmp( 10,  5,  2)=     0.71787411_dp
       Clmp( 10,  4,  2)=     0.71415150_dp
       Clmp( 10,  3,  2)=     0.71390396_dp
       Clmp( 10,  2,  2)=     0.71663243_dp
       Clmp( 10,  1,  2)=     0.72189075_dp
       Clmp( 10,  0,  1)=     0.65590489_dp
       Clmp( 10,  1,  1)=     0.72189069_dp
       Clmp( 10,  2,  1)=     0.71663243_dp
       Clmp( 10,  3,  1)=     0.71390414_dp
       Clmp( 10,  4,  1)=     0.71415144_dp
       Clmp( 10,  5,  1)=     0.71787423_dp
       Clmp( 10,  6,  1)=     0.72590035_dp
       Clmp( 10,  7,  1)=     0.73974365_dp
       Clmp( 10,  8,  1)=     0.76244813_dp
       Clmp( 10,  9,  1)=     0.80130792_dp
       Clmp( 10, 10,  1)=     0.88188988_dp
       i=p
       if( i < 1) i = 2

       Dlmp=Clmp(l,m,i)*Real_Spher_Harm_Ucvec(L,M,P,U)

       return
    End Function Real_Spher_HarmCharge_Ucvec

    !!--++
    !!--++ Function Start1(X,Mp) Result (Start)
    !!--++    real(kind=dp), intent(in) :: x        !  In -> Argument of Jn(x)
    !!--++    integer, intent(in)       :: mp       !  In -> Value of magnitude
    !!--++    integer                   :: start    ! Out -> Starting point
    !!--++
    !!--++    (PRIVATE)
    !!--++    Determine the starting point for backward
    !!--++    recurrence such that the magnitude of Jn(x) at that point is
    !!--++    about 10^(-MP).
    !!--++
    !!--++ Update: February - 2005
    !!
    Pure Function Start1(X,Mp) Result (Start)
       !---- Arguments ----!
       real(kind=dp), intent(in) :: x
       integer, intent(in)       :: mp
       integer                   :: start

       !---- Local variables ----!
       integer      :: n1,n0,nn, it
       real(kind=dp):: f,f0,f1,a0

       a0=abs(x)
       n0=int(1.1_dp*a0)+1
       f0=envj(n0,a0)-mp
       n1=n0+5
       f1=envj(n1,a0)-mp
       do it=1,20
          nn=n1-(n1-n0)/(1.0_dp-f0/f1)
          f=envj(nn,a0)-mp
          if (abs(nn-n1) < 1) exit
          n0=n1
          f0=f1
          n1=nn
          f1=f
       end do
       start=nn

       return
    End Function Start1

    !!--++
    !!--++ Function Start2(X,N,Mp) Result (Start)
    !!--++    real(kind=dp), intent(in) :: x      !  In -> Argument of Jn(x)
    !!--++    integer, intent(in)       :: n      !  In -> Order of Jn(x)
    !!--++    integer, intent(in)       :: mp     !  In -> Significant digit
    !!--++    integer                   :: start  ! Out -> Starting point
    !!--++
    !!--++    (PRIVATE)
    !!--++    Determine the starting point for backward
    !!--++    recurrence such that all Jn(x) has MP significants digits
    !!--++
    !!--++ Update: February - 2005
    !!
    Pure Function Start2(X,N,Mp) Result(Start)
       !---- Arguments ----!
       real(kind=dp), intent(in) :: x
       integer,       intent(in) :: n
       integer,       intent(in) :: mp
       integer                   :: start

       !---- Local variables ----!
       real(kind=dp) :: a0, hmp, ejn, obj,f,f0,f1
       integer       :: n0,n1,nn, it

       a0=abs(x)
       hmp=0.5_dp*mp
       ejn=envj(n,a0)
       if (ejn <= hmp) then
          obj=mp
          n0=int(1.1_dp*a0)+1  ! +1 was absent in the original version ... this solves the problem of
       else                    ! Intel, gfortran and g95 compilers ... Lahey was calculating well event if n0=0!
          obj=hmp+ejn
          n0=n
       end if
       f0=envj(n0,a0)-obj
       n1=n0+5
       f1=envj(n1,a0)-obj
       do it=1,20
          nn=n1-(n1-n0)/(1.0_dp-f0/f1)
          f=envj(nn,a0)-obj
          if (abs(nn-n1) < 1) exit
          n0=n1
          f0=f1
          n1=nn
          f1=f
       end do
       start=nn+10

       return
    End Function Start2

    !---------------------!
    !---- Subroutines ----!
    !---------------------!

    !!----
    !!---- Subroutine Init_Err_Spher()
    !!----
    !!----    Initialize the errors flags in CFML_Spherical_Harmonics
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Init_Err_Spher()

       ERR_Spher=.false.
       ERR_Spher_Mess=" "

       return
    End Subroutine Init_Err_Spher

    !!----
    !!---- Subroutine Pikout_Lj_Cubic(Group,Lj,Ncoef,Lun)
    !!----    character (len=*), intent (in)        :: group   !  In ->
    !!----    integer, dimension(2,11), intent(out) :: lj      ! Out ->
    !!----    integer, intent (out)                 :: ncoef   ! Out ->
    !!----    integer, optional, intent (in)        :: lun     !  In ->
    !!--<<
    !!----    Picking out rules for indices of cubic harmonics for the 5 cubic groups
    !!----    Only up to tenth order Given in Table 4 of reference M.Kara &
    !!----    K. Kurki-Suonio, Acta Cryt. A37, 201 (1981)
    !!-->>
    !!---- Update: February - 2005
    !!
    Subroutine Pikout_Lj_Cubic(Group,Lj,Ncoef,Lun)
       !---- Arguments ----!
       character (len=*),        intent (in) :: group
       integer, dimension(2,11), intent(out) :: lj
       integer, intent (out)                 :: ncoef
       integer, optional, intent (in)        :: lun

       !---- Local variables ----!
       integer :: k,lc

       lc=-1
       call Init_Err_Spher()

       if (present(lun)) lc=lun

       Select case (group)
          Case("23")
             ncoef=11
             lj(1,1) =0
             lj(2,1) =1
             lj(1,2) =3
             lj(2,2) =1
             lj(1,3) =4
             lj(2,3) =1
             lj(1,4) =6
             lj(2,4) =1
             lj(1,5) =6
             lj(2,5) =2
             lj(1,6) =7
             lj(2,6) =1
             lj(1,7) =8
             lj(2,7) =1
             lj(1,8) =9
             lj(2,8) =1
             lj(1,9) =9
             lj(2,9) =2
             lj(1,10)=10
             lj(2,10)=1
             lj(1,11)=10
             lj(2,11)=2
             if (lc > 0) then
                write(unit=lc,fmt="(a)") " => Allowed cubic harmonics for point group < T, 23 >"
                write(unit=lc,fmt="(a)") " => Eleven coefficients"
                write(unit=lc,fmt="(11(a,2i3,a))") ( "  (",lj(1,k),lj(2,k),")",  k=1,ncoef)
             end if
          Case("m3")
             ncoef=7
             lj(1,1) =0
             lj(2,1) =1
             lj(1,2) =4
             lj(2,2) =1
             lj(1,3) =6
             lj(2,3) =1
             lj(1,4) =6
             lj(2,4) =2
             lj(1,5) =8
             lj(2,5) =1
             lj(1,6) =10
             lj(2,6) =1
             lj(1,7) =10
             lj(2,7) =2
             if (lc > 0) then
                write(unit=lc,fmt="(a)") " => Allowed cubic harmonics for point group < Th, m3 >"
                write(unit=lc,fmt="(a)") " => Seven coefficients"
                write(unit=lc,fmt="(7(a,2i3,a))") ( "  (",lj(1,k),lj(2,k),")",  k=1,ncoef)
             end if
          Case("432")
             ncoef=6
             lj(1,1) =0
             lj(2,1) =1
             lj(1,2) =4
             lj(2,2) =1
             lj(1,3) =6
             lj(2,3) =1
             lj(1,4) =8
             lj(2,4) =1
             lj(1,5) =9
             lj(2,5) =2
             lj(1,6) =10
             lj(2,6) =1
             if (lc > 0) then
                write(unit=lc,fmt="(a)") " => Allowed cubic harmonics for point group < O, 432 >"
                write(unit=lc,fmt="(a)") " => Six coefficients"
                write(unit=lc,fmt="(6(a,2i3,a))") ( "  (",lj(1,k),lj(2,k),")",  k=1,ncoef)
             end if
          Case("-43m")
             ncoef=8
             lj(1,1) =0
             lj(2,1) =1
             lj(1,2) =3
             lj(2,2) =1
             lj(1,3) =4
             lj(2,3) =1
             lj(1,4) =6
             lj(2,4) =1
             lj(1,5) =7
             lj(2,5) =1
             lj(1,6) =8
             lj(2,6) =1
             lj(1,7) =9
             lj(2,7) =1
             lj(1,8) =10
             lj(2,8) =1
             if (lc > 0) then
                write(unit=lc,fmt="(a)") " => Allowed cubic harmonics for point group < Td, -43m >"
                write(unit=lc,fmt="(a)") " => Eight coefficients"
                write(unit=lc,fmt="(8(a,2i3,a))") ( "  (",lj(1,k),lj(2,k),")",  k=1,ncoef)
             end if
          case("m3m")
             ncoef=5
             lj(1,1) =0
             lj(2,1) =1
             lj(1,2) =4
             lj(2,2) =1
             lj(1,3) =6
             lj(2,3) =1
             lj(1,4) =8
             lj(2,4) =1
             lj(1,5) =10
             lj(2,5) =1
             if (lc > 0) then
                write(unit=lc,fmt="(a)") " => Allowed cubic harmonics for point group < Oh, m3m >"
                write(unit=lc,fmt="(a)") " => Five coefficients"
                write(unit=lc,fmt="(5(a,2i3,a))") ( "  (",lj(1,k),lj(2,k),")",  k=1,ncoef)
             end if
          case default

             err_spher=.true.
             ERR_Spher_Mess=" Wrong cubic point group passed to subroutine: pikout_lj_cubic "

       end select

       return
    End Subroutine Pikout_Lj_Cubic

    !!----
    !!---- Subroutine Sphjn(N,X,Nm,Jn,Djn)
    !!----    integer,       intent(in)                   :: n    !  In -> Order of jn(x) (n=0,1,2,3,...)
    !!----    real(kind=dp), intent(in)                   :: x    !  In -> Argument of jn(x)
    !!----    integer, intent(out)                        :: nm   ! Out -> Highest order computed
    !!----    real(kind=dp), dimension(0:n), intent(out)  :: jn   ! Out -> Array with spherical Bessel functions jn(x)
    !!----    real(kind=dp), dimension(0:n), intent(out)  :: djn  ! Out -> Array with derivatives jn'(x)
    !!----
    !!----    Compute spherical Bessel functions jn(x) and their derivatives
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Sphjn(n,x,nm,jn,djn)
       !---- Arguments ----!
       integer,                       intent(in)  :: n   !Order of jn(x) (n=0,1,2,3,...)
       real(kind=dp),                 intent(in)  :: x   !Argument of jn(x)
       integer,                       intent(out) :: nm  !Highest order computed
       real(kind=dp), dimension(0:n), intent(out) :: jn  !array with spherical Bessel functions jn(x)
       real(kind=dp), dimension(0:n), intent(out) :: djn !array with derivatives jn'(x)

       !---- Local variables ----!
       integer       :: k,m
       real(kind=dp) :: sa,sb, f,f1,f0, cs

       nm=n
       if (abs(x) <= 1.0e-30_dp) then
          do k=0,n
             jn(k) = 0.0_dp
             djn(k)= 0.0_dp
          end do
          jn(0)=1.0_dp
          djn(1)=1.0_dp/3.0_dp
          return
       end if

       jn(0)=sin(x)/x
       jn(1)=(jn(0)-cos(x))/x
       if (n >= 2) then
          sa=jn(0)
          sb=jn(1)
          m=start1(x,200)
          if (m < n) then
             nm=m
          else
             m=start2(x,n,15)
          end if
          f0=0.0_dp
          f1=1.0e-100_dp
          do k=m,0,-1
             f=(2.0_dp*k+3.0_dp)*f1/x-f0
             if (k <= nm) jn(k)=f
             f0=f1
             f1=f
          end do
          if (abs(sa) > abs(sb)) cs=sa/f
          if (abs(sa) <= abs(sb)) cs=sb/f0
          do k=0,nm
             jn(k)=cs*jn(k)
          end do
       end if

       djn(0)=(cos(x)-sin(x)/x)/x
       do k=1,nm
          djn(k)=jn(k-1)-(k+1.0_dp)*jn(k)/x
       end do

       return
    End Subroutine Sphjn

 End Module CFML_Spherical_Harmonics




