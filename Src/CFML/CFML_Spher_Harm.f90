!!----
!!---- Copyleft(C) 1999-2009,              Version: 4.0
!!---- Juan Rodriguez-Carvajal & Javier Gonzalez-Platas
!!----
!!---- MODULE: CFML_Spherical_Harmonics
!!----   INFO: Spherical Harmonics routines
!!----
!!---- HISTORY
!!----    Update: January - 2005
!!----            January - 2000    Based in public codes. Created by JRC
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
    logical, public    :: ERR_Spher

    !!----
    !!---- ERR_Spher_Mess
    !!----    character(len=150), public :: ERR_Spher_Mess
    !!----
    !!----    String containing information about the last error
    !!----
    !!---- Update: February - 2005
    !!
    character(len=150), public :: ERR_Spher_Mess

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

       y=0.5_dp*log10(6.28_dp*n)-n*log10(1.36_dp*x/n)

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
             Int_Slater_Besselv = (2*s/(s*s+z*z))**l*fact/(s*s+z*z)
             return
          end if
          if (l == n) then
             Int_Slater_Besselv = (2*l+2)*z* (2*s/(s*s+z*z))**l*fact/(s*s+z*z)/(s*s+z*z)
             return
          end if

          !---- We arrive here iif l>=n+1
          !---- We calculate the case n=l-1
          isb0 = (2*s/(s*s+z*z))**l*fact/(s*s+z*z)

          !---- We calculate the case n=l
          Int_Slater_Besselv = (2*l+2)*z* (2*s/(s*s+z*z))**l*fact/(s*s+z*z)/(s*s+z*z)

          !---- And apply the recursivity up to n
          do i=l+1,n
             tmp= Int_Slater_Besselv
             Int_Slater_Besselv = (2.0*(i+1)*z*Int_Slater_Besselv-&
                                  (i+l+1.0)*(i-l)*ISB0)/(s*s+z*z)
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
          somx2=sqrt((1.0-x)*(1.0+x))
          fact=1.0
          do i=1,m
             pmm=-pmm*fact*somx2
             fact=fact+2.0
          end do
       end if

       if (l == m) then
          plmx=pmm
       else
          pmmp1=x*(2*m+1)*pmm           !Compute P(m,m+1)
          if (l == m+1) then
             plmx=pmmp1                 !Compute P(m,L), L > m+1
          else
             do ll=m+2,l
                pll=(x*real(2*ll-1)*pmmp1-real(ll+m-1)*pmm)/real(ll-m)
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
    !!----    real(kind=sp) Spherical Harmonics: M.Kara & K. Kurki-Suonio, Acta Cryt. A37, 201 (1981)
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
       norm=real(2*l+1)

       do i=l-m+1,l+m
          norm=norm/real(i)
       end do
       if (m == 0) then
          norm=sqrt(norm/(4.0*pi))
       else
          norm=sqrt(norm/(2.0*pi))
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
    !!---- Update: February - 2005
    !!
    Pure Function Real_Spher_Harm_Ucvec(l,m,p,u) result(ylmp)
       !---- Arguments ----!
       integer,                    intent (in) :: l,m,p
       real(kind=cp),dimension(3), intent (in) :: u
       real(kind=cp)                           :: ylmp

       !---- Local Variables ----!
       real(kind=dp)             :: norm
       real(kind=cp)             :: pphi,x,ss
       integer                   :: i
       real(kind=cp),dimension(3):: v

       v=u
       ss=dot_product(v,v)
       if (abs(ss-1.0) > 2.0*eps) then  !Test the provided unit vector
          if(ss > eps) then
            v=v/sqrt(ss)             !Normalize the input vector
          else
            ylmp=0.0_cp
            return
          end if
       end if
       x=v(3)                    !costheta
       if (abs(x-1.0) < eps .or. abs(x+1.0) < eps) then
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
       norm=real(2*l+1)

       do i=l-m+1,l+m
          norm=norm/real(i)
       end do
       if (m == 0) then
          norm=sqrt(norm/(4.0*pi))
       else
          norm=sqrt(norm/(2.0*pi))
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
       real,dimension(0:10,0:10,2) :: Clmp
       if(l > 10 .or. m > 10 .or. m > l) then
          Dlmp=0.0
          return
       end if
       !Assign values to the normalization constants
       Clmp = 0.0
       Clmp(  0,  0,  1)=     0.28209481
       Clmp(  1,  1,  2)=     0.65147001
       Clmp(  1,  0,  1)=     0.65147001
       Clmp(  1,  1,  1)=     0.65146995
       Clmp(  2,  2,  2)=     0.68646848
       Clmp(  2,  1,  2)=     0.68646836
       Clmp(  2,  0,  1)=     0.65552908
       Clmp(  2,  1,  1)=     0.68646836
       Clmp(  2,  2,  1)=     0.68646848
       Clmp(  3,  3,  2)=     0.71929127
       Clmp(  3,  2,  2)=     0.69189525
       Clmp(  3,  1,  2)=     0.70087820
       Clmp(  3,  0,  1)=     0.65613419
       Clmp(  3,  1,  1)=     0.70087814
       Clmp(  3,  2,  1)=     0.69189525
       Clmp(  3,  3,  1)=     0.71929145
       Clmp(  4,  4,  2)=     0.74899846
       Clmp(  4,  3,  2)=     0.70616245
       Clmp(  4,  2,  2)=     0.69879586
       Clmp(  4,  1,  2)=     0.70847464
       Clmp(  4,  0,  1)=     0.65620995
       Clmp(  4,  1,  1)=     0.70847464
       Clmp(  4,  2,  1)=     0.69879586
       Clmp(  4,  3,  1)=     0.70616263
       Clmp(  4,  4,  1)=     0.74899834
       Clmp(  5,  5,  2)=     0.77591366
       Clmp(  5,  4,  2)=     0.72266090
       Clmp(  5,  3,  2)=     0.70547515
       Clmp(  5,  2,  2)=     0.70407319
       Clmp(  5,  1,  2)=     0.71306634
       Clmp(  5,  0,  1)=     0.65617186
       Clmp(  5,  1,  1)=     0.71306634
       Clmp(  5,  2,  1)=     0.70407319
       Clmp(  5,  3,  1)=     0.70547533
       Clmp(  5,  4,  1)=     0.72266084
       Clmp(  5,  5,  1)=     0.77591383
       Clmp(  6,  6,  2)=     0.80047959
       Clmp(  6,  5,  2)=     0.73945135
       Clmp(  6,  4,  2)=     0.71554828
       Clmp(  6,  3,  2)=     0.70703226
       Clmp(  6,  2,  2)=     0.70801049
       Clmp(  6,  1,  2)=     0.71609598
       Clmp(  6,  0,  1)=     0.65610999
       Clmp(  6,  1,  1)=     0.71609592
       Clmp(  6,  2,  1)=     0.70801049
       Clmp(  6,  3,  1)=     0.70703244
       Clmp(  6,  4,  1)=     0.71554822
       Clmp(  6,  5,  1)=     0.73945147
       Clmp(  6,  6,  1)=     0.80047959
       Clmp(  7,  7,  2)=     0.82308108
       Clmp(  7,  6,  2)=     0.75586903
       Clmp(  7,  5,  2)=     0.72696197
       Clmp(  7,  4,  2)=     0.71344930
       Clmp(  7,  3,  2)=     0.70896560
       Clmp(  7,  2,  2)=     0.71099448
       Clmp(  7,  1,  2)=     0.71822095
       Clmp(  7,  0,  1)=     0.65604913
       Clmp(  7,  1,  1)=     0.71822089
       Clmp(  7,  2,  1)=     0.71099448
       Clmp(  7,  3,  1)=     0.70896578
       Clmp(  7,  4,  1)=     0.71344924
       Clmp(  7,  5,  1)=     0.72696209
       Clmp(  7,  6,  1)=     0.75586903
       Clmp(  7,  7,  1)=     0.82308108
       Clmp(  8,  8,  2)=     0.84402764
       Clmp(  8,  7,  2)=     0.77168250
       Clmp(  8,  6,  2)=     0.73882592
       Clmp(  8,  5,  2)=     0.72154719
       Clmp(  8,  4,  2)=     0.71311980
       Clmp(  8,  3,  2)=     0.71081281
       Clmp(  8,  2,  2)=     0.71331161
       Clmp(  8,  1,  2)=     0.71977973
       Clmp(  8,  0,  1)=     0.65599412
       Clmp(  8,  1,  1)=     0.71977967
       Clmp(  8,  2,  1)=     0.71331161
       Clmp(  8,  3,  1)=     0.71081293
       Clmp(  8,  4,  1)=     0.71311975
       Clmp(  8,  5,  1)=     0.72154731
       Clmp(  8,  6,  1)=     0.73882592
       Clmp(  8,  7,  1)=     0.77168250
       Clmp(  8,  8,  1)=     0.84402764
       Clmp(  9,  9,  2)=     0.86356676
       Clmp(  9,  8,  2)=     0.78682709
       Clmp(  9,  7,  2)=     0.75072247
       Clmp(  9,  6,  2)=     0.73046678
       Clmp(  9,  5,  2)=     0.71902418
       Clmp(  9,  4,  2)=     0.71349597
       Clmp(  9,  3,  2)=     0.71247119
       Clmp(  9,  2,  2)=     0.71514851
       Clmp(  9,  1,  2)=     0.72096473
       Clmp(  9,  0,  1)=     0.65595061
       Clmp(  9,  1,  1)=     0.72096467
       Clmp(  9,  2,  1)=     0.71514851
       Clmp(  9,  3,  1)=     0.71247137
       Clmp(  9,  4,  1)=     0.71349585
       Clmp(  9,  5,  1)=     0.71902430
       Clmp(  9,  6,  1)=     0.73046678
       Clmp(  9,  7,  1)=     0.75072247
       Clmp(  9,  8,  1)=     0.78682709
       Clmp(  9,  9,  1)=     0.86356682
       Clmp( 10, 10,  2)=     0.88188988
       Clmp( 10,  9,  2)=     0.80130792
       Clmp( 10,  8,  2)=     0.76244813
       Clmp( 10,  7,  2)=     0.73974365
       Clmp( 10,  6,  2)=     0.72590035
       Clmp( 10,  5,  2)=     0.71787411
       Clmp( 10,  4,  2)=     0.71415150
       Clmp( 10,  3,  2)=     0.71390396
       Clmp( 10,  2,  2)=     0.71663243
       Clmp( 10,  1,  2)=     0.72189075
       Clmp( 10,  0,  1)=     0.65590489
       Clmp( 10,  1,  1)=     0.72189069
       Clmp( 10,  2,  1)=     0.71663243
       Clmp( 10,  3,  1)=     0.71390414
       Clmp( 10,  4,  1)=     0.71415144
       Clmp( 10,  5,  1)=     0.71787423
       Clmp( 10,  6,  1)=     0.72590035
       Clmp( 10,  7,  1)=     0.73974365
       Clmp( 10,  8,  1)=     0.76244813
       Clmp( 10,  9,  1)=     0.80130792
       Clmp( 10, 10,  1)=     0.88188988
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
       n0=int(1.1*a0)+1
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
          n0=int(1.1*a0)
       else
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
       if (abs(x) <= 1.0e-100_dp) then
          do k=0,n
             jn(k)=0.0_dp
             djn(k)=0.0_dp
          end do
          jn(0)=1.0_dp
          djn(1)=0.3333333333333333_dp
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




