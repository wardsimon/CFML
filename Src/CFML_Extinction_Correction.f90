!!-------------------------------------------------------
!!---- Crystallographic Fortran Modules Library (CrysFML)
!!-------------------------------------------------------
!!---- The CrysFML project is distributed under LGPL. In agreement with the
!!---- Intergovernmental Convention of the ILL, this software cannot be used
!!---- in military applications.
!!----
!!---- Copyright (C) 1999-2015  Institut Laue-Langevin (ILL), Grenoble, FRANCE
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
!!---- MODULE: CFML_Extinction_correction
!!----   INFO: Small module for extinction corrections
!!----
!!---- HISTORY
!!----    Update: 13/05/2015
!!----
!!----
!!---- DEPENDENCIES
!!----
!!--++     Use CFML_GlobalDeps, only :: cp
!!----
!!---- VARIABLES
!!----    ERR_EXTINC
!!----    ERR_EXTINC_MESS
!!----
!!---- PUBLIC PROCEDURES
!!----    Functions:
!!----       AG_THETA
!!----       AL_THETA
!!----       BG_THETA
!!----       BL_THETA
!!----
!!----    Subroutines:
!!----       BECKER_COPPENS
!!----       SHELX_EXTINCTION,
!!----
!!
  Module CFML_Extinction_correction
    Use CFML_GlobalDeps, only : cp
    Implicit none
    private
    ! Public Functions
    public :: AG_theta,AL_theta,BG_theta,BL_theta
    ! Public Subroutines
    public :: SHELX_Extinction, Becker_Coppens
    ! Global variables
    logical,            public :: err_extinc  = .false.
    character(len=120), public :: err_extinc_mess=" "

    Contains

       Function AG_theta(cos2t) result (ag)   !A(theta) Gaussian
         real(kind=cp), intent (in) :: cos2t
         real(kind=cp)              :: ag
          ag = 0.58 + (0.48 + 0.24 * cos2t) * cos2t
         return
       End Function AG_theta

       Function AL_theta(cos2t) result (al)       !A(theta) Lorentzian
         real(kind=cp), intent (in) :: cos2t
         real(kind=cp)              :: al
          al = 0.025 + 0.285 * cos2t
         return
       End Function AL_theta

       Function BG_theta(cos2t) result (bg)   !B(theta) Gaussian
         real(kind=cp), intent (in) :: cos2t
         real(kind=cp)              :: bg
          bg = 0.02 - 0.025 * cos2t
         return
       End Function BG_theta

       Function BL_theta(cos2t) result (bl)   !B(theta) Lorentzian
         real(kind=cp), intent (in) :: cos2t
         real(kind=cp)              :: bl
         if(cos2t > 0.0) then
          bl = 0.15 - 0.2 * (0.75-cos2t)**2
         else
          bl =  -0.45 * cos2t
         end if
         return
       End Function BL_theta

       !!----  Subroutine Becker_Coppens(iext,f2,cext,r,g,ys,dydr,dydg)
       !!----     integer,                     intent (in) :: iext !1: Gaussian, 2: Lorentzian
       !!----     real(kind=cp),               intent (in) :: f2   !Square of structure factor
       !!----     real(kind=cp), dimension(4), intent (in) :: cext !coefficients for Becker-Coppens correction
       !!----     real(kind=cp),               intent (in) :: r,g  !Radius of blocks, width of mosaic distribution
       !!----     real(kind=cp),               intent(out) :: ys   !Extinction correction factor: Icorr = I.ys
       !!----     real(kind=cp), optional,     intent(out) :: dydr !Derivative of "ys" w.r.t "r"
       !!----     real(kind=cp), optional,     intent(out) :: dydg !Derivative of "ys" w.r.t "g"
       !!----
       !!----   The secondary extinction correction is calculated according to the equations developped
       !!----   by Becker&Coppens in Acta Cryst A30, 129 (1974).
       !!----   The fortran code (name of variables) is inspired from the CCSL extinction calculations.
       !!----   The coefficients of the Becker-Coppens correction for secondary extinctions are the
       !!----   same as those used in CCSL (P.J. Brown and J.C. Matthewman):
       !!----   The following coefficients are provided externally to the subroutine
       !!----   using the functions of this module A(Theta) and B(Theta) defined above
       !!----
       !!----    cext(1) = !Tbar*1000*Lambda^3/(V^2 sin2t) (CW)
       !!----            = !Tbar*1000*Lambda^4/(V^2 sint^2) (TOF)
       !!----    cext(2) = Lambda/sin2t (CW) (Lambda/sint)**2 (TOF)
       !!----    cext(3) = A(Theta)
       !!----    cext(4) = B(Theta)
       !!----
       !!----    The free parameters of the model are "r" and "g"
       !!----

      Subroutine Becker_Coppens(iext,f2,cext,r,g,ys,dydr,dydg)
         integer,                     intent (in) :: iext
         real(kind=cp),               intent (in) :: f2
         real(kind=cp), dimension(4), intent (in) :: cext
         real(kind=cp),               intent (in) :: r,g
         real(kind=cp),               intent(out) :: ys
         real(kind=cp), optional,     intent(out) :: dydr
         real(kind=cp), optional,     intent(out) :: dydg
       ! Local variables
         real(kind=cp) :: a, b, c, d, h, factor, x, x2, xx, c4, yy, der

          a = 1.5*r/cext(2)     ! r = domr in Jane's notation
                                ! a = 3/2 * r * sin2t/L = <alpha>
                                ! <alpha>= l* sin2t/L =>   l = 3/2 * r
          b = cext(1)*f2/1.5    ! b = 2/3*1000*L^3/(sin2t* V^2)*F^2
          c = 1.5*g             ! g = amosc in Jane's notation
          h = 2.0*g*g

          if(iext == 1) then          !Gaussian
              d=1.0/sqrt(1.0+ a*a /h)
          else if(iext == 2) then     !Lorentzian
              d=1.0/(1.0+a/c)
          end if
          x=b*a*d
          x2=2.0*x
          xx=x*x
          c4=1.0+cext(4)*x
          yy=1.0/(1.0+x2+cext(3)*xx/c4)
          ys=sqrt(yy)
          if(iext == 1) then          !Gaussian
            factor=a*d/h
          else if(iext == 2) then     !Lorentzian
            factor=1.0/c
          end if
          if(present(dydr) .and. present(dydg)) then
            der=-yy*ys*(1.0+cext(3)*(1.0-0.5*xx*cext(4)/c4)/c4)
            dydr = der * x*(1.0-a*d*factor)/r
            dydg = der * 1.5*factor*x*d*a/c
          end if
        return
      End Subroutine Becker_Coppens

      Subroutine SHELX_Extinction(job,iext,Lambda,ssnn,hkl,f2,extc,ys,der)
        integer,                    intent (in) :: job     ! =0,2 for x-rays, =1,3 for neutrons
        integer,                    intent (in) :: iext    ! Extinction model (1: isotropic, 2:anisotropic)
        real(kind=cp),              intent (in) :: Lambda  ! Wavelength
        real(kind=cp),              intent (in) :: ssnn    !(SinTheta/Lambda)^2
        real(kind=cp), dimension(3),intent (in) :: hkl     ! Components of the scattering vectors
        real(kind=cp),              intent (in) :: f2      ! Square of the structure factor
        real(kind=cp), dimension(6),intent (in) :: extc    ! Extinction coefficients
        real(kind=cp),              intent(out) :: ys      ! Extinction correction factor: Icorr = I.ys
        real(kind=cp), dimension(6),optional,intent(out) :: der     ! Derivatives of ys w.r.t. extinction coefficients
       !---- Local variables ----!
        real(kind=cp):: a,b,c, g, r, sin2t, coefa, qq
        real(kind=cp), dimension(6) :: hkl_q

        coefa=0.001
        if(job == 0 .or. job == 2 ) coefa=0.000001
        if(present(der)) der=0.0
        ys=1.0
        a = lambda*lambda*lambda
        qq=SQRT(ssnn)*lambda
        if(abs(qq) > 1.0) then
           err_extinc=.true.
           write(err_extinc_mess,fmt="(a,3f8.2,a)") " ----> WARNING!: Extinction correction fixed to 1.0 for reflection: ",&
               hkl, "  Check cell parameters !!!"
           return
        else
           sin2t=SIN(2.0*ASIN(qq))
        end if

        Select Case(iext)
          Case(1)
               b = 0.001*f2*a/sin2t
               c = 1.0 + b * extc(1)
               ys=1.0/SQRT(c)
               if(present(der)) der(1)=-0.5*c**(-1.5)*b

           Case (2,3,4)   !Anisotropic Shelx-like extinction correction
               b = coefa*f2*a/sin2t
               b = b * 0.25 /ssnn
               hkl_q=(/ hkl(1)*hkl(1),hkl(2)*hkl(2),hkl(3)*hkl(3),hkl(1)*hkl(2), &
                      hkl(1)*hkl(3),hkl(2)*hkl(3) /)
               r=dot_product(extc,hkl_q)
               c = 1.0 + b * r
               ys=1.0/SQRT(c)
               if(present(der)) then
                 g=-0.5*c**(-1.5)*b !component of derivative
                 der= g * hkl_q
               end if
        End Select
        return
      End Subroutine SHELX_Extinction

  End Module CFML_Extinction_correction