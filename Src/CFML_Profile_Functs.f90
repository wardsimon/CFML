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
!!---- MODULE:  CFML_PowderProfiles_CW
!!----   INFO:  Calculation of peak profile functions
!!----
!!---- HISTORY
!!----    Update: 03/03/2011
!!----
!!---- DEPENDENCIES
!!----    CFML_GlobalDeps, only: cp
!!----
!!---- VARIABLES
!!----
!!---- PROCEDURES
!!----
!!----    Functions:
!!----      BACK_TO_BACK_EXP
!!----      EXPONENTIAL
!!----      GAUSSIAN
!!----      HAT
!!----      IKEDA_CARPENTER
!!----      LORENTZIAN
!!----      PSEUDOVOIGT
!!----      SPLIT_PSEUDOVOIGT
!!----      TCH_PVOIGT
!!----
!!----    Subroutines:
!!----      BACK_TO_BACK_EXP_DER
!!----      EXPONENTIAL_DER
!!----      GAUSSIAN_DER
!!----      HAT_DER
!!----      IKEDA_CARPENTER_DER
!!----      LORENTZIAN_DER
!!----      PSEUDOVOIGT_DER
!!----      SPLIT_PSEUDOVOIGT_DER
!!----      TCH_PVOIGT_DER
!!----
!!
 Module CFML_PowderProfiles_CW
    !---- Use Files ----!
    Use CFML_GlobalDeps, only: cp

    !---- Variables ----!
    implicit none

    !---- List of Public Functions ----!
    public :: Pseudovoigt,Lorentzian,Gaussian,Back_To_Back_Exp, Ikeda_Carpenter, Exponential, &
              Hat, Split_Pseudovoigt, TCH_pVoigt

    !---- List of Public Subroutines ----!
    public :: Pseudovoigt_Der,Lorentzian_Der,Gaussian_Der,Back_To_Back_Exp_Der, Ikeda_Carpenter_Der, &
              Exponential_Der, Hat_Der, Split_Pseudovoigt_Der, TCH_pVoigt_Der


 Contains
    !!----
    !!---- Pure Function Back_To_Back_Exp(X,Par)
    !!----    real(kind=cp),              intent(in) :: x
    !!----    real(kind=cp), dimension(:),intent(in) :: Par
    !!----
    !!----
    !!---- Update: October - 2005
    !!
    Pure Function Back_To_Back_Exp(X,Par) Result (Bb_Val)
       !---- Arguments ----!
       real(kind=cp),              intent(in) :: x
       real(kind=cp), dimension(:),intent(in) :: par
       real(kind=cp)                          :: bb_val

       !--- Local variables ---!
       real(kind=cp) :: alfa,beta,N

       alfa=par(1)
       beta=par(2)
       N= 0.5*alfa*beta/(alfa+beta)
       if ( x < 0.0) then
          bb_val =  N*exp(alfa*x)
       else
          bb_val =  N*exp(-beta*x)
       end if

       return
    End Function Back_To_Back_Exp

    !!----
    !!---- Pure Function Exponential(X,Par)
    !!----    real(kind=cp),              intent(in) :: x
    !!----    real(kind=cp), dimension(:),intent(in) :: Par
    !!----
    !!----
    !!---- Update: October - 2005
    !!
    Pure Function Exponential(X,Par) Result (Ex_Val)
       !---- Arguments ----!
       real(kind=cp),              intent(in) :: x
       real(kind=cp), dimension(:),intent(in) :: par
       real(kind=cp)                          :: ex_val

       !--- Local variables ---!
       real(kind=cp) :: alfa

       if (x < 0.0) then
          ex_val=0.0
       else
          alfa=par(1)
          ex_val = alfa*exp(-alfa*x)
       end if

       return
    End Function Exponential

    !!----
    !!---- Pure Function Gaussian(X,Par)
    !!----    real(kind=cp),              intent(in) :: x
    !!----    real(kind=cp), dimension(:),intent(in) :: Par
    !!----
    !!----
    !!---- Update: October - 2005
    !!
    Pure Function Gaussian(X,Par) Result (Gauss_Val)
       !---- Arguments ----!
       real(kind=cp),              intent(in) :: x
       real(kind=cp), dimension(:),intent(in) :: par
       real(kind=cp)                          :: gauss_val

       !--- Local variables ---!
       real(kind=cp) :: H,ag,bg

       H=par(1)
       ag= 0.93943727869965133377234032841018/H
       bg= 2.7725887222397812376689284858327/(H*H)
       gauss_val = ag* exp(-bg*x*x)

       return
    End Function Gaussian

    !!----
    !!---- Pure Function Hat(X,Par)
    !!----    real(kind=cp),              intent(in) :: x
    !!----    real(kind=cp), dimension(:),intent(in) :: Par
    !!----
    !!----
    !!---- Update: October - 2005
    !!
    Pure Function Hat(X,Par) Result (H_Val)
       !---- Arguments ----!
       real(kind=cp),              intent(in) :: x
       real(kind=cp), dimension(:),intent(in) :: par
       real(kind=cp)                          :: h_val

       !--- Local variables ---!
       real :: hw

       hw=par(1)*0.5

       if (x < -hw .or. x > hw) then
          h_val=0.0
       else
          h_val =1.0/par(1)
       end if

       return
    End Function Hat

    !!----
    !!---- Pure Function Ikeda_Carpenter(X,Par)
    !!----    real(kind=cp),              intent(in) :: x
    !!----    real(kind=cp), dimension(:),intent(in) :: Par
    !!----
    !!----
    !!---- Update: October - 2005
    !!
    Pure Function Ikeda_Carpenter(X,Par) Result (Ik_Val)
       !---- Arguments ---!
       real(kind=cp),              intent(in) :: x
       real(kind=cp), dimension(:),intent(in) :: par
       real(kind=cp)                          :: ik_val

       !--- Local variables ---!
       real :: alfa,beta,R,amb,a3,x2,exb,exa,ramb,poly

       if (x < 0.0) then
          ik_val=0.0
       else
          x2=x*x
          alfa=par(1)
          beta=par(2)
          R   =par(3)
          amb=alfa-beta
          ramb=1.0/(amb*amb*amb)
          a3=alfa*alfa*alfa
          exb=exp(-beta*x)
          exa=exp(-alfa*x)
          poly=1.0+(1.0+0.5*amb*x)*amb*x
          ik_val = 0.5*a3*((1.0-R)*x2*exa+2.0*R*beta*ramb*(exb-exa*poly))
       end if

       return
    End Function Ikeda_Carpenter

    !!----
    !!---- Pure Function Lorentzian(X,Par)
    !!----    real,              intent(in) :: x
    !!----    real, dimension(:),intent(in) :: Par
    !!----
    !!----
    !!---- Update: October - 2005
    !!
    Pure Function Lorentzian(X,Par) Result (Lor_Val)
       !---- Arguments ----!
       real(kind=cp),              intent(in) :: x
       real(kind=cp), dimension(:),intent(in) :: par
       real(kind=cp)                          :: lor_val

       !--- Local variables ---!
       real(kind=cp) :: H,al,bl

       H=par(1)
       al= 0.63661977236758134307553505349006/H
       bl= 4.0/(H*H)
       lor_val = al/(1.0+bl*x*x)

       return
    End Function Lorentzian

    !!----
    !!---- Pure Function Pseudovoigt(X,Par)
    !!----    real(kind=cp),              intent(in)  :: x
    !!----    real(kind=cp), dimension(:),intent(in) :: Par
    !!----
    !!----
    !!---- Update: October - 2005
    !!
    Pure Function Pseudovoigt(X,Par) Result (Pv_Val)
       !---- Arguments ----!
       real(kind=cp),              intent(in) :: x
       real(kind=cp), dimension(:),intent(in) :: par
       real(kind=cp)                          :: pv_val

       !--- Local variables ---!
       real(kind=cp) :: eta,H,x2,ag,bg,al,bl,lor,gauss

       H=par(1)
       eta=par(2)
       x2=x*x
       ag= 0.93943727869965133377234032841018/H
       bg= 2.7725887222397812376689284858327/(H*H)
       al= 0.63661977236758134307553505349006/H
       bl= 4.0/(H*H)
       gauss = ag* exp(-bg*x2)
       lor   = al/(1.0+bl*x2)
       pv_val = eta*lor + (1.0 - eta)*gauss

       return
    End Function Pseudovoigt

    !!----
    !!---- Pure Function Split_Pseudovoigt(X,Par)
    !!----    real(kind=cp),              intent(in)  :: x
    !!----    real(kind=cp), dimension(:),intent(in)  :: Par
    !!----
    !!----
    !!---- Update: October - 2005
    !!
    Pure Function Split_Pseudovoigt(X,Par) Result (Pv_Val)
       !---- Arguments ----!
       real(kind=cp),              intent(in) :: x
       real(kind=cp), dimension(:),intent(in) :: par
       real(kind=cp)                          :: pv_val

       !--- Local variables ---!
       real(kind=cp) :: eta1,eta2,eta,H1,H2,hsq,Norm,x2,bg,bl,lor,gauss

       H1=par(1)
       H2=par(2)
       eta1=par(3)
       eta2=par(4)
       Norm= 0.25*H1*(eta1*1.0126586+2.128934)+ &
             0.25*H2*(eta2*1.0126586+2.128934)
       if (x < 0.0) then
          eta=eta1
          hsq=h1*h1
       else
          eta=eta2
          hsq=h2*h2
       end if
       x2=x*x
       bg= 2.7725887222397812376689284858327/hsq
       bl= 4.0/hsq
       gauss =  exp(-bg*x2)
       lor   = 1.0/(1.0+bl*x2)
       pv_val = (eta*lor + (1.0 - eta)*gauss)/Norm

       return
    End Function Split_Pseudovoigt

    !!----
    !!---- Pure Function TCH_pVoigt(X,Par)
    !!----    real(kind=cp),              intent(in)  :: x
    !!----    real(kind=cp), dimension(:),intent(in) :: Par
    !!----
    !!----
    !!---- Update: April - 2009
    !!
    Pure Function TCH_pVoigt(X,Par) Result (Pv_Val)
       !---- Arguments ----!
       real(kind=cp),              intent(in) :: x
       real(kind=cp), dimension(:),intent(in) :: par
       real(kind=cp)                          :: pv_val

       !--- Local variables ---!
       real(kind=cp)            :: Hg,Hl,eta,H,x2,ag,bg,al,bl,lor,gauss,r
       real(kind=cp), parameter :: o1= 2.69269, o2=2.42843, o3=4.47163, o4= 0.07842
       real(kind=cp), parameter :: e1= 1.36603, e2=0.47719, e3=0.11116

       Hg=par(1)
       Hl=par(2)

       !Calculate H and eta from Tomson-Cox-Hasting formula
       H=hg**5+o1*hg**4*hl+o2*hg**3*hl**2+o3*hg**2*hl**3+o4*hg*hl**4+hl**5
       H=abs(H)**0.2_cp
       r = hl/H                       !HL/FWHM
       eta = max( 1.0e-06_cp, r*(e1 -(e2 + e3*r)*r) )  !eta
       x2=x*x
       ag= 0.93943727869965133377234032841018_cp/H
       bg= 2.7725887222397812376689284858327_cp/(H*H)
       al= 0.63661977236758134307553505349006_cp/H
       bl= 4.0_cp/(H*H)
       gauss = ag* exp(-bg*x2)
       lor   = al/(1.0_cp+bl*x2)
       pv_val = eta*lor + (1.0_cp - eta)*gauss

       return
    End Function TCH_pVoigt

    !!----
    !!---- Pure Subroutine Back_To_Back_Exp_Der(X,Par,Bb_Val,Dpar)
    !!----    real(kind=cp),                        intent(in)  :: x
    !!----    real(kind=cp),           dimension(:),intent(in)  :: par
    !!----    real(kind=cp),                        intent(out) :: bb_val
    !!----    real(kind=cp), optional, dimension(:),intent(out) :: dpar
    !!----
    !!---- Update: October - 2005
    !!
    Pure Subroutine Back_To_Back_Exp_Der(X,Par,Bb_Val,Dpar)
       !---- Arguments ----!
       real(kind=cp),                        intent(in)  :: x
       real(kind=cp),           dimension(:),intent(in)  :: par
       real(kind=cp),                        intent(out) :: bb_val
       real(kind=cp), optional, dimension(:),intent(out) :: dpar

       !--- Local variables ---!
       real(kind=cp) :: alfa,beta,N,derX,derAlf,derBet

       alfa=par(1)
       beta=par(2)
       N= 0.5*alfa*beta/(alfa+beta)   !dN/da= 2(N/bet)**2 , dN/db= 2(N/alf)**2

       if ( x < 0.0) then
          bb_val =  N*exp(alfa*x)
          if (present(dpar)) then
             derX=bb_val*alfa
             derAlf=bb_val*(2.0*N/beta/beta+x)    !df/da= dN/da*exp(a*x)+Nexp(a*x)*(x)= (2N/bet**2+x)* f
             derBet=bb_val* 2.0*N/alfa/alfa       !df/db= dN/db * exp(-a*x) = 2N/alf**2 *f
             dpar(1:3) = (/derX,derAlf,derBet/)
          end if
       else
          bb_val =  N*exp(-beta*x)
          if (present(dpar)) then
             derX=bb_val*beta
             derAlf=bb_val*2.0*N/beta/beta      !df/da= dN/da*exp(-b*x)= 2N/bet**2 * f
             derBet=bb_val*(2.0*N/alfa/alfa-x)  !df/db= dN/db * exp(-b*x)+Nexp(-b*x)*(-x) = (2N/alf**2 -x) *f
             dpar(1:3) = (/derX,derAlf,derBet/)
          end if
       end if

       return
    End Subroutine Back_To_Back_Exp_Der

    !!----
    !!---- Pure Subroutine Exponential_Der(X,Par,Ex_Val,Dpar)
    !!----    real(kind=cp),                        intent(in)  :: x
    !!----    real(kind=cp),           dimension(:),intent(in)  :: par
    !!----    real(kind=cp),                        intent(out) :: ex_val
    !!----    real(kind=cp), optional, dimension(:),intent(out) :: dpar
    !!----
    !!---- Update: October - 2005
    !!
    Pure Subroutine Exponential_Der(X,Par,Ex_Val,Dpar)
       !---- Arguments ----!
       real(kind=cp),                        intent(in) :: x
       real(kind=cp),           dimension(:),intent(in) :: par
       real(kind=cp),                        intent(out):: ex_val
       real(kind=cp), optional, dimension(:),intent(out):: dpar

       !--- Local variables ---!
       real(kind=cp) :: alfa

       if (x < 0.0) then
          ex_val=0.0
          if (present(dpar)) then
             dpar(1:2) = (/0.0,0.0/)
          end if
       else
          alfa=par(1)
          ex_val = alfa*exp(-alfa*x)
          if (present(dpar)) then
             dpar(1:2) = ex_val*(/-alfa, 1.0/alfa - x/)   !derX, derAlfa
          end if
       end if

       return
    End Subroutine Exponential_Der

    !!----
    !!---- Pure Subroutine Gaussian_Der(X,Par,Gauss_Val,Dpar)
    !!----    real(kind=cp),                        intent(in)  :: x
    !!----    real(kind=cp),           dimension(:),intent(in)  :: par
    !!----    real(kind=cp),                        intent(out) :: Gauss_val
    !!----    real(kind=cp), optional, dimension(:),intent(out) :: dpar
    !!----
    !!---- Update: October - 2005
    !!
    Pure Subroutine Gaussian_Der(X,Par,Gauss_Val,Dpar)
       !---- Arguments ----!
       real(kind=cp),                       intent(in) :: x
       real(kind=cp),          dimension(:),intent(in) :: par
       real(kind=cp),                       intent(out):: gauss_val
       real(kind=cp), optional,dimension(:),intent(out):: dpar

       !--- Local variables ---!
       real(kind=cp) :: H,ag,bg, gaussp,dgaussH,x2

       H=par(1)
       ag= 0.93943727869965133377234032841018/H
       bg= 2.7725887222397812376689284858327/(H*H)
       x2=x*x
       gauss_val = ag* exp(-bg*x2)
       if (present(dpar)) then
          gaussp  = -2.0 * gauss_val * bg * x
          dgaussH = (2.0*bg*x2-1.0)*gauss_val/H
          dpar(1:2)=(/gaussp,dgaussH/)
       end if

       return
    End Subroutine Gaussian_Der

    !!----
    !!---- Pure Subroutine Hat_Der(X,Par,H_Val,Dpar)
    !!----    real(kind=cp),                        intent(in)  :: x
    !!----    real(kind=cp),           dimension(:),intent(in)  :: par
    !!----    real(kind=cp),                        intent(out) :: H_val
    !!----    real(kind=cp), optional, dimension(:),intent(out) :: dpar
    !!----
    !!---- Update: October - 2005
    !!
    Pure Subroutine Hat_Der(X,Par,H_Val,Dpar)
       !---- Arguments ----!
       real(kind=cp),                        intent(in)  :: x
       real(kind=cp),           dimension(:),intent(in)  :: par
       real(kind=cp),                        intent(out) :: h_val
       real(kind=cp), optional, dimension(:),intent(out) :: dpar

       !--- Local variables ---!
       real(kind=cp) :: hw

       hw=par(1)*0.5

       if (x < -hw .or. x > hw) then
          h_val=0.0_cp
          if (present(dpar)) then
             dpar(1:2) = (/0.0_cp,0.0_cp/)
          end if
       else
          h_val =1.0/par(1)
          if (present(dpar)) then
             dpar(1:2) = (/0.0_cp,-1.0_cp/(par(1)*par(1))/)
          end if
       end if

       return
    End Subroutine Hat_Der

    !!----
    !!---- Pure Subroutine Ikeda_Carpenter_Der(X,Par,Ik_Val,Dpar)
    !!----    real(kind=cp),                        intent(in)  :: x
    !!----    real(kind=cp),           dimension(:),intent(in)  :: par
    !!----    real(kind=cp),                        intent(out) :: Ik_val
    !!----    real(kind=cp), optional, dimension(:),intent(out) :: dpar
    !!----
    !!---- Update: October - 2005
    !!
    Pure Subroutine Ikeda_Carpenter_Der(X,Par,Ik_Val,Dpar)
       !---- Arguments ----!
       real(kind=cp),                        intent(in)  :: x
       real(kind=cp),           dimension(:),intent(in)  :: par
       real(kind=cp),                        intent(out) :: ik_val
       real(kind=cp), optional, dimension(:),intent(out) :: dpar

       !--- Local variables ---!
       real(kind=cp) :: alfa,beta,R,amb,a3,x2,exb,exa,ramb,poly
       real(kind=cp) :: derX, derAlf, derBet, derR

       if (x < 0.0) then
          ik_val=0.0
          if (present(dpar)) then
             derX=0.0
             derAlf=0.0
             derBet=0.0
             derR=0.0
             dpar(1:4) = (/derX,derAlf,derBet,derR/)
          end if
       else
          x2=x*x
          alfa=par(1)
          beta=par(2)
          R   =par(3)
          amb=alfa-beta         !damb/da=1; damb/db=-1; d(b*amb^(-3))/db=amb^(-3)-3b amb^-3/amb= ramb*(1+3b/amb)
          ramb=1.0/(amb*amb*amb)!dramb/da=-3 amb^(-4);  dramb/db= 3 amb^(-4)
          a3=alfa*alfa*alfa
          exb=exp(-beta*x)
          exa=exp(-alfa*x)                 !dexa/da= exa * (-x)
          poly=1.0+(1.0+0.5*amb*x)*amb*x   !dpoly/dx= (1+amb*x)*amb ; dpoly/da= x+0.5*x*x*2*amb = (1+x*amb)*x
          ik_val = 0.5*a3*((1.0-R)*x2*exa+2.0*beta*ramb*(exb-exa*poly)*R)
          if (present(dpar)) then
             derX=0.5*a3*((1.0-R)*(2.0*x*exa-x2*exa*alfa)+2.0*beta*ramb*R*(-exb*beta+exa*alfa*poly-exa*(1.0+amb*x)*amb))
             derAlf=3.0*poly/a3 + 0.5*a3*(-(1.0-R)*x2*exa*x + 2.0*beta*R*ramb*( -3.0/amb *(exb-exa*poly) + &
                    exa*( x*poly-(1.0+x*amb)*x )))
             derBet=a3*R * ( ramb*(1.0+3.0*beta/amb)*(exb-exa*poly) + beta*ramb*(-exb*x + exa*(1.0+x*amb)*x))
             derR=0.5*a3*(-x2*exa+2.0*beta*ramb*(exb-exa*poly))
             dpar(1:4) = (/derX,derAlf,derBet,derR/)
          end if
       end if

       return
    End Subroutine Ikeda_Carpenter_Der

    !!----
    !!---- Pure Subroutine Lorentzian_Der(X,Par,Lor_Val,Dpar)
    !!----    real(kind=cp),                        intent(in)  :: x
    !!----    real(kind=cp),           dimension(:),intent(in)  :: par
    !!----    real(kind=cp),                        intent(out) :: Lor_val
    !!----    real(kind=cp), optional, dimension(:),intent(out) :: dpar
    !!----
    !!---- Update: October - 2005
    !!
    Pure Subroutine Lorentzian_Der(X,Par,Lor_Val,Dpar)
       !---- Arguments ----!
       real(kind=cp),                        intent(in) :: x
       real(kind=cp),           dimension(:),intent(in) :: par
       real(kind=cp),                        intent(out):: lor_val
       real(kind=cp), optional, dimension(:),intent(out):: dpar

       !--- Local variables ---!
       real(kind=cp) :: H,al,bl,lorp,dlorh,x2

       H=par(1)
       al= 0.63661977236758134307553505349006_cp/H
       bl= 4.0_cp/(H*H)
       x2=x*x
       lor_val = al/(1.0_cp+bl*x2)
       if (present(dpar)) then
          lorp  = -2.0_cp *lor_val*lor_val*bl*x/al
          dlorH = (2.0_cp*bl*lor_val*x2/al -1.0_cp)*lor_val/H
          dpar(1:2)=(/lorp,dlorH/)
       end if

       return
    End Subroutine Lorentzian_Der

    !!----
    !!---- Pure Subroutine Pseudovoigt_Der(X,Par,Pv_Val,Dpar)
    !!----    real(kind=cp),                        intent(in)  :: x
    !!----    real(kind=cp),           dimension(:),intent(in)  :: par
    !!----    real(kind=cp),                        intent(out) :: Pv_val
    !!----    real(kind=cp), optional, dimension(:),intent(out) :: dpar
    !!----
    !!----  Pv_Val is the value of the function for the argument x
    !!----  Par=(/H,eta/)  and Dpar(1:3)=(/derx,derH,derEta/)
    !!----
    !!---- Update: October - 2005
    !!
    Pure Subroutine Pseudovoigt_Der(X,Par,Pv_Val,Dpar)
       !---- Arguments ----!
       real(kind=cp),                       intent(in) :: x
       real(kind=cp), dimension(:),         intent(in) :: par
       real(kind=cp),                       intent(out):: pv_val
       real(kind=cp), optional,dimension(:),intent(out):: dpar

       !--- Local variables ---!
       real(kind=cp) :: eta,H,x2,ag,bg,al,bl,lor,gauss, invH,invH2
       real(kind=cp) :: derH,derEta,derx,dlorH,dgaussH,lorp,gaussp

       H=par(1)
       eta=par(2)
       x2=x*x
       invH=1.0_cp/H
       invH2=invH*invH
       ag= 0.93943727869965133377234032841018_cp*invH
       bg= 2.7725887222397812376689284858327_cp*invH2
       al= 0.63661977236758134307553505349006_cp*invH
       bl= 4.0_cp*invH2
       gauss = ag* exp(-bg*x2)
       lor   = al/(1.0_cp+bl*x2)
       pv_val = eta*lor + (1.0_cp - eta)*gauss

       if (present(dpar)) then
          derEta= lor-gauss  !Eta

          lorp = -2.0_cp *lor*lor*bl*x/al  !x
          gaussp = -2.0_cp * gauss * bg * x  !x
          derx=eta*lorp+(1.0_cp-eta)*gaussp  !x

          !dalH= -al*invH
          !dblH= -2.0*bl*invH
          !dlorH= lor/al *dalH - lor*lor/al *x2 * dblH
          dlorH= (2.0_cp*bl*lor*x2/al -1.0_cp)*invH*lor

          !dagH=-ag*invH
          !dbgH=-2.0*bg*invH
          !dgaussH= dagH*gauss/ag - gauss * x2*dbgH = -invH*gauss+2.0*bg*invH*gauss*x2
          dgaussH= (2.0_cp*bg*x2-1.0_cp)*invH*gauss

          derH=eta*dlorH + (1.0_cp-eta) * dgaussH
          dpar(1:3)=(/derx,derH,derEta/)
       end if

       return
    End Subroutine Pseudovoigt_Der

    !!----
    !!---- Pure Subroutine Split_Pseudovoigt_Der(X,Par,Pv_Val,Dpar)
    !!----    real(kind=cp),                        intent(in)  :: x
    !!----    real(kind=cp),           dimension(:),intent(in)  :: par
    !!----    real(kind=cp),                        intent(out) :: Pv_val
    !!----    real(kind=cp), optional, dimension(:),intent(out) :: dpar
    !!----
    !!---- Update: October - 2005
    !!
    Pure Subroutine Split_Pseudovoigt_Der(X,Par,Pv_Val,Dpar)
       !---- Arguments ----!
       real(kind=cp),                       intent(in) :: x
       real(kind=cp), dimension(:),         intent(in) :: par
       real(kind=cp),                       intent(out):: pv_val
       real(kind=cp), optional,dimension(:),intent(out):: dpar

       !--- Local variables ---!
       real(kind=cp) :: eta1,eta2,eta,H1,H2,hsq,Norm,x2,bg,bl,lor,gauss,  &
                        derH1,derH2,derEta1,derEta2,derx,dlorH,dgaussH,   &
                        lorp,gaussp, invH, Numer

       H1=par(1)
       H2=par(2)
       eta1=par(3)
       eta2=par(4)
       Norm= 0.25_cp*H1*(eta1*1.0126586_cp+2.128934_cp)+ &
             0.25_cp*H2*(eta2*1.0126586_cp+2.128934_cp)
       if (x < 0.0) then
          eta=eta1
          hsq=h1*h1
          invH=1.0_cp/h1
       else
          eta=eta2
          hsq=h2*h2
          invH=1.0_cp/h2
       end if
       x2=x*x
       bg= 2.7725887222397812376689284858327_cp/hsq
       bl= 4.0_cp/hsq
       gauss =  exp(-bg*x2)
       lor   = 1.0_cp/(1.0_cp+bl*x2)
       pv_val = (eta*lor + (1.0_cp - eta)*gauss)/Norm

       if (present(dpar)) then

            lorp = -2.0_cp *lor*lor*bl*x     !x
          gaussp = -2.0_cp * gauss * bg * x  !x
          derx=(eta*lorp+(1.0_cp-eta)*gaussp)/Norm  !x
          !
          !dblH= -2.0*bl*invH
          !dlorH=  - lor*lor *x2 * dblH = lor*lor *x2 * 2.0*bl*invH
          dlorH= 2.0_cp*bl*x2*invH*lor*lor

          !
          !dbgH=-2.0*bg*invH
          !dgaussH= - gauss * x2*dbgH = 2.0*bg*invH*gauss*x2
          dgaussH= 2.0_cp*bg*x2*invH*gauss

          Numer  = eta*dlorH + (1.0_cp-eta) * dgaussH

          if (x < 0.0_cp) then

             ! dNormEi= 0.25*Hi*1.0126586 = Hi*0.25316465
             derEta1= (lor-gauss- pv_val*H1*0.25316465_cp)/Norm   !Eta
             derEta2= -pv_val*h2*0.25316465_cp/Norm

             ! pv_val = (eta*lor + (1.0 - eta)*gauss)/Norm
             ! derH1 =  (eta*dLorH1+(1-eta)*dgaussH1)*Norm - pv_val*Norm*dNormH1/Norm**2
             ! derH1 =  (Numer - pv_val*dNormH1)/Norm
             ! dNormH= 0.25*(eta*1.0126586+2.128934)
             ! derH2 = pv_val*Norm (-1/Norm**2)  * dNormH2 = - pv_val*0.25*(eta2*1.0126586+2.128934)/Norm

             derH1  = (numer - pv_val* 0.25_cp*(eta*1.0126586_cp+2.128934_cp))/Norm
             derH2  = - pv_val*0.25_cp*(eta2*1.0126586_cp+2.128934_cp)/Norm

          else

             ! dNormEi= 0.25*Hi*1.0126586 = Hi*0.25316465
             derEta2= (lor-gauss- pv_val*H2*0.25316465_cp)/Norm   !Eta
             derEta1= -pv_val*h1*0.25316465_cp/Norm

             ! pv_val = (eta*lor + (1.0 - eta)*gauss)/Norm
             ! derH2 =  (eta*dLorH2+(1-eta)*dgaussH2)*Norm - pv_val*Norm*dNormH2/Norm**2
             ! derH2 =  (Numer - pv_val*dNormH1)/Norm
             ! dNormH= 0.25*(eta*1.0126586+2.128934)
             ! derH1 = pv_val*Norm (-1/Norm**2)  * dNormH1 = - pv_val*0.25*(eta1*1.0126586+2.128934)/Norm

             derH2  = (numer - pv_val* 0.25_cp*(eta*1.0126586_cp+2.128934_cp))/Norm
             derH1  = - pv_val*0.25_cp*(eta1*1.0126586_cp+2.128934_cp)/Norm

          end if

          dpar(1:5)=(/derx,derH1,derH2,derEta1,derEta2/)

       end if

       return
    End Subroutine Split_Pseudovoigt_Der
    !!----
    !!---- Pure Subroutine TCH_pVoigt_Der(X,Par,Pv_Val,dPar)
    !!----    real(kind=cp),              intent(in)  :: x
    !!----    real(kind=cp), dimension(:),intent(in) :: Par
    !!----
    !!----  Pv_Val is the value of the function for the argument x
    !!----  Par=(/HG,HL/)  and Dpar(1:3)=(/derx,derHG,derHL/)
    !!----
    !!---- Update: April - 2009
    !!
    Pure Subroutine TCH_pVoigt_Der(X,Par,Pv_Val,dPar)
       !---- Arguments ----!
       real(kind=cp),                       intent(in) :: x
       real(kind=cp), dimension(:),         intent(in) :: par
       real(kind=cp),                       intent(out):: pv_val
       real(kind=cp), optional,dimension(:),intent(out):: dpar

       !--- Local variables ---!
       real(kind=cp), parameter :: o1= 2.69269_cp, o2=2.42843_cp, o3=4.47163_cp, o4= 0.07842_cp
       real(kind=cp), parameter :: e1= 1.36603_cp, e2=0.47719_cp, e3=0.11116_cp
       real(kind=cp) :: Hg,Hl,eta,H,x2,ag,bg,al,bl,lor,gauss, invH,invH2,r
       real(kind=cp) :: derH,derEta,derx,dlorH,dgaussH,lorp,gaussp, &
                        dhdhg,dhdhl,deta,detag,detal,derHg,derHl

       Hg=par(1)
       Hl=par(2)
       !Calculate H and eta from Tomson-Cox-Hasting formula
       H=hg**5+o1*hg**4*hl+o2*hg**3*hl**2+o3*hg**2*hl**3+o4*hg*hl**4+hl**5
       H=abs(H)**0.2_cp
       r = hl/H                       !HL/FWHM
       eta = max( 1.0e-06_cp, r*(e1 -(e2 + e3*r)*r) )  !eta
       x2=x*x
       invH=1.0_cp/H
       invH2=invH*invH
       ag= 0.93943727869965133377234032841018_cp*invH
       bg= 2.7725887222397812376689284858327_cp*invH2
       al= 0.63661977236758134307553505349006_cp*invH
       bl= 4.0_cp*invH2
       gauss = ag* exp(-bg*x2)
       lor   = al/(1.0_cp+bl*x2)
       pv_val = eta*lor + (1.0_cp - eta)*gauss

       if(present(dpar)) then
          dhdhg = 0.2_cp/H**4*(5.0_cp*hg**4+ 4.0_cp*o1* hg**3*hl+  &
                  3.0_cp*o2*hg*hg*hl*hl + 2.0_cp*o3*hg*hl**3 + o4*hl**4)
          dhdhl = 0.2_cp/H**4*(o1*hg**4+ 2.0_cp*o2*hg**3*hl+  &
                  3.0_cp*o3*hg*hg*hl*hl + 4.0_cp*o4*hg*hl**3 + 5.0_cp*hl**4)
           deta = e1-(2.0_cp*e2 - 3.0_cp*e3*r)*r  !derivative of ETA w.r.t. r
          detag = -r*deta*dhdhg*invH
          detal = (1.0_cp-r*dhdhl)*deta*invH

          derEta= lor-gauss  !Eta
          lorp = -2.0_cp *lor*lor*bl*x/al  !x
          gaussp = -2.0_cp * gauss * bg * x  !x
          derx=eta*lorp+(1.0-eta)*gaussp  !x

          dlorH= (2.0_cp*bl*lor*x2/al -1.0_cp)*invH*lor
          dgaussH= (2.0_cp*bg*x2-1.0_cp)*invH*gauss
          derH=eta*dlorH + (1.0_cp-eta) * dgaussH

          derHG= derH * dhdhg + derEta * detag  !Chain rule
          derHL= derH * dhdhl + derEta * detal
          dpar(1:3)=(/derx,derHG,derHL/)
       end if

       return
    End Subroutine TCH_pVoigt_Der

 End Module CFML_PowderProfiles_CW
