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
!!---- MODULE: CFML_Extinction_correction
!!----   INFO: Small module for extinction corrections
!!----
!!----
!!
Module CFML_ExtinCorr
   !---- Use Modules ----!
   Use CFML_GlobalDeps, only: CP, Err_CFML, Clear_Error
   
   !---- Variables ----!
   Implicit none
   
   private
   
   !---- Public Functions ----!
   public :: AG_theta, AL_theta, BG_theta, BL_theta
    
   !---- Public Subroutines ----!
   public :: SHELX_Extinction, Becker_Coppens, Correct_FlippingRatios
   
   Interface
      Module Subroutine Becker_Coppens(iext,f2,cext,r,g,ys,dydr,dydg)
         !---- Arguments ----! 
         integer,                     intent (in) :: iext       
         real(kind=cp),               intent (in) :: f2         
         real(kind=cp), dimension(4), intent (in) :: cext       
         real(kind=cp),               intent (in) :: r,g        
         real(kind=cp),               intent(out) :: ys         
         real(kind=cp), optional,     intent(out) :: dydr       
         real(kind=cp), optional,     intent(out) :: dydg       
      End Subroutine Becker_Coppens
      
      Module Subroutine Correct_FlippingRatios(iext,Lambda,q,extc,ssnn,hkl,AN,BN,AM,BM,&
                                               yp,ym,ypm,dyp,dym,dypm,dymag)
         !---- Arguments ----!
         integer,                               intent(in)  :: iext
         real(kind=cp),                         intent(in)  :: lambda
         real(kind=cp),                         intent(in)  :: q 
         real(kind=cp), dimension(6),           intent(in)  :: extc
         real(kind=cp),                         intent(in)  :: ssnn
         real(kind=cp), dimension(3),           intent(in)  :: hkl
         real(kind=cp),                         intent(in)  :: AN
         real(kind=cp),                         intent(in)  :: BN
         real(kind=cp),                         intent(in)  :: AM
         real(kind=cp),                         intent(in)  :: BM
         real(kind=cp),                         intent(out) :: yp
         real(kind=cp),                         intent(out) :: ym
         real(kind=cp),                         intent(out) :: ypm
         real(kind=cp), dimension(:), optional, intent(out) :: dyp
         real(kind=cp), dimension(:), optional, intent(out) :: dym
         real(kind=cp), dimension(:), optional, intent(out) :: dypm
         real(kind=cp), dimension(:), optional, intent(out) :: dymag  
      End Subroutine Correct_FlippingRatios
      
      Module Subroutine SHELX_Extinction(job,iext,Lambda,ssnn,hkl,f2,extc,ys,der,derf2)
         !---- Arguments ----!
         integer,                    intent (in) :: job     
         integer,                    intent (in) :: iext    
         real(kind=cp),              intent (in) :: Lambda  
         real(kind=cp),              intent (in) :: ssnn    
         real(kind=cp), dimension(3),intent (in) :: hkl     
         real(kind=cp),              intent (in) :: f2      
         real(kind=cp), dimension(6),intent (in) :: extc    
         real(kind=cp),              intent(out) :: ys      
         real(kind=cp), dimension(6),optional,intent(out) :: der   
         real(kind=cp),              optional,intent(out) :: derf2 
      End Subroutine SHELX_Extinction                         
      
   End Interface

 Contains

    !!---- 
    !!---- AG_THETA
    !!---- 
    !!---- 23/04/2019
    !! 
    Pure Function AG_theta(cos2t) result (ag)   !A(theta) Gaussian
       !---- Arguments ----! 
       real(kind=cp), intent (in) :: cos2t
       real(kind=cp)              :: ag
       
       ag = 0.58_cp + (0.48_cp + 0.24_cp * cos2t) * cos2t
       
       return
    End Function AG_theta

    !!---- 
    !!---- AL_THETA
    !!---- 
    !!---- 23/04/2019
    !! 
    Pure Function AL_theta(cos2t) result (al)       !A(theta) Lorentzian
       !---- Arguments ----!
       real(kind=cp), intent (in) :: cos2t
       real(kind=cp)              :: al
       
       al = 0.025_cp + 0.285_cp * cos2t
       
       return
    End Function AL_theta

    !!---- 
    !!---- BG_THETA
    !!---- 
    !!---- 23/04/2019
    !! 
    Pure Function BG_theta(cos2t) result (bg)   !B(theta) Gaussian
       !---- Arguments ----!
       real(kind=cp), intent (in) :: cos2t
       real(kind=cp)              :: bg
       
       bg = 0.02_cp - 0.025_cp * cos2t
         
       return
    End Function BG_theta

    !!---- 
    !!---- BL_THETA
    !!---- 
    !! 
    Pure Function BL_theta(cos2t) result (bl)   !B(theta) Lorentzian
       !---- Arguments ----! 
       real(kind=cp), intent (in) :: cos2t
       real(kind=cp)              :: bl
       
       if (cos2t > 0.0_cp) then
          bl = 0.15_cp - 0.2_cp * (0.75_cp-cos2t)**2
       else
          bl =  -0.45_cp * cos2t
       end if
       
       return
    End Function BL_theta

 End Module CFML_ExtinCorr