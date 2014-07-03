
 subroutine calcul_new_profile(FWHM_inst, eta_inst, FWHM, eta, X)
  use cryscalc_module,        only : wavelength, pi 
  USE Pattern_profile_module, only : particle_size
  USE CFML_Math_General,      ONLY : cosd

  implicit none
   real, intent(in)        :: FWHM_inst
   real, intent(in)        :: eta_inst
   real, intent(inout)     :: FWHM
   real, intent(inout)     :: eta   
   real, intent(in)        :: X
   real                    :: HG, HL, HL_size
 

   
   call calcul_HGHL(FWHM_inst, eta_inst, HG, HL)  ! conversion FWHM, eta ==> HG, HL
   
   HL_size = 2.*wavelength/(pi*particle_size*cosd(X/2.)) * 180/pi
   HL = HL + HL_size
 
   call calcul_H_eta(HG, HL, FWHM, eta)  ! conversion HG, HL ==> FWHM, eta
   
  return 
 end subroutine calcul_new_profile  

 

subroutine calcul_HGHL(FWHM, eta, HG, HL) 
 ! calcul de HG et HL à partir de FWHM et eta
 use cryscalc_module, only : debug_proc
 implicit none
  real, intent(in)    :: FWHM
  real, intent(in)    :: eta
  real, intent(out)   :: HG
  real, intent(out)   :: HL
  real                :: ZBRENT
  real                :: ratio, tol
  external            :: FETA, FGAU
  real                :: z
  
 !z = 1-0.74417*eta - 0.24781*eta**2. - 0.00810*eta**3.
 ! if(z > 0.) then
 !  HG = FWHM*sqrt(z)
 ! else
 !  HG = 0.
 ! end if
 ! HL = FWHM*(0.72928*eta + 0.19289*eta**2. + 0.07783*eta**3.)
  
  
  if(debug_proc%level_3) call write_debug_proc_level(3, "calculate HG and HL from FWHM and eta")
   
  tol = 1.E-06
 
  ratio = ZBRENT(FETA, 0., 1., tol)
  HL = FWHM * ratio
  HG = ZBRENT(FGAU, 0., FWHM, tol)
 
 
 return
end subroutine calcul_HGHL

subroutine calcul_H_eta(HG, HL,FWHM, eta)
 use cryscalc_module, only : debug_proc
 implicit none
  real, intent(in)   :: HG
  real, intent(in)   :: HL
  real, intent(out)  :: FWHM
  real, intent(out)  :: eta
  real               :: ratio
  
  if(debug_proc%level_3) call write_debug_proc_level(3, "calculate FWHM and eta from HG and HL")
   
  FWHM = (HG**5. +2.69269*HG**4.*HL + 2.42843*HG**3.*HG**2. + 4.47163*HG**2.*HL**3. +   &
         0.07842*HG*HL**4. + HL**5.)**0.2 
  ratio = HL/FWHM
  eta = 1.36603*ratio - 0.47719*ratio**2. + 0.11116*ratio**3.
     
 return
end subroutine calcul_H_eta


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      function FETA(X)
      USE pattern_profile_module, ONLY  : eta
      implicit none
       REAL, INTENT(IN)    :: x
       REAL                :: feta

       feta= eta - (1.36603 - 0.47719*X + 0.11116*X*X)*X
       return
      end function
!
!------------------------------------------------------------------
      function Fgau(X)
      USE pattern_profile_module, ONLY : HG, HL, H
      implicit none
       REAL, INTENT(IN)    :: x
       REAL                :: fgau

        fgau = H - (X**5. +2.69269*X**4.*HL + 2.42843*X**3.*HL**2. + 4.47163*X**2.*HL**3. +   &
                   0.07842*X*HL**4. + HL**5.)**0.2
      return
      end function

!------------------------------------------------------------------
      FUNCTION zbrent(func,x1,x2,tol)
       USE IO_module
	   USE cryscalc_module, only : lecture_OK
 
       implicit none
       REAL               :: zbrent
       REAL               :: func
       REAL, INTENT(IN)   :: x1, x2, tol
       INTEGER, PARAMETER :: ITMAX = 100
       REAL, parameter    :: eps = 3.e-08
       integer            :: iter
       REAL               :: a,b,c,d,e, fa,fb, fc, xm,p,q,r,s
       REAL               :: tol1

      a=x1
      b=x2
      fa=func(a)
      fb=func(b)
	  
      if((fa > 0. .and. fb > 0.) .or. (fa < 0. .and. fb < 0.)) then
	   !call write_info('')
       !call write_info(' > Wrong values to apply T.C.H. formulae !!')
	   !call write_info('')
       lecture_ok = .false.
       !stop     ! ' => root must be bracketed for zbrent'
      else
       lecture_ok = .true.
      endif


      c=b
      fc=fb
      do iter=1,ITMAX
        if((fb > 0. .and. fc > 0.) .or. (fb < 0. .and. fc < 0.))then
          c=a
          fc=fa
          d=b-a
          e=d
        endif
        if(abs(fc) < abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
        endif
        tol1=2.*EPS*abs(b)+0.5*tol
        xm=.5*(c-b)
        !if(abs(xm) < tol1 .or. fb == 0.)then
        if(abs(xm) < tol1 .or. ABS(fb) < eps)then
          zbrent=b
          return
        endif
        if(abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
          s=fb/fa
          !if(a ==c) then
          IF(ABS(a-c) < eps) then
            p=2.*xm*s
            q=1.-s
          else
            q=fa/fc
            r=fb/fc
            p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
            q=(q-1.)*(r-1.)*(s-1.)
          endif
          if(p > 0.) q=-q
          p=abs(p)
          if(2.*p < min(3.*xm*q-abs(tol1*q),abs(e*q))) then
            e=d
            d=p/q
          else
            d=xm
            e=d
          endif
        else
          d=xm
          e=d
        endif
        a=b
        fa=fb
        if(abs(d) > tol1) then
          b=b+d
        else
          b=b+sign(tol1,xm)
        endif
        fb=func(b)
     END do

      zbrent=b
      return
      END function


