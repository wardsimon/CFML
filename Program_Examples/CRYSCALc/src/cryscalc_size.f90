
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


