!     Last change:  TR   14 Apr 2006   11:59 am
subroutine X_space_calculation(input_string)
 USE cryscal_module
 USE IO_module,      ONLY : write_info
 implicit none
  CHARACTER (LEN=*), INTENT(IN) :: input_string
  INTEGER                       :: i
  REAL                          :: X, X_dhkl, X_dstar, X_Qhkl, X_stl, X_theta
  REAL                          :: Z
  REAL                          :: eps
  eps=0.001

  select case (input_string)
      case ('STL')
       call write_info('')
       IF(keyword_wave) then
        call write_info('       STL(A-1)         Q(A-1)           d(A)     dstar(A-1)     theta(deg)   2theta(deg)')
       else
        call write_info('       STL(A-1)         Q(A-1)           d(A)     dstar(A-1)')
       endif
       call write_info('')
       do i = 1, nb_stl_value
        X = stl_value(i)
        IF(X < eps) cycle
        X_Qhkl  = 4.*pi* X
        X_dhkl  = 1./(2*X)
        X_dstar = 1./X_dhkl
        IF(keyword_WAVE) then
         Z = wavelength / (2 * X_dhkl)
         IF (Z**2 < 1.) then
          X_theta = ATAN(Z / (SQRT(1. - Z **2))) * 180. / PI
         else
          X_theta = 180.
         endif
         WRITE(message_text, '(6F15.5)') stl_value(i), X_Qhkl, X_dhkl, X_dstar, X_theta, 2.*X_theta
        else
         WRITE(message_text, '(4F15.5)') stl_value(i), X_Qhkl, X_dhkl, X_dstar
        endif
        call write_info(TRIM(message_text))
       end do
       call write_info('')

      case ('DHKL')
       IF(keyword_wave) then
        call write_info('           d(A)     dstar(A-1)       STL(A-1)         Q(A-1)     theta(deg)    2theta(deg)')
       else
        call write_info('           d(A)     dstar(A-1)       STL(A-1)         Q(A-1)')
       endif
       call write_info('')
       do i=1, nb_dhkl_value
        X = dhkl_value(i)
        IF(X < eps) cycle
        X_STL  = 1./(2.*X)
        X_Qhkl = 4.*pi* X_STL
        X_dstar = 1./X
        IF(keyword_WAVE) then
         Z = wavelength / (2 * X )
         IF (Z**2 < 1.) then
          X_theta = ATAN(Z / (SQRT(1. - Z **2))) * 180. / PI
         else
          X_theta = 180.
         endif
         WRITE(message_text, '(6F15.5)') dhkl_value(i), X_dstar, X_stl, X_Qhkl, X_theta, 2.*X_theta
        else
         WRITE(message_text, '(4F15.5)') dhkl_value(i), X_dstar, X_stl, X_Qhkl
        endif
        call write_info(TRIM(message_text))
       end do
       call write_info('')

      case ('DSTAR')
       IF(keyword_wave) then
        call write_info('        d*(A-1)           d(A)       STL(A-1)         Q(A-1)     theta(deg)    2theta(deg)')
       else
        call write_info('        d*(A-1)           d(A)       STL(A-1)         Q(A-1)')
       endif
       call write_info('')
       do i=1, nb_dstar_value
        IF(dstar_value(i) < eps) cycle
        X = 1./dstar_value(i)
        X_STL  = 1./(2.*X)
        X_Qhkl = 4.*pi* X_STL
        IF(keyword_WAVE) then
         Z = wavelength / (2 * X )
         IF (Z**2 < 1.) then
          X_theta = ATAN(Z / (SQRT(1. - Z **2))) * 180. / PI
         else
          X_theta = 180.
         endif
         WRITE(message_text, '(6F15.5)') dstar_value(i), X, X_stl, X_Qhkl, X_theta, 2.*X_theta
        else
         WRITE(message_text, '(4F15.5)') dstar_value(i), X, X_stl, X_Qhkl
        endif
        call write_info(TRIM(message_text))
       end do
       call write_info('')

      case ('QHKL')
       IF(keyword_wave) then
        call write_info('         Q(A-1)           d(A)     dstar(A-1)       STL(A-1)     theta(deg)    2theta(deg)')
       else
        call write_info('         Q(A-1)           d(A)     dstar(A-1)       STL(A-1)')
       endif
       call write_info('')
       do i=1, nb_Qhkl_value
        X = Qhkl_value(i)
        IF(X < eps) cycle
        X_STL  = X/(4.*pi)
        X_dhkl = 1./(2*X_STL)
        X_dstar = 1./X_dhkl
        IF(keyword_WAVE) then
         Z = wavelength / (2 * X_dhkl)
         IF (Z**2 < 1.) then
          X_theta = ATAN(Z / (SQRT(1. - Z **2))) * 180. / PI
         else
          X_theta = 180.
         endif
         WRITE(message_text, '(6F15.5)') Qhkl_value(i), X_dhkl, X_dstar, X_STL, X_theta, 2.*X_theta
        else
         WRITE(message_text, '(4F15.5)') Qhkl_value(i), X_dhkl, X_dstar, X_STL
        endif
        call write_info(TRIM(message_text))
       end do
       call write_info('')



      case ('THETA')
       IF(.NOT. keyword_wave) then
        call write_info('')
        call write_info('  !! WAVE keyword is a mandatory keyword !!')
        call write_info('')
        return
       endif
       call write_info('     Theta(deg)           d(A)     dstar(A-1)         Q(A-1)       STL(A-1)')
       call write_info('')

       do i=1, nb_theta_value
        X = theta_value(i)
        IF(X < eps) cycle
        X_STL  = SIN(X*pi/180)/ wavelength
        X_dhkl = 1/(2*X_STL)
        X_dstar = 1./X_dhkl
        X_Qhkl = 4.*pi*X_STL
        WRITE(message_text, '(5F15.5)') theta_value(i),  X_dhkl, X_dstar, X_Qhkl, X_stl
        call write_info(TRIM(message_text))
       end do
       call write_info('')

      case ('2THETA')
       IF(.NOT. keyword_wave) then
        call write_info('')
        call write_info('  !! WAVE keyword is a mandatory keyword !!')
        call write_info('')
        return
       endif
       call write_info('    2Theta(deg)           d(A)     dstar(A-1)         Q(A-1)       STL(A-1)')
       call write_info('')

       do i=1, nb_2theta_value
        X = two_theta_value(i)/2.
        IF(X < eps) cycle
        X_STL  = SIN(X*pi/180)/ wavelength
        X_dhkl = 1/(2*X_STL)
        X_dstar = 1./X_dhkl
        X_Qhkl = 4.*pi*X_STL
        WRITE(message_text, '(5F15.5)') TWO_theta_value(i), X_dhkl, X_dstar, X_Qhkl, X_stl
        call write_info(TRIM(message_text))
       end do
       call write_info('')


      case default
  end select

 RETURN
end subroutine X_space_calculation
