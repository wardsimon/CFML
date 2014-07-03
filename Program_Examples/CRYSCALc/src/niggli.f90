!     Last change:  TR    5 Dec 2006    6:11 pm
subroutine Niggli_CELL_TR()
! algorithm:
! Krivy and Gruber, Acta Cryst 1976, A32, 297

 USE cryscalc_module,   ONLY : unit_cell, message_text, debug_proc
 USE macros_module,     ONLY : signe
 USE IO_module,         ONLY : write_info
 USE CFML_Math_General, ONLY : cosd, acosd
  implicit none
   REAL, DIMENSION(6)            :: A
   REAL                          :: CosAlfa, CosBeta, CosGamma
   INTEGER                       :: i
   REAL                          :: tmp_value
   REAL, parameter               :: eps = 0.1

   if(debug_proc%level_2)  call write_debug_proc_level(2, "Niggli_cell_TR")
   
  a(1) = unit_cell%param(1) **2.
  a(2) = unit_cell%param(2) **2.
  a(3) = unit_cell%param(3) **2.

  CosAlfa  = COSd(unit_cell%param(4))
  Cosbeta  = COSd(unit_cell%param(5))
  Cosgamma = COSd(unit_cell%param(6))

  a(4) = 2*unit_cell%param(2) * unit_cell%param(3) * CosAlfa
  a(5) = 2*unit_cell%param(1) * unit_cell%param(3) * CosBeta
  a(6) = 2*unit_cell%param(1) * unit_cell%param(2) * Cosgamma


   ! N1

   do

   IF( (a(1) > a(2))  .or.    &
       (abs(a(1) - a(2)) < eps)     .and.  (ABS(a(4)) > ABS(a(5)))    ) then
     call change(a(1), a(2))
     call change(a(4), a(5))
   END if

  ! N2
    IF( (a(2) > a(3))  .or.    &
       (abs(a(2) - a(3)) < eps)     .and.  (ABS(a(5)) > ABS(a(6)))    ) then
     call change(a(2), a(3))
     call change(a(5), a(6))
     cycle
   END if

  ! N3
  IF( a(4) * a(5) * a(6) > 0.) then
   a(4) = ABS(a(4))
   a(5) = ABS(a(5))
   a(6) = ABS(a(6))
  ! N4
  else
   a(4) = -ABS(a(4))
   a(5) = -ABS(a(5))
   a(6) = -ABS(a(6))
  endif

  ! N5
  IF( (ABS(a(4)) > a(2)+eps)                               .OR.   &
      ((ABS(a(4)-a(2)) < eps) .AND. (2*a(5) < a(6)))       .OR.   &
      ((ABS(a(4)+a(2)) < eps) .AND. (a(6)<0.))           ) then
    a(3) = a(2) + a(3) - a(4)*signe(a(4))
    a(5) = a(5) - a(6)*signe(a(4))
    a(4) = a(4) - 2*a(2)*signe(a(4))
    cycle
  END if

  ! N6
  IF( (ABS(a(5)) > a(1)+eps)                          .OR.   &
      ((ABS(a(5)-a(1)) < eps) .AND. (2*a(4) < a(6)))  .OR.   &
      ((ABS(a(5)+a(1)) < eps) .AND. (a(6)<0.))           ) then
      a(3) = a(1) + a(3) - a(5)*signe(a(5))
      a(4) = a(4) - a(6)*signe(a(5))
      a(5) = a(5) - 2*a(1)*signe(a(5))
      cycle
  END if

  ! N7
  IF( (ABS(a(6)) > a(1)+eps)                          .OR.   &
      ((ABS(a(6)-a(1)) < eps) .AND. (2*a(4) < a(5)))  .OR.   &
      ((ABS(a(6)+a(1)) < eps) .AND. (a(5)<0.))           ) then
      a(2) = a(1) + a(2) - a(6)*signe(a(6))
      a(4) = a(4) - a(5)*signe(a(6))
      a(6) = a(6) - 2*a(1)*signe(a(6))
      cycle
  END if

  ! N8
  IF( (a(1)+a(2)+a(4)+a(5)+a(6) < 0.) .or.                                     &
     ((a(1)+a(2)+a(4)+a(5)+a(6) < eps) .AND. (2*(a(1)+a(5))+a(6)>eps)   )     ) then
      a(3) = SUM(a(1:6))
      a(4) = 2*a(2) + a(4) + a(5)
      a(5) = 2*a(1) + a(5) + a(6)
      cycle
  else
   exit
  END if

  END do

  unit_cell%niggli(1) = SQRT(a(1))
  unit_cell%niggli(2) = SQRT(a(2))
  unit_cell%niggli(3) = SQRT(a(3))
  unit_cell%niggli(4) = acosd(a(4) / (2*unit_cell%niggli(2)*unit_cell%niggli(3)))
  unit_cell%niggli(5) = acosd(a(5) / (2*unit_cell%niggli(1)*unit_cell%niggli(3)))
  unit_cell%niggli(6) = acosd(a(6) / (2*unit_cell%niggli(1)*unit_cell%niggli(2)))

   call write_info('            >>> Niggli cell parameters:  ')

   WRITE(message_text,'(a,F10.5,a)') '                 .     a = ', unit_cell%niggli(1),' (A)'
    call write_info(TRIM(message_text))
   WRITE(message_text,'(a,F10.5,a)') '                 .     b = ', unit_cell%niggli(2),' (A)'
    call write_info(TRIM(message_text))
   WRITE(message_text,'(a,F10.5,a)') '                 .     c = ', unit_cell%niggli(3),' (A)'
    call write_info(TRIM(message_text))
   WRITE(message_text,'(a,F10.5,a)') '                 .  alfa = ', unit_cell%niggli(4),' (deg)'
    call write_info(TRIM(message_text))
   WRITE(message_text,'(a,F10.5,a)') '                 .  beta = ', unit_cell%niggli(5),' (deg)'
    call write_info(TRIM(message_text))
   WRITE(message_text,'(a,F10.5,a)') '                 . gamma = ', unit_cell%niggli(6),' (deg)'
    call write_info(TRIM(message_text))


   call get_cell_sym(unit_cell%Niggli, unit_cell%lattice)
    call write_info('')
    WRITE(message_text,'(3a)') '             >>> ', TRIM(unit_cell%lattice), ' lattice.'
    call write_info(TRIM(message_text))
    call write_info('')

  RETURN
end subroutine Niggli_CELL_TR

!--------------------------------------------------------
subroutine change(a, b)
 implicit none
 REAL, INTENT(INOUT)  :: a, b
 REAL                 :: tmp_value

  tmp_value = a
  a = b
  b = tmp_value

end subroutine change
!--------------------------------------------------------

subroutine get_cell_sym(a, cell_sym)
 implicit none
 REAL, DIMENSION(6), INTENT(IN)   :: a
 CHARACTER (LEN=32), INTENT(OUT)  :: cell_sym
 REAL, parameter                  :: eps     = 0.05
 REAL, parameter                  :: eps_ang = 0.1

 if (  (a(4) > 90. - eps_ang)   .and.  (a(4) < 90. + eps_ang)   .AND.   &     ! alpha = 90.
       (a(5) > 90. - eps_ang)   .and.  (a(5) < 90. + eps_ang)   .AND.   &     ! beta  = 90.
       (a(6) > 90. - eps_ang)   .and.  (a(6) < 90. + eps_ang)  ) then         ! gamma = 90.

       IF( (ABS(a(1) - a(2)) < eps)   ) THEN   ! a = b
           if ((ABS(a(2) - a(3)) < eps)) then  ! b = c
            cell_sym = 'cubic'
           else
            cell_sym = 'tetragonal'
           end if
       else
        cell_sym = 'orthorhombic'
       endif

 ELSEIF((a(4) > 90. - eps_ang)   .and.  (a(4) < 90. + eps_ang)   .AND.    &     ! alpha = 90.
         (a(5) > 90. - eps_ang)   .and.  (a(5) < 90. + eps_ang)   .AND.   &     ! beta  = 90.
         (a(6) > 120. - eps_ang)  .and.  (a(6) < 120. + eps_ang)  ) then        ! gamma = 120.
       IF( (ABS(a(1) - a(2)) < eps)   ) THEN   ! a = b
        cell_sym = 'hexagonal'
       endif

 ELSEIF( (ABS(a(4) - a(5)) < eps) .AND. (ABS(a(5) - a(6)) < eps)  ) THEN   ! alpha = beta = gamma
       IF( (ABS(a(1) - a(2)) < eps) .AND.  (ABS(a(2) - a(3)) < eps)) then  ! a = b = c
        cell_sym = 'rhomboedral'
       endif

 elseif (  (a(4) > 90. - eps_ang)   .and.  (a(4) < 90. + eps_ang)   .AND.   &     ! alpha = 90.
           (a(6) > 90. - eps_ang)   .and.  (a(6) < 90. + eps_ang)  ) then         ! gamma = 90.*
        cell_sym = 'monoclinic'

 else
        cell_sym = 'triclinic'

 endif


end subroutine   get_cell_sym


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Niggli_cell_CFML
 USE CFML_Crystal_Metrics, only : Niggli_cell
 USE cryscalc_module,      ONLY : unit_cell, message_text, debug_proc
 USE IO_module,            ONLY : write_info
 implicit none
  real, parameter       :: eps_par = 0.01
  real, parameter       :: eps_ang = 0.05
  logical               :: original_Niggli
 
 if(debug_proc%level_2)  call write_debug_proc_level(2, "Niggli_cell_CFML")
 
    !!--.. Three non coplanar vectors {a,b,c} generates a lattice using integer linear combinations
    !!--.. There are an infinite number of primitive unit cells generating the same lattice L.
    !!--.. N={a,b,c} is a Buerger cell if and only if |a|+|b|+|c| is a minimal value for all primitive
    !!--.. cells of L.
    !!--.. N is a Niggli cell of L if  (i) it is as Buerger cell of L and
    !!--..                            (ii) |90-alpha| + |90-beta| + |90-gamma| -> maximum
    !!--..                  / a.a  b.b  c.c \       /  s11  s22  s33 \
    !!--..   Niggli matrix  |               |   =   |                |
    !!--..                  \ b.c  a.c  a.b /       \  s23  s13  s12 /
    !!--..


 unit_cell%tmp = unit_cell%param
 call Niggli_cell(unit_cell%param)          ! unit_cell : 
 
 unit_cell%niggli = unit_cell%param
 unit_cell%param  = unit_cell%tmp
 
   call write_info('            >>> Niggli cell parameters:  ')

   WRITE(message_text,'(a,F10.5,a)') '                 .     a = ', unit_cell%niggli(1),' (A)'
    call write_info(TRIM(message_text))
   WRITE(message_text,'(a,F10.5,a)') '                 .     b = ', unit_cell%niggli(2),' (A)'
    call write_info(TRIM(message_text))
   WRITE(message_text,'(a,F10.5,a)') '                 .     c = ', unit_cell%niggli(3),' (A)'
    call write_info(TRIM(message_text))
   WRITE(message_text,'(a,F10.5,a)') '                 .  alfa = ', unit_cell%niggli(4),' (deg)'
    call write_info(TRIM(message_text))
   WRITE(message_text,'(a,F10.5,a)') '                 .  beta = ', unit_cell%niggli(5),' (deg)'
    call write_info(TRIM(message_text))
   WRITE(message_text,'(a,F10.5,a)') '                 . gamma = ', unit_cell%niggli(6),' (deg)'
    call write_info(TRIM(message_text))

   original_Niggli = .false.	
   if(abs(unit_cell%param(1) - unit_cell%niggli(1)) < eps_par) then
    if(abs(unit_cell%param(2) - unit_cell%niggli(2)) < eps_par) then
	 if(abs(unit_cell%param(3) - unit_cell%niggli(3)) < eps_par) then
	  if(abs(unit_cell%param(4) - unit_cell%niggli(4)) < eps_ang) then
	   if(abs(unit_cell%param(5) - unit_cell%niggli(5)) < eps_ang) then
	    if(abs(unit_cell%param(6) - unit_cell%niggli(6)) < eps_ang) then
		 original_Niggli = .true.
        end if
       end if
      end if
     end if
    end if
   end if
   
   call write_info('')
   if(original_Niggli) then
    call write_info('    >> Original cell WAS the Niggli one !')
   else
    call write_info('    >> Original cell WAS NOT the Niggli one !')
   end if 

   !call Niggli_cell_TR ! juste pour comparer les resultats
   
 return
endsubroutine Niggli_cell_CFML