!     Last change:  TR   23 Oct 2006   11:01 am

module magnetic_module

!# Extrait de :
!#        C. Kittel
!#        Physique de l'etat solide
!#        Dunod universite
!#        5me edition, 1983
!#
!# ec:     electronic configuration
!# level:  niveau
!# L:      moment magnetique orbital
!# S:      moment magnetique de spin
!# J:      moment magnetique effectif
!# g:      coef. de Lande
!# p:
!# gJ:
!#
!# ion   ec      CN      SP      CR      IR
!#

 implicit none

  private

  PUBLIC ::  set_mag_3d, set_mag_4f


  ! definitions
  TYPE, PUBLIC :: mag_4f_type
   CHARACTER (LEN=4)               :: ion
   CHARACTER (LEN=4)               :: ec        ! config.
   CHARACTER (LEN=6)               :: level     ! niveau
   CHARACTER (LEN=1)               :: L         ! L
   CHARACTER (LEN=3)               :: S         ! S
   CHARACTER (LEN=4)               :: J
   character (LEN=5)               :: g         !
   character (LEN=6)               :: p         !
   character (LEN=5)               :: gJ         !
  END TYPE mag_4f_type
  TYPE (mag_4f_type), ALLOCATABLE, DIMENSION(:), PUBLIC  :: mag_4f

  INTEGER, PARAMETER, PUBLIC :: nb_mag_4f_lines = 13

  TYPE, PUBLIC :: mag_3d_type
   CHARACTER (LEN=4)               :: ion
   CHARACTER (LEN=3)               :: ec        ! config.
   CHARACTER (LEN=5)               :: level     ! niveau
   CHARACTER (LEN=5)               :: g         !
   CHARACTER (LEN=4)               :: p_calc_g  ! p_calc = g(J(J+1)**0.5
   CHARACTER (LEN=3)               :: S
   character (LEN=4)               :: p_calc_S  ! p_calc = 2[S(S+1)]**0.5
   character (LEN=3)               :: p_exp     !
  END TYPE mag_3d_type
  TYPE (mag_3d_type), ALLOCATABLE, DIMENSION(:), PUBLIC  :: mag_3d
  INTEGER, PARAMETER, PUBLIC :: nb_mag_3d_lines = 13


contains

 subroutine set_mag_3d
  if (.not. ALLOCATED( mag_3d)) ALLOCATE( mag_3d(nb_mag_3d_lines))
   !                          ion     ec       level       g    p_calc_g    S        p_calc_S        p_exp
   mag_3d(  1) = mag_3d_type('Ti3+', '3d1', '2D3/2',  '4/5  ', '1.55', '1/2', '1.73', '1.8')
   mag_3d(  2) = mag_3d_type('V4+ ', '3d1', '2D3/2',  '4/5  ', '1.55', '1/2', '1.73', '1.8')
   mag_3d(  3) = mag_3d_type('V3+ ', '3d2', '3F2  ',  '2/3  ', '1.63', '2/2', '2.83', '2.8')
   mag_3d(  4) = mag_3d_type('Cr3+', '3d3', '4F3/2',  '2/5  ', '0.77', '3/2', '3.87', '3.8')
   mag_3d(  5) = mag_3d_type('V2+ ', '3d3', '4F3/2',  '2/5  ', '0.77', '3/2', '3.87', '3.8')
   mag_3d(  6) = mag_3d_type('Mn3+', '3d4', '5D0  ',  '     ', '0   ', '4/2', '4.90', '4.9')
   mag_3d(  7) = mag_3d_type('Cr2+', '3d4', '5D0  ',  '     ', '0   ', '4/2', '4.90', '4.9')
   mag_3d(  8) = mag_3d_type('Fe3+', '3d5', '6S5/2',  '2    ', '5.92', '5/2', '5.92', '4.9')
   mag_3d(  9) = mag_3d_type('Mn2+', '3d5', '6S5/2',  '2    ', '5.92', '5/2', '5.92', '4.9')
   mag_3d( 10) = mag_3d_type('Fe2+', '3d6', '5D4  ',  '3/2  ', '6.70', '4/2', '4.90', '5.4')
   mag_3d( 11) = mag_3d_type('Co2+', '3d7', '4F9/2',  '4/3  ', '6.63', '3/2', '3.87', '4.8')
   mag_3d( 12) = mag_3d_type('Ni2+', '3d8', '3F4  ',  '5/4  ', '5.59', '2/2', '2.83', '3.2')
   mag_3d( 13) = mag_3d_type('Cu2+', '3d9', '2D5/2',  '13/10', '3.55', '1/2', '1.73', '1.9')

 end subroutine set_mag_3d


 subroutine set_mag_4f

  if (.not. ALLOCATED( mag_4f)) ALLOCATE( mag_4f(nb_mag_4f_lines))
   !                          ion     ec       level       L    S      J          g           p        gJ

   mag_4f(  1) = mag_4f_type('Ce3+', '4f1 ',  '2F5/2 ',   '3', '1/2',  '5/2 ',   '0.857' ,  ' 2.535', ' 2.14' )
   mag_4f(  2) = mag_4f_type('Pr3+', '4f2 ',  '3H4   ',   '5', '1  ',  '4   ',   '0.8  ' ,  ' 3.577', ' 3.2 ' )
   mag_4f(  3) = mag_4f_type('Nd3+', '4f3 ',  '4I9/2 ',   '6', '3/2',  '9/2 ',   '0.727' ,  ' 3.617', ' 3.27' )
   mag_4f(  4) = mag_4f_type('Pm3+', '4f4 ',  '5I4   ',   '6', '2  ',  '4   ',   '0.6  ' ,  ' 2.68 ', ' 2.4 ' )
   mag_4f(  5) = mag_4f_type('Sm3+', '4f5 ',  '6H5/2 ',   '5', '5/2',  '5/2 ',   '0.286' ,  ' 0.845', ' 0.71' )
   mag_4f(  6) = mag_4f_type('Eu3+', '4f6 ',  '7F0   ',   '3', '3  ',  '0   ',   '     ' ,  ' 0    ', '     ' )
   mag_4f(  7) = mag_4f_type('Gd3+', '4f7 ',  '8S7.2 ',   '0', '7/2',  '7/2 ',   '2    ' ,  ' 7.937', ' 7   ' )
   mag_4f(  8) = mag_4f_type('Tb3+', '4f8 ',  '7F6   ',   '3', '3  ',  '6   ',   '1.5  ' ,  ' 9.72 ', ' 9   ' )
   mag_4f(  9) = mag_4f_type('Dy3+', '4f9 ',  '6H15/2',   '5', '5/2',  '15/2',   '1.33 ' ,  '10.64 ', '10   ' )
   mag_4f( 10) = mag_4f_type('Ho3+', '4f10',  '5I8   ',   '6', '3  ',  '8   ',   '1.25 ' ,  '10.607', '10   ' )
   mag_4f( 11) = mag_4f_type('Er3+', '4f11',  '4I15/2',   '6', '3/2',  '15/2',   '1.20 ' ,  ' 9.58 ', ' 9   ' )
   mag_4f( 12) = mag_4f_type('Tm3+', '4f12',  '3H6   ',   '5', '1  ',  '6   ',   '1.166' ,  ' 7.56 ', ' 7   ' )
   mag_4f( 13) = mag_4f_type('Yb3+', '4f13',  '2F7/2 ',   '3', '1/2',  '7/2 ',   '1.143' ,  ' 4.536', ' 4   ' )

end subroutine set_mag_4f




end module magnetic_module

