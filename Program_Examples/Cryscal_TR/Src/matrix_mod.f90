!     Last change:  TR    6 Nov 2006    5:57 pm
module matrix_module

  implicit none

  PUBLIC                   :: get_matrix_complements_algebriques_33,  &
                              get_matrix_transposee_33,               &
                              matrix_determinant_33,                  &
                              MV_product,                             &
                              MATMUL_33,                              &
                              get_Mat_from_setting



  contains

!------------------------------------------------------------
! matrice des complements algebriques
 subroutine  get_matrix_complements_algebriques_33(Mi, Mo)
  implicit none
   REAL, DIMENSION(3,3), INTENT(IN)    :: Mi
   REAL, DIMENSION(3,3), INTENT(OUT)   :: Mo

   Mo(1,1) =   Mi(2,2)*Mi(3,3) - Mi(3,2)*Mi(2,3)
   Mo(1,2) = -(Mi(2,1)*Mi(3,3) - Mi(3,1)*Mi(2,3))
   Mo(1,3) =   Mi(2,1)*Mi(3,2) - Mi(3,1)*Mi(2,2)


   Mo(2,1) = -(Mi(1,2)*Mi(3,3) - Mi(3,2)*Mi(1,3))
   Mo(2,2) =   Mi(1,1)*Mi(3,3) - Mi(3,1)*Mi(1,3)
   Mo(2,3) = -(Mi(1,1)*Mi(3,2) - Mi(3,1)*Mi(1,2))

   Mo(3,1) =   Mi(1,2)*Mi(2,3) - Mi(2,2)*Mi(1,3)
   Mo(3,2) = -(Mi(1,1)*Mi(2,3) - Mi(2,1)*Mi(1,3))
   Mo(3,3) =   Mi(1,1)*Mi(2,2) - Mi(2,1)*Mi(1,2)

 end subroutine   get_matrix_complements_algebriques_33

!------------------------------------------------------------
! matrice transposee
 subroutine  get_matrix_transposee_33(Mi, Mout)
  implicit none
   REAL, DIMENSION(3,3), INTENT(IN)    :: Mi
   REAL, DIMENSION(3,3), INTENT(OUT)   :: Mout

   Mout(1,1) = Mi(1,1)
   Mout(1,2) = Mi(2,1)
   Mout(1,3) = Mi(3,1)
   Mout(2,1) = Mi(1,2)
   Mout(2,2) = Mi(2,2)
   Mout(2,3) = Mi(3,2)
   Mout(3,1) = Mi(1,3)
   Mout(3,2) = Mi(2,3)
   Mout(3,3) = Mi(3,3)

 end subroutine   get_matrix_transposee_33

!------------------------------------------------------------

   FUNCTION matrix_determinant_33(Mij)  RESULT(det)
   ! calcul du determinant d'une matrice 3*3 suivant la premiere ligne
   implicit none
    REAL, DIMENSION(3,3), INTENT(IN)          :: Mij
    REAL                                      :: det


    det = Mij(1,1) * (-1)**(1+1) * (Mij(2,2) * Mij(3,3) - Mij(3,2) * Mij(2,3))   +   &
          Mij(1,2) * (-1)**(1+2) * (Mij(2,1) * Mij(3,3) - Mij(3,1) * Mij(2,3))   +   &
          Mij(1,3) * (-1)**(1+3) * (Mij(2,1) * Mij(3,2) - Mij(3,1) * Mij(2,2))

   END  FUNCTION matrix_determinant_33

!------------------------------------------------------------

  FUNCTION MV_product(Vi,M)  RESULT (Vo)
  ! calcul produit Matrice_33 * vecteur colonne
  implicit none
   REAL, DIMENSION(3),   INTENT(IN)    :: Vi   ! vecteur initial
   REAL, DIMENSION(3,3), INTENT(IN)    :: M
   REAL, DIMENSION(3)                  :: Vo

   Vo(1) = M(1,1) * vi(1)  +  M(1,2) * vi(2) +  M(1,3) * vi(3)
   Vo(2) = M(2,1) * vi(1)  +  M(2,2) * vi(2) +  M(2,3) * vi(3)
   Vo(3) = M(3,1) * vi(1)  +  M(3,2) * vi(2) +  M(3,3) * vi(3)


 return
 end function MV_product
!------------------------------------------------------------

  FUNCTION MATMUL_33(A, B)  RESULT(C)
  ! produit de 2 matrices 3*3
  implicit none
   REAL, DIMENSION(3,3), INTENT(IN)   :: A
   REAL, DIMENSION(3,3), INTENT(IN)   :: B
   REAL, DIMENSION(3,3)               :: C

   ! Cij = SUM_1_3( A(i,k)*B(k,j))
   ! i j      i k    k j      i k    k j      i k    k j
   C(1,1) = A(1,1)*B(1,1) + A(1,2)*B(2,1) + A(1,3)*B(3,1)
   C(1,2) = A(1,1)*B(1,2) + A(1,2)*B(2,2) + A(1,3)*B(3,2)
   C(1,3) = A(1,1)*B(1,3) + A(1,2)*B(2,3) + A(1,3)*B(3,3)

   C(2,1) = A(2,1)*B(1,1) + A(2,2)*B(2,1) + A(2,3)*B(3,1)
   C(2,2) = A(2,1)*B(1,2) + A(2,2)*B(2,2) + A(2,3)*B(3,2)
   C(2,3) = A(2,1)*B(1,3) + A(2,2)*B(2,3) + A(2,3)*B(3,3)

   C(3,1) = A(3,1)*B(1,1) + A(3,2)*B(2,1) + A(3,3)*B(3,1)
   C(3,2) = A(3,1)*B(1,2) + A(3,2)*B(2,2) + A(3,3)*B(3,2)
   C(3,3) = A(3,1)*B(1,3) + A(3,2)*B(2,3) + A(3,3)*B(3,3)

   return
  END function MATMUL_33


!------------------------------------------------------------
! determination de la matrice associee a une chaine a b c
 subroutine get_Mat_from_setting(input_string, ligne, input_car, Matrice)

 implicit none
  CHARACTER(LEN=*),                   INTENT(IN)    :: input_string
  INTEGER,                            INTENT(IN)    :: ligne
  REAL,              DIMENSION(3,3),  INTENT(OUT)   :: Matrice
  CHARACTER(LEN=1),                   INTENT(IN)    :: input_car
  ! local variables
  integer                           :: col
  CHARACTER(LEN=1)                  :: string
  INTEGER                           :: k
  REAL                              :: sign_R


 k= INDEX(input_string,input_car)
 IF    (input_car =='A') then
  col = 1
 ELSEIF(input_car =='B') then
  col = 2
 ELSEIF(input_car =='C') then
  col = 3
 endif

 IF(k==0) then
  Matrice(ligne, col) = 0.
  sign_R              = 1.
 else

  IF(k==1) then
   Matrice(ligne, col ) = 1
   sign_R               = 1.
  ELSE

   READ(input_string(k-1:k-1),*) string
   IF(string == "+") then
    sign_R          = 1.
    Matrice(ligne, col) = 1.
   ELSEIF(string == '-') then
    sign_R          = -1.
    Matrice(ligne, col) = 1
   else
    READ(string, *) Matrice(ligne, col)
    READ(input_string(k-2:k-2),*) string
    IF(string == "+") then
     sign_R        = 1.
    ELSEIF(string == '-') then
     sign_R        = -1.
    endif
   endif

  endif
 endif

 Matrice(ligne, col) = sign_R * Matrice(ligne,col)
 return

 end subroutine  get_Mat_from_setting



end module matrix_module


