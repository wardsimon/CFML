!     Last change:  TR    1 Feb 2006    6:34 pm
module math_module

 public

 contains

 function transmission(mu, x) RESULT(T)
 ! calcul du coef. de transmission T = exp(-mu.x)
  REAL, INTENT(IN)              ::  mu
  REAL, INTENT(IN)              ::  x
  REAL                          ::  T

  T = EXP(-mu*x)

 end function transmission

!-----------------------------------------------------------------------
!------GET THE CRYST.-ORTHO. MATRIX E CORRESPONDING TO CELL EDGES A & B
 subroutine orthog(A, B, E)     ! extrait de CRYSTAL (JRC)
  USE CFML_Math_General, ONLY: cosd, sind
  implicit none
   REAL, DIMENSION(3),   INTENT(IN)  :: A
   REAL, DIMENSION(3),   INTENT(IN)  :: B
   REAL, DIMENSION(3,3), INTENT(OUT) :: E

   REAL :: COSGAS, SINGAS

    COSGAS = (COSD(B(1))*COSD(B(2))-COSD(B(3)))
    COSGAS = COSGAS/(SIND(B(1))*SIND(B(2)))
    SINGAS = SQRT(1.-COSGAS**2)
    E(1,1) = A(1)*SIND(B(2))*SINGAS
    E(1,2) = 0
    E(1,3) = 0
    E(2,1) = -A(1)*SIND(B(2))*COSGAS
    E(2,2) = A(2)*SIND(B(1))
    E(2,3) = 0
    E(3,1) = A(1)*COSD(B(2))
    E(3,2) = A(2)*COSD(B(1))
    E(3,3) = A(3)
  return
 end subroutine orthog


!----------------------------------------------------------------------
!------ GET THE METRIC TENSOR E CORRESPONDING TO CELL EDGES A & ANGLES B
      SUBROUTINE METRIC(A,B,E)
       USE CFML_Math_General, ONLY: cosd, sind
       implicit none
       REAL, DIMENSION(3),   INTENT(IN)    :: A
       REAL, DIMENSION(3),   INTENT(IN)    :: B
       REAL, DIMENSION(3,3), INTENT(OUT)   :: E
         integer :: i

         DO  I=1,3
          E(I,I)=A(I)*A(I)
         end do
         E(1,2)=A(1)*A(2)*COSD(B(3))
         E(2,1)=E(1,2)
         E(1,3)=A(1)*A(3)*COSD(B(2))
         E(3,1)=E(1,3)
         E(2,3)=A(2)*A(3)*COSD(B(1))
         E(3,2)=E(2,3)
      return
      end subroutine Metric


!----------------------------------------------------------------------
!------ INVERT 3X3 MATRIX A, PUT THE RESULT IN B
    SUBROUTINE MATINV(A,B,IFAIL)
     implicit NONE
      REAL, DIMENSION(3,3), INTENT(IN)  :: A
      REAL, DIMENSION(3,3), INTENT(OUT) :: B
      INTEGER,              INTENT(OUT) :: ifail
      REAL, DIMENSION(3,3)              :: E
      REAL                              :: DMAT
      INTEGER                           :: i, j

      IFAIL=0
      E(1,1) =   A(2,2)*A(3,3) - A(2,3)*A(3,2)
      E(2,1) = -(A(2,1)*A(3,3) - A(2,3)*A(3,1))
      E(3,1) =   A(2,1)*A(3,2) - A(2,2)*A(3,1)
      E(1,2) = -(A(1,2)*A(3,3) - A(1,3)*A(3,2))
      E(2,2) =   A(1,1)*A(3,3) - A(1,3)*A(3,1)
      E(3,2) = -(A(1,1)*A(3,2) - A(1,2)*A(3,1))
      E(1,3) =   A(1,2)*A(2,3) - A(1,3)*A(2,2)
      E(2,3) = -(A(1,1)*A(2,3) - A(1,3)*A(2,1))
      E(3,3) =   A(1,1)*A(2,2) - A(1,2)*A(2,1)
      DMAT = A(1,1)*E(1,1)+A(1,2)*E(2,1)+A(1,3)*E(3,1)
      IF(ABS(DMAT) < 1.E-08) THEN
        IFAIL=1
        RETURN
      ENDIF
        DO  I = 1,3
          DO J = 1,3
           B(I,J) = E(I,J)/DMAT
          END DO
        END DO
      RETURN

      END subroutine MATINV


!----------------------------------------------------------------------
!------LENGTH OF VECTOR B WHEN A IS THE CRYST. TO ORTHOGONAL MATRIX
      SUBROUTINE LENGTH(A,B,S)
      implicit none
       REAL, DIMENSION(3,3), INTENT(IN)  :: A   ! cryst. to orthogonal matrix
       REAL, DIMENSION(3),   INTENT(IN)  :: B   ! vecteur
       REAL,                 INTENT(OUT) :: S   ! longueud du vecteur
       INTEGER            :: i, j
       REAL, DIMENSION(3) :: V

        DO  I = 1,3
         V(I) = 0.
          DO J = 1,3
           V(I) = V(I)+A(I,J)*B(J)
          END DO
        END DO
        S = SQRT(V(1)**2+V(2)**2+V(3)**2)
      return
      end subroutine length


!----------------------------------------------------------------------
!------MULTIPLY MATRIX A BY VECTOR B
      SUBROUTINE MATVEC(A,B,V)
       implicit none
       REAL, DIMENSION(3,3), INTENT(IN)  :: A  ! matrix 3*3
       REAL, DIMENSION(3),   INTENT(IN)  :: B  ! vecteur B
       REAL, DIMENSION(3),   INTENT(out) :: V  ! V=A*B
       INTEGER        :: i, j


        DO  I = 1,3
         V(I) = 0.
          DO  J = 1,3
           V(I) = V(I)+A(I,J)*B(J)
          end do
        end do

        return
      end subroutine matvec

!----------------------------------------------------------------------
!------SCALAR PRODUCT OF VECTORS A AND B
      SUBROUTINE SCALPR(A,B,S)
       implicit none
        REAL, DIMENSION(3), INTENT(IN)  :: A  ! vecteur A
        REAL, DIMENSION(3), INTENT(IN)  :: B  ! vecteur B
        REAL,               INTENT(OUT) :: S  ! resultat du produit scalaire A.B
         INTEGER :: i

        S = 0.0
        DO  I = 1,3
         S = S+A(I)*B(I)
        end do

      return
      end subroutine scalpr

end module math_module

