!!----
!!---- SUBMODULE CFML_Math_General
!!----
!!----
!!
Submodule (CFML_Math_General) Determinant
 Contains
 
    !!--++ SUBROUTINE DETERMINANT_C
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the determinant of a real square matrix.
    !!--++    Calculates the pseudo-determinant of a complex square matrix.
    !!--++    The calculated value is only useful for linear dependency purposes.
    !!--++    It tell us if the complex matrix is singular or not.
    !!--++
    !!--++    P R O V I S I O N A L (The determinant of A is not calculated at present)
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Subroutine Determinant_C(A,n,determ)    
       !---- Arguments ----!
       complex, dimension(:,:), intent( in) :: A         !input square matrix (n,n)
       integer,                 intent( in) :: n         !actual dimension of A
       real(kind=cp),           intent(out) :: determ    !det(A) if real and det(AR)^2 + det(AI)^2 if complex

       !---- local variables ----!
       real(kind=cp),    dimension(2*n,2*n) :: AC   !real square matrix
       real(kind=cp)                        :: d
       integer                              :: i,nn
       logical                              :: singular

       nn=2*n
       AC(  1:n ,  1:n ) =  real(A(1:n ,1:n))
       AC(n+1:nn,  1:n ) = aimag(A(1:n ,1:n))
       AC(n+1:nn,n+1:nn) =    AC(  1:n ,1:n)
       AC(  1:n ,n+1:nn) =   -AC(n+1:nn,1:n)

       call lu_decomp(ac(1:nn,1:nn),d,singular)

       if (singular) then
          determ=0.0
       else
          determ=0.0
          do i=1,nn
             d=d*sign(1.0_cp,ac(i,i))
             determ=determ+ log(abs(ac(i,i)))
          end do
          determ=d*exp(determ)
       end if

       return
    End Subroutine Determinant_C
 
    !!--++ SUBROUTINE DETERMINANT_R
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the determinant of a real square matrix.
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Subroutine Determinant_R(A,n,determ)    
       !---- Arguments ----!
       real(kind=cp), dimension(:,:), intent( in) :: A      ! Input array NxN
       integer,                       intent( in) :: n      ! Dimension of A
       real(kind=cp),                 intent(out) :: determ ! Value

       !---- local variables ----!
       real(kind=cp),    dimension(n,n)  :: AC
       real(kind=cp)                     :: d
       integer                           :: i
       logical                           :: singular

       ac=A(1:n,1:n)
       call lu_decomp(ac,d,singular)

       if (singular) then
          determ=0.0
       else
          determ=0.0
          do i=1,n
             d=d*sign(1.0_cp,ac(i,i))
             determ=determ + log(abs(ac(i,i)))
          end do
          determ=d*exp(determ)
       end if

       return
    End Subroutine Determinant_R
 
   
End Submodule Determinant
