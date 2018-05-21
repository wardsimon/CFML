!!----
!!---- SUBMODULE CFML_Math_General
!!----
!!----
!!
Submodule (CFML_Math_General) CFML_MG_10
 Contains
 
    !!---- SUBROUTINE SMOOTHINGVEC
    !!----
    !!----    Procedure to smooth the vector values
    !!----
    !!---- Update: January - 2006
    !!
    Module Subroutine SmoothingVec(Y, N, Niter, Ys)    
       !---- Arguments ----!
       real(kind=cp),dimension(:),            intent(in out) :: Y         !  In Out-> Array to be smoothed
       integer,                               intent(in)     :: n         !  In -> Number of points
       integer,                               intent(in)     :: niter     !  In -> Number of iterations
       real(kind=cp),dimension(:), optional,  intent(out)    :: Ys        !  Out-> Array smoothed

       !---- Local Variables ----!
       integer                     :: n1, n2
       integer                     :: i, iter
       real(kind=cp), dimension (n):: datYs


       n1 = 4
       n2 = n-3

       do iter = 1 ,niter
          datYs(n1-1)=((Y(n1-2)+Y(n1))*10.0+(Y(n1-3)+Y(n1+1))*5.0+Y(n1+2))/31.0
          datYs(n1-2)=((Y(n1-3)+Y(n1-1))*10.0+Y(n1)*5.0+Y(n1+1))/26.0
          datYs(n1-3)=(Y(n1-2)*10.0+Y(n1-1)*5.0+Y(n1))/16.0

          do i=n1,n2
             datYs(i)=(Y(i-3)+Y(i+3)+5.0*(Y(i-2)+Y(i+2))+10.0*(Y(i-1)+Y(i+1)))/ 32.0
          end do

          datYs(n2+1)=((Y(n2+2)+Y(n2))*10.0+(Y(n2+3)+Y(n2-1))*5.0+Y(n2-2))/31.0
          datYs(n2+2)=((Y(n2+3)+Y(n2+1))*10.0+Y(n2)*5.0+Y(n2-1))/26.0
          datYs(n2+3)=(Y(n2+2)*10.0+Y(n2+1)*5.0+Y(n2))/16.0

          if(present(Ys)) then
             Ys(1:n) = datYs(1:n)
          else
             Y(1:n) = datYs(1:n)
          end if
       end do

       return
    End Subroutine SmoothingVec
 
End Submodule CFML_MG_10
