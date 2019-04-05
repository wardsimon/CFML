!!----
!!---- SUBMODULE CFML_Math_General
!!----
!!----
!!
Submodule (CFML_Maths) CFML_Math_009
 Contains

    !!----
    !!---- SPLINE
    !!----    Spline N points
    !!----
    !!---- 04/04/2019 
    !!
    Module Pure Subroutine Spline(x,y,n,yp1,ypn,ys)
       !---- Arguments ----!
       real(kind=cp), dimension(:), intent(in)  :: x               !  In -> Array X
       real(kind=cp), dimension(:), intent(in)  :: y               !  In -> Array Yi=F(Xi)
       integer ,                    intent(in)  :: n               !  In -> Dimension of X, Y
       real(kind=cp),               intent(in)  :: yp1             !  In -> Derivate of Point 1
       real(kind=cp),               intent(in)  :: ypn             !  In -> Derivate of Point N
       real(kind=cp), dimension(:), intent(out) :: ys              ! Out -> array containing second derivatives

       !---- Local Variables ----!
       integer                     :: i, k
       real(kind=cp), dimension(n) :: u
       real(kind=cp)               :: sig, p, qn, un

       if (yp1 > huge(1.0_cp)) then
          ys(1)=0.0
          u(1)=0.0
       else
          ys(1)=-0.5
          u(1)=(3.0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
       end if

       do i=2,n-1
          sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
          p=sig*ys(i-1)+2.0
          ys(i)=(sig-1.0)/p
          u(i)=(6.0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))  &
               /(x(i+1)-x(i-1))-sig*u(i-1))/p
       end do
       if (ypn > huge(1.0_cp)) then
          qn=0.0
          un=0.0
       else
          qn=0.5
          un=(3.0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
       end if
       ys(n)=(un-qn*u(n-1))/(qn*ys(n-1)+1.0)
       do k=n-1,1,-1
          ys(k)=ys(k)*ys(k+1)+u(k)
       end do

       return
    End Subroutine Spline

End Submodule CFML_Math_009
