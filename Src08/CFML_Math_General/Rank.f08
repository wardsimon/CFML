!!----
!!---- SUBMODULE CFML_Math_General
!!----
!!----
!!
Submodule (CFML_Math_General) Rank
 Contains
 
    !!--++ SUBROUTINE RANK_DP
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Computes the rank (in algebraic sense) of the rectangular matrix A.
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Subroutine Rank_dp(a,tol,r)    
       !---- Arguments ----!
       real(kind=dp), dimension(:,:),intent( in)      :: a     ! Input array
       real(kind=dp),                intent( in)      :: tol   ! Tolerance
       integer,                      intent(out)      :: r

       !---- Arguments ----!
       real(kind=dp), dimension(size(a,1),size(a,2))  :: u
       real(kind=dp), dimension(size(a,2))            :: w
       real(kind=dp), dimension(size(a,2),size(a,2))  :: v
       integer                                        :: i

       u=a
       call svdcmp(u,w,v)
       if (ERR_CFML) then
          r=0
       else
          r=0
          do i=1,size(a,2)
             if(w(i) > tol) r=r+1
          end do
       end if

       return
    End Subroutine Rank_dp
 
    !!--++ SUBROUTINE RANK_SP
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Computes the rank (in algebraic sense) of the rectangular matrix A.
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Subroutine Rank_sp(a,tol,r)    
       !---- Arguments ----!
       real(kind=sp), dimension(:,:),intent( in)      :: a
       real(kind=sp),                intent( in)      :: tol
       integer,                      intent(out)      :: r

       !---- Local variables ----!
       real(kind=sp), dimension(size(a,1),size(a,2))  :: u
       real(kind=sp), dimension(size(a,2))            :: w
       real(kind=sp), dimension(size(a,2),size(a,2))  :: v
       integer :: i

       u=a
       call svdcmp(u,w,v)
       if (ERR_CFML) then
          r=0
       else
          r=0
          do i=1,size(a,2)
             if(w(i) > tol) r=r+1
          end do
       end if

       return
    End Subroutine Rank_sp
 
   
End Submodule Rank
