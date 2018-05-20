!!----
!!---- SUBMODULE CFML_Math_General
!!----
!!----
!!
Submodule (CFML_Math_General) In_Limits
 Contains
 
    !!--++ FUNCTION IN_LIMITS_DP
    !!--++
    !!--++   Logical function that is true if all the components of the vector vect
    !!--++   are within the limits:   limits(1,i)  <= vect(i) <=  limits(2,i), for all i.
    !!--++
    !!--++   Updated: March - 2013
    !!
    Module Function in_limits_dp(v,limits,n) result(ok)    
       !---- Arguments ----!
       real(kind=dp), dimension(:),             intent(in) :: v        ! Input vector
       real(kind=dp), dimension(:,:),           intent(in) :: limits   ! Normally (2,n)
       integer,                       optional, intent(in) :: n        ! Dimension of Vect
       logical                                             :: ok

       !---- Local Variables ----!
       integer :: i,ndim,ndim2

       !> Init
       ndim=size(v)
       ndim2=size(limits,2)
       if (ndim /= ndim2) then
          ok=.false.
          return
       end if

       if (present(n)) then
          ndim=min(n,ndim)
       end if
       ndim=max(ndim,1)

       ok=.true.
       do i=1,ndim
          if (v(i) >= limits(1,i) .and. v(i) <= limits(2,i)) cycle
          ok=.false.
          exit
       end do

       return
    End Function in_limits_dp
 
    !!
    Module Function in_limits_int(v,limits,n) result(ok)    
       !---- Arguments ----!
       integer, dimension(:),             intent(in) :: v        ! Input Vector
       integer, dimension(:,:),           intent(in) :: limits   ! Normally (2,n)
       integer,                 optional, intent(in) :: n        ! Dimension of vect
       logical                                       :: ok

       !---- Local arguments ----!
       integer :: i,ndim,ndim2

       !> Init
       ndim=size(v)
       ndim2=size(limits,2)
       if (ndim /= ndim2) then
          ok=.false.
          return
       end if

       if (present(n)) then
          ndim=min(n,ndim)
       end if
       ndim=max(ndim,1)

       ok=.true.
       do i=1,ndim
          if (v(i) >= limits(1,i) .and. v(i) <= limits(2,i)) cycle
          ok=.false.
          exit
       end do

       return
    End Function in_limits_int
 
    !!--++ FUNCTION IN_LIMITS_SP
    !!--++
    !!--++   Logical function that is true if all the components of the vector vect
    !!--++   are within the limits:   limits(1,i)  <= vect(i) <=  limits(2,i), for all i.
    !!--++
    !!--++   Updated: March - 2013
    !!
    Module Function in_limits_sp(v,limits,n) result(ok)    
       !---- Arguments ----!
       real(kind=sp), dimension(:),             intent(in) :: v        ! Input vector
       real(kind=sp), dimension(:,:),           intent(in) :: limits   ! Normally (2,n)
       integer,                       optional, intent(in) :: n        ! Dimension of the vector V
       logical                                             :: ok

       !---- Local Variables ----!
       integer :: i,ndim,ndim2

       !> Init
       ndim=size(v)
       ndim2=size(limits,2)
       if (ndim /= ndim2) then
          ok=.false.
          return
       end if

       if (present(n)) then
          ndim=min(n,ndim)
       end if
       ndim=max(ndim,1)

       ok=.true.
       do i=1,ndim
          if (v(i) >= limits(1,i) .and. v(i) <= limits(2,i)) cycle
          ok=.false.
          exit
       end do

       return
    End Function in_limits_sp
 
   
End Submodule In_Limits
