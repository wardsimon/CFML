!!----
!!---- SUBMODULE CFML_Math_General
!!----
!!----
!!
Submodule (CFML_Math_General) Co_Linear
 Contains
 
    !!--++ FUNCTION CO_LINEAR_C
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if two complex vectors are co-linear
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Function Co_linear_C(a,b,n) Result(co_linear)    
       !---- Argument ----!
       complex, dimension(:), intent(in)           :: a,b    ! Complex vectors
       integer,               optional, intent(in) :: n      ! Dimension of the vectors
       logical                                     :: co_linear

       !---- Local variables ----!
       integer :: i,ia,ib,ndim
       complex :: c

       !> Init
       ndim=size(a)
       if (present(n)) ndim=n

       co_linear=.true.

       do i=1,ndim
          if (abs(a(i)) > epss) then
             ia=i
             exit
          end if
       end do

       do i=1,ndim
          if (abs(b(i)) > epss) then
             ib=i
             exit
          end if
       end do

       if (ia /= ib) then
          co_linear=.false.
          return
       else
          c=a(ia)/b(ib)
          do i=1,ndim
             if (abs(a(i)-c*b(i)) > epss) then
                co_linear=.false.
                return
             end if
          end do
       end if

       return
    End Function Co_linear_C
 
    !!--++ FUNCTION CO_LINEAR_I
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if two integer vectors are co-linear
    !!--++
    !!--++ Update: October - 2008
    !!
    Module Function Co_linear_I(a,b,n) Result(co_linear)    
       !---- Argument ----!
       integer, dimension(:),           intent(in) :: a,b        ! Input vectors
       integer,               optional, intent(in) :: n          ! Dimension of the vector
       logical                                     :: co_linear

       !---- Local variables ----!
       integer       :: i,ia,ib,ndim
       real(kind=cp) :: c

       !> init
       ndim=size(a)
       if (present(n)) ndim=n

       co_linear=.true.

       do i=1,ndim
          if (abs(a(i)) > 0) then
             ia=i
             exit
          end if
       end do

       do i=1,ndim
          if (abs(b(i)) > 0) then
             ib=i
             exit
          end if
       end do

       if (ia /= ib) then
          co_linear=.false.
          return
       else
          c=real(a(ia))/real(b(ib))
          do i=1,ndim
             if (abs( real(a(i))-c*real(b(i)) ) > epss) then
                co_linear=.false.
                return
             end if
          end do
       end if

       return
    End Function Co_linear_I
 
    !!--++ FUNCTION CO_LINEAR_R
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if two real vectors are co-linear
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Function Co_linear_R(a,b,n) Result(co_linear)    
       !---- Argument ----!
       real(kind=cp), dimension(:),           intent(in) :: a,b        ! Input real vectors
       integer,                     optional, intent(in) :: n          ! Dimension of the vectors
       logical                                           :: co_linear

       !---- Local variables ----!
       integer       :: i,ia,ib,ndim
       real(kind=cp) :: c

       !> init
       ndim=size(a)
       if (present(n)) ndim=n

       co_linear=.true.

       do i=1,ndim
          if (abs(a(i)) > epss) then
             ia=i
             exit
          end if
       end do

       do i=1,ndim
          if (abs(b(i)) > epss) then
             ib=i
             exit
          end if
       end do
       if (ia /= ib) then
          co_linear=.false.
          return
       else
          c=a(ia)/b(ib)
          do i=1,ndim
             if (abs(a(i)-c*b(i)) > epss) then
                co_linear=.false.
                return
             end if
          end do
       end if

       return
    End Function Co_linear_R
 
End Submodule Co_Linear
