!!----
!!---- SUBMODULE CFML_Math_General
!!----
!!----
!!
Submodule (CFML_Math_General) Equal_Vector
 Contains
 
    !!--++ FUNCTION EQUAL_VECTOR_I
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if two integer vectors are equal in N
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Function Equal_Vector_I(a,b,n) result(info)    
       !---- Argument ----!
       integer, dimension(:),           intent(in) :: a,b    ! Input vectors
       integer,               optional, intent(in) :: n      ! Dimension of the vectors
       logical                                     :: info

       !---- Local variables ----!
       integer :: i,ndim,ndim2

       !> Init
       info=.false.

       !> Check
       ndim=size(a)
       ndim2=size(b)
       if (ndim /= ndim2) return

       if (present(n)) then
          ndim=min(ndim,n)
       end if
       ndim=max(ndim,1)

       do i=1,ndim
          if (a(i) /= b(i)) return
       end do
       info=.true.

       return
    End Function Equal_Vector_I
 
    !!--++ FUNCTION EQUAL_VECTOR_R
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if two real(kind=sp) vectors are equal in N
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Function Equal_Vector_R(a,b,n) result(info)    
       !---- Argument ----!
       real(kind=cp), dimension(:),           intent(in) :: a,b      ! Input vectors
       integer,                     optional, intent(in) :: n        ! Dimension of the vector
       logical                                           :: info

       !---- Local variables ----!
       integer :: i,ndim,ndim2

       !> init
       info=.false.

       !> Check
       ndim=size(a)
       ndim2=size(b)
       if (ndim /= ndim2) return

       if (present(n)) then
          ndim=min(ndim,n)
       end if
       ndim=max(ndim,1)

       do i=1,ndim
          if (abs(a(i) - b(i)) > epss ) return
       end do
       info=.true.

       return
    End Function Equal_Vector_R
 
End Submodule Equal_Vector
