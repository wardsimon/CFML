!!----
!!---- SUBMODULE CFML_Math_General
!!----
!!----
!!
Submodule (CFML_Math_General) Equal_Matrix
 Contains
 
    !!--++ FUNCTION EQUAL_MATRIX_I
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if two integer arrays are equal in NxN
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Function Equal_Matrix_I(a,b,n) result(info)    
       !---- Argument ----!
       integer, dimension(:,:),          intent(in) :: a,b     ! Input arrays (NxN)
       integer,                optional, intent(in) :: n       ! Dimension of Arrays
       logical                                      :: info

       !---- Local variables ----!
       integer :: i,j,ndim,ndim2

       !> Init
       info=.false.

       !> Check on a
       ndim=size(a,dim=1)
       ndim2=size(a,dim=2)
       if (ndim /= ndim2) return

       !> Check on b
       ndim=size(b,dim=1)
       ndim2=size(b,dim=2)
       if (ndim /= ndim2) return

       !> Same shape
       ndim=size(a,dim=1)
       ndim2=size(b,dim=1)
       if (ndim /= ndim2) return

       if (present(n)) then
          if (n > ndim) return
          ndim=n
       end if

       do i=1,ndim
          do j=1,ndim
             if (a(i,j) /= b(i,j)) return
          end do
       end do
       info=.true.

       return
    End Function Equal_Matrix_I
 
    !!--++ FUNCTION EQUAL_MATRIX_R
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if two integer arrays are equal in NxN
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Function Equal_Matrix_R(a,b,n) result(info)    
       !---- Argument ----!
       real(kind=cp), dimension(:,:),          intent(in) :: a,b      ! Input arrays NxN
       integer,                      optional, intent(in) :: n        ! Dimensions N
       logical                                            :: info

       !---- Local variables ----!
       integer :: i,j,ndim,ndim2

       !> init
       info=.false.

       !> Check on a
       ndim=size(a,dim=1)
       ndim2=size(a,dim=2)
       if (ndim /= ndim2) return

       !> Check on b
       ndim=size(b,dim=1)
       ndim2=size(b,dim=2)
       if (ndim /= ndim2) return

       !> Same shape
       ndim=size(a,dim=1)
       ndim2=size(b,dim=1)
       if (ndim /= ndim2) return

       if (present(n)) then
          if (n > ndim) return
          ndim=n
       end if

       do i=1,ndim
          do j=1,ndim
             if (abs(a(i,j) - b(i,j)) > epss ) return
          end do
       end do
       info=.true.

       return
    End Function Equal_Matrix_R
 
End Submodule Equal_Matrix
