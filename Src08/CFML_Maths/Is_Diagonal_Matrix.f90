!!----
!!---- SUBMODULE CFML_Maths
!!----
!!----
!!
Submodule (CFML_Maths) Is_Diagonal_Matrix
 Contains
    !!----
    !!---- IS_DIAGONAL_MATRIX_I
    !!----    Determine if the matrix is diagonal
    !!----
    !!---- 28/03/2019 
    !!
    Module Pure Function Is_Diagonal_Matrix_I(A) Result(info)
       !---- Arguments ----! 
       integer, dimension(:,:), intent(in)  :: A
       logical                              :: info

       !---- Local Variables ----!
       integer :: i,j

       !> Init
       info = .true.
       
       do j=1, size(A,2)
          do i=1, size(A,1)
             if (j /= i .and. A(i,j) /= 0) then
                info = .false.
                return
             end if
          end do
       end do

       return
    End Function Is_Diagonal_Matrix_I

    !!----
    !!---- IS_DIAGONAL_MATRIX_R
    !!----    Determine if the matrix is diagonal
    !!----
    !!---- 28/03/2019 
    !!
    Module Pure Function Is_Diagonal_Matrix_R(A) Result(info)
       !---- Arguments ----!
       real(kind=cp), dimension(:,:), intent(in)  :: A
       logical                                    :: info

       !---- Local Variables ----!
       integer :: i,j

       !> Init
       info = .true.
        
       do j=1, size(A,2)
          do i=1 , size(A,1)
             if (j /= i .and. abs(A(i,j)) > epsilon(1.0_sp)) then
                info = .false.
                return
             end if
          end do
       end do
       
       return
    End Function Is_Diagonal_Matrix_R

End Submodule Is_Diagonal_Matrix
