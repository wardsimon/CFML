!!----
!!---- SUBMODULE CFML_Maths
!!----
!!----
!!
Submodule (CFML_Maths) Is_Null_Vector
 Contains
    !!----
    !!---- IS_NULL_VECTOR_I
    !!----    Determine if a vector is null
    !!----
    !!---- 28/03/2019 
    !!
    Module Pure Function Is_Null_Vector_I(V) Result(info)
       !---- Arguments ----!
       integer,  dimension(:), intent(in)  :: V
       logical                             :: Info
       
       !---- Local Variables ----!
       integer :: i
       
       info= .true.
       do i=1, size(v)
          if (v(i) /= 0) then
             info= .false.
             exit
          end if
       end do
       
       return
    End Function Is_Null_Vector_I

    !!----
    !!---- IS_NULL_VECTOR_R
    !!----    Determine if a vector is null
    !!----
    !!---- 28/03/2019 
    !!
    Module Pure Function Is_Null_Vector_R(V) result(info)
       !---- Arguments ----!
       real(kind=cp), dimension(:), intent(in)  :: v
       logical                                  :: info
       
       !---- Local Variables ----!
       integer :: i

       info= .true.
       do i=1, size(v)
          if (abs(v(i)) > epsilon(1.0_sp)) then
             info= .false.
             exit
          end if
       end do

       return
    End Function Is_Null_Vector_R
    
End Submodule Is_Null_Vector
