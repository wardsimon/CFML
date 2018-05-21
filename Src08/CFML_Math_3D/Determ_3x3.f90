!!----
!!---- SUBMODULE CFML_Math_3D
!!----
!!----
!!
Submodule (CFML_Math_3D) Determ_3x3
 Contains
 
    !!--++ FUNCTION DETERM_A_I
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the determinant of an integer 3x3 matrix
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Function Determ_A_I(Array) Result(determ)    
       !---- Argument ----!
       integer, dimension(3,3), intent(in) :: array
       integer                             :: determ

       determ=array(1,1)*array(2,2)*array(3,3)+ &
              array(2,1)*array(3,2)*array(1,3)+ &
              array(1,2)*array(2,3)*array(3,1)- &
              array(1,3)*array(2,2)*array(3,1)- &
              array(1,1)*array(3,2)*array(2,3)- &
              array(1,2)*array(2,1)*array(3,3)

       return
    End Function Determ_A_I
    
    !!--++ FUNCTION DETERM_A_R
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the determinant of a real 3x3 matrix
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Function Determ_A_R(Array) Result (determ)    
       !---- Argument ----!
       real(kind=cp), dimension(3,3), intent(in) :: array
       real(kind=cp)                             :: determ

       determ=array(1,1)*array(2,2)*array(3,3)+ &
              array(2,1)*array(3,2)*array(1,3)+ &
              array(1,2)*array(2,3)*array(3,1)- &
              array(1,3)*array(2,2)*array(3,1)- &
              array(1,1)*array(3,2)*array(2,3)- &
              array(1,2)*array(2,1)*array(3,3)

       return
    End Function Determ_A_R
 
End Submodule Determ_3x3
