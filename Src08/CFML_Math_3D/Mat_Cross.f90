!!----
!!---- SUBMODULE CFML_Math_3D
!!----
!!----
!!
Submodule (CFML_Math_3D) Mat_Cross
 Contains
 
    !!--++ FUNCTION MAT_CROSS_CMPL_DP
    !!--++
    !!--++ OVERLOADED
    !!--++    Calculates the matrix corresponding to the operator u x
    !!--++    Antisymmetric matrix of the form:
    !!--++                /  0   -u(3)  u(2)\
    !!--++    M=[u]cross=|  u(3)   0   -u(1) |
    !!--++                \-u(2)  u(1)   0  /
    !!--++
    !!--++  Updated: June - 2012
    !!
    Module Function Mat_Cross_cmpl_dp(Vec) Result(M)    
       !---- Argument ----!
       complex(kind=dp), dimension(3), intent( in) :: Vec
       complex(kind=dp), dimension(3,3)            :: M

       M = reshape( (/  (0.0_dp,0.0_dp),   -vec(3),         vec(2),  &
                            vec(3),   (0.0_dp,0.0_dp),     -vec(1),  &
                           -vec(2),           vec(1),   (0.0_dp,0.0_dp)/),(/3,3/))
       return
    End Function Mat_Cross_cmpl_dp
 
    !!--++ FUNCTION MAT_CROSS_CMPL_SP
    !!--++
    !!--++ OVERLOADED
    !!--++    Calculates the matrix corresponding to the operator u x
    !!--++    Antisymmetric matrix of the form:
    !!--++                /  0   -u(3)  u(2)\
    !!--++    M=[u]cross=|  u(3)   0   -u(1) |
    !!--++                \-u(2)  u(1)   0  /
    !!--++
    !!--++  Updated: June - 2012
    !!
    Module Function Mat_Cross_cmpl_sp(Vec) Result(M)    
       !---- Argument ----!
       complex(kind=sp), dimension(3), intent( in) :: Vec
       complex(kind=sp), dimension(3,3)            :: M

       M = reshape( (/  (0.0_sp,0.0_sp),   -vec(3),         vec(2),  &
                            vec(3),   (0.0_sp,0.0_sp),     -vec(1),  &
                           -vec(2),           vec(1),   (0.0_sp,0.0_sp)/),(/3,3/))
       return
    End Function Mat_Cross_cmpl_sp
 
    !!--++ FUNCTION MAT_CROSS_DP
    !!--++
    !!--++ OVERLOADED
    !!--++    Calculates the matrix corresponding to the operator u x
    !!--++    Antisymmetric matrix of the form:
    !!--++                /  0   -u(3)  u(2)\
    !!--++    M=[u]cross=|  u(3)   0   -u(1) |
    !!--++                \-u(2)  u(1)   0  /
    !!--++
    !!--++  Updated: June - 2012
    !!
    Module Function Mat_Cross_dp(Vec) Result(M)    
       !---- Argument ----!
       real(kind=dp), dimension(3), intent( in) :: Vec
       real(kind=dp), dimension(3,3)            :: M

       M = reshape( (/ 0.0_dp,   -vec(3),     vec(2),  &
                        vec(3),    0.0_dp,   -vec(1),  &
                       -vec(2),     vec(1),    0.0_dp/),(/3,3/))
       return
    End Function Mat_Cross_dp
 
    !!--++ FUNCTION MAT_CROSS_IN
    !!--++
    !!--++ OVERLOADED
    !!--++    Calculates the matrix corresponding to the operator u x
    !!--++    Antisymmetric matrix of the form:
    !!--++                /  0   -u(3)  u(2)\
    !!--++    M=[u]cross=|  u(3)   0   -u(1) |
    !!--++                \-u(2)  u(1)   0  /
    !!--++
    !!--++  Updated: June - 2012
    !!
    Module Function Mat_Cross_in(Vec) Result(M)    
       !---- Argument ----!
       integer, dimension(3), intent( in) :: Vec
       integer, dimension(3,3)            :: M

       M = reshape( (/   0,    -vec(3),    vec(2),  &
                        vec(3),    0,     -vec(1),  &
                       -vec(2),   vec(1),     0 /),(/3,3/))
       return
    End Function Mat_Cross_in
 
    !!--++ FUNCTION MAT_CROSS_SP
    !!--++
    !!--++ OVERLOADED
    !!--++    Calculates the matrix corresponding to the operator u x
    !!--++    Antisymmetric matrix of the form:
    !!--++                /  0   -u(3)  u(2)\
    !!--++    M=[u]cross=|  u(3)   0   -u(1) |
    !!--++                \-u(2)  u(1)   0  /
    !!--++
    !!--++  Updated: June - 2012
    !!
    Module Function Mat_Cross_sp(Vec) Result(M)    
       !---- Argument ----!
       real(kind=sp), dimension(3), intent( in) :: Vec
       real(kind=sp), dimension(3,3)            :: M

       M = reshape( (/ 0.0_sp, -vec(3),    vec(2),  &
                        vec(3),  0.0_sp,  -vec(1),  &
                       -vec(2),   vec(1),   0.0_sp/),(/3,3/))
       return
    End Function Mat_Cross_sp
 
   
End Submodule Mat_Cross
