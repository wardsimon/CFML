!!----
!!---- SUBMODULE CFML_Math_3D
!!----
!!----
!!
Submodule (CFML_Maths) CFML_Math_010
 Contains
    !!----
    !!---- VEC_LENGTH
    !!----    Length of vector B when A is the Crystallographic
    !!----    to orthogonal matrix length=c
    !!----
    !!---- Update: February - 2005
    !!
    Pure Module Function Vec_Length(G,Vec) Result(c)
       !---- Arguments ----!
       real(kind=cp), intent(in)  , dimension(3,3)       :: G      ! Metric array
       real(kind=cp), intent(in)  , dimension(3  )       :: Vec    ! Vector
       real(kind=cp)                                     :: c      ! Length of Vector

       !---- Local variables ----!
       integer                     :: i,j
       real(kind=cp), dimension(3) :: v

       v=0.0_cp
       do i = 1,3
          do j = 1,3
             v(i) = v(i)+G(i,j)*Vec(j)
          end do
       end do

       c = sqrt(v(1)**2+v(2)**2+v(3)**2)

       return
    End Function Vec_Length

End Submodule CFML_Math_010
