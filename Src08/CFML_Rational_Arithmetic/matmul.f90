Submodule (CFML_Rational_Arithmetic) Matmul

 Contains
    !!----
    !!---- FUNCTION RATIONAL_MATMUL_MATVEC
    !!----
    !!----
    !!
    Module Function Rational_Matmul_Matvec(Mat,Vec) Result(Vec_Out)
       !---- Arguments ----!
       type(rational), dimension(:,:), intent (in) :: mat
       type(rational), dimension(:),   intent (in) :: vec
       type(rational), dimension(size(vec))        :: vec_out
       
       !---- Local Variables ----!
       integer :: n1,n2,n3,i
 
       n1=size(mat,dim=1); n2=size(mat,dim=2); n3=size(vec)
       if (n1 == n2 .and. n2 == n3) then
          do i=1,n3
             vec_out(i) = rational_simplify(dot_product(mat(i,:),vec))
          end do
       else
          Err_CFML=.true.
          Err_CFML_Flag=0
          write(unit=Err_CFML_Msg,fmt="(a,3i4)") &
               "Error in MATMUL: the dimensions of the arguments are non-conformable",n1,n2,n3
       end if
       
       return
    End Function Rational_Matmul_Matvec

    !!----
    !!---- FUNCTION RATIONAL_MATMUL_MATMAT
    !!----
    !!----
    !!
    Module Function Rational_Matmul_Matmat(Mat1,Mat2) Result(Mat_Out)
       !---- Arguments ----! 
       type(rational), dimension(:,:), intent (in)                 :: mat1
       type(rational), dimension(:,:), intent (in)                 :: mat2
       type(rational),dimension(size(mat1,dim=1),size(mat2,dim=2)) :: mat_out

       !---- Local Variables ----!
       integer :: n1,n2,n3,n4,i,j
      
       n1=size(mat1,dim=1); n2=size(mat1,dim=2); n3=size(mat2,dim=1); n4=size(mat2,dim=2)
       if (n2 == n3) then
          do j=1,n4
             do i=1,n1
                mat_out(i,j) = rational_simplify(dot_product(mat1(i,:),mat2(:,j)))
             end do
          end do
       else
          Err_CFML=.true.
          Err_CFML_Flag=2
          write(unit=Err_CFML_Msg,fmt="(a,4i4)") &
               "Error in MATMUL: the dimensions of the arguments are non-conformable",n1,n2,n3,n4
       end if
       
       return
    End Function Rational_Matmul_Matmat
    
End Submodule Matmul