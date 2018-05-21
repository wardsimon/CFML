!!----
!!---- SUBMODULE CFML_Math_General
!!----
!!----
!!
Submodule (CFML_Math_General) Norm
 Contains
 
    !!----  FUNCTION EUCLIDEAN_NORM
    !!----
    !!----  This function calculates safely the Euclidean norm of a vector.
    !!----  Intermediate overflows are avoided using this function. The original
    !!----  name "enorm" from MINPACK has been changed and the subroutine has
    !!----  been translated to Fortran 90.
    !!----
    !!----
    !!--..  Original documentation (from MINPACK):
    !!--..
    !!--..  Function enorm
    !!--..
    !!--..  Given an n-vector x, this function calculates the euclidean norm of x.
    !!--..
    !!--..  The euclidean norm is computed by accumulating the sum of squares in
    !!--..  three different sums.  The sums of squares for the small and large
    !!--..  components are scaled so that no overflows occur.  Non-destructive
    !!--..  underflows are permitted.  Underflows and overflows do not occur in the
    !!--..  computation of the unscaled sum of squares for the intermediate
    !!--..  components.  The definitions of small, intermediate and large components
    !!--..  depend on two constants, rdwarf and rgiant.  The main restrictions on
    !!--..  these constants are that rdwarf**2 not underflow and rgiant**2 not
    !!--..  overflow.  The constants given here are suitable for every known computer.
    !!--..
    !!--..  The function statement is
    !!--..
    !!--..    REAL (kind=cp) function enorm(n,x)
    !!--..
    !!--..  where
    !!--..
    !!--..    n is a positive integer input variable.
    !!--..
    !!--..    x is an input array of length n.
    !!--..
    !!--..  Subprograms called
    !!--..
    !!--..    Fortran-supplied ... ABS,SQRT
    !!--..
    !!--..  Argonne National Laboratory. MINPACK project. march 1980.
    !!--..  Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
    !!----
    !!----  Update: August - 2009
    !!----
    Module Function Euclidean_Norm(x,n) Result(Fn_Val)    
       !---- Arguments ----!
       Real (Kind=cp), Dimension(:),           intent(In)  :: x      ! Input vector
       Integer,                      optional, intent(In)  :: n      ! Dimension of the vector
       Real (Kind=cp)                                      :: Fn_Val ! Return value

       !--- Local Variables ---!
       Integer                   :: i,ndim
       Real (Kind=cp)            :: agiant, floatn, s1, s2, s3, xabs, x1max, x3max
       Real (Kind=cp), Parameter :: one = 1.0_cp, zero = 0.0_cp, rdwarf = 3.834e-20_cp,  &
                                    rgiant = 1.304e+19_cp

       !> check
       ndim=size(x)
       if (present(n)) ndim=n

       s1 = zero
       s2 = zero
       s3 = zero
       x1max = zero
       x3max = zero
       floatn = ndim
       agiant = rgiant/floatn

       do i = 1, ndim
          xabs = Abs(x(i))
          if (.Not. (xabs > rdwarf .AND. xabs < agiant)) then
             ! sum for large components.
             if (xabs > rdwarf) then
                if (xabs > x1max) then
                   s1 = one + s1*(x1max/xabs)**2
                   x1max = xabs
                   cycle
                end if
                s1 = s1 + (xabs/x1max)**2
                cycle
             End If

             ! sum for small components.
             If (xabs > x3max) Then
                s3 = one + s3*(x3max/xabs)**2
                x3max = xabs
                Cycle
             End If

             If (xabs /= zero) s3 = s3 + (xabs/x3max)**2
             Cycle
          End if

          !  sum for intermediate components.
          s2 = s2 + xabs**2
       End Do

       ! calculation of norm.
       If (s1 /= zero) Then
          Fn_Val = x1max*Sqrt(s1 + (s2/x1max)/x1max)
          Return
       End If

       If (s2 /= zero) Then
          If (s2 >= x3max) Fn_Val = Sqrt(s2*(one + (x3max/s2)*(x3max*s3)))
          If (s2 < x3max) Fn_Val = Sqrt(x3max*((s2/x3max) + (x3max*s3)))
          Return
       End If

       Fn_Val = x3max*Sqrt(s3)

       Return
    End Function Euclidean_Norm
 
    !!--++ FUNCTION NORM_I
    !!--++
    !!--++    Calculate the Norm of a vector
    !!--++
    !!--++ Update: April - 2009
    !!
    Module Function Norm_I(X,G) Result(R)    
       !---- Arguments ----!
       integer,       dimension(:),   intent(in) :: x    ! Input vector
       real(kind=cp), dimension(:,:), intent(in) :: g    ! Metric array
       real(kind=cp)                             :: r    ! Norm of the input vector

       if (size(x)*size(x) /= size(g)) then
          r=tiny(0.0)
       else
          r=sqrt(dot_product(real(x), matmul(g,real(x))))
       end if

       return
    End Function Norm_I
 
    !!--++ FUNCTION NORM_R
    !!--++
    !!--++    Calculate the Norm of a vector
    !!--++
    !!--++ Update: April - 2009
    !!
    Module Function Norm_R(X,G) Result(R)    
       !---- Arguments ----!
       real(kind=cp), dimension(:),   intent(in) :: x   ! Input vector
       real(kind=cp), dimension(:,:), intent(in) :: g   ! Metrics
       real(kind=cp)                             :: r   ! Norm of the vector

       if (size(x)*size(x) /= size(g)) then
          r=tiny(0.0)
       else
          r=sqrt(dot_product(x, matmul(g,x)))
       end if

       return
    End Function Norm_R
 
   
End Submodule Norm
