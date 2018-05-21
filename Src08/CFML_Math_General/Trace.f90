!!----
!!---- SUBMODULE CFML_Math_General
!!----
!!----
!!
Submodule (CFML_Math_General) Trace
 Contains
 
    !!--++ FUNCTION TRACE_C
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the trace of a complex nxn array
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Function Trace_C(a) Result(b)    
       !---- Argument ----!
       complex, dimension(:,:), intent(in) :: a
       complex                             :: b

       !---- Local variables ----!
       integer :: i,imax

       b=(0.0,0.0)
       imax=min(size(a,1),size(a,2))
       do i=1,imax
          b=b+a(i,i)
       end do

       return
    End Function Trace_C
 
    !!--++ FUNCTION TRACE_I
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the trace of an integer 3x3 array
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Function Trace_I(a) Result(b)    
       !---- Argument ----!
       integer, dimension(:,:), intent(in) :: a
       integer                             :: b

       !---- Local variables ----!
       integer :: i,imax

       b=0
       imax=min(size(a,1),size(a,2))
       do i=1,imax
          b=b+a(i,i)
       end do

       return
    End Function Trace_I
 
    !!--++ FUNCTION TRACE_R
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Calculates the trace of a real 3x3 array
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Function Trace_R(a) Result(b)    
       !---- Argument ----!
       real(kind=cp), dimension(:,:), intent(in) :: a
       real(kind=cp)                             :: b

       !---- Local variables ----!
       integer :: i,imax

       b=0.0
       imax=min(size(a,1),size(a,2))
       do i=1,imax
          b=b+a(i,i)
       end do

       return
    End Function Trace_R
 
   
End Submodule Trace
