!!----
!!---- SUBMODULE CFML_Math_General
!!----
!!----
!!
Submodule (CFML_Math_General) Zbelong
 Contains
 
    !!--++ FUNCTION ZBELONGM
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if a real array is an Integer matrix
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Function ZbelongM(v) Result(belong)    
       !---- Argument ----!
       real(kind=cp),   dimension(:,:), intent( in) :: v        ! Input array
       logical                                      :: belong

       !---- Local variables ----!
       real(kind=cp),   dimension(size(v,1),size(v,2)) :: vec

       vec= abs(real(nint (v))-v)
       belong=.not. ANY(vec > epss)

       return
    End Function ZbelongM
 
    !!--++ FUNCTION ZBELONGN
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if a real number is an Integer
    !!--++
    !!--++ Update: February - 2005
   !!
    Module Function ZbelongN(a) Result(belong)    
       !---- Argument ----!
       real(kind=cp), intent( in) :: a              ! Input number
       logical                    :: belong

       belong=.false.
       if (abs(real(nint (a))-a) > epss) return
       belong=.true.

       return
    End Function ZbelongN
 
    !!--++ FUNCTION ZBELONGV
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Determines if a real vector is an Integer vector
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Function ZbelongV(v) Result(belong)    
       !---- Argument ----!
       real(kind=cp),   dimension(:), intent( in) :: v      ! Input vector
       logical                                    :: belong

       !---- Local variables ----!
       integer                             :: i
       real(kind=cp),   dimension(size(v)) :: vec

       belong=.false.
       vec= abs(real(nint (v))-v)
       do i=1,size(v)
          if (vec(i) > epss) return
       end do
       belong=.true.

       return
    End Function ZbelongV
 
   
End Submodule Zbelong
