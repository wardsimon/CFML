!!----
!!---- SUBMODULE CFML_Math_General
!!----
!!----
!!
Submodule (CFML_Math_General) Swap
 Contains
 
    !!--++ SUBROUTINE SWAP_C
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Elemental Subroutine Swap_C(a,b)    
       !---- Arguments ----!
       complex, intent(in out) :: a
       complex, intent(in out) :: b

       !---- Local variables ----!
       complex :: dum

       dum=a
       a=b
       b=dum

       return
    End Subroutine Swap_C
 
    !!--++  SUBROUTINE SWAP_I
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Elemental Subroutine Swap_I(A,B)    
       !---- Arguments ----!
       integer , intent(in out) :: a
       integer , intent(in out) :: b

       !---- Local variables ----!
       integer  :: dum

       dum=a
       a=b
       b=dum

       return
    End Subroutine Swap_I
 
    !!--++  SUBROUTINE SWAP_R
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Elemental Subroutine Swap_R(A,B)    
       !---- Arguments ----!
       real(kind=cp), intent(in out) :: a
       real(kind=cp), intent(in out) :: b

       !---- Local variables ----!
       real(kind=cp) :: dum

       dum=a
       a=b
       b=dum

       return
    End Subroutine Swap_R
 
    !!--++  SUBROUTINE SWAP_MASKED_R
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Swap the contents of a and b if mask=.true.
    !!--++
    !!--++ Update: February - 2005
    !!
    Module Elemental Subroutine Swap_Masked_R(A,B,Mask)    
       !---- Arguments ----!
       real(kind=cp), intent(in out) :: a
       real(kind=cp), intent(in out) :: b
       logical,           intent(in) :: mask

       !---- Local Variables ----!
       real(kind=cp) :: swp

       if (mask) then
          swp=a
          a=b
          b=swp
       end if

       return
    End Subroutine Swap_Masked_R
 
   
End Submodule Swap
