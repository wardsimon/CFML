!!----
!!---- SUBMODULE CFML_Math_General
!!----
!!----
!!
Submodule (CFML_Math_General) CFML_MG_09
 Contains
 
    !!---- FUNCTION PPCM
    !!----
    !!----    Function calculating the minimum common multiple of two integers
    !!----
    !!---- Update: February - 2005
    !!
    Module Elemental Function Ppcm(a,b) result(mcm)    
       !---- Arguments ----!
       integer, intent(in) :: a,b    ! Integers
       integer             :: mcm    ! Minimum common multiple

       !---- Local variables ----!
       integer :: u,v,w,i

       u=max(a,b)
       v=min(a,b)
       mcm=1
       if (v <= 1) then
          mcm=u
          return
       end if
       w=int(sqrt(real(u)))+1
       do i=2,w
          do
             if(.not. ((mod(u,i)==0) .or. (mod(v,i)==0)) ) exit
             mcm=mcm*i
             if (modulo(u,i) == 0) u=u/i
             if (modulo(v,i) == 0) v=v/i
          end do
       end do

       return
    End Function Ppcm
 
End Submodule CFML_MG_09
