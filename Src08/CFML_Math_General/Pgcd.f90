!!----
!!---- SUBMODULE CFML_Math_General
!!----
!!----
!!
Submodule (CFML_Math_General) CFML_MG_06
 Contains
 
    !!---- FUNCTION PGCD
    !!----
    !!----    Function calculating the maximum common divisor of two integers
    !!----
    !!---- Update: February - 2005
    !!
    Module Elemental Function Pgcd(a,b) Result(mcd)    
       !---- Arguments ----!
       integer, intent(in) :: a,b   ! Integers
       integer             :: mcd   ! Maximum common divisor

       !---- Local variables ----!
       integer  :: u,v,m

       u=max(a,b)
       v=min(a,b)
       m=0
       do
          if (m == 1) exit
          m=mod(u,v)
          u=v
          v=m
       end do
       mcd=u

       return
    End Function Pgcd
 
End Submodule CFML_MG_06
