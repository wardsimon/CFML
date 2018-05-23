Submodule (Test01) Test04
   implicit none
   
 Contains
 
    Module Function prueba10(a1, i2) result(v)
       real,             intent(in) :: a1
       integer,          intent(in) :: i2
       real                         :: v
       
       v= a1*rmax+real(i2)*nmax
       
       return
    End Function prueba10
    
End Submodule Test04