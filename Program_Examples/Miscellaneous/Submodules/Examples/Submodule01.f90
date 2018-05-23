Submodule (Test01) Prueba
   implicit none
   
 Contains
 
    Module Subroutine prueba00(c1) 
       character(len=*), intent(in) :: c1

       print*," ===> "//trim(c1)        
       
       return
    End Subroutine prueba00
    
    Module Subroutine prueba01(i1) 
       integer, intent(in) :: i1

       print*," ===> ",i1        
       
       return
    End Subroutine prueba01
    
    Module Subroutine prueba02(r1) 
       real, intent(in) :: r1

       print*," ===> ",r1        
       
       return
    End Subroutine prueba02
    
End Submodule Prueba