Submodule (Test01) Test03
   implicit none
   
 Contains
 
    Module Subroutine prueba11(a1, i2, c3)
       real,             intent(in) :: a1
       integer,          intent(in) :: i2
       character(len=*), intent(in) :: c3
       
       !---- Local Variables ----!
       real :: v
       
       v=(a1*real(i2)+nmax)*rmax
       
       print*,"The value for NMax is ", nmax
       print*,"The valor for RMax is ",rmax
       print*," "
       print*,"Final value: ",v
       print*,"********"
       call prueba(c3)
       call prueba(nmax)
       print*,"********"
       
       return
    End Subroutine prueba11  
    
End Submodule Test03