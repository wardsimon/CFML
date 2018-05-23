Program EjSubmod
   Use Test01
   
   Implicit none
   
   real :: valor
   
   valor=prueba10(0.2341,6)
   
   print*,' Calling the function PRUEBA10... -> ',valor
   print*,' '
   print*,' Executing Subroutine PRUEBA11...'
   call prueba11(2.111,3,"Incredible!!!")
   
End Program EjSubmod