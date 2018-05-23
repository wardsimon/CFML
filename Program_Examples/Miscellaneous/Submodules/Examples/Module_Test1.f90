Module Test01
   implicit none
   
   integer, parameter :: NMAX=10
   real,    parameter :: RMAX=1.2345
   
   Interface  prueba
       Module Procedure prueba00
       Module Procedure prueba01
       Module Procedure prueba02
   End Interface
    
   Interface
      module subroutine prueba00(c1)
         character(len=*), intent(in) :: c1
      end subroutine prueba00
      
      module subroutine prueba01(i1)
         integer, intent(in) :: i1
      end subroutine prueba01
      
      module subroutine prueba02(r1)
         real, intent(in) :: r1
      end subroutine prueba02
      
      module function prueba10(a1, i2) result(v)
         real,             intent(in) :: a1
         integer,          intent(in) :: i2
         real                         :: v
      end function prueba10
      
      module subroutine prueba11(a1, i2, c3)
         real,             intent(in) :: a1
         integer,          intent(in) :: i2
         character(len=*), intent(in) :: c3
      end subroutine prueba11
   End Interface
   
   !Interface  prueba
   !    Module Procedure prueba00
   !    Module Procedure prueba01
   !    Module Procedure prueba02
   !End Interface
   
End Module Test01