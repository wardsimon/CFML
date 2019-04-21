SubModule (CFML_Groups) CFML_Grp_019
   Contains
   
   !!----
   !!---- ALLOCATE_OPERATOR
   !!----
   !!---- 19/04/19
   !!
   Module Subroutine Allocate_Operator(d, Op)
      !---- Arguments ----!
      integer,              intent(in)     :: d
      type(Symm_Oper_Type), intent(in out) :: Op
      
      
      if (allocated(Op%Mat)) deallocate(Op%Mat)
      allocate(Op%Mat(d,d))
       
      !>Inititalize to identity matrix
      call Set_Identity_Matrix(d)
      
      Op%Mat=Identity_Matrix
      Op%time_inv=1
      Op%dt=1
      
      return
   End Subroutine Allocate_Operator
   
   !!----
   !!---- ALLOCATE_OPERATORS
   !!----    integer,              intent(in)     :: d,multip
   !!----    type(Symm_Oper_Type), dimension(:), allocatable, intent(in out) :: Op
   !!----
   !!----  Multip is the expected maximum number of operators
   !!----
   Module Subroutine Allocate_Operators(D,Multip,Op)
      !---- Arguments ----!
      integer,                                         intent(in)     :: d       ! Dimension
      integer,                                         intent(in)     :: multip  ! is the expected maximum number of operators
      type(Symm_Oper_Type), dimension(:), allocatable, intent(in out) :: Op
      
      !---- Local Variables ----!
      integer :: i
      
      !> Init
      if(allocated(Op)) deallocate(Op)
      allocate(Op(multip))
      
      do i=1,multip
         call Allocate_Operator(d,Op(i))
      end do
      
      return
   End Subroutine Allocate_Operators
   
End SubModule CFML_Grp_019 