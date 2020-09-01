!!----
!!----
!!----
!!
SubModule (CFML_gSpaceGroups) SPG_Sort_Oper
   implicit none
   Contains

   !!----
   !!---- SORT_OP
   !!----
   !!---- 20/04/19
   !!
   Module Subroutine Sort_Oper(N, Op, Cod)
      !---- Arguments ----!
      integer,                            intent(in)     :: n
      type(Symm_Oper_Type) ,dimension(n), intent(in out) :: Op
      character(len=*),                   intent(in)     :: cod

      !---- Local Variables ----!
      type(Symm_Oper_Type)  :: Ops
      integer               :: i, j,opso
      integer, dimension(n) :: option

      !> Init
      if (cod == "tim") then
         option=Op(:)%time_inv
      else
         option=Op(:)%dt
      end if

      do i=2, n
         Ops = Op(i)
         opso= Ops%dt
         if (cod == "tim") opso= Ops%time_inv
         j = i - 1
         do while (j >= 1)
            if (option(j) >= opso) exit
            Op(j + 1) = Op(j)
            option(j + 1) =  option(j)
            j = j - 1
         end do
         Op(j + 1) = Ops
         option(j+1) = opso
      end do
   End Subroutine Sort_Oper

End SubModule SPG_Sort_Oper

