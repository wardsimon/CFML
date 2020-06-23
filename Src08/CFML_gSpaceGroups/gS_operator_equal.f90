SubModule (CFML_gSpaceGroups) SPG_Oper_Equal
   implicit none
   Contains

   !!----
   !!---- EQUAL_SYMM_OPER
   !!----
   !!---- 19/04/2019
   !!
   Module Function Equal_Symm_Oper(Op1, Op2) Result(info)
      !---- Arguments ----!
      type(Symm_Oper_Type), intent(in) :: Op1,Op2
      logical                          :: info

      !> Init
      info=.false.

      if (Op1%time_inv == Op2%time_inv) then
         if (Rational_Equal(Op1%Mat,Op2%Mat)) info=.true.
      end if
   End Function Equal_Symm_Oper

   !!----
   !!---- EQUAL_GROUP
   !!----
   !!---- 19/04/2019
   !!
   Module Function Equal_Group(Gr1, Gr2) Result(info)
      !---- Arguments ----!
      class(Spg_Type), intent(in) :: Gr1
      class(Spg_Type), intent(in) :: Gr2
      logical                    :: info

      !---- Local Vatiables ----!
      integer :: i,j
      logical :: esta

      !> Init
      info=.false.

      if (Gr1%multip /= Gr2%multip) return
      do i=2,Gr1%multip
         esta=.false.
         do j=2,Gr2%multip
            if (Gr1%Op(i) == Gr2%Op(j)) then
               esta=.true.
               exit
            end if
         end do
         if (.not. esta) return
      end do
      info=.true.
   End Function Equal_Group

End SubModule SPG_Oper_Equal