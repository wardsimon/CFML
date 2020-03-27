!!----
!!----
!!----
!!
SubModule (CFML_Reflections) RFL_001
   Contains
   !!----
   !!---- H_EQUAL
   !!----     .True. if both reflections are equal
   !!----
   !!---- 20/06/2019
   !!
   Module Function H_Equal(H,K) Result (Info)
      !---- Arguments ----!
      integer, dimension(:), intent(in) :: h
      integer, dimension(:), intent(in) :: k
      logical                           :: info

      !---- Local Variables ----!
      integer :: i

      !> Init
      info=.true.
      do i=1,size(H)
         if (h(i) /= k(i)) then
            info=.false.
            return
         end if
      end do
   End Function H_Equal

End SubModule RFL_001