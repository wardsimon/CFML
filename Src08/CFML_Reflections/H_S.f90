!!----
!!----
!!----
SubModule (CFML_Reflections) RFL_005
   Contains
   !!----
   !!---- H_S
   !!----   Calculates: sin_theta/lamda = 1/(2d)
   !!----
   !!---- 20/06/2019 
   !!
   Module Function H_S(H, Cell, Nk, Kv) Result(S)
      !---- Arguments ----!
      integer, dimension(:),                  intent(in) :: h
      class(Cell_G_Type),                     intent(in) :: Cell
      integer, optional,                      intent(in) :: Nk
      real(kind=cp),dimension(:,:), optional, intent(in) :: Kv
      real(kind=cp)                                      :: S
      
      !--- Local variables ---!
      integer                     :: i
      real(kind=cp), dimension(3) :: hkl

      hkl=h(1:3)
      if (present(nk) .and. present(kv)) then
         do i=1,nk
            hkl=hkl+h(3+i)*kv(:,i)
         end do
      end if
      S= 0.5*sqrt( hkl(1)*hkl(1)*Cell%GR(1,1) +     hkl(2)*hkl(2)*Cell%GR(2,2) + &
                   hkl(3)*hkl(3)*Cell%GR(3,3) + 2.0*hkl(1)*hkl(2)*Cell%GR(1,2) + &
               2.0*hkl(1)*hkl(3)*Cell%GR(1,3) + 2.0*hkl(2)*hkl(3)*Cell%GR(2,3) )
   End Function H_S
End SubModule RFL_005 