!!----
!!----
!!----
SubModule (CFML_Reflections) RFL_012
   Contains
   !!----
   !!---- H_Equiv_List
   !!----    Calculate the multiplicity of the reflection and the list of all
   !!----    equivalent reflections. Friedel law assumed if Friedel=.true.
   !!----
   !!---- 24/06/2019
   !!
   Module Subroutine H_Equiv_List(H, SpG, Friedel, Mult, H_List)
      !---- Arguments ----!
      integer, dimension(:),    intent (in) :: H
      class(SpG_Type),          intent (in) :: SpG
      Logical,                  intent (in) :: Friedel
      integer,                  intent(out) :: mult
      integer, dimension(:,:),  intent(out) :: h_list

      !---- Local Variables ----!
      integer, dimension(size(h))  :: k
      integer                      :: i,j,ng
      logical                      :: esta
      integer,dimension(3,3)       :: Op

      !> Init
      h_list = 0
      mult=1
      h_list(:,1)=h(:)
      
      ng=SpG%multip
      do i=2,ng
         Op=SpG%Op(i)%Mat(1:3,1:3)
         k = matmul(h,Op)
         
         esta=.false.
         do j=1,mult
            if (h_equal(k,h_list(:,j)) .or. (h_equal(-k,h_list(:,j)) .and. Friedel)) then
               esta=.true.
               exit
            end if
         end do
         if (esta) cycle
         
         mult=mult+1
         h_list(:,mult) = k
      end do

      if (Friedel .or. SpG%centred == 2) then
         j=mult
         mult=mult*2
         do i=j+1,mult
            h_list(:,i)=-h_list(:,i-j)
         end do
      end if

   End Subroutine H_Equiv_List
   
End SubModule RFL_012   