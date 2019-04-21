!!----
!!----
!!----
SubModule (CFML_Groups) CFML_GRP_002
   Contains
   
   !!----
   !!---- REDUCED_TRANSLATION
   !!----
   !!---- 19/04/19
   !!
   Module Subroutine Reduced_Translation(Mat)
      !---- Arguments ----!
   	type(rational), dimension(:,:), intent(inout) :: Mat
   	
   	!---- Local Variables ----!
   	integer :: d,i,n,m
   	
   	!> Init
   	n=size(Mat,dim=1)
   	
   	d=n-1
      do i=1,d
         if (Rational_Is_Integer(Mat(i,n))) then
            Mat(i,n) = (0//1)
         else
            m=Mat(i,n)%numerator / Mat(i,n)%denominator
            if (Mat(i,n) > 0_LI) then
               Mat(i,n) = Mat(i,n) - (m//1)
            else
               Mat(i,n) = Mat(i,n) - (m//1) + (1//1)
            end if
         end if
      end do
     
      return
   End Subroutine Reduced_Translation
   
End SubModule CFML_GRP_002   