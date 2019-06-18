!!----
SubModule (CFML_gSpaceGroups) Spg_102
   Contains
   !!----
   !!---- GET_STABILIZER
   !!----    Subroutine to obtain the list of symmetry operator of a space group that leaves
   !!----    invariant an atomic position. This subroutine provides a pointer to the symmetry
   !!----    operators of the site point group and the additional translation with respect to
   !!----    the canonical representant.
   !!----
   !!---- 13/06/2019
   !!
   Module Subroutine Get_Stabilizer(X, Spg,Order,Ptr,Atr)
      !---- Arguments ----!
      real(kind=cp), dimension(3),  intent (in)  :: x     ! real space position (fractional coordinates)
      class(Spg_Type),              intent (in)  :: Spg   ! Space group
      integer,                      intent(out)  :: order ! Number of sym.op. keeping invariant the position x
      integer, dimension(:),        intent(out)  :: ptr   ! Array pointing to the symmetry operators numbers
                                                          ! of the stabilizer of x
      real(kind=cp), dimension(:,:),intent(out)  :: atr   ! Associated additional translation to the symmetry operator
      
      !---- Local variables ----!
      real(kind=cp), dimension(3)    :: xx, tr
      integer                        :: j,n1,n2,n3

      !> Init
      order = 1              ! Identity belongs always to the stabilizer
      ptr   = 0; ptr(1)= 1
      atr   = 0.0_cp

      do n1=-1,1
         do n2=-1,1
            do n3=-1,1
               tr=real([n1, n2, n3])
               do j=2,Spg%multip
                  xx=Apply_OP(Spg%Op(j),x) + tr - x
                  if (sum(abs(xx)) > 2.0 * EPSS) cycle
                  order=order+1
                  ptr(order)=j
                  atr(:,order)=tr
               end do
            end do
         end do
      end do

   End Subroutine Get_Stabilizer
    
End SubModule Spg_102   