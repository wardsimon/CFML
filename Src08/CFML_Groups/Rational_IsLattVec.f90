!!----
!!----
!!----
SubModule (CFML_Groups) CFML_GRP_001
   Contains
   
   !!----
   !!---- IS_LATTICE_VEC
   !!----
   !!---- 19/04/19
   !!
   Module Function Is_Lattice_Vec(V,Ltr,Nlat) Result(Lattice)
      !---- Argument ----!
      type(rational), dimension(:),   intent( in) :: v
      type(rational), dimension(:,:), intent( in) :: Ltr
      integer,                        intent( in) :: nlat
      logical                                     :: Lattice

      !---- Local variables ----!
      type(rational), dimension(size(v)) :: vec
      integer                            :: i

      !> Init 
      Lattice=.false.

      !> if V is an integral vector then V is a lattice vector 
      if (Rational_Is_Integer(v)) then       
         Lattice=.true.
         return
      end if

         
      do i=1, nlat
         vec=Ltr(:,i)-v
         if (Rational_Is_Integer(vec)) then
            Lattice=.true.
            return
         end if
      end do
      
      return
   End Function Is_Lattice_Vec
    
End SubModule CFML_GRP_001  