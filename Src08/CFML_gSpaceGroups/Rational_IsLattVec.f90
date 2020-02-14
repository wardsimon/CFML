!!----
!!----
!!----
SubModule (CFML_gSpaceGroups) SPG_003
   Contains

   !!----
   !!---- IS_LATTICE_VEC
   !!----
   !!---- 19/04/19
   !!
   Module Function Is_Lattice_Vec_rat(V,Ltr) Result(Lattice)
      !---- Argument ----!
      type(rational), dimension(:),   intent( in) :: v
      type(rational), dimension(:,:), intent( in) :: Ltr
      logical                                     :: Lattice

      !---- Local variables ----!
      type(rational), dimension(size(v)) :: vec
      integer                            :: i, NLat

      !> Init
      Lattice=.false.
      nlat=size(Ltr,dim=2)
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
   End Function Is_Lattice_Vec_rat

   Module Function is_Lattice_vec_real(v,Ltr) Result(Lattice)
      !---- Arguments ----!
      real(kind=cp), dimension(:),   intent(in) :: v
      type(rational),dimension(:,:), intent(in) :: Ltr
      logical                                   :: Lattice
      !---- Local Variables ----!
      integer :: i,nlat
      real(kind=cp), dimension(size(v)) :: vec,vl
      !> Init
      Lattice_Transl=.false.
      nlat=size(Ltr,dim=2)
       if (Zbelong(v)) then       ! if v is an integral vector =>  v is a lattice vector
          Lattice=.true.
       else                       ! if not look for lattice type
          do i=1,nlat
            vl=Ltr(:,i)
            vec=vl-v
            if (Zbelong(vec)) then
              Lattice=.true.
              exit
            end if
          end do
       end if
   End Function is_Lattice_vec_real


End SubModule SPG_003