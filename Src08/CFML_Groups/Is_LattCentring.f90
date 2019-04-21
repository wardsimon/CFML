!!----
!!----
!!----
!!
SubModule (CFML_Groups) CFML_GRP_011
   Contains
   
   !!----
   !!---- IS_LATTICE_CENTRING
   !!----
   !!---- 19/04/19
   !!
   Module Function Is_Lattice_Centring(Op) Result(Info)
      !---- Arguments ----!
      type(Symm_Oper_Type), intent(in) :: Op
      logical                          :: info
      
      !---- Local Variables ----!
      integer :: dr
      
      !> Init
      info=.false.
      if (.not. allocated(identity_matrix)) return
      
      dr=size(Op%Mat(:,1))-1
      if (Rational_Equal(identity_matrix(1:dr,1:dr),Op%Mat(1:dr,1:dr)) .and. Op%time_inv == 1) info=.true.
      
      return
   End Function Is_Lattice_Centring

End SubModule CFML_GRP_011   
   
