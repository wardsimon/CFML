!!----
!!----
!!----
!!
SubModule (CFML_gSpaceGroups) SPG_008
   Contains

   !!----
   !!---- Allocate_SpaceGroup
   !!----
   !!---- 19/04/19
   !!
   Module Subroutine Allocate_SpaceGroup(D, Multip, Grp)
      !---- Arguments ----!
      integer,         intent(in)     :: d,multip
      class(Spg_Type), intent(in out) :: Grp

      !---- Local Variables ----!
      integer :: i

      !> Init
      if (allocated(Grp%Op)) deallocate(Grp%Op)
      allocate(Grp%Op(multip))

      Grp%d=d
      Grp%multip=multip
      do i=1,multip
         call Allocate_Op(d,Grp%Op(i))
      end do

      if(allocated(Grp%Symb_Op)) deallocate(Grp%Symb_Op)
      allocate(Grp%Symb_Op(multip))

      if(allocated(Grp%centre_coord)) deallocate(Grp%centre_coord)
      allocate(Grp%centre_coord(d-1))
      Grp%centre_coord=0_LI

      if(allocated(Grp%anticentre_coord)) deallocate(Grp%anticentre_coord)
      allocate(Grp%anticentre_coord(d-1))
      Grp%centre_coord=0_LI

   End Subroutine Allocate_SpaceGroup

   Module Subroutine Allocate_KVector(nk, nq, Kvec)
      !---- Arguments ----!
      integer,               intent(in)    :: nk  ! Number of independent k-vectors
      integer,               intent(in)    :: nq  !number of effective set of Q_coeff > nk
      type(Kvect_Info_Type), intent(inout) :: Kvec

      !> Init
      if (nk == 0) then
         Kvec%nk=0
         Kvec%nq=0
         if (allocated(Kvec%kv)) deallocate(kvec%kv)
         if (allocated(Kvec%sintlim)) deallocate(kvec%sintlim)
         if (allocated(Kvec%nharm)) deallocate(kvec%nharm)
         if (allocated(Kvec%Q_Coeff)) deallocate(kvec%Q_Coeff)
         return
      end if

      Kvec%nk=nk
      allocate(Kvec%kv(3,nk))
      allocate(Kvec%sintlim(nk))
      allocate(Kvec%nharm(nk))
      Kvec%kv=0.0_cp
      Kvec%sintlim=0.0_cp
      Kvec%nharm=0

      Kvec%nq=0
      if (nq < nk) return

      Kvec%nq=nq
      allocate(Kvec%q_coeff(nk,nq))
      Kvec%q_coeff=0

   End Subroutine Allocate_KVector

End SubModule SPG_008

