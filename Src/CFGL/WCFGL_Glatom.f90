module WCFGL_glatom
!------------------------------------------------
! Written by Laurent C.Chapon
! July 2004.
! Updated :
!------------------------------------------------
  use CFML_Math_General,              only : modulo_lat
  use CFML_Crystallographic_Symmetry, only : space_group_type, ApplySO

  implicit none

  real, parameter, private   :: eps=0.0002

!-------------------------------------------------------------------------------
  type :: gl_atom
    character(len=10)                            :: label
    character(len=2)                             :: symbol         ! Chemical Symbol.
    real                                         :: radius         ! Radius.
    real                                         :: Biso           ! Biso.
    real                                         :: color(4)       ! Color.
    real                                         :: xf(3)          ! Pos in crystallographic unit cell.
    integer                                      :: multip         ! Multiplicity of the position
    logical                                      :: colorfromtable ! If color from table is requested
    real         , dimension(:,:), pointer       :: xf_eq          ! Equivalent positions in crystallographic unit cell.
  end type gl_atom
!-------------------------------------------------------------------------------
  !If molecule=.true. the symmetry operators of the space group are not applied as well as
  !the modulo lattice translations
  logical, save   :: molecule

  contains

!-------------------------------------------------------------------------------
  subroutine new_gl_atom(label, symbol, xf, Biso, radius, color, SpaceG, atom_out)
    character(len=10),      intent(in)                 :: label
    character(len=2),       intent(in)                 :: symbol
    real         ,          intent(in)                 :: xf(3)
    real         ,          intent(in)                 :: Biso
    real         ,          intent(in)                 :: radius
    real         ,          intent(in)                 :: color(4)
    type(space_group_type), intent(in)                 :: SpaceG
    type(gl_atom),          intent(out)                :: atom_out

      atom_out%xf_eq => null()

      atom_out%label = label

      atom_out%symbol=symbol

      atom_out%radius=radius

      atom_out%color=color

      if(molecule) then
        atom_out%xf=xf
      else
        atom_out%xf=modulo_lat(xf)      ! Make sure the position is in the unit cell.
      end if

      atom_out%Biso=Biso

      call calculate_eq_positions(atom_out,SpaceG)

     return

  end subroutine new_gl_atom
!-------------------------------------------------------------------------------
  subroutine calculate_eq_positions(atom,spaceG)
    type(space_group_type), intent(in)      :: spaceG
    type(gl_atom),          intent(inout)   :: atom
    ! Local variables---------------------
    real, dimension(:,:), allocatable :: xtemp
    real                              :: xx(3), v(3)
    integer                           :: L, j, nt

    if (associated(atom%xf_eq)) deallocate(atom%xf_eq)

    atom%xf_eq => null()

    allocate(xtemp(3,spaceG%multip))
    L=1
    xtemp(:,1)=atom%xf
    if(.not. molecule) then
      do_eq:do j=2, spaceG%multip
        xx=ApplySO(spaceG%Symop(j),xtemp(:,1))
        xx=modulo_lat(xx)

        do nt=1,L
          v=xtemp(:,nt)-xx(:)
          if (sum(abs((v))) < eps ) cycle do_eq
        end do

        L=L+1
        xtemp(:,L)=xx(:)
      end do do_eq
    end if
    allocate(atom%xf_eq(3,L))

    atom%xf_eq(:,:)=xtemp(:,1:L)
    atom%multip=L

    if (allocated(xtemp)) deallocate(xtemp)

    return

  end subroutine calculate_eq_positions
!-------------------------------------------------------------------------------
end module WCFGL_glatom