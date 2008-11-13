module WCFGL_glpoly
!------------------------------------------------------------------------------

  use WCFGL_geometry
  use WCFGL_glatom
  use WCFGL_glbond
  use WCFGL_atom_tree
  use WCFGL_bond_tree
  use WCFGL_chull3D
  use winteracter
  implicit none

  type ::  gl_poly
    real, dimension(4)              :: poly_color
    real, dimension(3)              :: edge_color
    real                            :: edge_thickness
    integer                         :: coord
    real, dimension(:,:),  pointer  :: x_center ! Contains positions of center atoms
    real, dimension(:,:,:), pointer :: xf       ! First dimension is xf_eq,Second is 3, third is number of connections
    integer, dimension(:),  pointer :: dlist    ! Store xf_eq display lists with the polyhedra
    logical                         :: show_edges
  end type gl_poly


  real, parameter, private                       :: eps=0.00001

  contains

!------------------------------------------------------------------------------
  subroutine new_gl_polys(label,poly_out,success,show_edges,ecolor,eradius,color)
    character(len=10), intent(in)                            :: label
    type(gl_poly)    , intent(out),dimension(:), allocatable :: poly_out
    logical          , intent(out)                           :: success
    logical          , intent(in)                            :: show_edges
    real             , intent(in)                            :: ecolor(3)
    real             , intent(in)                            :: eradius
    real             , intent(in), optional                  :: color(4)
    ! Local variables
    type(gl_atom), dimension(:), allocatable :: atomlist
    type(gl_bond), dimension(:), allocatable :: bondlist1, bondlist2, bondlisttemp, bondlist
    real, dimension(:,:), allocatable :: polybuffer
    integer :: i, j, k, l, nat, neq, nb, nbeq, count, countxf, newcon, aa
    logical :: success1, success2, successcon, success3

    success=  .false.
    success1= .false.
    success2= .false.
    success3= .false.

    call search_atom_list_by_label(label,atomlist,success1)

    if (.not.(success1)) return

    call search_bond_list_by_label(label,bondlist1,success2)

    if (success2) then
    allocate(bondlisttemp(size(bondlist1)))
    bondlisttemp=bondlist1
    end if

    newcon=0

    do aa=1,size(atomlist)
      if (allocated(bondlist2)) deallocate(bondlist2)
      call search_bond_list_by_label(atomlist(aa)%symbol//"        ",bondlist2,success3)
      if (success3) then
        newcon=newcon+1

        successcon=.true.
        if (newcon .eq. 1 .and. success2) then
          allocate(bondlist(size(bondlist1)+size(bondlist2)))
          bondlist(1:size(bondlist1))=bondlisttemp
          bondlist(size(bondlist1)+1:size(bondlist1)+size(bondlist2))=bondlist2
        else if (newcon .eq. 1 .and. .not.(success2)) then
          allocate(bondlist(size(bondlist2)))
          bondlist(1:size(bondlist2))=bondlist2
        else
        if (allocated(bondlisttemp)) deallocate(bondlisttemp)
          allocate(bondlisttemp(size(bondlist)))
          bondlisttemp=bondlist
          deallocate(bondlist)
          allocate(bondlist(size(bondlisttemp)+size(bondlist2)))
          bondlist(1:size(bondlisttemp))=bondlisttemp
          bondlist(size(bondlisttemp)+1:size(bondlisttemp)+size(bondlist2))=bondlist2
          deallocate(bondlisttemp)
        end if
      end if
    end do

    if (.not.(success2) .and. .not.(successcon) ) return

    if (allocated(bondlist1)) deallocate(bondlist1)
    if (allocated(bondlist2)) deallocate(bondlist2)
    if (allocated(bondlisttemp)) deallocate(bondlisttemp)

    allocate(poly_out(size(atomlist))) ! Allocate poly_out with the number of
    allocate(polybuffer(3,1000))       ! This is a buffer, first dimension is

    nat=size(atomlist)
    nb=size(bondlist)

    do i=1,nat
      neq=size(atomlist(i)%xf_eq, DIM=2)
      allocate(poly_out(i)%x_center(3,neq))
      allocate(poly_out(i)%dlist(neq))
      poly_out(i)%show_edges=show_edges
      poly_out(i)%edge_color=ecolor
      poly_out(i)%edge_thickness=eradius
      poly_out(i)%x_center=atomlist(i)%xf_eq
      do j=1, neq
        count=0
        do k=1, nb
          nbeq=size(bondlist(k)%xf, DIM=2)
          do l=1, nbeq
            if (.norm.(atomlist(i)%xf_eq(:,j)-bondlist(k)%xf(1:3,l))<eps) then
              count=count+1 ! Count the number of neighbors
              polybuffer(:,count)=bondlist(k)%xf(4:6,l)
            end if
          end do
        end do
      if (j==1) then
        allocate(poly_out(i)%xf(neq,3,count))
        poly_out(i)%coord=count
        if (present(color)) then
          poly_out(i)%poly_color=color
        else
          poly_out(i)%poly_color=atomlist(i)%color
        end if
      end if
      poly_out(i)%xf(j,:,:)=polybuffer(:,1:count)
      poly_out(i)%dlist(j)=glGenLists(1)
      call glNewList(poly_out(i)%dlist(j),GL_COMPILE)
      call construct_poly(f2c_array(poly_out(i)%xf(j,:,:),poly_out(i)%coord),poly_out(i)%poly_color,show_edges,ecolor,eradius)
      call glEndList()
      end do
    end do

    if (count>0) success=.true.

    if (allocated(polybuffer)) deallocate(polybuffer)
    if (allocated(atomlist)) deallocate(atomlist)
    if (allocated(bondlist)) deallocate(bondlist)
    if (allocated(bondlist2)) deallocate(bondlist2)
    return

  end subroutine new_gl_polys
!------------------------------------------------------------------------------
end module WCFGL_glpoly
