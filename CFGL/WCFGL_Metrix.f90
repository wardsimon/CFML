module WCFGL_metrix
!------------------------------------------------
! Written by Laurent C.Chapon
! July 2004.
! Updated History :
! 8/07/2004 Problem with current_cell creator fixed
!
!------------------------------------------------
  use OPENGL
  use WCFGL_constant,             only : deg2rad
  use WCFGL_geometry
  use crystallographic_symmetry, only : space_group_type, set_spacegroup, Err_symm
  use WCFGL_objects_definition
    
  implicit none

!------------------------------------------------------------------------------
  type :: gl_axes
    real,dimension(3) :: pos
    real              :: lenght
    logical           :: dead
    integer           :: dlist
  end type
!------------------------------------------------------------------------------
  type :: gl_cell
    real,dimension(3)                :: a,angle
    real                             :: volume
    real, dimension(3,3)             :: t_mtx
    real, dimension(3,3)             :: t_mtx_unit
    real, dimension(4)               :: color
    real                             :: width
    logical                          :: dead
    logical                          :: multiple
    integer                          :: dlist
  end type
!------------------------------------------------------------------------------

  type(gl_cell), pointer, public,  save    :: current_cell => null()
  type(gl_axes), pointer, public,  save    :: current_axes => null()

  type(space_group_type), pointer, public, save  :: current_space_group => null()
  real,    dimension(:),  pointer, public, save  :: current_box         => null()

  logical, public, save                          :: new_cell=.false., new_box=.false., &
                                                    new_spacegroup=.false. 
                                              
  real(kind=glfloat)   ,  dimension(:), pointer, public, save :: current_baricentre => null() 

  real(kind=gldouble)                 , pointer, public ,save :: current_maxsize    => null()
  interface f2c
    module procedure f2c_r
    module procedure f2c_i
  end interface

  contains
!------------------------------------------------------------------------------
  subroutine define_box(xmin,xmax,ymin,ymax,zmin,zmax)
    real, intent(in) :: xmin, xmax, ymin, ymax, zmin,zmax

    if (.not.(associated(current_box))) then
      allocate(current_box(6))
      current_box=(/xmin, xmax, ymin, ymax, zmin,zmax/)
    else
      current_box=(/xmin, xmax, ymin, ymax, zmin,zmax/)
    end if
    
    new_box=.true.
    
    if (.not.(associated(current_baricentre))) allocate(current_baricentre(3))
    current_baricentre=f2c((current_box(2:6:2)+current_box(1:5:2))/2.0)

    if (.not.(associated(current_maxsize))) allocate(current_maxsize)
      
    if (associated(current_cell)) then
      current_maxsize=maxval(current_cell%a*(current_box(2:6:2)-current_box(1:5:2)))
    else
      current_maxsize=maxval((current_box(2:6:2)-current_box(1:5:2)))
    end if


    return

  end subroutine define_box
!------------------------------------------------------------------------------
  subroutine init_box()

    if (.not.(associated(current_box))) allocate(current_box(6))

    current_box(1:5:2)=0.0
    current_box(2:6:2)=1.0
    
    if (.not.(associated(current_baricentre))) allocate(current_baricentre(3))
    current_baricentre=f2c((current_box(2:6:2)+current_box(1:5:2))/2.0)
    
    if (.not.(associated(current_maxsize))) allocate(current_maxsize)
      
    if (associated(current_cell)) then
      current_maxsize=maxval(current_cell%a*(current_box(2:6:2)-current_box(1:5:2)))
    else
      current_maxsize=maxval((current_box(2:6:2)-current_box(1:5:2)))
    end if
    
    return

  end subroutine init_box
!------------------------------------------------------------------------------
  subroutine define_cell(a,angle)
    real, intent(in) :: a(3), angle(3)
    ! Local variables---------------------
	real, dimension(3) :: c100, c010, c001, c110, c101, c011, c111, cspher

    if (.not.(associated(current_box))) call init_box()

	if (angle(1)<=0 .or. angle(2)<=0 .or. angle(3)<=0 &
        .or. a(1)<=0 .or. a(2)<=0 .or. a(3)<=0) return

    if(.not.(associated(current_cell))) then   ! New cell is created, load parameters
      allocate(current_cell)                   ! per default and create a display list.
      current_cell%dead=.false.
      current_cell%multiple=.false.
      current_cell%color=(/0.5,0.5,0.5,1.0/)
      current_cell%width=1.45
      current_cell%dlist = glgenlists(1)
      allocate(current_axes)                   ! Same think for axes.
      current_axes%pos   = (/-0.3,-0.3,-0.3/)
      current_axes%lenght= 2.0
      current_axes%dead  =.false.
      current_axes%dlist = glgenlists(1)
    end if

    current_cell%a(1)=a(1)
    current_cell%a(2)=a(2)
    current_cell%a(3)=a(3)

    current_cell%angle(1)=angle(1)
    current_cell%angle(2)=angle(2)
    current_cell%angle(3)=angle(3)

    current_cell%volume=a(1)*a(2)*a(3)* &
    sqrt(1.0-(cos(deg2rad*angle(1)))**2- &
    (cos(deg2rad*angle(2)))**2- &
    (cos(deg2rad*angle(3)))**2+ &
    2.0*cos(deg2rad*angle(1))*cos(deg2rad*angle(2))*cos(deg2rad*angle(3)))

    current_cell%t_mtx(1,1)=a(1)
    current_cell%t_mtx(2,1)=0.0
    current_cell%t_mtx(3,1)=0.0
    current_cell%t_mtx(1,2)=a(2)*cos(deg2rad*angle(3))
    current_cell%t_mtx(2,2)=a(2)*sin(deg2rad*angle(3))
    current_cell%t_mtx(3,2)=0.0
    current_cell%t_mtx(1,3)=a(3)*cos(deg2rad*angle(2))
    current_cell%t_mtx(2,3)=a(3)/sin(deg2rad*angle(3))*(cos(deg2rad*angle(1))-cos(deg2rad*angle(2))*cos(deg2rad*angle(3)))
    current_cell%t_mtx(3,3)=current_cell%volume/(a(1)*a(2)*sin(deg2rad*angle(3)))

    current_cell%t_mtx_unit(1,1)=1.0
    current_cell%t_mtx_unit(2,1)=0.0
    current_cell%t_mtx_unit(3,1)=0.0
    current_cell%t_mtx_unit(1,2)=cos(deg2rad*angle(3))
    current_cell%t_mtx_unit(2,2)=sin(deg2rad*angle(3))
    current_cell%t_mtx_unit(3,2)=0.0
    current_cell%t_mtx_unit(1,3)=cos(deg2rad*angle(2))
    current_cell%t_mtx_unit(2,3)=1.0/sin(deg2rad*angle(3))*(cos(deg2rad*angle(1))-cos(deg2rad*angle(2))*cos(deg2rad*angle(3)))
    current_cell%t_mtx_unit(3,3)=current_cell%volume/(a(1)*a(2)*a(3)*sin(deg2rad*angle(3)))

    c100=f2c(p100)
    c010=f2c(p010)
    c001=f2c(p001)
    c110=f2c(p110)
    c101=f2c(p101)
    c011=f2c(p011)
    c111=f2c(p111)

    call glnewlist(current_cell%dlist, GL_COMPILE) ! This is the object cell to plot
       call glbegin(GL_LINES)
          call glvertex3fv(p000)
          call glvertex3fv(c100)
		  call glvertex3fv(p000)
		  call glvertex3fv(c010)
		  call glvertex3fv(p000)
		  call glvertex3fv(c001)
		  call glvertex3fv(c100)
		  call glvertex3fv(c110)
		  call glvertex3fv(c100)
		  call glvertex3fv(c101)
		  call glvertex3fv(c010)
		  call glvertex3fv(c110)
		  call glvertex3fv(c010)
		  call glvertex3fv(c011)
		  call glvertex3fv(c001)
	      call glvertex3fv(c101)
		  call glvertex3fv(c001)
		  call glvertex3fv(c011)
		  call glvertex3fv(c110)
		  call glvertex3fv(c111)
		  call glvertex3fv(c101)
		  call glvertex3fv(c111)
		  call glvertex3fv(c011)
		  call glvertex3fv(c111)
	    call glend()
	call glendlist()

    c100=current_axes%lenght*(p100)
    c010=current_axes%lenght*(p010)
    c001=current_axes%lenght*(p001)
    
   if (.not.(associated(axes_definition))) call init_axes()

    call glnewlist(current_axes%dlist,GL_COMPILE) ! This is the object axes to plot
        call glpushmatrix()
        cspher=cart2spher(matmul(current_cell%t_mtx_unit,c100))
        call glrotatef(cspher(2),0.0,0.0,1.0)
        call glrotatef(-cspher(3),0.0,1.0,0.0)
        call glscalef(cspher(1),0.7,0.7)
        call gltranslatef(0.5,0.0,0.0)
        call glmaterialfv(GL_FRONT,GL_DIFFUSE,(/1.0,0.0,0.0,1.0/))
        call glcalllist(axes_definition)
        call glpopmatrix()
        call glpushmatrix()
        call gltranslatef(c100(1),c100(2),c100(3))
        call glcolor3f(1.0,0.0,0.0)
        call glscalef(1.5,1.5,1.5)
        call wgltextstring('a')
        call glpopmatrix()
        call glpushmatrix()
        cspher=cart2spher(matmul(current_cell%t_mtx_unit,c010))
        call glrotatef(cspher(2),0.0,0.0,1.0)
        call glrotatef(-cspher(3),0.0,1.0,0.0)
        call glscalef(cspher(1),0.7,0.7)
        call gltranslatef(0.5,0.0,0.0)
        call glmaterialfv(GL_FRONT,GL_DIFFUSE,(/0.0,1.0,0.0,1.0/))
        call glcalllist(axes_definition)
        call glpopmatrix()
        call glpushmatrix()
        call gltranslatef(c010(1),c010(2),c010(3))
        call glcolor3f(0.0,1.0,0.0)
        call glrotatef(cspher(2),0.0,0.0,1.0)
        call glrotatef(-cspher(3),0.0,1.0,0.0)
        call glscalef(1.5,1.5,1.5)
        call wgltextstring('b')
        call glpopmatrix()
        call glpushmatrix()
        cspher=cart2spher(matmul(current_cell%t_mtx_unit,c001))
        call glrotatef(cspher(2),0.0,0.0,1.0)
        call glrotatef(-cspher(3),0.0,1.0,0.0)
        call glscalef(cspher(1),0.7,0.7)
        call gltranslatef(0.5,0.0,0.0)
        call glmaterialfv(GL_FRONT,GL_DIFFUSE,(/0.0,0.0,1.0,1.0/))
        call glcalllist(axes_definition)
        call glpopmatrix()
        call glpushmatrix()
        call gltranslatef(c001(1),c001(2),c001(3))
        call glcolor3f(0.0,0.0,1.0)
        call glrotatef(cspher(2),0.0,0.0,1.0)
        call glrotatef(-cspher(3),0.0,1.0,0.0)
        call glscalef(1.5,1.5,1.5)
        call wgltextstring('c')
        call glpopmatrix()
    call glendlist()

    new_cell=.true.
    
    if (.not.(associated(current_baricentre))) allocate(current_baricentre(3))
    current_baricentre=f2c((current_box(2:6:2)+current_box(1:5:2))/2.0)
    
    if (.not.(associated(current_maxsize))) allocate(current_maxsize)
    current_maxsize=maxval(current_cell%a*(current_box(2:6:2)-current_box(1:5:2)))

    return

  end subroutine define_cell
!------------------------------------------------------------------------------
  subroutine init_spacegroup()

    if (.not.(associated(current_space_group))) allocate(current_space_group)

    call set_spacegroup("P 1",current_space_group)

    new_spacegroup=.true.

    return

  end subroutine init_spacegroup
!------------------------------------------------------------------------------
  subroutine define_spacegroup(symbol,genr,ng)
    character(len=*), intent(in) :: symbol
    character(len=*), dimension(:),optional,intent(in) :: genr
    integer,optional, intent(in) :: ng
    

    if (.not.(associated(current_space_group))) allocate(current_space_group)
      

    if(present(genr) .and. present(ng)) then
        call set_spacegroup(symbol,current_space_group,gen=genr,ngen=ng,mode="GEN")
    else
        call set_spacegroup(symbol,current_space_group)
    end if
    
    if (Err_Symm) call init_spacegroup()
       
    new_spacegroup=.true.

    return

  end subroutine define_spacegroup
!------------------------------------------------------------------------------
  function f2c_r(v_in) result(v_out)
    real(glfloat),dimension(3), intent(in) :: v_in
    real(glfloat),dimension(3)             :: v_out

    if(associated(current_cell)) then
      v_out=matmul(current_cell%t_mtx,v_in)
    else
      v_out=v_in
    end if

    return

  end function f2c_r
  !------------------------------------------------------------------------------
  function f2c_array(a_in,n) result(a_out) 
    real(glfloat),dimension(:,:), intent(in)  :: a_in
    integer                     , intent(in)  :: n
    real(glfloat),dimension(3,n)              :: a_out
    integer :: j 
    
    if(associated(current_cell)) then
      do j=1, n
      a_out(:,j)=matmul(current_cell%t_mtx,a_in(:,j))
      end do
    else
      a_out=a_in
    end if

    return

  end function f2c_array
  !------------------------------------------------------------------------------
  function f2c_i(v_in) result(v_out)
    integer,dimension(3), intent(in) :: v_in
    real,dimension(3)                :: v_out
    real,dimension(3)                :: temp

    temp=real(v_in)

    if(associated(current_cell)) then
      v_out=matmul(current_cell%t_mtx,temp)
    else
      v_out=temp
    end if

    return

  end function f2c_i
!------------------------------------------------------------------------------
  subroutine draw_cell()
  integer           :: i,j,k
  real,dimension(3) :: trans

  	if (.not.(associated(current_cell))) return

    if (current_cell%dead) return

    call gldisable(GL_LIGHTING)
    call gllinewidth(current_cell%width)
    call glcolor4fv(current_cell%color)

    if (current_cell%multiple) then
      do i=int(current_box(1)), int(current_box(2))-1
        do j=int(current_box(3)), int(current_box(4))-1
          do k=int(current_box(5)), int(current_box(6))-1
            trans=f2c((/i,j,k/))
            call glpushmatrix()
            call gltranslatef(trans(1),trans(2),trans(3))
	        call glcalllist(current_cell%dlist)
            call glpopmatrix()
          end do
        end do
      end do
    else
      trans=f2c((/current_box(1),current_box(2),current_box(3)/))
      call glcalllist(current_cell%dlist)
    end if
    
    call glenable(GL_LIGHTING)

    return

  end subroutine draw_cell
!------------------------------------------------------------------------------
  subroutine draw_axes()
 ! Local variables---------------------
  real,dimension(3)  :: cpos

  	if (.not.(associated(current_axes))) return

    if (current_axes%dead) return

    cpos=f2c(current_axes%pos)
    call glpushmatrix()
    call gltranslatef(cpos(1),cpos(2),cpos(3))
    call glcalllist(current_axes%dlist)
    call glpopmatrix()
    
    return

  end subroutine draw_axes
!------------------------------------------------------------------------------
  function in_limit(vector) result(isinlimit)
  real,dimension(3), intent(in) :: vector
  logical          :: isinlimit

  if (.not.(associated(current_box))) call init_box()

  if (vector(1)>=current_box(1).and. &
      vector(1)<=current_box(2).and. &
      vector(2)>=current_box(3).and. &
      vector(2)<=current_box(4).and. &
      vector(3)>=current_box(5).and. &
      vector(3)<=current_box(6)) then
    isinlimit=.true.
  else
    isinlimit=.false.
  end if

  return

  end function in_limit
!------------------------------------------------------------------------------
  function intersect_volume(point,vector) result(points_i)
  real,dimension(3), intent(in) :: point, vector
  real,dimension(3,6)           :: points_i

  if (.not.(associated(current_box))) call init_box()

     points_i=0.0

     if (vector(1) /= 0.0) then
       points_i(2:3,1)=point(2:3)
       points_i(2:3,2)=point(2:3)
       points_i(1,1)=current_box(1)
       points_i(1,2)=current_box(2)
     end if

     if (vector(2) /= 0.0) then
       points_i(1:3:2,3)=point(1:3:2)
       points_i(1:3:2,4)=point(1:3:2)
       points_i(2,3)=current_box(3)
       points_i(2,4)=current_box(4)
     end if

     if (vector(3) /= 0.0) then
       points_i(1:2,5)=point(1:2)
       points_i(1:2,6)=point(1:2)
       points_i(3,5)=current_box(5)
       points_i(3,6)=current_box(6)
     end if

  return

  end function intersect_volume
!------------------------------------------------------------------------------
end module WCFGL_metrix