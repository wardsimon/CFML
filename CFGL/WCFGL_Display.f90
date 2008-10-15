module WCFGL_display

  use OPENGL
  use WCFGL_metrix
  use WCFGL_trackball
  use WCFGL_atom_tree,  only : draw_atoms
  use WCFGL_matom_tree, only : draw_matoms
  use WCFGL_bond_tree,  only : draw_bonds
  use WCFGL_poly_tree,  only : draw_polys
  use winteracter

  implicit none

  real(kind=glfloat),               public, pointer  :: current_zoomfactor    => null()
  integer,            dimension(:), public, pointer  :: current_centroffset   => null()
  integer                         , public, pointer  :: current_screen_width  => null(), &
                                                        current_screen_height => null()

  real, dimension(4), public, save :: background_color=(/0.9,0.9,0.9,1.0/),&
                                      lamp_pos =(/0.0,0.0,50.0,0.0/),&
                                      lamp_diff=(/0.9,0.9,0.9,1.0/),&
                                      mat_specular=(/1.0,1.0,1.0,1.0/)

  integer, public, save            :: mat_shininess=80

  contains
!-----------------------------------------------------------------------------
  subroutine display_init()

    call glReadBuffer(GL_BACK)
    call glDrawBuffer(GL_BACK)

    call glEnable (GL_LIGHTING)
    call glEnable(GL_LIGHT0)

    if (.not.(associated(current_centroffset))) then
      allocate(current_centroffset(2))
      current_centroffset= 0
    end if

    if (.not.(associated(current_zoomfactor))) then
      allocate(current_zoomfactor)
      current_zoomfactor = 1.0
    end if

    if (.not.(associated(current_maxsize))) then
      allocate(current_maxsize)
      current_maxsize    = 0.0
    end if

    if (.not.(associated(current_baricentre))) then
      allocate(current_baricentre(3))
      current_baricentre = 0.0
    end if

    if (.not.(associated(current_screen_width))) allocate(current_screen_width)
    if (.not.(associated(current_screen_height))) allocate(current_screen_height)

    call init_rot_matrix()

    call glenable(GL_DEPTH_TEST)
    call glShadeModel(GL_SMOOTH)
    call glEnable (GL_LIGHTING)
    call glEnable(GL_LIGHT0)
    call glenable(GL_LINE_SMOOTH)
    call glenable(GL_BLEND)
    call glenable(GL_NORMALIZE)
    call glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)

    call wgltextfont(IFAMILY=FFHelvetica,ZEXTR=0.10)


    return


  end subroutine display_init
!-----------------------------------------------------------------------------
  subroutine display()

    call glClearColor(background_color(1),background_color(2),background_color(3),background_color(4))

    call glmaterialfv(GL_FRONT, GL_SPECULAR, mat_specular )
    call glMateriali(GL_FRONT, GL_SHININESS, mat_shininess)

    call glviewport(current_centroffset(1),current_centroffset(2)+(min(current_screen_width,current_screen_height)-&
                       max(current_screen_width,current_screen_height))/2,&
                       max(current_screen_width,current_screen_height),&
                       max(current_screen_width,current_screen_height))

    call glmatrixmode(GL_PROJECTION)
    call glloadidentity()
    call glortho(-current_zoomfactor*current_maxsize,&
                  current_zoomfactor*current_maxsize,&
                 -current_zoomfactor*current_maxsize,&
                  current_zoomfactor*current_maxsize,&
                 -2*current_zoomfactor*current_maxsize,&
                  2*current_zoomfactor*current_maxsize)
    call glmatrixmode(GL_MODELVIEW)
    call glloadidentity()

    call glclear(GL_COLOR_BUFFER_BIT+GL_DEPTH_BUFFER_BIT)
    call gllightfv(GL_LIGHT0, GL_AMBIENT, (/0.05,0.05,0.05,1.0/))
    call gllightfv(GL_LIGHT0, GL_DIFFUSE, lamp_diff )
    call gllightfv(GL_LIGHT0, GL_POSITION, lamp_pos)

    call gltranslatef(-current_baricentre(1),-current_baricentre(2),-current_baricentre(3))

    call glpushmatrix()
      call gltranslatef(current_baricentre(1),current_baricentre(2),current_baricentre(3))
      call glmultmatrixf(current_rot_matrix)
      call gltranslatef(-current_baricentre(1),-current_baricentre(2),-current_baricentre(3))

      call draw_cell()
      call draw_axes()
      call draw_atoms()
      call draw_bonds()
      call draw_matoms()
      call draw_polys()
    call glpopmatrix()

    call Wglswapbuffers()

    new_box=        .false.
    new_cell=       .false.
    new_spacegroup =.false.

    return

  end subroutine display
!-----------------------------------------------------------------------------
  subroutine calculate_offset(x1,y1,x2,y2)
    integer, intent(in) :: x1,y1,x2,y2

    current_centroffset(1)=current_centroffset(1)+x2-x1
    current_centroffset(2)=current_centroffset(2)+y1-y2

    return
  end subroutine calculate_offset
!-----------------------------------------------------------------------------
end module WCFGL_display
