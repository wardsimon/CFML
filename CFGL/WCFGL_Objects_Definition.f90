module WCFGL_objects_definition
!------------------------------------------------
! Written by Laurent C.Chapon
! June 2004.
! Updated :
!------------------------------------------------
  use OPENGL

  implicit none

  type(GLUquadricObj), pointer   ,  private            :: object
  integer,                          public, save       :: highres_moment_subdiv(2)=(/20,1/), &
                                                          highres_bond_subdiv(2)=(/20,1/), &
                                                          highres_atom_subdiv=20, moment_subdiv(2)=(/20,1/), &
                                                          bond_subdiv(2)=(/20,1/), atom_subdiv=20

  real(gldouble),                   public ,save       :: moment_thickness=0.07, &
                                                          cone_radius=0.27, &
                                                          bond_thickness=0.050, &
                                                          atom_scale=0.3, moment_length=0.7
  real,                             public, save       :: arrow_pos = -0.5

  integer            , public,  pointer :: atom_definition   => null(), &
                                           moment_definition => null(), &
                                           bond_definition   => null(), &
                                           axes_definition   => null()
  logical,             private, save    :: highres=.true.

  contains
!-----------------------------------------------------------------------------
  subroutine change_moment_length(value)
   integer, intent(in) :: value


   if (value>95 .or. value<5) return

   moment_length=(100-value)/100.0

   call init_moment()

  end subroutine change_moment_length
!-----------------------------------------------------------------------------
  subroutine change_moment_width(value)
   integer, intent(in) :: value

   if (value<5.or.(value/100.0<moment_thickness)) return

   cone_radius=value/100.0

   call init_moment()

  end subroutine change_moment_width
!-----------------------------------------------------------------------------
  subroutine change_moment_thickness(value)
   integer, intent(in) :: value

   if (value<5) return

   moment_thickness=value/100.0

   call init_moment()

  end subroutine change_moment_thickness
!-----------------------------------------------------------------------------
  subroutine change_res_objects(quality)
    integer, intent(in) :: quality

    if (quality<6 .or. quality>60) return

    highres_moment_subdiv(1) = quality
    highres_bond_subdiv(1)   = quality
    highres_atom_subdiv      = quality

    call switch_highres_objects()

    return
  end subroutine change_res_objects
!-----------------------------------------------------------------------------
  subroutine switch_lowres_objects()
    highres=.false.
    moment_subdiv=(/6,1/)
    bond_subdiv=(/2,1/)
    atom_subdiv=6
    call init_atom()
    call init_moment()
    call init_bond()
    return
  end subroutine switch_lowres_objects
!-----------------------------------------------------------------------------
  subroutine switch_highres_objects()
    highres=.true.
    moment_subdiv=highres_moment_subdiv
    bond_subdiv  =highres_bond_subdiv
    atom_subdiv  =highres_atom_subdiv
    call init_atom()
    call init_moment()
    call init_bond()
    return
  end subroutine switch_highres_objects
!-----------------------------------------------------------------------------
  subroutine init_atom()
    if(.not.(associated(atom_definition))) then
      allocate(atom_definition)
      atom_definition = glgenlists(1)
    end if

    allocate(object)

    object = glunewquadric()
    if (highres) then
      call gluquadricnormals(object,GLU_SMOOTH)
      call gluquadricdrawstyle(object,GLU_FILL)
    else
      call gluquadricnormals(object,GLU_SMOOTH)
      call gluquadricdrawstyle(object,GLU_FILL)
    end if
    call glnewlist(atom_definition,GL_COMPILE)
    call gluquadricorientation(object,GLU_OUTSIDE)
    call glusphere(object,atom_scale,atom_subdiv,atom_subdiv)
    call glendlist()

    deallocate(object)

    return

  end subroutine init_atom
!-----------------------------------------------------------------------------
  subroutine init_moment()
    if(.not.(associated(moment_definition))) then
      allocate(moment_definition)
      moment_definition = glgenlists(1)
    end if

    allocate(object)

    object = glunewquadric()

    if (highres) then
      call gluquadricnormals(object,GLU_SMOOTH)
      call gluquadricdrawstyle(object,GLU_FILL)
    else
      call gluquadricnormals(object,GLU_SMOOTH)
      call gluquadricdrawstyle(object,GLU_FILL)
    end if
    call glnewlist(moment_definition,GL_COMPILE)
    call glpushmatrix()
    !call gltranslatef(-0.5,0.0,0.0)
    call gltranslatef(arrow_pos,0.0,0.0)
    call glpushmatrix()
    call glrotatef(90.0,0.0,1.0,0.0)
    call gluquadricorientation(object,GLU_INSIDE)
    call gludisk(object,0.0_gldouble,moment_thickness,moment_subdiv(1),moment_subdiv(2))
    call gluquadricorientation(object,GLU_OUTSIDE)
    call glucylinder(object,moment_thickness,moment_thickness,&
                     moment_length,moment_subdiv(1),moment_subdiv(2))
    call gltranslatef(0.0,0.0,real(moment_length))
    call gluquadricorientation(object,GLU_INSIDE)
    call gludisk(object,0.0_gldouble,cone_radius,moment_subdiv(1),moment_subdiv(2))
    call gluquadricorientation(object,GLU_OUTSIDE)
    call glucylinder(object,cone_radius,0.005_gldouble,1.0_gldouble-moment_length,&
                     moment_subdiv(1),moment_subdiv(2))
    call glpopmatrix()
    call glpopmatrix()
    call glendlist()

    deallocate(object)

    return

  end subroutine init_moment
!-----------------------------------------------------------------------------
  subroutine init_bond()
    if(.not.(associated(bond_definition))) then
      allocate(bond_definition)
      bond_definition = glgenlists(1)
    end if

    allocate(object)

    object = glunewquadric()

    if (highres) then
      call gluquadricnormals(object,GLU_SMOOTH)
      call gluquadricdrawstyle(object,GLU_FILL)
    else
      call gluquadricnormals(object,GLU_SMOOTH)
      call gluquadricdrawstyle(object,GLU_FILL)
    end if
    call glnewlist(bond_definition,GL_COMPILE)
    call glpushmatrix()
    call glrotatef(90.0,0.0,1.0,0.0)
    call gluquadricorientation(object,GLU_INSIDE)
    call gludisk(object,0.0_gldouble,bond_thickness,bond_subdiv(1),bond_subdiv(2))
    call gluquadricorientation(object,GLU_OUTSIDE)
    call glucylinder(object,bond_thickness,bond_thickness,&
                     1.0_gldouble,bond_subdiv(1),bond_subdiv(2))
    call gltranslatef(0.0,0.0,1.0)
    call gluquadricorientation(object,GLU_INSIDE)
    call gludisk(object,0.0_gldouble,bond_thickness,bond_subdiv(1),bond_subdiv(2))
    call glpopmatrix()
    call glendlist()

    deallocate(object)

    return

  end subroutine init_bond
!-----------------------------------------------------------------------------
  subroutine init_axes()
    if(.not.(associated(axes_definition))) then
      allocate(axes_definition)
      axes_definition = glgenlists(1)
    end if

    allocate(object)

    object = glunewquadric()

    call gluquadricnormals(object,GLU_SMOOTH)
    call gluquadricdrawstyle(object,GLU_FILL)
    call glnewlist(axes_definition,GL_COMPILE)
    call glpushmatrix()
    call gltranslatef(-0.5,0.0,0.0)
    call glpushmatrix()
    call glrotatef(90.0,0.0,1.0,0.0)
    call gluquadricorientation(object,GLU_INSIDE)
    call gludisk(object,0.0_gldouble,0.05_gldouble,15,1)
    call gluquadricorientation(object,GLU_OUTSIDE)
    call glucylinder(object,0.07_gldouble,0.07_gldouble,0.7_gldouble,15,1)
    call gltranslatef(0.0,0.0,0.7)
    call gluquadricorientation(object,GLU_INSIDE)
    call gludisk(object,0.0_gldouble,0.27_gldouble,15,1)
    call gluquadricorientation(object,GLU_OUTSIDE)
    call glucylinder(object,0.27_gldouble,0.005_gldouble,0.3_gldouble,15,1)
    call glpopmatrix()
    call glpopmatrix()
    call glendlist()

    deallocate(object)

    return

  end subroutine init_axes
!-----------------------------------------------------------------------------
end module WCFGL_objects_definition