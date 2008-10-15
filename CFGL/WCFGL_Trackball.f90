module WCFGL_trackball
!------------------------------------------------------------------------------
! Written by LCC April 2004.
! This module contains the subroutine to
! handle rotations via a Virtual
! Trackball system.
!------------------------------------------------------------------------------
  use OPENGL
  use WCFGL_constant, only : PI, deg2rad
  use WCFGL_geometry
  use WCFGL_quaternion

  implicit none
  
  real(kind=glfloat), dimension(:,:), pointer, public, save :: current_rot_matrix => null()
  real(kind=glfloat), dimension(:,:), pointer, public, save :: current_invrot_matrix => null()

  contains
!------------------------------------------------------------------------------
  subroutine init_rot_matrix()
  
    if (.not.(associated(current_rot_matrix))) allocate(current_rot_matrix(4,4))
    if (.not.(associated(current_invrot_matrix))) allocate(current_invrot_matrix(4,4))  
    current_rot_matrix=0.0
    current_rot_matrix(1,1)=1.0
    current_rot_matrix(2,2)=1.0
    current_rot_matrix(3,3)=1.0
    current_rot_matrix(4,4)=1.0
    
    current_invrot_matrix=0.0
    current_invrot_matrix(1,1)=1.0
    current_invrot_matrix(2,2)=1.0
    current_invrot_matrix(3,3)=1.0
    current_invrot_matrix(4,4)=1.0

    return
    
  end subroutine init_rot_matrix
!------------------------------------------------------------------------------
! Given a (x,y) point on the viewport,
! this returns a point projected on the hemisphere
! for mouse trackball system.
! Assume that the z-axis is up and we are using
! a right handed coordinate system
  subroutine viewport2hemi(xview,yview,width,height,v_out)
    integer,           intent(in)  :: xview, yview, width, height
    real,dimension(3), intent(out) :: v_out
    !
    real,dimension(2)          :: v_circle
    real                       :: v_circle_norm

    if (xview < 0 ) then
      v_circle(1)=-1.0
    else if (xview > width) then
      v_circle(1)=1.0
    else
      v_circle(1)=(2.0*xview-width)/width
    end if

   if (yview < 0) then
     v_circle(2)=1.0
   else if (yview > height) then
     v_circle(2)=-1.0
   else
     v_circle(2)=(height-2.0*yview)/height
   end if

   v_circle_norm=sqrt(sum(v_circle*v_circle))

   if (v_circle_norm > 1.0) then
     v_circle=v_circle/v_circle_norm
     v_out(1:2)=v_circle
     v_out(3)=0.0
   else
     v_out(1:2)=v_circle
     v_out(3)=sqrt(1.0-(v_circle_norm)**2)
   end if

   return

end subroutine viewport2hemi
!------------------------------------------------------------------------------
! Given two points on the hemisphere,
! this calculate the angle of rotation in degrees and
! axis of rotation (unitary vector) and create a quaternion
  subroutine calc_rot_param(v_in1, v_in2, q_out)
    real,dimension(3), intent(in) :: v_in1, v_in2
    type(quat),       intent(out) :: q_out
    !Local variables
    real, dimension(3)            :: v_rot
    real                          :: angle ! In deg
    real                          :: speed=1.25

    angle=90.0*speed*(.norm.(v_in1-v_in2))
    v_rot=.unit.(v_in1 .vect. v_in2)
    
    if (.norm.(v_in1-v_in2) <= 0.00000001) then ! Unity quaternion if the two vectors are very close.
      q_out=quat(1.0,(/0.0,0.0,0.0/))! This is to avoid underflow problems
      return
    endif
    
    q_out=quat(cos(angle*deg2rad/2.0),sin(angle*deg2rad/2.0)*v_rot)
    
    return

  end subroutine calc_rot_param
!------------------------------------------------------------------------------
end module WCFGL_trackball
