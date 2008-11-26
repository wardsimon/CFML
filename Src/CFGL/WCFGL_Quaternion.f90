module WCFGL_quaternion
!------------------------------------------------
! Written by Laurent C.Chapon
! June 2004.
! Updated :
!------------------------------------------------
  use WCFGL_constant
  use WCFGL_geometry
  use CFML_Math_General, only: sind,cosd

  implicit none

  public :: operator(*)                !Operators
  public :: q_from_rot, q_rot_vect     !Functions
  public :: quat2mat,rot_from_q        !Subroutines
  private:: q_mult_q

  type,public :: quat ! Define a quaternion as a real part and a 3D vector
    real               :: s
    real, dimension(3) :: v
  end type

  interface operator (*)
    module procedure q_mult_q
  end interface

  real, parameter, private :: eps=0.000001

  contains
!------------------------------------------------------------------------------
  function q_mult_q(q1,q2) result(q_out) ! Multiply two quaternions
    type(quat), intent(in) :: q1, q2
    type(quat) :: q_out

    q_out%s=(q1%s*q2%s)-(q1%v .scal. q2%v)
    q_out%v=(q1%v .vect. q2%v)+ (q1%s*q2%v)+(q2%s*q1%v)

    return

  end function q_mult_q
!------------------------------------------------------------------------------
!This function constructs the quaternion corresponding to a rotation given
!by the angle (in degrees) and the rotation axis (it may be provided non-normalized!)
  function q_from_rot(ang,axis) result(q)
    real,             intent(in) :: ang
    real,dimension(3),intent(in) :: axis
    type(quat)                   :: q
    !--- Local variables ---!
    real              :: norm

    norm = .norm. axis
    if(norm >= 0.00000001) then
      q=quat(cosd(0.5*ang), axis*sind(0.5*ang)/norm)
    else
      q=quat(1.0,(/0.0,0.0,0.0/))
    end if
    return

  end function q_from_rot
!---------------------------------------------------------------------------------
!This function applies a rotation (given as a quaternion) to a vector
  function q_rot_vect(q,v) result(u)
    type(quat),       intent(in) :: q
    real,dimension(3),intent(in) :: v
    real,dimension(3)            :: u
    !--- Local variables ---!
    type(quat) :: qv,qc, qu

    qv = quat(0.0,   v)
    qc = quat(q%s,-q%v)
    qu = (q  * qv) * qc
    u  = qu%v
    return

  end function q_rot_vect
!------------------------------------------------------------------------------
!Gives the rotation angle in degrees and the rotation axis
!The axis vector is normalized if the rotation angle is greater than zero
  subroutine rot_from_q(q,ang,axis)
    type(quat),       intent( in) :: q
    real,             intent(out) :: ang
    real,dimension(3),intent(out) :: axis

    !--- Local variables ---!
    real              :: s

    ang=2.0*acosd(q%s)
    s=sind(ang*0.5)
    if (abs(s) <= 1.0e-10) then
      axis=q%v
    else
      axis=q%v/s
    end if
    return
  end subroutine rot_from_q
!------------------------------------------------------------------------------
  subroutine quat2mat(q_in,mat_out) ! Quaternion to rotation matrix
    type(quat), intent(in)  :: q_in
    real      , intent(out) :: mat_out(4,4)

    mat_out(1,1)=1.0 - 2.0*(q_in%v(2))**2 - 2.0*(q_in%v(3))**2
    mat_out(2,1)=2.0*q_in%v(1)*q_in%v(2) + 2.0*q_in%v(3)*q_in%s
    mat_out(3,1)=2.0*q_in%v(1)*q_in%v(3) - 2.0*q_in%v(2)*q_in%s
    mat_out(4,1)=0.0
    mat_out(1,2)=2.0*q_in%v(1)*q_in%v(2) - 2.0*q_in%v(3)*q_in%s
    mat_out(2,2)=1.0 - 2.0*(q_in%v(1))**2 - 2.0*(q_in%v(3))**2
    mat_out(3,2)=2.0*q_in%v(2)*q_in%v(3) + 2.0*q_in%v(1)*q_in%s
    mat_out(4,2)=0.0
    mat_out(1,3)=2.0*q_in%v(1)*q_in%v(3) + 2.0*q_in%v(2)*q_in%s
    mat_out(2,3)=2.0*q_in%v(2)*q_in%v(3) - 2.0*q_in%v(1)*q_in%s
    mat_out(3,3)=1.0 - 2.0*(q_in%v(1))**2 - 2.0*(q_in%v(2))**2
    mat_out(4,3)=0.0
    mat_out(1,4)=0.0
    mat_out(2,4)=0.0
    mat_out(3,4)=0.0
    mat_out(4,4)=1.0

    return

  end subroutine quat2mat
!------------------------------------------------------------------------------
end module WCFGL_quaternion
