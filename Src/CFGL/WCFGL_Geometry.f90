 module WCFGL_geometry
   !------------------------------------------------------------------------------
   ! Written by LCC March 2004
   ! Updated : July 2004. Simplified and add transformations such as spheric to cartesian...
   !------------------------------------------------------------------------------
     use WCFGL_constant, only : rad2deg

     implicit none

     private

     public  :: operator(.unit.), operator (.norm.), operator (.scal.), operator (.vect.)
     public  :: cart2spher
     private :: vect_unit_3D,vect_norm_3D, vect3D_scal_vect3D, vect3D_vect_vect3D
   !------------------------------------------------------------------------------
     real, dimension(3), public , parameter ::  p000=(/0.0,0.0,0.0/), &
                                                p100=(/1.0,0.0,0.0/), &
                                                p010=(/0.0,1.0,0.0/), &
                                                p001=(/0.0,0.0,1.0/), &
                                                p110=(/1.0,1.0,0.0/), &
                                                p101=(/1.0,0.0,1.0/), &
                                                p011=(/0.0,1.0,1.0/), &
                                                p111=(/1.0,1.0,1.0/)
   !------------------------------------------------------------------------------
       interface operator (.unit.)
         module procedure vect_unit_3D
       end interface

       interface operator (.norm.)
         module procedure vect_norm_3D
       end interface

       interface operator (.scal.)
          module procedure vect3D_scal_vect3D
       end interface

       interface operator (.vect.)
          module procedure vect3D_vect_vect3D
       end interface

       contains
     !------------------------------------------------------------------------------
       function vect_unit_3D(v_in) result(v_out) ! Normalise a vector
         real,dimension(3), intent(in) :: v_in
         real,dimension(3)             :: v_out

         v_out=v_in/(.norm.v_in)

         return

       end function vect_unit_3D
     !------------------------------------------------------------------------------
       function vect_norm_3D(v_in) result(v_norm) ! Return the norm of a vector
         real,dimension(3), intent(in) :: v_in
         real                          :: v_norm

         v_norm=sqrt(dot_product(v_in,v_in))
         return

       end function vect_norm_3D
     !------------------------------------------------------------------------------
     function vect3D_scal_vect3D(v_in1,v_in2) result(scalpro) ! Return scalar product
         real,dimension(3), intent(in) :: v_in1, v_in2
         real             :: scalpro

        scalpro=dot_product(v_in1,v_in2)

        return

     end function vect3D_scal_vect3D
     !------------------------------------------------------------------------------
     function vect3D_vect_vect3D(v_in1,v_in2) result(v_out) ! Return vectorial product
         real,dimension(3), intent(in) :: v_in1, v_in2
         real,dimension(3)             :: v_out

        v_out(1)=v_in1(2)*v_in2(3)-v_in1(3)*v_in2(2)
        v_out(2)=v_in1(3)*v_in2(1)-v_in1(1)*v_in2(3)
        v_out(3)=v_in1(1)*v_in2(2)-v_in1(2)*v_in2(1)

        return

     end function vect3D_vect_vect3D
     !------------------------------------------------------------------------------
     function cart2spher(v_in) result(v_out)
       real,dimension(3) , intent(in) :: v_in  !(x,y,z)
       real,dimension(3)              :: v_out !(r,phi,theta) angles in radians

       v_out(1)=.norm.(v_in)

       if (v_in(1) == 0.0 .and. v_in(2) == 0.0) then
        v_out(2) = 0.0
      else
        v_out(2) = atan2(v_in(2),v_in(1))
      end if

      if (v_out(1) == 0.0) then
        v_out(3) = 0.0
      else
        v_out(3) = acos(v_in(3)/v_out(1))
      end if

      v_out(2:3)=rad2deg*(v_out(2:3))

      v_out(3)=90.0-v_out(3)

      return

    end function cart2spher
    !------------------------------------------------------------------------------
  end module WCFGL_geometry
