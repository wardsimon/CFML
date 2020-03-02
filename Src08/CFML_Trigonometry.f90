  !> @brief
  !! Module CFML_Trigonometry
  !! This module implements elemental trigonometric functions working with
  !! angular variable given in degrees or resulting values in degrees
  !!
  !! The major part of the functions are already intrinsically implemented in Intel Fortran but
  !! not in GFortran. Only the function Rtan is not implemented in Ifort.
  !!
  !! @author
  !!
  Module CFML_Trigonometry

    use CFML_GlobalDeps, only: sp, dp, PI, TO_RAD, to_DEG
    use CFML_Maths,      only: epss


    implicit none

    Private

    public:: Acosd, Asind, Atan2d, Atand, Cosd, Sind, Tand, Rtan

    Interface  Acosd
       Module Procedure Acosd_dp
       Module Procedure Acosd_sp
    End Interface

    Interface  Asind
       Module Procedure Asind_dp
       Module Procedure Asind_sp
    End Interface

    Interface  Atan2d
       Module Procedure Atan2d_dp
       Module Procedure Atan2d_sp
    End Interface

    Interface  Atand
       Module Procedure Atand_dp
       Module Procedure Atand_sp
    End Interface

    Interface  Cosd
       Module Procedure Cosd_dp
       Module Procedure Cosd_sp
    End Interface

    Interface  Sind
       Module Procedure Sind_dp
       Module Procedure Sind_sp
    End Interface

    Interface  Rtan
       Module Procedure Rtan_dp
       Module Procedure Rtan_sp
    End Interface

    Interface  Tand
       Module Procedure Tand_dp
       Module Procedure Tand_sp
    End Interface


  contains
    !!----
    !!---- Elemental Function Acosd(x) Result(arc_cos)
    !!----    real(kind=sp/dp), intent(in) :: x
    !!----    real(kind=sp/dp)             :: arc_cos
    !!----
    !!----    Inverse cosine function -> output in Degrees
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Elemental Function Acosd_dp(x) Result(arc_cos)
    !!--++    real(kind=dp), intent(in) :: x
    !!--++    real(kind=dp)             :: arc_cos
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Inverse cosine function -> output in Degrees
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Acosd_dp(x) Result(arc_cos)
       !---- Argument ----!
       real(kind=dp), intent(in) :: x
       real(kind=dp)             :: arc_cos

       if (abs(x) > 1.0_dp ) then
          if (x > 0.0_dp)  then
             arc_cos=0.0_dp
          else
             arc_cos=180.0_dp
          end if
       else
          arc_cos=acos(x)*to_DEG
       end if

       return
    End Function Acosd_dp

    !!--++
    !!--++ Elemental Function Acosd_sp(x) Result(arc_cos)
    !!--++    real(kind=sp), intent(in) :: x
    !!--++    real(kind=sp)             :: arc_cos
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Inverse cosine function -> output in Degrees
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Acosd_sp(x) Result(arc_cos)
       !---- Argument ----!
       real(kind=sp), intent(in) :: x
       real(kind=sp)             :: arc_cos

       if (abs(x) > 1.0_sp ) then
          if (x > 0.0_sp)  then
             arc_cos=0.0_sp
          else
             arc_cos=180.0_sp
          end if
       else
          arc_cos=acos(x)*to_DEG
       end if

       return
    End Function Acosd_sp

    !!----
    !!---- Function Asind(x) Result(arc_sin)
    !!----    real(kind=sp/dp), intent(in) :: x
    !!----    real(kind=sp/dp)             :: arc_sin
    !!----
    !!----    Inverse sine function -> output in Degrees
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Elemental Function Asind_dp(x) result(arc_sin)
    !!--++    real(kind=dp), intent(in) :: x
    !!--++    real(kind=dp)             :: arc_sin
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Inverse sine function -> output in Degrees
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Asind_dp(x) Result(arc_sin)
       !---- Argument ----!
       real(kind=dp), intent(in) :: x
       real(kind=dp)             :: arc_sin

       if (abs(x) > 1.0_dp ) then
          if (x > 0.0_dp) then
             arc_sin=90.0_dp
          else
             arc_sin=-90.0_dp
          end if
       else
          arc_sin=asin(x)*to_DEG
       end if

       return
    End Function Asind_dp

    !!--++
    !!--++ Elemental Function Asind_sp(x) result(arc_sin)
    !!--++    real(kind=sp), intent(in) :: x
    !!--++    real(kind=sp)             :: arc_sin
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Inverse sine function -> output in Degrees
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Asind_sp(x) Result(arc_sin)
       !---- Argument ----!
       real(kind=sp), intent(in) :: x
       real(kind=sp)             :: arc_sin

       if (abs(x) > 1.0_sp ) then
          if (x > 0.0_sp) then
             arc_sin=90.0_sp
          else
             arc_sin=-90.0_sp
          end if
       else
          arc_sin=asin(x)*to_DEG
       end if

       return
    End Function Asind_sp

    !!----
    !!---- Elemental Function Atan2d(y,x) Result(atande)
    !!----    real(kind=sp/dp), intent(in) :: y,x
    !!----    real(kind=sp/dp)             :: atande
    !!----
    !!----    Inverse tangent function of y/x
    !!----    y,x have the same units -> output in Degrees
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Elemental Function Atan2d_dp(y,x) Result(atande)
    !!--++    real(kind=dp), intent(in) :: y,x
    !!--++    real(kind=dp)             :: atande
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Inverse tangent function of y/x
    !!--++    y,x have the same units -> output in Degrees
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Atan2d_dp(y,x) Result(atand)
       !---- Argument ----!
       real(kind=dp), intent(in) :: y,x
       real(kind=dp)             :: atand

       atand=atan2(y,x)*to_DEG

       return
    End Function Atan2d_dp

    !!--++
    !!--++ Elemental Function Atan2d_sp(y,x) Result(atande)
    !!--++    real(kind=sp), intent(in) :: y,x
    !!--++    real(kind=sp)             :: atande
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Inverse tangent function of y/x
    !!--++    y,x have the same units -> output in Degrees
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Atan2d_sp(y,x) Result(atande)
       !---- Argument ----!
       real(kind=sp), intent(in) :: y,x
       real(kind=sp)             :: atande

       atande=atan2(y,x)*to_DEG

       return
    End Function Atan2d_sp

    !!----
    !!---- Elemental Function Atand(x) Result(atande)
    !!----    real(kind=sp/dp), intent(in) :: x
    !!----    real(kind=sp/dp)             :: atande
    !!----
    !!----    Inverse tangent function, X no units -> output in Degrees
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Elemental Function Atand_dp(x) result(atande)
    !!--++    real(kind=dp), intent(in) :: x
    !!--++    real(kind=dp)             :: atande
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Inverse tangent function, X no units -> output in Degrees
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Atand_dp(x) Result(atand)
       !---- Argument ----!
       real(kind=dp), intent(in) :: x
       real(kind=dp)             :: atand

       atand=atan(x)*to_DEG

       return
    End Function Atand_dp

    !!--++
    !!--++ Function Atand_sp(x) result(atande)
    !!--++    real(kind=sp), intent(in) :: x
    !!--++    real(kind=sp)             :: atande
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Inverse tangent function, X no units -> output in Degrees
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Atand_sp(x) Result(atande)
       !---- Argument ----!
       real(kind=sp), intent(in) :: x
       real(kind=sp)             :: atande

       atande=atan(x)*to_DEG

       return
    End Function Atand_sp



    !!----
    !!---- Elemental Function Cosd(x) Result(cosine)
    !!----    real(kind=sp/dp), intent(in) :: x
    !!----    real(kind=sp/dp)             :: cosine
    !!----
    !!----    Cosine function, X in degrees
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Elemental Function Cosd_dp(x) Result(cosine)
    !!--++    real(kind=dp), intent(in) :: x
    !!--++    real(kind=dp)             :: cosine
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Cosine function, X in degrees
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Cosd_dp(x) Result(cosine)
       !---- Argument ----!
       real(kind=dp), intent(in) :: x
       real(kind=dp)             :: cosine

       cosine=cos(to_RAD*x)

       return
    End Function Cosd_dp

    !!--++
    !!--++ Elemental Function Cosd_sp(x) Result(cosine)
    !!--++    real(kind=sp), intent(in) :: x
    !!--++    real(kind=sp)             :: cosine
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Cosine function, X in degrees
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Cosd_sp(x) Result(cosine)
       !---- Argument ----!
       real(kind=sp), intent(in) :: x
       real(kind=sp)             :: cosine

       cosine=cos(to_RAD*x)

       return
    End Function Cosd_sp

    !!----
    !!---- Elemental Function Sind(x) Result(sine)
    !!----    real(kind=sp/dp), intent(in) :: x
    !!----    real(kind=sp/dp)             :: sine
    !!----
    !!----    Sine function, X in degrees
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Elemental Function Sind_dp(x) Result(sine)
    !!--++    real(kind=dp), intent(in) :: x
    !!--++    real(kind=dp)             :: sine
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Sine function, X in degrees
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Sind_dp(x) Result(sine)
       !---- Argument ----!
       real(kind=dp), intent(in) :: x
       real(kind=dp)             :: sine

       sine=sin(to_RAD*x)

       return
    End Function Sind_dp

    !!--++
    !!--++ Elemental Function Sind_sp(x) Result(sine)
    !!--++    real(kind=sp), intent(in) :: x
    !!--++    real(kind=sp)             :: sine
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Sine function, X in degrees
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Sind_sp(x) Result(sine)
       !---- Argument ----!
       real(kind=sp), intent(in) :: x
       real(kind=sp)             :: sine

       sine=sin(to_RAD*x)

       return
    End Function Sind_sp

    !!----
    !!---- Subroutine Rtan(y,x,ang,deg)
    !!----    real(sp/dp),               intent( in) :: x,y
    !!----    real(sp/dp),               intent(out) :: ang
    !!----    character(len=*),optional, intent( in) :: deg
    !!----
    !!----    Returns ang=arctan(y/x) in the quadrant where the signs sin(ang) and
    !!----    cos(ang) are those of y and x. If deg is present, return ang in degrees.
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Subroutine Rtan_dp(y,x,ang,deg)
    !!--++    real(dp),                  intent( in) :: x,y
    !!--++    real(dp),                  intent(out) :: ang
    !!--++    character(len=*),optional, intent( in) :: deg
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Returns ang=arctan(y/x) in the quadrant where the signs sin(ang) and
    !!--++    cos(ang) are those of y and x. If deg is present, return ang in degrees.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Rtan_dp(y,x,ang,deg)
       !---- Arguments ----!
       real(kind=dp),              Intent( In)   :: x,y
       real(kind=dp),              Intent(Out)   :: ang
       character(len=*), optional, Intent( In)   :: deg

       !---- Local variables ----!
       real(kind=dp):: abx,aby

       abx=abs(x)
       aby=abs(y)
       if ((abx < epss) .and. (aby < epss)) then
          ang = 0.0_dp
          return
       else if(abx < epss) then
          ang = pi/2.0_dp
       else if(aby < abx) then
          ang = atan(aby/abx)
          if(x < 0.0_dp) ang = pi-ang
       else
          ang = pi/2.0_dp - atan(abx/aby)
          if(x < 0.0_dp) ang = pi-ang
       end if
       if (y < 0.0_dp) ang = -ang
       if (present(deg)) ang = ang*to_deg

       return
    End Subroutine Rtan_dp

    !!--++
    !!--++ Subroutine Rtan_sp(x,y,ang,deg)
    !!--++    real(sp),                  intent( in) :: x,y
    !!--++    real(sp),                  intent(out) :: ang
    !!--++    character(len=*),optional, intent( in) :: deg
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Returns ang=arctan(y/x) in the quadrant where the signs sin(ang) and
    !!--++    cos(ang) are those of y and x. If deg is present, return ang in degrees.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Rtan_sp(y,x,ang,deg)
       !---- Arguments ----!
       real(kind=sp),              Intent( In)   :: x,y
       real(kind=sp),              Intent(Out)   :: ang
       character(len=*), optional, Intent( In)   :: deg

       !---- local variables ----!
       real(kind=sp):: abx,aby

       abx=abs(x)
       aby=abs(y)
       if ((abx < epss) .and. (aby < epss)) then
          ang = 0.0_sp
          return
       else if(abx < epss) then
          ang = pi/2.0_sp
       else if(aby < abx) then
          ang = atan(aby/abx)
          if(x < 0.0_sp) ang = pi-ang
       else
          ang = pi/2.0_sp - atan(abx/aby)
          if(x < 0.0_sp) ang = pi-ang
       end if
       if(y < 0.0_sp) ang = -ang
       if (present(deg)) ang = ang*to_deg

       return
    End Subroutine Rtan_sp

    !!----
    !!---- Elemental Function Tand(x) Result(tande)
    !!----    real(kind=sp/dp), intent(in) :: x
    !!----    real(kind=sp/dp)             :: tande
    !!----
    !!----    Tangent function, X in degrees
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Elemental Function Tand_dp(x) Result(tande)
    !!--++    real(kind=dp), intent(in) :: x
    !!--++    real(kind=dp)             :: tande
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Tangent function, X in degrees
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Tand_dp(x) Result(tand)
       !---- Argument ----!
       real(kind=dp), intent(in) :: x
       real(kind=dp)             :: tand

       tand=tan(to_RAD*x)

       return
    End Function Tand_dp

    !!--++
    !!--++ Elemental Function Tand_sp(x) Result(tande)
    !!--++    real(kind=sp), intent(in) :: x
    !!--++    real(kind=sp)             :: tande
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Tangent function, X in degrees
    !!--++
    !!--++ Update: February - 2005
    !!
    Elemental Function Tand_sp(x) Result(tande)
       !---- Argument ----!
       real(kind=sp), intent(in) :: x
       real(kind=sp)             :: tande

       tande=tan(to_RAD*x)

       return
    End Function Tand_sp


  End Module CFML_Trigonometry