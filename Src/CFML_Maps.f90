!!-------------------------------------------------------
!!---- Crystallographic Fortran Modules Library (CrysFML)
!!-------------------------------------------------------
!!---- The CrysFML project is distributed under LGPL. In agreement with the
!!---- Intergovernmental Convention of the ILL, this software cannot be used
!!---- in military applications.
!!----
!!---- Copyright (C) 1999-2012  Institut Laue-Langevin (ILL), Grenoble, FRANCE
!!----                          Universidad de La Laguna (ULL), Tenerife, SPAIN
!!----                          Laboratoire Leon Brillouin(LLB), Saclay, FRANCE
!!----
!!---- Authors: Juan Rodriguez-Carvajal (ILL)
!!----          Javier Gonzalez-Platas  (ULL)
!!----
!!---- Contributors: Laurent Chapon     (ILL)
!!----               Marc Janoschek     (Los Alamos National Laboratory, USA)
!!----               Oksana Zaharko     (Paul Scherrer Institute, Switzerland)
!!----               Tierry Roisnel     (CDIFX,Rennes France)
!!----               Eric Pellegrini    (ILL)
!!----
!!---- This library is free software; you can redistribute it and/or
!!---- modify it under the terms of the GNU Lesser General Public
!!---- License as published by the Free Software Foundation; either
!!---- version 3.0 of the License, or (at your option) any later version.
!!----
!!---- This library is distributed in the hope that it will be useful,
!!---- but WITHOUT ANY WARRANTY; without even the implied warranty of
!!---- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!!---- Lesser General Public License for more details.
!!----
!!---- You should have received a copy of the GNU Lesser General Public
!!---- License along with this library; if not, see <http://www.gnu.org/licenses/>.
!!----
!!----
!!---- MODULE: CFML_Maps_Calculations
!!----   INFO: Subroutines related to operations on the array's map
!!----
!!---- HISTORY
!!----    Update: 07/03/2011
!!----
!!---- DEPENDENCIES
!!--++    Use CFML_GlobalDeps,                only: cp
!!--++    Use CFML_Crystallographic_Symmetry, only: Space_Group_Type
!!--++    Use CFML_Crystal_Metrics,           only: Crystal_Cell_Type
!!--++    Use CFML_Geometry_Calc,             only: Distance
!!----
!!---- VARIABLES
!!----    CUBE_INFO_TYPE
!!----    CUBE_INFO
!!----    ERR_MAPS
!!----    ERR_MAPS_MESS
!!----    MAX_POINTS
!!----
!!---- PROCEDURES
!!----    Functions:
!!----       INDEX_CUBE
!!--++       PEAK_POSITION                    [Private]
!!----       VERTICE_POINT
!!--++       VERTICE_POINT_CAL                [Overloaded]
!!--++       VERTICE_POINT_FIX                [Overloaded]
!!----       VERTICES_CUBE
!!----       VPOINT_IN_CUBE
!!----       VPOINT_IN_LINE
!!----       VPOINT_IN_SQUARE
!!--++       XY_SECT                          [Private]
!!----
!!----    Subroutines:
!!----       CALCULATE_CONTOUR2D
!!----       CALCULATE_MESH
!!----       INIT_ERR_MAPS
!!----       LOAD_EXTENDEDMAP
!!----       LOAD_SECTION
!!--++       PEAK_LIST           [Private]
!!----       SEARCH_PEAKS
!!----       SET_CUBE_INFO
!!----       STATISTIC_MAP
!!----
!!
 Module CFML_Maps_Calculations

    !---- Use Modules ----!
    use CFML_GlobalDeps,                only: cp, eps
    use CFML_Crystallographic_Symmetry, only: Space_Group_Type, ApplySO
    use CFML_Crystal_Metrics,           only: Crystal_Cell_Type
    use CFML_Geometry_Calc,             only: Distance

    implicit none

    private
    !---- List of public functions ----!
    public :: Index_Cube, Vertice_Point, Vertices_Cube, Vpoint_in_Cube, &
              Vpoint_in_Square, Vpoint_in_Line


    !---- List of public subroutines ----!
    public :: Init_Err_Maps, Load_ExtendedMap, Load_Section, Search_Peaks, Statistic_Map, &
              Calculate_Contour2D, Calculate_Mesh, Set_Cube_Info

    !---- List of private functions ----!
    private :: Peak_Position, xy_sect

    !---- List of private subroutines ----!
    private :: Peak_List, Vertice_point_cal, Vertice_point_fix


    !---- Definitions ----!

    !!----
    !!----  TYPE :: CUBE_INFO_TYPE
    !!--..
    !!----  Type :: Cube_Info_Type
    !!----     integer                :: NElem     ! Number of Elemens
    !!----     integer                :: Code      ! Code of Elements
    !!----                                           1:               4:
    !!----                                           2:               5:
    !!----                                           3:               6:
    !!----     integer, dimension(12) :: Edges     ! Code for Edge connections
    !!----  End Type Cube_Info_Type
    !!----
    !!----  Update: February - 2005
    !!
    Type, public :: Cube_Info_Type
       integer                :: NElem
       integer                :: Code
       integer, dimension(12) :: Edges
    End Type Cube_Info_Type

    !!----
    !!---- CUBE_INFO
    !!----     Type (Cube_Info_Type), dimension(0:255) :: Cube_Info
    !!----
    !!----     Information of Mesh in a cube
    !!----
    !!---- Update: February - 2005
    !!
    Type (Cube_Info_Type), dimension(0:255), public :: Cube_Info

    !!----
    !!---- ERR_MAPS
    !!----    logical, public  :: err_maps
    !!----
    !!----    Logical Variable indicating an error in CFML_Maps_Calculations module
    !!----
    !!---- Update: February - 2005
    !!
    logical, public  :: err_maps

    !!----
    !!---- ERR_MAPS_MESS
    !!----    character(len=150), public :: ERR_Maps_Mess
    !!----
    !!----    String containing information about the last error
    !!----
    !!---- Update: February - 2005
    !!
    character(len=150), public :: ERR_Maps_Mess

    !!----
    !!---- MAX_POINTS
    !!----    integer, parameter, public ::  Max_Points
    !!----
    !!----    Number of maximum points permitted
    !!----
    !!---- Update: February - 2005
    !!
    integer, public, parameter ::  Max_Points    = 500000


    !---- Interfaces - Overlapp ----!
    Interface  Vertice_Point
       Module Procedure Vertice_Point_Fix
       Module Procedure Vertice_Point_Cal
    End Interface

 Contains

    !-------------------!
    !---- Functions ----!
    !-------------------!

    !!----
    !!---- Function Index_Cube(Iv,Mc) Result (Value)
    !!----    integer, dimension(8),intent(in) :: iv   ! In -> Vertices state On/Off
    !!----    logical,              intent(in) :: mc   ! Mc = .true. Only triangles (128-255 Code)
    !!----                                                    .false.               (  0-127 Code)
    !!----
    !!----    Return the index for Marching cubes algorithm
    !!----
    !!---- Update: February - 2005
    !!
    Function Index_Cube(iv, mc) Result(value)
       !---- Arguments ----!
       integer, dimension(8), intent(in) :: iv
       logical,               intent(in) :: mc
       integer                           :: value

       !---- local variables ----!
       integer :: i

       i = iv(1) + iv(2)*2 + iv(3)*4 + iv(4)*8 + iv(5)*16 + &
           iv(6)*32 + iv(7)*64 + iv(8)*128

       value=i
       if (mc) then
          if (value <= 127) value=255-i
       else
          if (value >= 128) value=255-i
       end if

       return
    End Function Index_Cube

    !!--++
    !!--++ Function Peak_Position(Rho,i,j,k) Result(Pto)
    !!--++    integer,dimension(:,:,:), intent(in) :: Rho     ! Density map scaled as integer values
    !!--++    integer,                  intent(in) :: i
    !!--++    integer,                  intent(in) :: j       ! (i,j,k) is the central point
    !!--++    integer,                  intent(in) :: k
    !!--++
    !!--++    (Private)
    !!--++    Return the position of the peak
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Peak_Position(nr3d,i,j,k) Result(Pto)
       !---- Arguments ----!
       integer, dimension(:,:,:),intent(in) :: nr3d
       integer,                  intent(in) :: i
       integer,                  intent(in) :: j
       integer,                  intent(in) :: k
       real(kind=cp), dimension(3)          :: pto

       !---- Local variables ----!
       integer       :: ntx,nty,ntz
       integer       :: i1,i2,i3,j1,j2,j3,k1,k2,k3
       real(kind=cp) :: a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19
       real(kind=cp) :: b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,b19
       real(kind=cp) :: b,c,d,h,kk,l,e,f,g,det,deltax,deltay,deltaz,x,y,z,dx,dy,dz

       !---- Calculation of the peak position ----!
       ntx=size(nr3d,1)
       nty=size(nr3d,2)
       ntz=size(nr3d,3)
       dx=1.0/real(ntx)
       dy=1.0/real(nty)
       dz=1.0/real(ntz)
       pto=0.0

       i2=i
       i1=i-1
       i3=i+1
       if (i1 <= 0)  i1=ntx+i1
       if (i3 > ntx) i3=i3-ntx

       j2=j
       j1=j-1
       j3=j+1
       if (j1 <= 0) j1=nty+j1
       if (j3 > nty) j3=j3-nty

       k2=k
       k1=k-1
       k3=k+1
       if (k1 <= 0) k1=ntz+k1
       if (k3 > ntz) k3=k3-ntz

        a1=real( MAX( nr3d(i2,j2,k2),1 ) )
        a2=real( MAX( nr3d(i1,j2,k2),1 ) )
        a3=real( MAX( nr3d(i3,j2,k2),1 ) )
        a4=real( MAX( nr3d(i3,j1,k2),1 ) )
        a5=real( MAX( nr3d(i2,j1,k2),1 ) )
        a6=real( MAX( nr3d(i1,j1,k2),1 ) )
        a7=real( MAX( nr3d(i1,j3,k2),1 ) )
        a8=real( MAX( nr3d(i2,j3,k2),1 ) )
        a9=real( MAX( nr3d(i3,j3,k2),1 ) )
       a10=real( MAX( nr3d(i2,j2,k1),1 ) )
       a11=real( MAX( nr3d(i1,j2,k1),1 ) )
       a12=real( MAX( nr3d(i3,j2,k1),1 ) )
       a13=real( MAX( nr3d(i2,j1,k1),1 ) )
       a14=real( MAX( nr3d(i2,j3,k1),1 ) )
       a15=real( MAX( nr3d(i2,j2,k3),1 ) )
       a16=real( MAX( nr3d(i1,j2,k3),1 ) )
       a17=real( MAX( nr3d(i3,j2,k3),1 ) )
       a18=real( MAX( nr3d(i2,j1,k3),1 ) )
       a19=real( MAX( nr3d(i2,j3,k3),1 ) )

        b1=LOG(a1)
        b2=LOG(a2)
        b3=LOG(a3)
        b4=LOG(a4)
        b5=LOG(a5)
        b6=LOG(a6)
        b7=LOG(a7)
        b8=LOG(a8)
        b9=LOG(a9)
       b10=LOG(a10)
       b11=LOG(a11)
       b12=LOG(a12)
       b13=LOG(a13)
       b14=LOG(a14)
       b15=LOG(a15)
       b16=LOG(a16)
       b17=LOG(a17)
       b18=LOG(a18)
       b19=LOG(a19)

       b=(b3+b4+b9+b12+b17-b2-b7-b6-b11-b16)/10.0
       c=(b7+b8+b9+b14+b19-b4-b5-b6-b13-b18)/10.0
       d=(b15+b16+b17+b18+b19-b10-b11-b13-b12-b14)/10.0
       h=(b13+b19-b18-b14)/4.0
       kk=(b11+b17-b16-b12)/4.0
       l=(b6+b9-b4-b7)/4.0
       e=(-10*b1+b2+b3+5*(b4+b6+b7+b9+b11+b12+b16+b17) &
                                   -6*(b5+b8+b10+b15)-2*(b13+b14+b18+b19))/21.0
       f=(-10*b1+b5+b8-6*(b2+b3+b10+b15)+5*(b4+b6+b7+b9+b13+b14+b18+b19) &
                                                     -2*(b11+b12+b16+b17))/21.0
       g=(-10*b1+b10+b15-6*(b2+b3+b5+b8)+5*(b11+b12+b13+b14+b16+b17+b18+b19) &
                                                         -2*(b4+b6+b7+b9))/21.0

       det=e*f*g+2*h*kk*l-kk*kk*f-h*h*e-l*l*g

       deltax=(-b*g*f-h*l*d-h*c*kk+kk*f*d+h*h*b+c*l*g)/det
       if (abs(deltax)-1.0 <= 0.0) then
          deltay=(-e*c*g-b*h*kk-l*d*kk+kk*kk*c+d*h*e+b*l*g)/det
          if (abs(deltay)-1.0 <= 0.0) then
             deltaz=(-e*f*d-l*c*kk-b*l*h+b*f*kk+l*l*d+h*c*e)/det
             if (abs(deltaz)-1.0 <= 0.0) then
                deltax=deltax*dx
                deltay=deltay*dy
                deltaz=deltaz*dz
             else
                deltax=0.0
                deltay=0.0
                deltaz=0.0
             end if
          else
            deltax=0.0
             deltay=0.0
             deltaz=0.0
          end if
       else
          deltax=0.0
          deltay=0.0
          deltaz=0.0
       end if

       x=(i-1)*dx
       y=(j-1)*dy
       z=(k-1)*dz

       x=x+deltax
       y=y+deltay
       z=z+deltaz

       x=mod(x+10.0_cp,1.0_cp)
       y=mod(y+10.0_cp,1.0_cp)
       z=mod(z+10.0_cp,1.0_cp)

       if (abs(x) <= 0.001_cp) x=0.0_cp
       if (abs(y) <= 0.001_cp) y=0.0_cp
       if (abs(z) <= 0.001_cp) z=0.0_cp
       if (abs(1.0-x) <= 0.001_cp ) x=0.0_cp
       if (abs(1.0-y) <= 0.001_cp ) y=0.0_cp
       if (abs(1.0-z) <= 0.001_cp ) z=0.0_cp

       pto(1)=x
       pto(2)=y
       pto(3)=z

       return
    End Function Peak_Position

    !!----
    !!---- Function Vertice_Point(Code_Edge,D0,D1,D2,D3,D4,D5,D6,D7,D8,D9) Result(V)
    !!----     integer, intent(in) :: code_edge
    !!----     or
    !!----     integer, intent(in) :: code_edge
    !!----     real(kind=cp),    intent(in) :: d0, d1,d2,d3,d4,d5,d6,d7,d8,d9
    !!----     real(kind=cp), dimension(3)  :: v
    !!----
    !!----     Return the relative position point from (i,j,k) of V1
    !!----     Given a binary dataset, linear interpolation is not needed to
    !!----     extract isosurfaces, When a cell edge in a binary dataset has
    !!----     both on and off corners, the midpoint of the edge is the
    !!----     intersection being looked for.
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Function Vertice_Point_Cal(Code_Edge,D0,D1,D2,D3,D4,D5,D6,D7,D8,D9) Result(V)
    !!--++     integer, intent(in) :: code_edge                         !  In ->
    !!--++     real(kind=cp),    intent(in) :: d0, d1,d2,d3,d4,d5,d6,d7,d8,d9    !  In ->
    !!--++     real(kind=cp), dimension(3)  :: v                                 ! Out ->
    !!--++
    !!--++     (OVERLOADED)
    !!--++     Return the relative position point from (i,j,k) of V1
    !!--++     Given a binary dataset, linear interpolation is not needed to
    !!--++     extract isosurfaces, When a cell edge in a binary dataset has
    !!--++     both on and off corners, the midpoint of the edge is the
    !!--++     intersection being looked for.
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Vertice_Point_Cal(Code_Edge,D0,D1,D2,D3,D4,D5,D6,D7,D8) Result(V)
       !---- Arguments ----!
       integer,          intent(in) :: code_edge
       real(kind=cp),    intent(in) :: d0,d1,d2,d3,d4,d5,d6,d7,d8
       real(kind=cp), dimension(3)  :: v

       !---- Local Variables ----!
       real(kind=cp)               :: dmin,dmax,dd

       v=0.0
       select case (code_edge)
          case ( 1)
             v(1)=0.5
             dmin=min(d2,d1)
             dmax=max(d2,d1)
             dd=dmax-dmin
             if (abs(dd) > eps ) v(1)=(d0-dmin)/dd

          case ( 2)
             v(1)=1.0
             v(2)=0.5
             dmin=min(d3,d2)
             dmax=max(d3,d2)
             dd=dmax-dmin
             if (abs(dd) > eps ) v(2)=(d0-dmin)/dd

          case ( 3)
             v(1)=0.5
             v(2)=1.0
             dmin=min(d4,d3)
             dmax=max(d4,d3)
             dd=dmax-dmin
             if (abs(dd) > eps ) v(1)=(d0-dmin)/dd

          case ( 4)
             v(2)=0.5
             dmin=min(d4,d1)
             dmax=max(d4,d1)
             dd=dmax-dmin
             if (abs(dd) > eps ) v(2)=(d0-dmin)/dd

          case ( 5)
             v(1)=0.5
             v(3)=1.0
             dmin=min(d6,d5)
             dmax=max(d6,d5)
             dd=dmax-dmin
             if (abs(dd) > eps ) v(1)=(d0-dmin)/dd

          case ( 6)
             v(1)=1.0
             v(2)=0.5
             v(3)=1.0
             dmin=min(d7,d6)
             dmax=max(d7,d6)
             dd=dmax-dmin
             if (abs(dd) > eps ) v(2)=(d0-dmin)/dd

          case ( 7)
             v(1)=0.5
             v(2)=1.0
             v(3)=1.0
             dmin=min(d8,d7)
             dmax=max(d8,d7)
             dd=dmax-dmin
             if (abs(dd) > eps ) v(1)=(d0-dmin)/dd

          case ( 8)
             v(2)=0.5
             v(3)=1.0
             dmin=min(d8,d5)
             dmax=max(d8,d5)
             dd=dmax-dmin
             if (abs(dd) > eps ) v(2)=(d0-dmin)/dd

          case ( 9)
             v(3)=0.5
             dmin=min(d5,d1)
             dmax=max(d5,d1)
             dd=dmax-dmin
             if (abs(dd) > eps ) v(3)=(d0-dmin)/dd

          case (10)
             v(1)=1.0
             v(3)=0.5
             dmin=min(d6,d2)
             dmax=max(d6,d2)
             dd=dmax-dmin
             if (abs(dd) > eps ) v(3)=(d0-dmin)/dd

          case (11)
             v(2)=1.0
             v(3)=0.5
             dmin=min(d8,d4)
             dmax=max(d8,d4)
             dd=dmax-dmin
             if (abs(dd) > eps ) v(3)=(d0-dmin)/dd

          case (12)
             v(1)=1.0
             v(2)=1.0
             v(3)=0.5
             dmin=min(d7,d3)
             dmax=max(d7,d3)
             dd=dmax-dmin
             if (abs(dd) > eps ) v(3)=(d0-dmin)/dd

       end select

       return
    End Function Vertice_Point_Cal

    !!--++
    !!--++ Function Vertice_Point_Fix(Code_Edge) Result(V)
    !!--++     integer,          intent(in) :: code_edge
    !!--++     real(kind=cp), dimension(3)  :: v
    !!--++
    !!--++     (OVERLOADED)
    !!--++     Return the relative position point from (i,j,k) of V1
    !!--++     Given a binary dataset.
    !!--++
    !!--++ Update: February - 2005
    !!
    Function Vertice_Point_Fix(Code_Edge) Result(V)
       !---- Arguments ----!
       integer,          intent(in) :: code_edge
       real(kind=cp), dimension(3)  :: v

       !---- Local Variables ----!

       v=0.0
       select case (code_edge)
          case ( 1)
             v(1)=0.5

          case ( 2)
             v(1)=1.0
             v(2)=0.5

          case ( 3)
             v(1)=0.5
             v(2)=1.0

          case ( 4)
             v(2)=0.5

          case ( 5)
             v(1)=0.5
             v(3)=1.0

          case ( 6)
             v   =1.0
             v(2)=0.5

          case ( 7)
             v   =1.0
             v(1)=0.5

          case ( 8)
             v(2)=0.5
             v(3)=1.0

          case ( 9)
             v(3)=0.5

          case (10)
             v(1)=1.0
             v(3)=0.5

          case (11)
             v(2)=1.0
             v(3)=0.5

          case (12)
             v   =1.0
             v(3)=0.5

       end select

       return
    End Function Vertice_Point_Fix

    !!----
    !!---- Function Vertices_Cube(Index_Cube) Result(Iv)
    !!----    integer, intent(in)   :: index_cube     !  In -> Index cibe
    !!----    integer, dimension(8) :: iv             ! Out -> Vertices state
    !!----
    !!----    Return the state of the 8 vertices of the cube in Marching
    !!----    cubes algorithm
    !!----
    !!---- Update: February - 2005
    !!
    Function Vertices_Cube(Index_Cube) Result(Iv)
       !---- Arguments ----!
       integer,intent(in)    :: index_cube
       integer, dimension(8) :: iv

       !---- Local Variables ----!
       integer               :: k

       k=index_cube
       iv=0

       iv(8)=k/128
       k=k-iv(8)*128

       iv(7)=k/64
       k=k-iv(7)*64

       iv(6)=k/32
       k=k-iv(6)*32

       iv(5)=k/16
       k=k-iv(5)*16

       iv(4)=k/8
       k=k-iv(4)*8

       iv(3)=k/4
       k=k-iv(3)*4

       iv(2)=k/2
       k=k-iv(2)*2

       iv(1)=k

       return
    End Function Vertices_Cube

    !!----
    !!---- Function VPoint_in_Cube(R,S,T,X000,X001,X010,X011,X100,X101 X110,X111) Result(X)
    !!----    real(kind=cp), intent(in) :: R
    !!----    real(kind=cp), intent(in) :: S
    !!----    real(kind=cp), intent(in) :: T
    !!----    real(kind=cp), intent(in) :: X000
    !!----    real(kind=cp), intent(in) :: X001
    !!----    real(kind=cp), intent(in) :: X010
    !!----    real(kind=cp), intent(in) :: X011
    !!----    real(kind=cp), intent(in) :: X100
    !!----    real(kind=cp), intent(in) :: X101
    !!----    real(kind=cp), intent(in) :: X110
    !!----    real(kind=cp), intent(in) :: X111
    !!----    real(kind=cp)             :: X
    !!----
    !!----    Function that interpolate the value into a cube
    !!----
    !!--<<   Diagram:
    !!----
    !!----     011--------------111
    !!----       |               |
    !!----       |               |
    !!----       |               |
    !!----       |               |
    !!----       |               |
    !!----     001--------------101
    !!----
    !!----
    !!----       *---------------*
    !!----       |               |
    !!----       |               |
    !!----       |      rst      |
    !!----       |               |
    !!----       |               |
    !!----       *---------------*
    !!----
    !!----
    !!----     010--------------110
    !!----       |               |
    !!----       |               |
    !!----       |               |
    !!----       |               |
    !!----       |               |
    !!----     000--------------100
    !!----
    !!----
    !!----   Formula:
    !!----
    !!----     Written as a polynomial in R, S and T, the interpolation map has the
    !!----     form:
    !!----
    !!----       X(R,S,T) =
    !!----         1         * ( + x000 )
    !!----       + r         * ( - x000 + x100 )
    !!----       +     s     * ( - x000        + x010 )
    !!----       +         t * ( - x000               + x001 )
    !!----       + r * s     * ( + x000 - x100 - x010                       + x110 )
    !!----       + r     * t * ( + x000 - x100        - x001        + x101 )
    !!----       +     s * t * ( + x000        - x010 - x001 + x011 )
    !!----       + r * s * t * ( - x000 + x100 + x010 + x001 - x011 - x101 - x110 + x111 )
    !!-->>
    !!--..   Reference:
    !!--..
    !!--..     William Gordon,
    !!--..     Blending-Function Methods of Bivariate and Multivariate Interpolation
    !!--..       and Approximation,
    !!--..     SIAM Journal on Numerical Analysis,
    !!--..     Volume 8, Number 1, March 1971, pages 158-177.
    !!--..
    !!--..     William Gordon and Charles Hall,
    !!--..     Transfinite Element Methods: Blending-Function Interpolation over
    !!--..       Arbitrary Curved Element Domains,
    !!--..     Numerische Mathematik,
    !!--..     Volume 21, Number 1, 1973, pages 109-129.
    !!--..
    !!--..     William Gordon and Charles Hall,
    !!--..     Construction of Curvilinear Coordinate Systems and Application to
    !!--..       Mesh Generation,
    !!--..     International Journal of Numerical Methods in Engineering,
    !!--..     Volume 7, 1973, pages 461-477.
    !!--..
    !!--..     Joe Thompson, Bharat Soni, Nigel Weatherill,
    !!--..     Handbook of Grid Generation,
    !!--..     CRC Press, 1999.
    !!----
    !!---- Update: April 2008
    !!
    Function VPoint_in_Cube(R,S,T,X000,X001,X010,X011,X100,X101,X110,X111) Result(X)
       !---- Arguments ----!
       real(kind=cp), intent(in) :: R
       real(kind=cp), intent(in) :: S
       real(kind=cp), intent(in) :: T
       real(kind=cp), intent(in) :: X000
       real(kind=cp), intent(in) :: X001
       real(kind=cp), intent(in) :: X010
       real(kind=cp), intent(in) :: X011
       real(kind=cp), intent(in) :: X100
       real(kind=cp), intent(in) :: X101
       real(kind=cp), intent(in) :: X110
       real(kind=cp), intent(in) :: X111
       real(kind=cp)             :: X

       x = &
               1.0E+00     * ( + x000 ) &
               + r         * ( - x000 + x100 ) &
               +     s     * ( - x000        + x010 ) &
               +         t * ( - x000               + x001 ) &
               + r * s     * ( + x000 - x100 - x010                      + x110 ) &
               + r     * t * ( + x000 - x100        - x001        + x101 ) &
               +     s * t * ( + x000        - x010 - x001 + x011 ) &
               + r * s * t * ( - x000 + x100 + x010 + x001 - x011 - x101 - x110 + x111 )

       return
    End Function VPoint_in_Cube

    !!----
    !!---- Function VPoint_in_Line(R, X0, X1) Result (X)
    !!----    real(kind=cp), intent(in) :: R       ! R is distance between the ends points
    !!----    real(kind=cp), intent(in) :: X0      ! Value of the Point 0
    !!----    real(kind=cp), intent(in) :: X1      ! Value of the Point 1
    !!----
    !!----    Function that interpolate the value
    !!----
    !!--<<    Diagram:
    !!----
    !!----      0-----r-----1
    !!-->>
    !!--..    Reference:
    !!--..
    !!--..      William Gordon,
    !!--..      Blending-Function Methods of Bivariate and Multivariate Interpolation
    !!--..        and Approximation,
    !!--..      SIAM Journal on Numerical Analysis,
    !!--..      Volume 8, Number 1, March 1971, pages 158-177.
    !!--..
    !!--..      William Gordon and Charles Hall,
    !!--..      Transfinite Element Methods: Blending-Function Interpolation over
    !!--..        Arbitrary Curved Element Domains,
    !!--..      Numerische Mathematik,
    !!--..      Volume 21, Number 1, 1973, pages 109-129.
    !!--..
    !!--..      William Gordon and Charles Hall,
    !!--..      Construction of Curvilinear Coordinate Systems and Application to
    !!--..        Mesh Generation,
    !!--..      International Journal of Numerical Methods in Engineering,
    !!--..      Volume 7, 1973, pages 461-477.
    !!--..
    !!--..      Joe Thompson, Bharat Soni, Nigel Weatherill,
    !!--..      Handbook of Grid Generation,
    !!--..      CRC Press, 1999.
    !!----
    !!---- Update: April 2008
    !!
    Function VPoint_in_Line(R, X0, X1) Result (X)
       !---- Arguments ----!
       real(kind=cp), intent(in) :: R
       real(kind=cp), intent(in) :: X0
       real(kind=cp), intent(in) :: X1
       real(kind=cp)             :: X

       x = ( 1.0 - r ) * x0 + r * x1

       return
    End Function VPoint_in_Line

    !!----
    !!---- Function VPoint_in_Square(R, S, X00, X01, X10, X11) Result(X)
    !!----    real(kind=cp), intent(in) :: r     !
    !!----    real(kind=cp), intent(in) :: s     !
    !!----    real(kind=cp), intent(in) :: x00   !
    !!----    real(kind=cp), intent(in) :: x01   !
    !!----    real(kind=cp), intent(in) :: x10   !
    !!----    real(kind=cp), intent(in) :: x11   !
    !!----    real(kind=cp)             :: x     !
    !!----
    !!----    Function that interpolate the value on square
    !!----
    !!--<<    Diagram:
    !!----
    !!----      01------------11
    !!----       |      .      |
    !!----       |      .      |
    !!----       |.....rs......|
    !!----       |      .      |
    !!----       |      .      |
    !!----      00------------10
    !!----
    !!----    Formula:
    !!----
    !!----      Written in terms of R and S, the map has the form:
    !!----
    !!----        X(R,S) =
    !!----                 1     * ( + x00 )
    !!----               + r     * ( - x00 + x10 )
    !!----               + s     * ( - x00       + x01 )
    !!----               + r * s * ( + x00 - x10 - x01 + x11 )
    !!----
    !!----      Written in terms of the coefficients, the map has the form:
    !!----
    !!----        X(R,S) = x00 * ( 1 - r - s + r * s )
    !!----               + x01 * (         s - r * s )
    !!----               + x10 * (     r     - r * s )
    !!----               + x11 * (             r * s )
    !!----
    !!----               = x00 * ( 1 - r ) * ( 1 - s )
    !!----               + x01 * ( 1 - r ) *       s
    !!----               + x10 *       r   * ( 1 - s )
    !!----               + x11 *       r           s
    !!----
    !!----      The nonlinear term ( r * s ) has an important role:
    !!----
    !!----        If ( x01 + x10 - x00 - x11 ) is zero, then the input data lies in
    !!----        a plane, and the mapping is affine.  All the interpolated data
    !!----        will lie on the plane defined by the four corner values.  In
    !!----        particular, on any line through the square, data values at
    !!----        intermediate points will lie between the values at the endpoints.
    !!----
    !!----        If ( x01 + x10 - x00 - x11 ) is not zero, then the input data does
    !!----        not lie in a plane, and the interpolation map is nonlinear.  On
    !!----        any line through the square, data values at intermediate points
    !!----        may lie above or below the data values at the endpoints.  The
    !!----        size of the coefficient of r * s will determine how severe this
    !!-->>        effect is.
    !!----
    !!--..    Reference:
    !!--..
    !!--..      William Gordon,
    !!--..      Blending-Function Methods of Bivariate and Multivariate Interpolation
    !!--..        and Approximation,
    !!--..      SIAM Journal on Numerical Analysis,
    !!--..      Volume 8, Number 1, March 1971, pages 158-177.
    !!--..
    !!--..      William Gordon and Charles Hall,
    !!--..      Transfinite Element Methods: Blending-Function Interpolation over
    !!--..        Arbitrary Curved Element Domains,
    !!--..      Numerische Mathematik,
    !!--..      Volume 21, Number 1, 1973, pages 109-129.
    !!--..
    !!--..      William Gordon and Charles Hall,
    !!--..      Construction of Curvilinear Coordinate Systems and Application to
    !!--..        Mesh Generation,
    !!--..      International Journal of Numerical Methods in Engineering,
    !!--..      Volume 7, 1973, pages 461-477.
    !!--..
    !!--..      Joe Thompson, Bharat Soni, Nigel Weatherill,
    !!--..      Handbook of Grid Generation,
    !!--..      CRC Press, 1999.
    !!--..
    !!----
    !!---- Update: April 2008
    !!
    Function VPoint_in_Square(R, S, X00, X01, X10, X11) Result(x)
       !---- Arguments ----!
       real(kind=cp), intent(in) :: r
       real(kind=cp), intent(in) :: s
       real(kind=cp), intent(in) :: x00
       real(kind=cp), intent(in) :: x01
       real(kind=cp), intent(in) :: x10
       real(kind=cp), intent(in) :: x11
       real(kind=cp)             :: x

       x =             + x00 &
           + r *     ( - x00 + x10 ) &
           + s *     ( - x00       + x01 ) &
           + r * s * ( + x00 - x10 - x01 + x11 )

       return
    End Function VPoint_in_Square

    !!--++
    !!--++ Pure Function xy_sect(p1,p2,h,xh) result(sect)
    !!--++    integer,                       intent(in) :: p1,p2
    !!--++    real(kind=cp), dimension(0:4), intent(in) :: h,xh
    !!--++    real(kind=cp)                             :: sect
    !!--++
    !!--++    (Private)
    !!--++    Calculates the intersection of lines
    !!--++    Used internally by Calculate_Contour2D
    !!--++
    !!--++ Update: August - 2005
    !!
    Pure Function xy_sect(p1,p2,h,xy) result(sect)
      integer,                       intent(in) :: p1,p2
      real(kind=cp), dimension(0:4), intent(in) :: h,xy
      real(kind=cp)                             :: sect
      sect= (h(p2)*xy(p1)-h(p1)*xy(p2))/(h(p2)-h(p1))
      return
    End Function xy_sect

    !---------------------!
    !---- Subroutines ----!
    !---------------------!

    !!----
    !!---- Subroutine Calculate_Contour2D(d,ilb,iub,jlb,jub,x,y,z,nlv,ntp,xyz)
    !!----    integer,                                    intent(in)     :: ilb,iub    ! Lower and upper limits on the first dimension
    !!----    integer,                                    intent(in)     :: jlb,jub    ! Lower and upper limits for the second dimension
    !!----    real(kind=cp), dimension (ilb:iub,jlb:jub), intent(in)     :: d          ! Section 2D
    !!----    real(kind=cp), dimension (ilb:iub),         intent(in)     :: x          ! Limits values on X
    !!----    real(kind=cp), dimension (jlb:jub),         intent(in)     :: y          ! Limits values on Y
    !!----    real(kind=cp), dimension (:),               intent(in)     :: z          ! Level values
    !!----    integer,                                    intent(in)     :: nlv        ! Number of levels
    !!----    integer,                                    intent(in out) :: ntp        ! Number of points
    !!----    real(kind=cp), dimension (:,:),             intent(out)    :: xyz        ! XY Points
    !!----
    !!----     Subroutine for Contour 2D
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Calculate_Contour2D(d,ilb,iub,jlb,jub,x,y,z,nlv,npt,xyz)
       !---- Arguments ----!
       integer,                                    intent(in)     :: ilb,iub
       integer,                                    intent(in)     :: jlb,jub
       real(kind=cp), dimension (ilb:iub,jlb:jub), intent(in)     :: d
       real(kind=cp), dimension (ilb:iub),         intent(in)     :: x
       real(kind=cp), dimension (jlb:jub),         intent(in)     :: y
       real(kind=cp), dimension (:),               intent(in)     :: z
       integer,                                    intent(in)     :: nlv ! numeri de liveli
       integer,                                    intent(in out) :: npt ! punti
       real(kind=cp), dimension (:,:),             intent(out)    :: xyz

       !---- Local variables ----!
       integer                             :: j,i,k,m,m1,m2,m3
       integer                             :: cases
       integer, dimension (0:4)            :: sh
       integer, dimension (1:4)            :: im=(/0,1,1,0/), jm=(/0,0,1,1/)
       integer, dimension (-1:1,-1:1,-1:1) :: castab
       real(kind=cp), dimension (0:4)      :: h, xh, yh
       real(kind=cp)                       :: dmin,dmax,x1,y1,x2,y2

       !---- Use statement functions for the line intersections => replaced by a pure private function
       ! xsect(p1,p2)=(h(p2)*xh(p1)-h(p1)*xh(p2))/(h(p2)-h(p1))
       ! ysect(p1,p2)=(h(p2)*yh(p1)-h(p1)*yh(p2))/(h(p2)-h(p1))

       !---- Init ----!
       castab= reshape ( (/0,0,9,0,1,5,7,4,8,0,3,6,2,3,2,6,3,0,8,4,7,5,1,0,9,0,0/),(/3,3,3/) )

       !---- Scan the arrays, top down, left to right within rows
       do j=jub-1,jlb,-1
          do i=ilb,iub-1
             dmin=min(d(i,j),d(i,j+1),d(i+1,j),d(i+1,j+1))
             dmax=max(d(i,j),d(i,j+1),d(i+1,j),d(i+1,j+1))

             if (dmax >= z(1) .and. dmin <= z(nlv)) then
                do k=1,nlv
                   if (z(k) >= dmin .and. z(k) <= dmax) then
                      do m=4,0,-1
                         if (m > 0) then
                            h(m)=d(i+im(m),j+jm(m)) -z(k)
                            xh(m)=x(i+im(m))
                            yh(m)=y(j+jm(m))
                         else
                            h(0)=0.25*(h(1)+h(2)+h(3)+h(4))
                            xh(0)=0.5*(x(i)+x(i+1))
                            yh(0)=0.5*(y(j)+y(j+1))
                         end if
                         if ( h(m) > 0.0) then
                            sh(m)=+1
                         else if ( h(m) < 0.0) then
                            sh(m)=-1
                         else
                            sh(m)=0
                         end if
                      end do

                      !--- Scan each triangle in the box
                      do m=1,4
                         m1=m
                         m2=0
                         if (m /= 4) then
                            m3=m+1
                         else
                            m3=1
                         end if
                         cases=castab(sh(m1),sh(m2),sh(m3))
                         if (cases /= 0) then
                            select case (cases)
                               case (1)
                                  ! Case 1 - Line between vertices 1 and 2
                                  x1=xh(m1)
                                  y1=yh(m1)
                                  x2=xh(m2)
                                  y2=yh(m2)

                               case (2)
                                  ! Case 2 - Line between vertices 2 and 3
                                  x1=xh(m2)
                                  y1=yh(m2)
                                  x2=xh(m3)
                                  y2=yh(m3)

                               case (3)
                                  ! Case 3 - Line between vertices 3 and 1
                                  x1=xh(m3)
                                  y1=yh(m3)
                                  x2=xh(m1)
                                  y2=yh(m1)

                               case (4)
                                  ! Case 4 - Line between vertices 1 and side 2-3
                                  x1=xh(m1)
                                  y1=yh(m1)
                                  x2= xy_sect(m2,m3,h,xh) !xsect(m2,m3)
                                  y2= xy_sect(m2,m3,h,yh) !ysect(m2,m3)

                               case (5)
                                  ! Case 5 - Line between vertices 2 and side 3-1
                                  x1=xh(m2)
                                  y1=yh(m2)
                                  x2= xy_sect(m3,m1,h,xh) !xsect(m3,m1)
                                  y2= xy_sect(m3,m1,h,yh) !ysect(m3,m1)

                               case (6)
                                  ! Case 6 - Line between vertices 3 and side 1-2
                                  x1=xh(m3)
                                  y1=yh(m3)
                                  x2= xy_sect(m1,m2,h,xh) !xsect(m1,m2)
                                  y2= xy_sect(m1,m2,h,yh) !ysect(m1,m2)

                               case (7)
                                  ! Case 7 - Line between sides 1-2 and  2-3
                                  x1= xy_sect(m1,m2,h,xh) !xsect(m1,m2)
                                  y1= xy_sect(m1,m2,h,yh) !ysect(m1,m2)
                                  x2= xy_sect(m2,m3,h,xh) !xsect(m2,m3)
                                  y2= xy_sect(m2,m3,h,yh) !ysect(m2,m3)

                               case (8)
                                  ! Case 8 - Line between sides 2-3 and  3-1
                                  x1= xy_sect(m2,m3,h,xh) !xsect(m2,m3)
                                  y1= xy_sect(m2,m3,h,yh) !ysect(m2,m3)
                                  x2= xy_sect(m3,m1,h,xh) !xsect(m3,m1)
                                  y2= xy_sect(m3,m1,h,yh) !ysect(m3,m1)

                               case (9)
                                  ! Case 9 - Line between sides 3-1 and  1-2
                                  x1= xy_sect(m3,m1,h,xh) !xsect(m3,m1)
                                  y1= xy_sect(m3,m1,h,yh) !ysect(m3,m1)
                                  x2= xy_sect(m1,m2,h,xh) !xsect(m1,m2)
                                  y2= xy_sect(m1,m2,h,yh) !ysect(m1,m2)
                            end select

                            if (npt+2 <= Max_Points) then
                               npt=npt+1
                               xyz(1,npt)=x1
                               xyz(2,npt)=y1
                               xyz(3,npt)=real(k)

                               npt=npt+1
                               xyz(1,npt)=x2
                               xyz(2,npt)=y2
                               xyz(3,npt)=real(k)
                            end if

                         end if
                      end do
                   end if
                end do
             end if
          end do
       end do

       return
    End Subroutine Calculate_Contour2D

    !!----
    !!---- Subroutine Calculate_Mesh(Rho,Ngrid,Nlevel,Levels,MC_Method,Npoints,Xyz,Limits,Step)
    !!----    real(kind=cp),    dimension(:,:,:),         intent(in)     :: Rho
    !!----    integer, dimension(3),                      intent(in)     :: Ngrid
    !!----    integer,                                    intent(in)     :: Nlevel
    !!----    real(kind=cp),    dimension(nlevel),        intent(in)     :: Levels
    !!----    character(len=*),                           intent(in)     :: MC_Method
    !!----    integer, dimension(nlevel),                 intent(out)    :: Npoints
    !!----    real(kind=cp),    dimension(:,:),           intent(out)    :: Xyz
    !!----    real(kind=cp),    dimension(2,3), optional, intent(in)     :: Limits
    !!----    integer, dimension(3), optional,            intent(in)     :: Step
    !!----
    !!----    Calculate the 3D Contour
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Calculate_Mesh(Rho,Ngrid,Nlevel,Levels,MC_Method,NPoints,Xyz,Limits,Step)
       !---- Arguments ----!
       real(kind=cp),    dimension(:,:,:),         intent(in)           :: Rho
       integer, dimension(3),                      intent(in)           :: ngrid
       integer,                                    intent(in)           :: nlevel
       real(kind=cp),    dimension(nlevel),        intent(in)           :: Levels
       character(len=*),                           intent(in)           :: MC_Method
       integer, dimension(nlevel),                 intent(out)          :: NPoints
       real(kind=cp),    dimension(:,:),           intent(out)          :: Xyz
       real(kind=cp),    dimension(2,3), optional, intent(in)           :: Limits
       integer, dimension(3), optional,            intent(in)           :: Step

       !---- Local Variables ----!
       character(len=2)               :: mc_char
       integer                        :: i,ii,j,k,n,lv,ntr
       integer                        :: ind,nelem,ncase
       integer                        :: nx_ini,nx_fin,ny_ini,ny_fin,nz_ini,nz_fin
       integer                        :: i1,j1,k1,i2,j2,k2,i3,j3,k3,i4,j4,k4
       integer                        :: i5,j5,k5,i6,j6,k6,i7,j7,k7,i8,j8,k8
       integer, dimension(8)          :: iv
       integer, dimension(10)         :: icode
       integer, dimension(3)          :: istep
       real(kind=cp)                  :: dx,dy,dz,dxs,dys,dzs
       real(kind=cp)                  :: den1,den2,den3,den4,den5,den6,den7,den8
       real(kind=cp)                  :: xmax,xmin,ymax,ymin,zmin,zmax
       real(kind=cp),dimension(3,10)  :: xp

       logical               :: mc

       !---- Init variables ----!
       call set_cube_info()

       dx=1.0/real(ngrid(1))
       dy=1.0/real(ngrid(2))
       dz=1.0/real(ngrid(3))
       dxs=dx
       dys=dy
       dzs=dz
       istep=1

       if (present(limits)) then
          xmin=minval(limits(:,1))
          xmax=maxval(limits(:,1))
          ymin=minval(limits(:,2))
          ymax=maxval(limits(:,2))
          zmin=minval(limits(:,3))
          zmax=maxval(limits(:,3))
       else
          xmin=0.0
          xmax=1.0
          ymin=0.0
          ymax=1.0
          zmin=0.0
          zmax=1.0
       end if

       nx_ini=nint(xmin/dx)-1
       ny_ini=nint(ymin/dy)-1
       nz_ini=nint(zmin/dz)-1
       nx_fin=nint(xmax/dx)+1
       ny_fin=nint(ymax/dy)+1
       nz_fin=nint(zmax/dz)+1

       if (present(step)) then
          istep=step
          dx=istep(1)*dx
          dy=istep(2)*dy
          dz=istep(3)*dz
       end if
       mc=.false.
       mc_char=adjustl(mc_method)
       if (mc_char=="TR" .or. mc_char=="tr") mc=.true.

       npoints=0
       ntr=0

       do lv=1,nlevel
          do i=nx_ini,nx_fin,istep(1)
             do j=ny_ini,ny_fin,istep(2)
                do k=nz_ini,nz_fin,istep(3)

                   i1=i
                   j1=j
                   k1=k
                   i1=mod(i1+4*ngrid(1),ngrid(1))
                   if (i1 == 0) i1=1
                   j1=mod(j1+4*ngrid(2),ngrid(2))
                   if (j1 == 0) j1=1
                   k1=mod(k1+4*ngrid(3),ngrid(3))
                   if (k1 == 0) k1=1

                   i2=i+istep(1)
                   j2=j
                   k2=k
                   i2=mod(i2+4*ngrid(1),ngrid(1))
                   if (i2 == 0) i2=1
                   j2=mod(j2+4*ngrid(2),ngrid(2))
                   if (j2 == 0) j2=1
                   k2=mod(k2+4*ngrid(3),ngrid(3))
                   if (k2 == 0) k2=1

                   i3=i+istep(1)
                   j3=j+istep(2)
                   k3=k
                   i3=mod(i3+4*ngrid(1),ngrid(1))
                   if (i3 == 0) i3=1
                   j3=mod(j3+4*ngrid(2),ngrid(2))
                   if (j3 == 0) j3=1
                   k3=mod(k3+4*ngrid(3),ngrid(3))
                   if (k3 == 0) k3=1

                   i4=i
                   j4=j+istep(2)
                   k4=k
                   i4=mod(i4+4*ngrid(1),ngrid(1))
                   if (i4 == 0) i4=1
                   j4=mod(j4+4*ngrid(2),ngrid(2))
                   if (j4 == 0) j4=1
                   k4=mod(k4+4*ngrid(3),ngrid(3))
                   if (k4 == 0) k4=1

                   i5=i
                   j5=j
                   k5=k+istep(3)
                   i5=mod(i5+4*ngrid(1),ngrid(1))
                   if (i5 == 0) i5=1
                   j5=mod(j5+4*ngrid(2),ngrid(2))
                   if (j5 == 0) j5=1
                   k5=mod(k5+4*ngrid(3),ngrid(3))
                   if (k5 == 0) k5=1

                   i6=i+istep(1)
                   j6=j
                   k6=k+istep(3)
                   i6=mod(i6+4*ngrid(1),ngrid(1))
                   if (i6 == 0) i6=1
                   j6=mod(j6+4*ngrid(2),ngrid(2))
                   if (j6 == 0) j6=1
                   k6=mod(k6+4*ngrid(3),ngrid(3))
                   if (k6 == 0) k6=1

                   i7=i+istep(1)
                   j7=j+istep(2)
                   k7=k+istep(3)
                   i7=mod(i7+4*ngrid(1),ngrid(1))
                   if (i7 == 0) i7=1
                   j7=mod(j7+4*ngrid(2),ngrid(2))
                   if (j7 == 0) j7=1
                   k7=mod(k7+4*ngrid(3),ngrid(3))
                   if (k7 == 0) k7=1

                   i8=i
                   j8=j+istep(2)
                   k8=k+istep(3)
                   i8=mod(i8+4*ngrid(1),ngrid(1))
                   if (i8 == 0) i8=1
                   j8=mod(j8+4*ngrid(2),ngrid(2))
                   if (j8 == 0) j8=1
                   k8=mod(k8+4*ngrid(3),ngrid(3))
                   if (k8 == 0) k8=1

                   den1=Rho(i1,j1,k1)
                   den2=Rho(i2,j2,k2)
                   den3=Rho(i3,j3,k3)
                   den4=Rho(i4,j4,k4)
                   den5=Rho(i5,j5,k5)
                   den6=Rho(i6,j6,k6)
                   den7=Rho(i7,j7,k7)
                   den8=Rho(i8,j8,k8)

                   iv=0
                   if (den1 >= levels(lv)) iv(1)=1
                   if (den2 >= levels(lv)) iv(2)=1
                   if (den3 >= levels(lv)) iv(3)=1
                   if (den4 >= levels(lv)) iv(4)=1
                   if (den5 >= levels(lv)) iv(5)=1
                   if (den6 >= levels(lv)) iv(6)=1
                   if (den7 >= levels(lv)) iv(7)=1
                   if (den8 >= levels(lv)) iv(8)=1

                   ind=index_cube(iv,mc)

                   nelem=cube_info(ind)%nelem
                   if (nelem == 0) cycle
                   ncase=cube_info(ind)%code
                   icode=0
                   xp=0.0
                   select case (ncase)
                      case (1) ! Triangle
                         do n=1,nelem
                            if (ntr+3 > Max_Points) then
                              ERR_Maps=.true.
                              ERR_Maps_Mess="The number of points (ntr+3) is higher than Max_Points: 500000!"
                              exit
                            end if
                            do ii=1,3
                               icode(ii)=cube_info(ind)%edges(3*(n-1)+ii)
                               xp(:,ii)=vertice_point(icode(ii))

                               npoints(lv)=npoints(lv)+1
                               ntr=ntr+1
                               xyz(1,ntr)=(i-1)*dxs + xp(1,ii)*dx
                               xyz(2,ntr)=(j-1)*dys + xp(2,ii)*dy
                               xyz(3,ntr)=(k-1)*dzs + xp(3,ii)*dz
                               xyz(4,ntr)=real(ncase)
                            end do
                         end do

                      case (2) ! Trapezoide
                         do n=1,nelem
                            if (ntr+4 > Max_Points) then
                              ERR_Maps=.true.
                              ERR_Maps_Mess="The number of points (ntr+4) is higher than Max_Points: 500000!"
                              exit
                            end if
                            do ii=1,4
                               icode(ii)=cube_info(ind)%edges(4*(n-1)+ii)
                               xp(:,ii)=vertice_point(icode(ii))

                               npoints(lv)=npoints(lv)+1
                               ntr=ntr+1
                               xyz(1,ntr)=(i-1)*dxs + xp(1,ii)*dx
                               xyz(2,ntr)=(j-1)*dys + xp(2,ii)*dy
                               xyz(3,ntr)=(k-1)*dzs + xp(3,ii)*dz
                               xyz(4,ntr)=real(ncase)
                            end do
                         end do

                      case (3) ! Triangle + Trapezoide
                         do n=1,nelem
                            if (ntr+7 > Max_Points) then
                              ERR_Maps=.true.
                              ERR_Maps_Mess="The number of points (ntr+7) is higher than Max_Points: 500000!"
                              exit
                            end if
                            do ii=1,7
                               icode(ii)=cube_info(ind)%edges(7*(n-1)+ii)
                               xp(:,ii)=vertice_point(icode(ii))

                               npoints(lv)=npoints(lv)+1
                               ntr=ntr+1
                               xyz(1,ntr)=(i-1)*dxs + xp(1,ii)*dx
                               xyz(2,ntr)=(j-1)*dys + xp(2,ii)*dy
                               xyz(3,ntr)=(k-1)*dzs + xp(3,ii)*dz
                               xyz(4,ntr)=real(ncase)
                            end do
                         end do

                      case (4) ! Triangle + Triangle + Trapezoide
                         do n=1,nelem
                            if (ntr+10 > Max_Points) then
                              ERR_Maps=.true.
                              ERR_Maps_Mess="The number of points (ntr+10) is higher than Max_Points: 500000!"
                              exit
                            end if
                            do ii=1,10
                               icode(ii)=cube_info(ind)%edges(10*(n-1)+ii)
                               xp(:,ii)=vertice_point(icode(ii))

                               npoints(lv)=npoints(lv)+1
                               ntr=ntr+1
                               xyz(1,ntr)=(i-1)*dxs + xp(1,ii)*dx
                               xyz(2,ntr)=(j-1)*dys + xp(2,ii)*dy
                               xyz(3,ntr)=(k-1)*dzs + xp(3,ii)*dz
                               xyz(4,ntr)=real(ncase)
                            end do
                         end do

                      case (5) ! Triangle + Line
                         do n=1,nelem
                            if (ntr+4 > Max_Points) then
                              ERR_Maps=.true.
                              ERR_Maps_Mess="The number of points (ntr+4) is higher than Max_Points: 500000!"
                              exit
                            end if
                            do ii=1,4
                               icode(ii)=cube_info(ind)%edges(4*(n-1)+ii)
                               xp(:,ii)=vertice_point(icode(ii))

                               npoints(lv)=npoints(lv)+1
                               ntr=ntr+1
                               xyz(1,ntr)=(i-1)*dxs + xp(1,ii)*dx
                               xyz(2,ntr)=(j-1)*dys + xp(2,ii)*dy
                               xyz(3,ntr)=(k-1)*dzs + xp(3,ii)*dz
                               xyz(4,ntr)=real(ncase)
                            end do
                         end do

                      case (6) ! Triangle + Triangle + Line
                         do n=1,nelem
                            if (ntr+7 > Max_Points) then
                              ERR_Maps=.true.
                              ERR_Maps_Mess="The number of points (ntr+7) is higher than Max_Points: 500000!"
                              exit
                            end if
                            do ii=1,7
                               icode(ii)=cube_info(ind)%edges(7*(n-1)+ii)
                               xp(:,ii)=vertice_point(icode(ii))

                               npoints(lv)=npoints(lv)+1
                               ntr=ntr+1
                               xyz(1,ntr)=(i-1)*dxs + xp(1,ii)*dx
                               xyz(2,ntr)=(j-1)*dys + xp(2,ii)*dy
                               xyz(3,ntr)=(k-1)*dzs + xp(3,ii)*dz
                               xyz(4,ntr)=real(ncase)
                            end do
                         end do

                   end select

                end do
             end do
          end do

       end do ! lv

       return
    End Subroutine Calculate_Mesh

    !!----
    !!---- Subroutine Init_Err_Maps( )
    !!----
    !!----    Initialize the errors flags in CFML_Maps_Calculations
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Init_Err_Maps()

       ERR_Maps=.false.
       ERR_Maps_Mess=" "

       return
    End Subroutine Init_Err_Maps

    !!----
    !!---- Subroutine Load_ExtentedMap(Rho,Ngrid,Limits,Rhonew)
    !!----    real(kind=cp), dimension(:,:,:), intent(in) :: rho
    !!----    integer, dimension(3),           intent(in) :: ngrid
    !!----    real(kind=cp),dimension(2,3),    intent(in) :: limits
    !!----    real(kind=cp), dimension(:,:,:), intent(out):: rhonew
    !!----
    !!----    Rhonew has one dimension more in each dimension that Rho
    !!----    This routine is useful for 2D representation.
    !!----        Rho(nx,ny,nz) -> Rhonew(nx+1,ny+1,nz+1)
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Load_ExtendedMap(Rho,Ngrid,Limits,Rhonew)
       !---- Arguments ----!
       real(kind=cp), dimension(:,:,:), intent(in) :: rho
       integer, dimension(3),           intent(in) :: ngrid
       real(kind=cp),dimension(2,3),    intent(in) :: limits
       real(kind=cp), dimension(:,:,:), intent(out):: rhonew

       !---- Local Variables ----!
       integer                     :: nx,ny,nz,nx1,ny1,nz1
       integer                     :: i,j,k !,ii,jj,kk
       integer                     :: ii1,ii2,jj1,jj2,kk1,kk2
       real(kind=cp)               :: dx, dy, dz !, xinc, yinc, zinc
       real(kind=cp)               :: xval,yval,zval
       real(kind=cp)               :: v1,v2,v3,v4,v5,v6,v7,v8,r,s,t,valor
       real(kind=cp)               :: x1,y1,z1,xx1,xx2,yy1,yy2,zz1,zz2
       real(kind=cp), dimension(3) :: rinc

       nx=ngrid(1)
       ny=ngrid(2)
       nz=ngrid(3)
       nx1=nx+1
       ny1=ny+1
       nz1=nz+1

       dx=1.0/real(nx)
       dy=1.0/real(ny)
       dz=1.0/real(nz)

       rinc(1)=abs(limits(2,1)-limits(1,1))*dx
       rinc(2)=abs(limits(2,2)-limits(1,2))*dy
       rinc(3)=abs(limits(2,3)-limits(1,3))*dz

       !---- Loading RhoNew from Rho ----!
       if (abs(limits(1,1)-0.0) <= eps .and. abs(limits(2,1)-1.0) <= eps .and. &
           abs(limits(1,2)-0.0) <= eps .and. abs(limits(2,2)-1.0) <= eps .and. &
           abs(limits(1,3)-0.0) <= eps .and. abs(limits(2,3)-1.0) <= eps ) then

           do i=1,nx1
               do j=1,ny1
                  do k=1,nz1

                     ii1=mod(i,nx1)
                     if (ii1 == 0) ii1=1
                     jj1=mod(j,ny1)
                     if (jj1 == 0) jj1=1
                     kk1=mod(k,nz1)
                     if (kk1 == 0) kk1=1

                     rhonew(i,j,k)=rho(ii1,jj1,kk1)

                  end do
               end do
            end do

       else

          do i=1,nx1
             do j=1,ny1
                do k=1,nz1
                   !---- Pto equivalente en (0,1) ----!
                   xval=limits(1,1) + (i-1)*rinc(1)
                   yval=limits(1,2) + (j-1)*rinc(2)
                   zval=limits(1,3) + (k-1)*rinc(3)

                   if (abs(xval) <= eps) then
                      xval=0.0
                   elseif (abs(xval-1.0) <= eps) then
                      xval=1.0
                   else
                      if (xval < 0.0) xval=xval+1.0
                      if (xval > 1.0) xval=xval-1.0
                   end if

                   if (abs(yval) <= eps) then
                      yval=0.0
                   elseif (abs(yval-1.0) <= eps) then
                      yval=1.0
                   else
                      if (yval < 0.0) yval=yval+1.0
                      if (yval > 1.0) yval=yval-1.0
                   end if

                   if (abs(zval) <= eps) then
                      zval=0.0
                   elseif (abs(zval-1.0) <= eps) then
                      zval=1.0
                   else
                      if (zval < 0.0) zval=zval+1.0
                      if (zval > 1.0) zval=zval-1.0
                   end if

                   !---- Entre que planos de X esta el Pto ----!
                   do ii1=1,nx
                      ii2=ii1+1
                      xx1=(ii1-1)*dx
                      xx2=ii1*dx
                      if (abs(xx1-xval) <= eps) then
                         exit
                      elseif (xval > xx1 .and. xval <= xx2) then
                         exit
                      end if
                   end do

                   !---- Entre que planos de Y esta el Pto ----!
                   do jj1=1,ny
                      jj2=jj1+1
                      yy1=(jj1-1)*dy
                      yy2=jj1*dy
                      if (abs(yy1-yval) <= eps) then
                         exit
                      else if (yval > yy1 .and. yval <= yy2) then
                         exit
                      end if
                   end do

                   !---- Entre que planos de Z esta el Pto ----!
                   do kk1=1,nz
                      kk2=kk1+1
                      zz1=(kk1-1)*dz
                      zz2=kk1*dz
                      if (abs(zz1-zval) <= eps) then
                         exit
                      else if (zval > zz1 .and. zval <= zz2) then
                         exit
                      end if
                   end do

                   ii1=mod(ii1,nx+1)
                   if (ii1 == 0) ii1=1
                   ii2=mod(ii2,nx+1)
                   if (ii2 == 0) ii2=1

                   jj1=mod(jj1,ny+1)
                   if (jj1 == 0) jj1=1
                   jj2=mod(jj2,ny+1)
                   if (jj2 == 0) jj2=1

                   kk1=mod(kk1,nz+1)
                   if (kk1 == 0) kk1=1
                   kk2=mod(kk2,nz+1)
                   if (kk2 == 0) kk2=1

                   v1=rho(ii1,jj1,kk1)
                   v2=rho(ii2,jj1,kk1)
                   v3=rho(ii2,jj2,kk1)
                   v4=rho(ii1,jj2,kk1)
                   v5=rho(ii1,jj1,kk2)
                   v6=rho(ii2,jj1,kk2)
                   v7=rho(ii2,jj2,kk2)
                   v8=rho(ii1,jj2,kk2)

                   !---- Defino puntos del cubo ----!
                   x1=(ii1-1)*dx
                   y1=(jj1-1)*dy
                   z1=(kk1-1)*dz

                   r=(xval-x1)/dx
                   s=(yval-y1)/dy
                   t=(zval-z1)/dz

                   if (abs(r) <= eps) r=0.0
                   if (abs(s) <= eps) s=0.0
                   if (abs(t) <= eps) t=0.0

                   !---- Interpolacion tridimensional a 8 puntos ----!
                   valor=(1.0-r)*(1.0-s)*(1.0-t)* v1 + &
                              r *(1.0-s)*(1.0-t)* v2 + &
                              r *     s *(1.0-t)* v3 + &
                         (1.0-r)*     s *(1.0-t)* v4 + &
                         (1.0-r)*(1.0-s)*     t * v5 + &
                              r *(1.0-s)*     t * v6 + &
                              r *     s *     t * v7 + &
                         (1.0-r)*     s *     t * v8

                   rhonew(i,j,k)=valor

                end do
             end do
          end do

       end if

       return
    End Subroutine Load_ExtendedMap

    !!----
    !!---- SUBROUTINE Load_Section(Rho,ngrid,imap,section,limits,ngrid2,dmap)
    !!----    real(kind=cp), dimension(:,:,:), intent(in) :: rho
    !!----    integer, dimension(3),           intent(in) :: ngrid
    !!----    integer,                         intent(in) :: imap
    !!----    integer,                         intent(in) :: section
    !!----    real(kind=cp), dimension(2,2),   intent(in) :: limits
    !!----    integer, dimension(2),           intent(in) :: ngrid2
    !!----    real(kind=cp), dimension(:,:),   intent(out):: dmap
    !!----
    !!----    This routine only works with fractional coordinates
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Load_Section(Rho,ngrid,imap,section,limits,ngrid2,dmap)
       !---- Arguments ----!
       real(kind=cp), dimension(:,:,:), intent(in) :: rho
       integer, dimension(3),           intent(in) :: ngrid
       integer,                         intent(in) :: imap
       integer,                         intent(in) :: section
       real(kind=cp), dimension(2,2),   intent(in) :: limits
       integer, dimension(2),           intent(in) :: ngrid2
       real(kind=cp), dimension(:,:),   intent(out):: dmap

       !---- Local Variables ----!
       integer :: ndimx, ndimy
       integer :: i,j,ii1,ii2,jj1,jj2,kk

       real(kind=cp)            :: uinc,vinc
       real(kind=cp)            :: uu,vv,x1,y1
       real(kind=cp)            :: dx,dy,f1,f2,f3,z1,z2,z3,z4
       real(kind=cp)            :: xx1,xx2,yy1,yy2

       !---- Contorno ----!
       dmap = 0.0

       select case (imap)
          case (1)
             ndimx=ngrid(1)+1
             ndimy=ngrid(2)+1
             dx=1.0/real(ngrid(1))
             dy=1.0/real(ngrid(2))

             kk=mod(section,ngrid(3)+1)
             if (kk == 0) kk=1

          case (2)
             ndimx=ngrid(2)+1
             ndimy=ngrid(3)+1
             dx=1.0/real(ngrid(2))
             dy=1.0/real(ngrid(3))

             kk=mod(section,ngrid(1)+1)
             if (kk == 0) kk=1

          case (3)
             ndimx=ngrid(3)+1
             ndimy=ngrid(1)+1
             dx=1.0/real(ngrid(3))
             dy=1.0/real(ngrid(1))

             kk=mod(section,ngrid(2)+1)
             if (kk == 0) kk=1
       end select

       uinc=(limits(1,2) - limits(1,1))/real(ngrid2(1)-1)
       vinc=(limits(2,2) - limits(2,1))/real(ngrid2(2)-1)

       do i=1,ngrid2(1)
          uu=limits(1,1)+uinc*(i-1)

          do j=1,ngrid2(2)
             vv=limits(2,1)+vinc*(j-1)

             x1=mod(uu+10.0_cp,1.0_cp)
             y1=mod(vv+10.0_cp,1.0_cp)

             !---- Entre que nodos X ----!
             do ii1=1,ndimx-1
                xx1=(ii1-1)*dx
                xx2=xx1+dx
                if (abs(xx1-x1) <= 0.0001) then
                   exit
                else if (x1 > xx1 .and. x1 <= xx2) then
                     exit
                end if
             end do
             ii2=ii1+1
             if (ii2 == ndimx) ii2=1

             !---- Entre que nodos Y ----!
             do jj1=1,ndimy-1
                yy1=(jj1-1)*dy
                yy2=yy1+dy
                if (abs(yy1-y1) <= 0.0001) then
                   exit
                else if (y1 > yy1 .and. y1 <= yy2) then
                   exit
                end if
             end do
             jj2=jj1+1
             if (jj2 == ndimy) jj2=1

             select case (imap)
                case (1)
                   z1=rho(ii1,jj1,kk)
                   z2=rho(ii2,jj1,kk)
                   z3=rho(ii1,jj2,kk)
                   z4=rho(ii2,jj2,kk)
                case (2)
                   z1=rho(kk,ii1,jj1)
                   z2=rho(kk,ii2,jj1)
                   z3=rho(kk,ii1,jj2)
                   z4=rho(kk,ii2,jj2)
                case (3)
                   z1=rho(jj1,kk,ii1)
                   z2=rho(jj1,kk,ii2)
                   z3=rho(jj2,kk,ii1)
                   z4=rho(jj2,kk,ii2)
             end select

             f1=( (x1-xx1)/(xx2-xx1) )*(z2-z1) + z1
             f2=( (y1-yy1)/(yy2-yy1) )*(z3-z1)
             f3=( (x1-xx1)/(xx2-xx1) )*( (y1-yy1)/(yy2-yy1) )*(z1-z2-z3+z4)

             dmap(i,j)=dmap(i,j) + f1+f2+f3

          end do
       end do

       return
    End Subroutine Load_Section

    !!--++
    !!--++ Subroutine Peak_List(Pto,Grp,Cell,Npeaks,Peaks)
    !!--++    real(kind=cp), dimension(4),   intent(in)     :: Pto      ! New position to add on the List
    !!--++    type(space_group_type),        intent(in)     :: Grp      ! Space Group
    !!--++    type(crystal_cell_type),       intent(in)     :: Cell     ! Cell
    !!--++    integer,                       intent(in out) :: NPeaks   ! Number of peaks on the list
    !!--++    real(kind=cp), dimension(:,:), intent(in out) :: Peaks    ! List of Peaks
    !!--++
    !!--++    (Private)
    !!--++    Add a new peak position on the list if there is no close peak (< 0.25).
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Peak_List(Pto,Grp,Cell,Npks,Pks)
       !---- Arguments ----!
       real(kind=cp), dimension(4),   intent(in)     :: Pto
       type(space_group_type),        intent(in)     :: Grp
       type(crystal_cell_type),       intent(in)     :: Cell
       integer,                       intent(in out) :: NPks
       real(kind=cp), dimension(:,:), intent(in out) :: Pks

       !---- Local variables ----!
       integer                      :: i,j
       real(kind=cp), dimension (3) :: pto1,pto2

       !---- Pto in asymmetric unit ----!
       do i=1,grp%multip
          pto1=ApplySO(grp%Symop(i),pto(1:3))
          pto1=mod(pto1+10.0_cp,1.0_cp)
          if (pto1(1) < grp%R_Asym_Unit(1,1)     .or. &
              pto1(1) > grp%R_Asym_Unit(1,2)+eps .or. &
              pto1(2) < grp%R_Asym_Unit(2,1)     .or. &
              pto1(2) > grp%R_Asym_Unit(2,2)+eps .or. &
              pto1(3) < grp%R_Asym_Unit(3,1)     .or. &
              pto1(3) > grp%R_Asym_Unit(3,2)+eps ) cycle
          exit
       end do

       if (npks == 0) then
          !---- First Peak ---!
          npks=1
          pks(1:3,1)=pto1
          pks(4,1)=pto(4)
       else
          !---- Searching if the peak in in the list ----!
          do j=1,npks
             do i=1,grp%multip
                pto2=ApplySO(grp%Symop(i),pks(1:3,j))
                pto2=mod(pto2+10.0_cp,1.0_cp)
                if (distance(pto1,pto2,cell) <= 0.25) return
             end do
          end do
          npks=npks+1
          pks(1:3,npks)=pto1
          pks(4,npks)=pto(4)
       end if

       return
    End Subroutine Peak_List

    !!----
    !!---- Subroutine Search_Peaks(Rho,Grp,Cell,Npeaks_to_Found,Peaks,Abs_Code)
    !!----    real(kind=cp), dimension(:,:,:),    intent(in)      :: Rho         ! The Map
    !!----    type(space_group_type),             intent(in)      :: Grp         ! Space Group
    !!----    type(crystal_cell_type),            intent(in)      :: Celda       ! Cell
    !!----    integer,                            intent(in out)  :: NPFound     ! Number of peaks to found
    !!----    real(kind=cp), dimension(4,NPfound),intent(out)     :: Peaks       ! Peak List
    !!----    logical, optional,                  intent(in)      :: Abs_Code    ! logical to use absolute value on Rho
    !!----
    !!----    General procedure to search peaks on Rho
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Search_Peaks(Rho,Grp,Cell,NPFound,Pks,Abs_Code)
       !---- Arguments ----!
       real(kind=cp), dimension(:,:,:),    intent(in)      :: Rho
       type(space_group_type),             intent(in)      :: Grp
       type(crystal_cell_type),            intent(in)      :: Cell
       integer,                            intent(in out)  :: NPFound
       real(kind=cp), dimension(:,:),      intent(out)     :: Pks
       logical, optional,                  intent(in)      :: Abs_Code

       !---- Local Variables ----!
       logical                                    :: mode_abs
       integer                                    :: nscan,npks
       integer                                    :: nu,nv,nw
       integer                                    :: ii_ini,jj_ini,kk_ini,ii_fin,jj_fin,kk_fin
       integer                                    :: ia,ja,ka
       integer                                    :: i,j,k,ji,ji11,ji21,i11,i21
       integer                                    :: izt1t,izt2t,izt3t
       real(kind=cp), parameter                   :: fac_max=1.0e3
       real(kind=cp)                              :: denabsmax,fac_scale
       real(kind=cp), dimension(4)                :: pto
       integer, dimension (size(rho,1),size(rho,2),size(rho,3)):: nr3d   !automatic array

       !---- Init ----!
       call init_err_maps()
       mode_abs=.false.
       if (present(abs_code)) mode_abs=abs_code
       npks=0
       pks =0.0

       !---- Memory for NR3D ----!
       nu=size(rho,1)
       nv=size(rho,2)
       nw=size(rho,3)

       nr3d=0

       !---- Loading Rho on Nr3D ----!
       denabsmax=max(maxval(rho),abs(minval(rho)))
       fac_scale=fac_max/denabsmax
       nr3d=nint(rho*fac_scale)
       if (mode_abs) nr3d=abs(nr3d)

       !---- Searching Zone ----!
       if (grp%R_Asym_Unit(2,1) < 0.0 ) then
          ii_ini=nint(grp%R_Asym_Unit(2,1)*nv)
       else
          ii_ini=nint(grp%R_Asym_Unit(2,1)*nv) + 1
       end if
       if (grp%R_Asym_Unit(1,1) < 0.0 ) then
          jj_ini=nint(grp%R_Asym_Unit(1,1)*nu)
       else
          jj_ini=nint(grp%R_Asym_Unit(1,1)*nu) + 1
       end if
       if (grp%R_Asym_Unit(3,1) < 0.0 ) then
          kk_ini=nint(grp%R_Asym_Unit(3,1)*nw)
       else
          kk_ini=nint(grp%R_Asym_Unit(3,1)*nw) + 1
       end if
       ii_fin=nint(grp%R_Asym_Unit(2,2)*nv) + 1
       jj_fin=nint(grp%R_Asym_Unit(1,2)*nu) + 1
       kk_fin=nint(grp%R_Asym_Unit(3,2)*nw) + 1
       ii_fin=min(ii_fin,nv)
       jj_fin=min(jj_fin,nu)
       kk_fin=min(kk_fin,nw)

       !---- Searching Procedure ----!
       search:do nscan=nint(fac_max),0,-10
          do ka=kk_ini,kk_fin
             if (ka <= 0) then
                k=nw+ka
             elseif (ka > nw) then
                k=ka-nw
             else
                k=ka
             end if
             do ia=ii_ini,ii_fin
                if (ia <= 0) then
                   i=nv+ia
                elseif (ia > nv) then
                   i=ia-nv
                else
                   i=ia
                end if
                do ja=jj_ini,jj_fin
                  if (ja <= 0) then
                      j=nu+ja
                   elseif(ja > nu) then
                      j=ja-nu
                   else
                      j=ja
                   end if

                   if (nr3d(j,i,k)-nscan <= 0) cycle
                   ji=j
                   ji11=j-1
                   ji21=j+1
                   i11=i-1
                   i21=i+1
                   izt2t=k
                   izt1t=k-1
                   izt3t=k+1

                   if (ji11  <= 0) ji11=nu+ji11
                   if (i11   <= 0) i11=nv+i11
                   if (izt1t <= 0) izt1t=nw+izt1t

                   if (ji21  > nu) ji21=ji21-nu
                   if (i21   > nv) i21=i21-nv
                   if (izt3t > nw) izt3t=izt3t-nw

                   if (nr3d(ji,i,izt2t)-nr3d(ji11,  i,izt2t) <= 0) cycle   !Punto 122 >
                   if (nr3d(ji,i,izt2t)-nr3d(ji21,  i,izt2t) <  0) cycle   !      322 >=
                   if (nr3d(ji,i,izt2t)-nr3d(ji21,i11,izt2t) <= 0) cycle   !      312 >
                   if (nr3d(ji,i,izt2t)-nr3d(ji  ,i11,izt2t) <= 0) cycle   !      212 >
                   if (nr3d(ji,i,izt2t)-nr3d(ji11,i11,izt2t) <= 0) cycle   !      112 >
                   if (nr3d(ji,i,izt2t)-nr3d(ji11,i21,izt2t) <  0) cycle   !      132 >=
                   if (nr3d(ji,i,izt2t)-nr3d(ji  ,i21,izt2t) <  0) cycle   !      232 >=
                   if (nr3d(ji,i,izt2t)-nr3d(ji21,i21,izt2t) <  0) cycle   !      332 >=
                   if (nr3d(ji,i,izt2t)-nr3d(ji  ,  i,izt1t) <= 0) cycle   !      221 >
                   if (nr3d(ji,i,izt2t)-nr3d(ji11,  i,izt1t) <= 0) cycle   !      121 >
                   if (nr3d(ji,i,izt2t)-nr3d(ji21,  i,izt1t) <= 0) cycle   !      321 >
                   if (nr3d(ji,i,izt2t)-nr3d(ji  ,i11,izt1t) <= 0) cycle   !      211 >
                   if (nr3d(ji,i,izt2t)-nr3d(ji  ,i21,izt1t) <= 0) cycle   !      231 >
                   if (nr3d(ji,i,izt2t)-nr3d(ji  ,  i,izt3t) <  0) cycle   !      223 >=
                   if (nr3d(ji,i,izt2t)-nr3d(ji11,  i,izt3t) <  0) cycle   !      123 >=
                   if (nr3d(ji,i,izt2t)-nr3d(ji21,  i,izt3t) <  0) cycle   !      323 >=
                   if (nr3d(ji,i,izt2t)-nr3d(ji  ,i11,izt3t) <  0) cycle   !      213 >=
                   if (nr3d(ji,i,izt2t)-nr3d(ji  ,i21,izt3t) <  0) cycle   !      233 >=

                   !---- Position of the peak ----!
                   pto(1:3)=peak_position(nr3d,j,i,k)
                   pto(4)=rho(j,i,k)

                   !---- Good peak? ----!
                   call peak_list(pto,Grp,Cell,npks,pks)
                   if (npks == npfound) exit search
                end do
             end do
          end do
       end do search! nscan

       !if (allocated(nr3d)) deallocate(nr3d)
       npfound=npks

       return
    End Subroutine Search_Peaks

    !!----
    !!---- Subroutine Set_Cube_Info()
    !!----
    !!----    Set values for Cube_Info.
    !!----    From 0 to 127 the code is defined according the next table.
    !!----    From 128 to 255 all is defined using triangles.
    !!----
    !!--<<   Code  Figure               Process
    !!----    =====================================================================
    !!----      1   Triangle             Pto1 -> Pto2 -> Pto3 -> Pto1
    !!----
    !!----      2   Trapezoide           Pto1 -> Pto2 -> Pto3 -> Pto4 -> Pto1
    !!----
    !!----      3   Triangle             Pto1 -> Pto2 -> Pto3 -> Pto1
    !!----          Trapezoide           Pto4 -> Pto5 -> Pto6 -> Pto7 -> Pto4
    !!----
    !!----      4   Triangle             Pto1 -> Pto2 -> Pto3 -> Pto1
    !!----          Triangle             Pto4 -> Pto5 -> Pto6 -> Pto4
    !!----          Trapezoide           Pto7 -> Pto8 -> Pto9 -> Pto10 -> Pto7
    !!----
    !!----      5   Triangle + Line      Pto1 -> Pto2 -> Pto3 -> Pto1 :  Pto1 -> Pto4
    !!----
    !!----      6   Triangle             Pto1 -> Pto2 -> Pto3 -> Pto1
    !!----          Triangle + Line      Pto4 -> Pto5 -> Pto6 -> Pto4 :  Pto4 -> Pto7
    !!-->>
    !!---- Update: February - 2005
    !!
    Subroutine Set_Cube_Info()
       !---- This is a polyline configuration ----!
       cube_info(  0) = cube_info_type (0,0,(/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/))
       cube_info(  1) = cube_info_type (1,1,(/ 1, 4, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0/))
       cube_info(  2) = cube_info_type (1,1,(/ 1, 2,10, 0, 0, 0, 0, 0, 0, 0, 0, 0/))
       cube_info(  3) = cube_info_type (1,2,(/ 2, 4, 9,10, 0, 0, 0, 0, 0, 0, 0, 0/))
       cube_info(  4) = cube_info_type (1,1,(/ 2, 3,12, 0, 0, 0, 0, 0, 0, 0, 0, 0/))
       cube_info(  5) = cube_info_type (2,1,(/ 1, 4, 9, 2, 3,12, 0, 0, 0, 0, 0, 0/))
       cube_info(  6) = cube_info_type (1,2,(/ 1, 3,12,10, 0, 0, 0, 0, 0, 0, 0, 0/))
       cube_info(  7) = cube_info_type (1,3,(/ 9,10,12, 3, 4, 9,12, 0, 0, 0, 0, 0/))
       cube_info(  8) = cube_info_type (1,1,(/ 3, 4,11, 0, 0, 0, 0, 0, 0, 0, 0, 0/))
       cube_info(  9) = cube_info_type (1,2,(/ 1, 3,11, 9, 0, 0, 0, 0, 0, 0, 0, 0/))
       cube_info( 10) = cube_info_type (2,1,(/ 1, 2,10, 3, 4,11, 0, 0, 0, 0, 0, 0/))
       cube_info( 11) = cube_info_type (1,3,(/ 9,10,11, 2, 3,11,10, 0, 0, 0, 0, 0/))
       cube_info( 12) = cube_info_type (1,2,(/ 2, 4,11,12, 0, 0, 0, 0, 0, 0, 0, 0/))
       cube_info( 13) = cube_info_type (1,3,(/ 9,11,12, 1, 2,12, 9, 0, 0, 0, 0, 0/))
       cube_info( 14) = cube_info_type (1,3,(/10,11,12, 1, 4,11,10, 0, 0, 0, 0, 0/))
       cube_info( 15) = cube_info_type (1,2,(/ 9,10,12,11, 0, 0, 0, 0, 0, 0, 0, 0/))
       cube_info( 16) = cube_info_type (1,1,(/ 5, 8, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0/))
       cube_info( 17) = cube_info_type (1,2,(/ 1, 4, 8, 5, 0, 0, 0, 0, 0, 0, 0, 0/))
       cube_info( 18) = cube_info_type (2,1,(/ 1, 2,10, 5, 8, 9, 0, 0, 0, 0, 0, 0/))
       cube_info( 19) = cube_info_type (1,3,(/ 2, 4, 8, 2, 8, 5,10, 0, 0, 0, 0, 0/))
       cube_info( 20) = cube_info_type (2,1,(/ 2, 3,12, 5, 8, 9, 0, 0, 0, 0, 0, 0/))
       cube_info( 21) = cube_info_type (1,3,(/ 2, 3,12, 1, 4, 8, 5, 0, 0, 0, 0, 0/))
       cube_info( 22) = cube_info_type (1,3,(/ 5, 8, 9, 1, 3,12,10, 0, 0, 0, 0, 0/))
       cube_info( 23) = cube_info_type (3,1,(/ 3, 4,12, 4, 5, 8, 5,10,12, 0, 0, 0/))
       cube_info( 24) = cube_info_type (2,1,(/ 3, 4,11, 5, 8, 9, 0, 0, 0, 0, 0, 0/))
       cube_info( 25) = cube_info_type (1,3,(/ 1, 3, 5, 3,11, 8, 5, 8, 0, 0, 0, 0/))
       cube_info( 26) = cube_info_type (3,1,(/ 1, 2,10, 3, 4,11, 5, 8, 9, 0, 0, 0/))
       cube_info( 27) = cube_info_type (2,5,(/ 2, 5,10, 3, 8,11, 3, 5, 0, 0, 0, 0/))
       cube_info( 28) = cube_info_type (1,3,(/ 5, 8, 9, 2, 4,11,12, 0, 0, 0, 0, 0/))
       cube_info( 29) = cube_info_type (3,1,(/ 1, 2, 5, 5, 8,11, 2,11,12, 0, 0, 0/))
       cube_info( 30) = cube_info_type (1,4,(/ 5, 8, 9,10,11,12, 1,10,11, 4, 0, 0/))
       cube_info( 31) = cube_info_type (1,6,(/ 5,10, 8,11,12,10, 8, 0, 0, 0, 0, 0/))
       cube_info( 32) = cube_info_type (1,1,(/ 5, 6,10, 0, 0, 0, 0, 0, 0, 0, 0, 0/))
       cube_info( 33) = cube_info_type (2,1,(/ 1, 4, 9, 5, 6,10, 0, 0, 0, 0, 0, 0/))
       cube_info( 34) = cube_info_type (1,2,(/ 1, 2, 6, 5, 0, 0, 0, 0, 0, 0, 0, 0/))
       cube_info( 35) = cube_info_type (1,3,(/ 2, 4, 6, 5, 6, 4, 9, 0, 0, 0, 0, 0/))
       cube_info( 36) = cube_info_type (2,1,(/ 2, 3,12, 5, 6,10, 0, 0, 0, 0, 0, 0/))
       cube_info( 37) = cube_info_type (3,1,(/ 1, 4, 9, 2, 3,12, 5, 6,10, 0, 0, 0/))
       cube_info( 38) = cube_info_type (1,3,(/ 1, 3, 5, 3, 5, 6,12, 0, 0, 0, 0, 0/))
       cube_info( 39) = cube_info_type (2,5,(/ 4, 5, 9, 3, 6, 3,12, 5, 0, 0, 0, 0/))
       cube_info( 40) = cube_info_type (2,1,(/ 3, 4,11, 5, 6,10, 0, 0, 0, 0, 0, 0/))
       cube_info( 41) = cube_info_type (1,3,(/ 5, 6,10, 1, 3,11, 9, 0, 0, 0, 0, 0/))
       cube_info( 42) = cube_info_type (1,3,(/ 3, 4,11, 1, 2, 6, 5, 0, 0, 0, 0, 0/))
       cube_info( 43) = cube_info_type (3,1,(/ 2, 3, 6, 5, 6, 9, 3, 9,11, 0, 0, 0/))
       cube_info( 44) = cube_info_type (1,3,(/ 5, 6,10, 2, 4,11,12, 0, 0, 0, 0, 0/))
       cube_info( 45) = cube_info_type (1,4,(/ 1, 2,10, 9,11,12, 5, 6,12, 9, 0, 0/))
       cube_info( 46) = cube_info_type (3,1,(/ 1, 4, 5, 4,11,12, 5, 6,12, 0, 0, 0/))
       cube_info( 47) = cube_info_type (1,3,(/ 9,11,12, 5, 6,12, 9, 0, 0, 0, 0, 0/))
       cube_info( 48) = cube_info_type (1,2,(/ 6, 8, 9,10, 0, 0, 0, 0, 0, 0, 0, 0/))
       cube_info( 49) = cube_info_type (1,3,(/ 4, 6, 8, 1, 4, 6,10, 0, 0, 0, 0, 0/))
       cube_info( 50) = cube_info_type (1,3,(/ 2, 6, 8, 1, 2, 8, 9, 0, 0, 0, 0, 0/))
       cube_info( 51) = cube_info_type (1,2,(/ 2, 4, 8, 6, 0, 0, 0, 0, 0, 0, 0, 0/))
       cube_info( 52) = cube_info_type (1,3,(/ 2, 3,12, 6, 8, 9,10, 0, 0, 0, 0, 0/))
       cube_info( 53) = cube_info_type (1,4,(/ 2, 3,12, 4, 6, 8, 1, 4, 6,10, 0, 0/))
       cube_info( 54) = cube_info_type (3,1,(/ 1, 3, 9, 6, 8, 9, 3, 6,12, 0, 0, 0/))
       cube_info( 55) = cube_info_type (1,6,(/ 4, 8, 6,12, 4, 3, 6, 0, 0, 0, 0, 0/))
       cube_info( 56) = cube_info_type (1,3,(/ 3, 4,11, 6, 8, 9,10, 0, 0, 0, 0, 0/))
       cube_info( 57) = cube_info_type (3,1,(/ 1, 3,10, 6, 8,10, 3, 8,11, 0, 0, 0/))
       cube_info( 58) = cube_info_type (1,4,(/ 3, 4,11, 2, 6, 8, 1, 2, 8, 9, 0, 0/))
       cube_info( 59) = cube_info_type (1,6,(/ 2, 6, 8,11, 3, 2, 8, 0, 0, 0, 0, 0/))
       cube_info( 60) = cube_info_type (2,2,(/ 2, 4,11,12, 6, 8, 9,10, 0, 0, 0, 0/))
       cube_info( 61) = cube_info_type (1,3,(/ 1, 2,10, 6, 8,11,12, 0, 0, 0, 0, 0/))
       cube_info( 62) = cube_info_type (1,3,(/ 1, 4, 9, 6, 8,11,12, 0, 0, 0, 0, 0/))
       cube_info( 63) = cube_info_type (1,2,(/ 6, 8,11,12, 0, 0, 0, 0, 0, 0, 0, 0/))
       cube_info( 64) = cube_info_type (1,1,(/ 6, 7,12, 0, 0, 0, 0, 0, 0, 0, 0, 0/))
       cube_info( 65) = cube_info_type (2,1,(/ 1, 4, 9, 6, 7,12, 0, 0, 0, 0, 0, 0/))
       cube_info( 66) = cube_info_type (2,1,(/ 1, 2,10, 6, 7,12, 0, 0, 0, 0, 0, 0/))
       cube_info( 67) = cube_info_type (1,3,(/ 6, 7,12, 2, 4, 9,10, 0, 0, 0, 0, 0/))
       cube_info( 68) = cube_info_type (1,2,(/ 2, 3, 7, 6, 0, 0, 0, 0, 0, 0, 0, 0/))
       cube_info( 69) = cube_info_type (1,3,(/ 1, 4, 9, 2, 3, 7, 6, 0, 0, 0, 0, 0/))
       cube_info( 70) = cube_info_type (1,3,(/ 1, 3, 7, 1, 7, 6,10, 0, 0, 0, 0, 0/))
       cube_info( 71) = cube_info_type (3,1,(/ 3, 4, 7, 4, 9,10, 6, 7,10, 0, 0, 0/))
       cube_info( 72) = cube_info_type (2,1,(/ 3, 4,11, 6, 7,12, 0, 0, 0, 0, 0, 0/))
       cube_info( 73) = cube_info_type (1,3,(/ 6, 7,12, 1, 3,11, 9, 0, 0, 0, 0, 0/))
       cube_info( 74) = cube_info_type (3,1,(/ 1, 2,10, 3, 4,11, 6, 7,12, 0, 0, 0/))
       cube_info( 75) = cube_info_type (1,4,(/ 6, 7,12, 9,10,11, 2, 3,11,10, 0, 0/))
       cube_info( 76) = cube_info_type (1,3,(/ 2, 4, 6, 4, 6, 7,11, 0, 0, 0, 0, 0/))
       cube_info( 77) = cube_info_type (3,1,(/ 1, 2, 6, 1, 9,11, 6, 7,11, 0, 0, 0/))
       cube_info( 78) = cube_info_type (2,5,(/ 4, 7,11, 1, 6,10, 1, 7, 0, 0, 0, 0/))
       cube_info( 79) = cube_info_type (1,3,(/ 9,10,11, 6, 7,11,10, 0, 0, 0, 0, 0/))
       cube_info( 80) = cube_info_type (2,1,(/ 5, 8, 9, 6, 7,12, 0, 0, 0, 0, 0, 0/))
       cube_info( 81) = cube_info_type (1,3,(/ 6, 7,12, 1, 4, 8, 5, 7, 0, 0, 0, 0/))
       cube_info( 82) = cube_info_type (3,1,(/ 1, 2,10, 5, 8, 9, 6, 7,12, 0, 0, 0/))
       cube_info( 83) = cube_info_type (1,4,(/ 6, 7,12, 2, 4, 8, 2, 8, 5,10, 0, 0/))
       cube_info( 84) = cube_info_type (1,3,(/ 5, 8, 9, 2, 3, 7, 6, 0, 0, 0, 0, 0/))
       cube_info( 85) = cube_info_type (2,2,(/ 1, 4, 8, 5, 2, 3, 7, 6, 0, 0, 0, 0/))
       cube_info( 86) = cube_info_type (1,4,(/ 5, 8, 9, 1, 3, 7, 1, 7, 6,10, 0, 0/))
       cube_info( 87) = cube_info_type (1,3,(/ 5, 6,10, 3, 4, 8, 7, 0, 0, 0, 0, 0/))
       cube_info( 88) = cube_info_type (3,1,(/ 3, 4,11, 5, 8, 9, 6, 7,12, 0, 0, 0/))
       cube_info( 89) = cube_info_type (1,4,(/ 6, 7,12, 1, 3, 5, 3, 5, 8,11, 0, 0/))
       cube_info( 90) = cube_info_type (4,1,(/ 1, 4, 9, 2, 3,12, 5, 6,10, 7, 8,11/))
       cube_info( 91) = cube_info_type (3,1,(/ 2, 3,12, 5, 6,10, 7, 8,11, 0, 0, 0/))
       cube_info( 92) = cube_info_type (1,4,(/ 5, 8, 9, 2, 4, 6, 4, 6, 7,11, 0, 0/))
       cube_info( 93) = cube_info_type (1,3,(/ 7, 8,11, 1, 2, 6, 5, 0, 0, 0, 0, 0/))
       cube_info( 94) = cube_info_type (3,1,(/ 1, 4, 9, 5, 6,10, 7, 8,11, 0, 0, 0/))
       cube_info( 95) = cube_info_type (2,1,(/ 5, 6,10, 7, 8,11, 0, 0, 0, 0, 0, 0/))
       cube_info( 96) = cube_info_type (1,2,(/ 5, 7,12,10, 0, 0, 0, 0, 0, 0, 0, 0/))
       cube_info( 97) = cube_info_type (1,3,(/ 1, 4, 9, 5, 7,12,10, 0, 0, 0, 0, 0/))
       cube_info( 98) = cube_info_type (1,3,(/ 1, 5, 7, 1, 2,12, 7, 0, 0, 0, 0, 0/))
       cube_info( 99) = cube_info_type (3,1,(/ 2, 4, 9, 5, 7, 9, 2, 7,12, 0, 0, 0/))
       cube_info(100) = cube_info_type (1,3,(/ 3, 5, 7, 2, 3, 5,10, 0, 0, 0, 0, 0/))
       cube_info(101) = cube_info_type (1,4,(/ 1, 4, 9, 3, 5, 7, 2, 3, 5,10, 0, 0/))
       cube_info(102) = cube_info_type (1,2,(/ 1, 3, 7, 5, 0, 0, 0, 0, 0, 0, 0, 0/))
       cube_info(103) = cube_info_type (1,3,(/ 3, 5, 7, 3, 4, 9, 5, 0, 0, 0, 0, 0/))
       cube_info(104) = cube_info_type (1,3,(/ 3, 4,11, 5, 7,12,10, 0, 0, 0, 0, 0/))
       cube_info(105) = cube_info_type (2,2,(/ 1, 3,11, 9, 5, 7,12,10, 0, 0, 0, 0/))
       cube_info(106) = cube_info_type (1,4,(/ 3, 4,11, 1, 5, 7, 1, 2,12, 7, 0, 0/))
       cube_info(107) = cube_info_type (1,3,(/ 2, 3,12, 5, 7,11, 9, 0, 0, 0, 0, 0/))
       cube_info(108) = cube_info_type (3,1,(/ 2, 4,10, 5, 7,10, 4, 7,11, 0, 0, 0/))
       cube_info(109) = cube_info_type (1,3,(/ 1, 2,10, 5, 7,11, 9, 0, 0, 0, 0, 0/))
       cube_info(110) = cube_info_type (1,3,(/ 1, 5, 7, 1, 4,11, 7, 0, 0, 0, 0, 0/))
       cube_info(111) = cube_info_type (1,2,(/ 5, 7,11, 9, 0, 0, 0, 0, 0, 0, 0, 0/))
       cube_info(112) = cube_info_type (1,3,(/ 9,10,12, 7, 8, 9,12, 0, 0, 0, 0, 0/))
       cube_info(113) = cube_info_type (3,1,(/ 1,10,12, 1, 4, 8, 7, 8,12, 0, 0, 0/))
       cube_info(114) = cube_info_type (2,5,(/ 2, 7,12, 1, 8, 9, 1, 7, 0, 0, 0, 0/))
       cube_info(115) = cube_info_type (1,3,(/ 2, 4, 8, 2, 8, 7,12, 0, 0, 0, 0, 0/))
       cube_info(116) = cube_info_type (3,1,(/ 2, 3, 7, 2, 9,10, 7, 8, 9, 0, 0, 0/))
       cube_info(117) = cube_info_type (1,3,(/ 1, 2,10, 3, 4, 8, 7, 0, 0, 0, 0, 0/))
       cube_info(118) = cube_info_type (1,3,(/ 1, 3, 7, 1, 7, 8, 9, 0, 0, 0, 0, 0/))
       cube_info(119) = cube_info_type (1,2,(/ 3, 4, 8, 7, 0, 0, 0, 0, 0, 0, 0, 0/))
       cube_info(120) = cube_info_type (1,4,(/ 3, 4,11, 9,10,12, 7, 8, 9,12, 0, 0/))
       cube_info(121) = cube_info_type (1,3,(/ 7, 8,11, 1, 3,12,10, 0, 0, 0, 0, 0/))
       cube_info(122) = cube_info_type (3,1,(/ 1, 4, 9, 2, 3,12, 7, 8,11, 0, 0, 0/))
       cube_info(123) = cube_info_type (2,1,(/ 2, 3,12, 7, 8,11, 0, 0, 0, 0, 0, 0/))
       cube_info(124) = cube_info_type (1,3,(/ 7, 8,11, 2, 4, 9,10, 0, 0, 0, 0, 0/))
       cube_info(125) = cube_info_type (2,1,(/ 1, 2,10, 7, 8,11, 0, 0, 0, 0, 0, 0/))
       cube_info(126) = cube_info_type (2,1,(/ 1, 4, 9, 7, 8,11, 0, 0, 0, 0, 0, 0/))
       cube_info(127) = cube_info_type (1,1,(/ 7, 8,11, 0, 0, 0, 0, 0, 0, 0, 0, 0/))

       !---- This is a triangle configuration ----!
       cube_info(128) = cube_info_type (1,1,(/ 7, 8,11, 0, 0, 0, 0, 0, 0, 0, 0, 0/))
       cube_info(129) = cube_info_type (2,1,(/ 1, 4, 9, 7, 8,11, 0, 0, 0, 0, 0, 0/))
       cube_info(130) = cube_info_type (2,1,(/ 1, 2,10, 7, 8,11, 0, 0, 0, 0, 0, 0/))
       cube_info(131) = cube_info_type (3,1,(/ 2, 4,10, 4, 9,10, 7, 8,11, 0, 0, 0/))
       cube_info(132) = cube_info_type (2,1,(/ 2, 3,12, 7, 8,11, 0, 0, 0, 0, 0, 0/))
       cube_info(133) = cube_info_type (3,1,(/ 1, 4, 9, 2, 3,12, 7, 8,11, 0, 0, 0/))
       cube_info(134) = cube_info_type (3,1,(/ 1, 3,10, 3,10,12, 7, 8,11, 0, 0, 0/))
       cube_info(135) = cube_info_type (4,1,(/ 3, 4,11, 7, 8, 9, 7, 9,12, 9,10,12/))
       cube_info(136) = cube_info_type (2,1,(/ 3, 4, 7, 4, 7, 8, 0, 0, 0, 0, 0, 0/))
       cube_info(137) = cube_info_type (3,1,(/ 1, 3, 7, 1, 7, 9, 7, 8, 9, 0, 0, 0/))
       cube_info(138) = cube_info_type (3,1,(/ 1, 2,10, 3, 4, 7, 4, 7, 8, 0, 0, 0/))
       cube_info(139) = cube_info_type (3,1,(/ 2, 3, 7, 2, 9,10, 7, 8, 9, 0, 0, 0/))
       cube_info(140) = cube_info_type (3,1,(/ 2, 4, 8, 2, 7, 8, 2, 7,12, 0, 0, 0/))
       cube_info(141) = cube_info_type (4,1,(/ 1, 2, 7, 2, 7, 8, 1, 8, 9, 2, 7,12/))
       cube_info(142) = cube_info_type (3,1,(/ 1,10,12, 1, 4, 8, 7, 8,12, 0, 0, 0/))
       cube_info(143) = cube_info_type (3,1,(/ 7, 8, 9, 7, 9,12, 9,10,12, 0, 0, 0/))
       cube_info(144) = cube_info_type (2,1,(/ 5, 7, 9, 7, 9,11, 0, 0, 0, 0, 0, 0/))
       cube_info(145) = cube_info_type (3,1,(/ 1, 4,11, 1, 5, 7, 1, 7,11, 0, 0, 0/))
       cube_info(146) = cube_info_type (3,1,(/ 1, 2,10, 5, 7, 9, 7, 9,11, 0, 0, 0/))
       cube_info(147) = cube_info_type (3,1,(/ 2, 4,10, 5, 7,10, 4, 7,11, 0, 0, 0/))
       cube_info(148) = cube_info_type (3,1,(/ 2, 3,12, 5, 7, 9, 7, 9,11, 0, 0, 0/))
       cube_info(149) = cube_info_type (4,1,(/ 1, 2,12, 1, 5, 7, 1, 7,12, 3, 4,11/))
       cube_info(150) = cube_info_type (4,1,(/ 1, 3, 9, 3, 9,11, 5, 7,10, 7,10,12/))
       cube_info(151) = cube_info_type (3,1,(/ 3, 4,11, 5, 7,10, 7,10,12, 0, 0, 0/))
       cube_info(152) = cube_info_type (3,1,(/ 3, 4, 5, 3, 5, 7, 4, 5, 9, 0, 0, 0/))
       cube_info(153) = cube_info_type (2,1,(/ 1, 3, 5, 3, 5, 7, 0, 0, 0, 0, 0, 0/))
       cube_info(154) = cube_info_type (4,1,(/ 1, 4, 9, 2, 3,10, 3, 5, 7, 3, 5,10/))
       cube_info(155) = cube_info_type (3,1,(/ 2, 3,10, 3, 5, 7, 3, 5,10, 0, 0, 0/))
       cube_info(156) = cube_info_type (3,1,(/ 2, 4, 9, 5, 7, 9, 2, 7,12, 0, 0, 0/))
       cube_info(157) = cube_info_type (3,1,(/ 1, 2, 7, 1, 5, 7, 2, 7,12, 0, 0, 0/))
       cube_info(158) = cube_info_type (3,1,(/ 1, 4, 9, 5, 7,10, 7,10,12, 0, 0, 0/))
       cube_info(159) = cube_info_type (2,1,(/ 5, 7,10, 7,10,12, 0, 0, 0, 0, 0, 0/))
       cube_info(160) = cube_info_type (2,1,(/ 5, 6,10, 7, 8,11, 0, 0, 0, 0, 0, 0/))
       cube_info(161) = cube_info_type (3,1,(/ 1, 4, 9, 5, 6,10, 7, 8,11, 0, 0, 0/))
       cube_info(162) = cube_info_type (3,1,(/ 1, 2, 5, 2, 5, 6, 7, 8,11, 0, 0, 0/))
       cube_info(163) = cube_info_type (4,1,(/ 2, 4, 6, 4, 6,11, 5, 8, 9, 6, 7,11/))
       cube_info(164) = cube_info_type (3,1,(/ 2, 3,12, 5, 6,10, 7, 8,11, 0, 0, 0/))
       cube_info(165) = cube_info_type (4,1,(/ 1, 4, 9, 2, 3,12, 5, 6,10, 7, 8,11/))
       cube_info(166) = cube_info_type (4,1,(/ 1, 3, 5, 3, 8, 5, 3, 8,11, 6, 7,12/))
       cube_info(167) = cube_info_type (3,1,(/ 3, 4,11, 5, 8, 9, 6, 7,12, 0, 0, 0/))
       cube_info(168) = cube_info_type (3,1,(/ 3, 4, 7, 4, 7, 8, 5, 6,10, 0, 0, 0/))
       cube_info(169) = cube_info_type (4,1,(/ 1, 3, 7, 1, 6,10, 1, 6, 7, 5, 8, 9/))
       cube_info(170) = cube_info_type (4,1,(/ 1, 4, 5, 4, 5, 8, 2, 3, 6, 3, 6, 7/))
       cube_info(171) = cube_info_type (3,1,(/ 2, 3, 6, 3, 6, 7, 5, 8, 9, 0, 0, 0/))
       cube_info(172) = cube_info_type (4,1,(/ 2, 4, 8, 2, 8,10, 5, 8,10, 6, 7,12/))
       cube_info(173) = cube_info_type (3,1,(/ 1, 2,10, 5, 8, 9, 6, 7,12, 0, 0, 0/))
       cube_info(174) = cube_info_type (3,1,(/ 1, 4, 5, 4, 5, 8, 6, 7,12, 0, 0, 0/))
       cube_info(175) = cube_info_type (2,1,(/ 5, 8, 9, 6, 7,12, 0, 0, 0, 0, 0, 0/))
       cube_info(176) = cube_info_type (3,1,(/ 6, 7,11, 6,10,11, 9,10,11, 0, 0, 0/))
       cube_info(177) = cube_info_type (4,1,(/ 1, 4, 6, 4, 6, 7, 1, 6,10, 4, 7,11/))
       cube_info(178) = cube_info_type (3,1,(/ 1, 2, 6, 1, 9,11, 6, 7,11, 0, 0, 0/))
       cube_info(179) = cube_info_type (3,1,(/ 2, 4, 6, 4, 6,11, 6, 7,11, 0, 0, 0/))
       cube_info(180) = cube_info_type (4,1,(/ 2, 3,10, 3,10,11, 6, 7,12, 9,10,11/))
       cube_info(181) = cube_info_type (3,1,(/ 1, 2,10, 3, 4,11, 6, 7,12, 0, 0, 0/))
       cube_info(182) = cube_info_type (3,1,(/ 1, 3, 9, 3, 9,11, 6, 7,12, 0, 0, 0/))
       cube_info(183) = cube_info_type (2,1,(/ 3, 4,11, 6, 7,12, 0, 0, 0, 0, 0, 0/))
       cube_info(184) = cube_info_type (3,1,(/ 3, 4, 7, 4, 9,10, 6, 7,10, 0, 0, 0/))
       cube_info(185) = cube_info_type (3,1,(/ 1, 3, 7, 1, 6, 7, 1, 6,10, 0, 0, 0/))
       cube_info(186) = cube_info_type (3,1,(/ 1, 4, 9, 2, 3, 6, 3, 6, 7, 0, 0, 0/))
       cube_info(187) = cube_info_type (2,1,(/ 2, 3, 6, 3, 6, 7, 0, 0, 0, 0, 0, 0/))
       cube_info(188) = cube_info_type (3,1,(/ 2, 4,10, 4, 9,10, 6, 7,12, 0, 0, 0/))
       cube_info(189) = cube_info_type (2,1,(/ 1, 2,10, 6, 7,12, 0, 0, 0, 0, 0, 0/))
       cube_info(190) = cube_info_type (2,1,(/ 1, 4, 9, 6, 7,12, 0, 0, 0, 0, 0, 0/))
       cube_info(191) = cube_info_type (1,1,(/ 6, 7,12, 0, 0, 0, 0, 0, 0, 0, 0, 0/))
       cube_info(192) = cube_info_type (2,1,(/ 6, 8,11, 6,11,12, 0, 0, 0, 0, 0, 0/))
       cube_info(193) = cube_info_type (3,1,(/ 1, 4, 9, 6, 8,12, 8,11,12, 0, 0, 0/))
       cube_info(194) = cube_info_type (3,1,(/ 1, 2,10, 6, 8,12, 8,11,12, 0, 0, 0/))
       cube_info(195) = cube_info_type (4,1,(/ 2, 4,12, 4,11,12, 6, 9,10, 6, 8, 9/))
       cube_info(196) = cube_info_type (3,1,(/ 2, 3,11, 2, 6, 8, 2, 8,11, 0, 0, 0/))
       cube_info(197) = cube_info_type (4,1,(/ 1, 2, 9, 2, 6, 8, 2, 8, 9, 3, 4,11/))
       cube_info(198) = cube_info_type (3,1,(/ 1, 3,10, 6, 8,10, 3, 8,11, 0, 0, 0/))
       cube_info(199) = cube_info_type (3,1,(/ 3, 4,11, 6, 8,10, 8, 9,10, 0, 0, 0/))
       cube_info(200) = cube_info_type (3,1,(/ 3, 4,12, 4, 6, 8, 4, 6,12, 0, 0, 0/))
       cube_info(201) = cube_info_type (3,1,(/ 1, 3, 9, 6, 8, 9, 3, 6,12, 0, 0, 0/))
       cube_info(202) = cube_info_type (4,1,(/ 1, 4, 6, 1, 6,10, 2, 3,12, 4, 6, 8/))
       cube_info(203) = cube_info_type (3,1,(/ 2, 3,12, 6, 8,10, 8, 9,10, 0, 0, 0/))
       cube_info(204) = cube_info_type (2,1,(/ 2, 4, 6, 4, 6, 8, 0, 0, 0, 0, 0, 0/))
       cube_info(205) = cube_info_type (3,1,(/ 1, 2, 9, 2, 6, 8, 2, 8, 9, 0, 0, 0/))
       cube_info(206) = cube_info_type (3,1,(/ 1, 4, 6, 1, 6,10, 4, 6, 8, 0, 0, 0/))
       cube_info(207) = cube_info_type (2,1,(/ 6, 8, 9, 6, 9,10, 0, 0, 0, 0, 0, 0/))
       cube_info(208) = cube_info_type (3,1,(/ 5, 6,12, 5, 9,12, 9,11,12, 0, 0, 0/))
       cube_info(209) = cube_info_type (3,1,(/ 1, 4, 5, 4,11,12, 5, 6,12, 0, 0, 0/))
       cube_info(210) = cube_info_type (4,1,(/ 1, 2,10, 5, 9,12, 5, 6,12, 9,11,12/))
       cube_info(211) = cube_info_type (3,1,(/ 2, 4,12, 4,11,12, 5, 6,10, 0, 0, 0/))
       cube_info(212) = cube_info_type (3,1,(/ 2, 3, 6, 5, 6, 9, 3, 9,11, 0, 0, 0/))
       cube_info(213) = cube_info_type (3,1,(/ 1, 2, 5, 2, 5, 6, 3, 4,11, 0, 0, 0/))
       cube_info(214) = cube_info_type (3,1,(/ 1, 3, 9, 3, 9,11, 5, 6,10, 0, 0, 0/))
       cube_info(215) = cube_info_type (2,1,(/ 3, 4,11, 5, 6,10, 0, 0, 0, 0, 0, 0/))
       cube_info(216) = cube_info_type (4,1,(/ 3, 4, 5, 4, 5, 6, 3, 6,12, 4, 5, 9/))
       cube_info(217) = cube_info_type (3,1,(/ 1, 3, 5, 3, 5,12, 5, 6,12, 0, 0, 0/))
       cube_info(218) = cube_info_type (3,1,(/ 1, 4, 9, 2, 3,12, 5, 6,10, 0, 0, 0/))
       cube_info(219) = cube_info_type (2,1,(/ 2, 3,12, 5, 6,10, 0, 0, 0, 0, 0, 0/))
       cube_info(220) = cube_info_type (3,1,(/ 2, 4, 6, 4, 5, 6, 4, 5, 9, 0, 0, 0/))
       cube_info(221) = cube_info_type (2,1,(/ 1, 2, 5, 2, 5, 6, 0, 0, 0, 0, 0, 0/))
       cube_info(222) = cube_info_type (2,1,(/ 1, 4, 9, 5, 6,10, 0, 0, 0, 0, 0, 0/))
       cube_info(223) = cube_info_type (1,1,(/ 5, 6,10, 0, 0, 0, 0, 0, 0, 0, 0, 0/))
       cube_info(224) = cube_info_type (3,1,(/ 5, 8,10, 8,10,11,10,11,12, 0, 0, 0/))
       cube_info(225) = cube_info_type (4,1,(/ 5, 8, 9,10,11,12, 1,10, 4, 4,10,11/))
       cube_info(226) = cube_info_type (3,1,(/ 1, 2, 5, 5, 8,11, 2,11,12, 0, 0, 0/))
       cube_info(227) = cube_info_type (3,1,(/ 2, 4,12, 4,11,12, 5, 8, 9, 0, 0, 0/))
       cube_info(228) = cube_info_type (4,1,(/ 2, 3, 5, 2, 5, 8, 2, 5,10, 3, 8,11/))
       cube_info(229) = cube_info_type (3,1,(/ 1, 2,10, 3, 4,11, 5, 8, 9, 0, 0, 0/))
       cube_info(230) = cube_info_type (3,1,(/ 1, 3, 5, 3, 5, 8, 3, 8,11, 0, 0, 0/))
       cube_info(231) = cube_info_type (2,1,(/ 3, 4,11, 5, 8, 9, 0, 0, 0, 0, 0, 0/))
       cube_info(232) = cube_info_type (3,1,(/ 3, 4,12, 4, 5, 8, 5,10,12, 0, 0, 0/))
       cube_info(233) = cube_info_type (3,1,(/ 1, 3,10, 3,10,12, 5, 8, 9, 0, 0, 0/))
       cube_info(234) = cube_info_type (3,1,(/ 1, 4, 5, 4, 5, 8, 2, 3,12, 0, 0, 0/))
       cube_info(235) = cube_info_type (2,1,(/ 2, 3,12, 5, 8, 9, 0, 0, 0, 0, 0, 0/))
       cube_info(236) = cube_info_type (3,1,(/ 2, 4, 8, 2, 8,10, 5, 8,10, 0, 0, 0/))
       cube_info(237) = cube_info_type (2,1,(/ 1, 2,10, 5, 8, 9, 0, 0, 0, 0, 0, 0/))
       cube_info(238) = cube_info_type (2,1,(/ 1, 4, 5, 4, 5, 8, 0, 0, 0, 0, 0, 0/))
       cube_info(239) = cube_info_type (1,1,(/ 5, 8, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0/))
       cube_info(240) = cube_info_type (2,1,(/ 9,10,11,10,11,12, 0, 0, 0, 0, 0, 0/))
       cube_info(241) = cube_info_type (3,1,(/ 1, 4,11, 1,10,11,10,11,12, 0, 0, 0/))
       cube_info(242) = cube_info_type (3,1,(/ 1, 2, 9, 2, 9,12, 9,11,12, 0, 0, 0/))
       cube_info(243) = cube_info_type (2,1,(/ 2, 4,11, 2,11,12, 0, 0, 0, 0, 0, 0/))
       cube_info(244) = cube_info_type (3,1,(/ 2, 3,10, 3,10,11, 9,10,11, 0, 0, 0/))
       cube_info(245) = cube_info_type (2,1,(/ 1, 2,10, 3, 4,11, 0, 0, 0, 0, 0, 0/))
       cube_info(246) = cube_info_type (2,1,(/ 1, 3, 9, 3, 9,11, 0, 0, 0, 0, 0, 0/))
       cube_info(247) = cube_info_type (1,1,(/ 3, 4,11, 0, 0, 0, 0, 0, 0, 0, 0, 0/))
       cube_info(248) = cube_info_type (3,1,(/ 3, 4,12, 4, 9,12, 9,10,12, 0, 0, 0/))
       cube_info(249) = cube_info_type (2,1,(/ 1, 3,10, 3,10,12, 0, 0, 0, 0, 0, 0/))
       cube_info(250) = cube_info_type (2,1,(/ 1, 4, 9, 2, 3,12, 0, 0, 0, 0, 0, 0/))
       cube_info(251) = cube_info_type (1,1,(/ 2, 3,12, 0, 0, 0, 0, 0, 0, 0, 0, 0/))
       cube_info(252) = cube_info_type (2,1,(/ 2, 4, 9, 2, 9,10, 0, 0, 0, 0, 0, 0/))
       cube_info(253) = cube_info_type (1,1,(/ 1, 2,10, 0, 0, 0, 0, 0, 0, 0, 0, 0/))
       cube_info(254) = cube_info_type (1,1,(/ 1, 4, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0/))
       cube_info(255) = cube_info_type (0,0,(/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/))

       return
    End Subroutine Set_Cube_Info

    !!----
    !!---- Subroutine Statistic_Map(Rho,MaxV,MinV,AveV,SigmaV)
    !!----    real(kind=cp), dimension(:,:,:), intent(in) :: Rho
    !!----    real(kind=cp),                   intent(out):: MaxV      ! Maximum value of Rho
    !!----    real(kind=cp),                   intent(out):: MinV      ! Minimum value of Rho
    !!----    real(kind=cp),                   intent(out):: AveV      ! Average value of Rho
    !!----    real(kind=cp),                   intent(out):: SigmaV    ! Sigma value of Rho
    !!----
    !!----    Some statistic parameters of the map
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Statistic_Map(Rho,MaxV,MinV,AveV,SigmaV)
       !---- Arguments ----!
       real(kind=cp), dimension(:,:,:), intent(in) :: Rho
       real(kind=cp),                   intent(out):: MaxV
       real(kind=cp),                   intent(out):: MinV
       real(kind=cp),                   intent(out):: AveV
       real(kind=cp),                   intent(out):: SigmaV

       !---- Local Variables ----!
       integer :: i,j,k
       integer :: nu,nv,nw

       call init_err_maps()

       nu=size(rho,1)
       nv=size(rho,2)
       nw=size(rho,3)
       if (nu*nv*nw == 0) then
          err_maps=.true.
          ERR_Maps_Mess="Some dimension on Rho is zero"
          return
       end if

       MaxV  =maxval(rho)
       MinV  =minval(rho)
       AveV  = 0.0
       SigmaV= 0.0

       do i=1,nu
         do j=1,nv
             do k=1,nw
                avev=avev + rho(i,j,k)
                sigmav= sigmav + rho(i,j,k)*rho(i,j,k)
             end do
          end do
       end do
       avev  = avev/real(nu*nv*nw)
       sigmav= sigmav/real(nu*nv*nw) - avev*avev
       if (sigmav > 0.0001) then
          sigmav=sqrt(sigmav)
       else
          sigmav=0.0
       end if

       return
    End Subroutine Statistic_Map


 End Module CFML_Maps_Calculations
