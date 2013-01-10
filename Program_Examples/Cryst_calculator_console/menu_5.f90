!!----
!!---- Menu: 2
!!---- Reflections
!!----
!!
 Module Menu_5
   !---- Use File ----!
   use CFML_Crystal_Metrics,  only: Crystal_Cell_Type, Set_Crystal_Cell, Write_Crystal_Cell, &
                                    Zone_Axis_Type, Cart_Vector, Change_Setting_Cell, Get_Basis_From_UVW
   use CFML_GlobalDeps,       only: cp, Pi, Eps
   use CFML_Math_3D,          only: Determ_V,invert_a, determ_a, Cross_Product
   use CFML_Math_General,     only: cosd, norm, scalar,co_prime,co_Linear,equal_vector, acosd, sind, asind
   use CFML_String_Utilities, only: L_Case, pack_string, Get_Mat_From_Symb, Err_String, ERR_String_Mess


   !---- Variables ----!
   implicit none

   type (Crystal_Cell_type) :: Celda
   Logical                  :: cell_Given=.false.
   integer                  :: ierr
   Character(len=1)         :: Cart_type="A"

   Type, public :: Zone_Planes_Type
     integer                            :: nplanes
     real,   dimension(  :),allocatable :: mh
     integer,dimension(:,:),allocatable :: h     !co-prime indices
     real,   dimension(:,:),allocatable :: hc
   End Type Zone_Planes_Type

 Contains

    !!----
    !!---- Subroutine Menu_Princ5
    !!----
    !!
    Subroutine Menu_Princ5()
       !---- Local Variables ----!
       character(len=2)  :: car

       do
          call system('cls')

          write(unit=*,fmt="(a)") "     GENERAL CRYSTALLOGRAPHY CALCULATOR "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "       Geometrical Calculations Menu "
          write(unit=*,fmt="(a)") "   ======================================"
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "  [0] Back..."
          write(unit=*,fmt="(a)") "  "
          write(unit=*,fmt="(a)") "  [1] Enter unit cell parameters (or show current cell)"
          write(unit=*,fmt="(a)") "  [2] Zone axis common to two planes"
          write(unit=*,fmt="(a)") "  [3] Angle between two directions in direct space"
          write(unit=*,fmt="(a)") "  [4] Angle between two directions in reciprocal space"
          write(unit=*,fmt="(a)") "  [5] Angle between two directions one reciprocal space and the other in direct space"
          write(unit=*,fmt="(a)") "  [6] Direction in direct space => unitary direction in reciprocal space"
          write(unit=*,fmt="(a)") "  [7] Direction in reciprocal space => unitary direction in direct space"
          write(unit=*,fmt="(a)") "  [8] Basis change matrix => get new unit cell parameters"
          write(unit=*,fmt="(a)") "  [9] Zone axis: list of zone planes and angles,... special angles"
          write(unit=*,fmt="(a)") " [10] Indexing edges of a trapezoidal plane"
          write(unit=*,fmt="(a)") " [11] List of planes intersecting a given one at a particular Zone Axis"
          write(unit=*,fmt="(a)") " [12] Enter Euler Angles Chi, Phi, Theta of frame [u,v,w]"
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)",advance="no") " OPTION: "
          read(*,'(a)') car
          if (len_trim(car) == 0) exit
          car=adjustl(car)

          select case (car)
             case ('0 ')
                exit

             case ('1 ')
                call Menu_Geom_1()

             case ('2 ')
                call Menu_Geom_2()

             case ('3 ')
                call Menu_Geom_3()

             case ('4 ')
                call Menu_Geom_4()

             case ('5 ')
                call Menu_Geom_5()

             case ('6 ')
                call Menu_Geom_6()

             case ('7 ')
                call Menu_Geom_7()

             case ('8 ')
                call Menu_Geom_8()

             case ('9 ')
                call Menu_Geom_9()

             case ('10')
                call Menu_Geom_10()

             case ('11')
                call Menu_Geom_11()
             case ('12')

                call Menu_Geom_12()

          end select
       end do

    End Subroutine Menu_Princ5


    !!----
    !!---- Subroutine Menu_Geom_1
    !!----
    !!
    Subroutine Menu_Geom_1()
       !---- Local Variables ----!
       character(len=80)        :: line
       character(len=1)         :: car
       real, dimension(3)       :: Cellv,Angl

       do
          call system('cls')
          write(unit=*,fmt="(a)") "     GENERAL CRYSTALLOGRAPHY CALCULATOR "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " [0] Back..."
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " [1] Enter unit cell parameters "
          write(unit=*,fmt="(a)") " [2] Show current Cell type"
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)",advance="no") " OPTION: "
          read(unit=*,fmt='(a)') car
          if (len_trim(car) == 0) exit
          car=adjustl(car)
          Select Case (car)
            Case("0")
               Exit
            Case("1")
               write(unit=*,fmt="(a)",advance="no") " => Cartesian System choice (A-> x//a, C-> z//c), enter A or C: "
               read(unit=*,fmt="(a)") Cart_Type

               Select Case(Cart_Type)
                 Case("A","a","C","c")
                 Case Default
                   Cart_Type="C"
               End Select
               ierr=0
               write(unit=*,fmt="(a)",advance="no") " => Enter a,b,c,alpha,beta,gamma: "
               read(*,'(a)') line
               if (len_trim(line)==0) exit
               read(unit=line,fmt=*,iostat=ierr) Cellv,Angl
               if(ierr /= 0) then
                  write(unit=*,fmt="(a)") " -> Error in the unit cell parameters ... retry!"
                  call system("pause ")
                  cycle
               end if
               call Set_Crystal_Cell(Cellv,Angl,Celda,Cartype=Cart_Type)
               call Write_Crystal_Cell(Celda)
               cell_given=.true.
            Case("2")
               if(cell_given) then
                 call Write_Crystal_Cell(Celda)
               else
                  write(unit=*,fmt="(a)") " -> You have to provide first the cell parameters!"
               end if
          End Select
          write(unit=*,fmt=*) " "
          call system("pause ")
       end do
    End Subroutine Menu_Geom_1
    !!----
    !!---- Subroutine Menu_Geom_2
    !!----
    !!
    Subroutine Menu_Geom_2()
       !---- Local Variables ----!
       character(len=80)      :: line
       integer, dimension(3)  :: h1,h2,u

       do
          call system('cls')
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "        GENERAL CRYSTALLOGRAPHY CALCULATOR "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "          Zone axis common to two planes"
          write(unit=*,fmt="(a)") "    ==========================================="
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " "

          write(unit=*,fmt="(a)",advance="no") " => Enter the indices of the  first plane (<cr> for exit): "
          read(*,'(a)') line
          if (len_trim(line)==0) exit
          line=adjustl(line)
          read(unit=line,fmt=*,iostat=ierr) h1
          if(ierr /= 0) then
              write(unit=*,fmt="(a)") " -> Error in the indices ... retry!"
              call system("pause ")
              cycle
          end if
          write(unit=*,fmt="(a)",advance="no") " => Enter the indices of the second plane (<cr> for exit): "
          read(*,'(a)') line
          if (len_trim(line)==0) exit
          line=adjustl(line)
          read(unit=line,fmt=*,iostat=ierr) h2
          if(ierr /= 0) then
              write(unit=*,fmt="(a)") " -> Error in the indices ... retry!"
              call system("pause ")
              cycle
          end if
          u=Cross_Product(h1,h2)
          write(unit=*,fmt="(a)")  "    Indices Plane1  Indices Plane2     Zone Axis    "
          write(unit=*,fmt="(a)")  "    ==============  ==============    ============ "
          write(unit=*,fmt="(3(tr4,3i4))") h1,h2,u
          write(unit=*,fmt=*) " "
          call system("pause ")
       end do

    End Subroutine Menu_Geom_2
    !!----
    !!---- Subroutine Menu_Geom_3
    !!----
    !!
    Subroutine Menu_Geom_3()
       !---- Local Variables ----!
       character(len=80)      :: line
       real, dimension(3)     :: u,v,uc,vc,ur,vr
       real                   :: angle,mu,mv

       do
          call system('cls')
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "        GENERAL CRYSTALLOGRAPHY CALCULATOR "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "     Angle between two directions in direct space"
          write(unit=*,fmt="(a)") " ===================================================="
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " "

          if (Cell_Given) then
             write(unit=*,fmt="(a)",advance="no") " => Enter the indices of the first direction (<cr> for exit): "
             read(*,'(a)') line
             if (len_trim(line)==0) exit
             line=adjustl(line)
             read(unit=line,fmt=*,iostat=ierr) u
             if(ierr /= 0) then
                 write(unit=*,fmt="(a)") " -> Error in the indices ... retry!"
                 call system("pause ")
                 cycle
             end if
             write(unit=*,fmt="(a)",advance="no") " => Enter the indices of the second direction (<cr> for exit): "
             read(*,'(a)') line
             if (len_trim(line)==0) exit
             line=adjustl(line)
             read(unit=line,fmt=*,iostat=ierr) v
             if(ierr /= 0) then
                 write(unit=*,fmt="(a)") " -> Error in the indices ... retry!"
                 call system("pause ")
                 cycle
             end if
             uc=Cart_Vector("D",u,Celda)
             ur=Matmul(Celda%GD,u)
             vc=Cart_Vector("D",v,Celda)
             vr=Matmul(Celda%GD,v)
             mu=sqrt(dot_Product(uc,uc))
             mv=sqrt(dot_Product(vc,vc))
             if(mu < eps .or. mv < eps) then
                 write(unit=*,fmt="(a)") " -> One of the directions is [0 0 0] ... retry!"
                 call system("pause ")
                 cycle
             end if
             angle=dot_Product(uc,vc)/mu/mv
             if(angle > 1.0) angle=1.0
             if(angle < -1.0) angle=-1.0
             angle=acosd(angle)
             write(unit=*,fmt="(a)")  "               Direct space[A]         Cartesian components[e]    Reciprocal Components[A*]"
             write(unit=*,fmt="(a,2(3f8.4,tr2),3f9.3)") "     First:",u,uc,ur
             write(unit=*,fmt="(a,2(3f8.4,tr2),3f9.3)") "    Second:",v,vc,vr
             write(unit=*,fmt="(2(a, f10.5))")     "     Angle:",angle, " degrees  or  180-Angle:",180.0-angle
          else
             write(unit=*,fmt="(a)") " -> The unit cell must be given first !"
             call system("pause ")
             exit
          end if
          write(unit=*,fmt=*) " "
          call system("pause ")
       end do

    End Subroutine Menu_Geom_3

    !!----
    !!---- Subroutine Menu_Geom_4
    !!----
    !!
    Subroutine Menu_Geom_4()
       !---- Local Variables ----!
       character(len=80)      :: line
       real, dimension(3)     :: h,k,hc,kc,hd,kd
       real                   :: angle,mh,mk

       do
          call system('cls')
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "        GENERAL CRYSTALLOGRAPHY CALCULATOR "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "     Angle between two directions in reciprocal space"
          write(unit=*,fmt="(a)") " ======================================================="
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " "

          if (Cell_Given) then
             write(unit=*,fmt="(a)",advance="no") " => Enter the indices of the first direction (<cr> for exit): "
             read(*,'(a)') line
             if (len_trim(line)==0) exit
             line=adjustl(line)
             read(unit=line,fmt=*,iostat=ierr) h
             if(ierr /= 0) then
                 write(unit=*,fmt="(a)") " -> Error in the indices ... retry!"
                 call system("pause ")
                 cycle
             end if
             write(unit=*,fmt="(a)",advance="no") " => Enter the indices of the second direction (<cr> for exit): "
             read(*,'(a)') line
             if (len_trim(line)==0) exit
             line=adjustl(line)
             read(unit=line,fmt=*,iostat=ierr) k
             if(ierr /= 0) then
                 write(unit=*,fmt="(a)") " -> Error in the indices ... retry!"
                 call system("pause ")
                 cycle
             end if

             hc=Cart_Vector("R",h,Celda)
             !hd=Matmul(Celda%BL_M,h)  !identical for z//c -> BL
             !write(*,"(3f8.4,tr4,3f8.4)")  hc,hd
             hd=Matmul(Celda%GR,h)
             kc=Cart_Vector("R",k,Celda)
             kd=Matmul(Celda%GR,k)
             mh=sqrt(dot_Product(hc,hc))
             mk=sqrt(dot_Product(kc,kc))
             if(mh < eps .or. mk < eps) then
                 write(unit=*,fmt="(a)") " -> One of the directions is [0 0 0] ... retry!"
                 call system("pause ")
                 cycle
             end if
             angle=dot_Product(hc,kc)/mh/mk
             if(angle > 1.0) angle=1.0
             if(angle < -1.0) angle=-1.0
             angle=acosd(angle)
             write(unit=*,fmt="(a)")   "             Reciprocal space[A*]      Cartesian components[e]     Direct Components[A]"
             write(unit=*,fmt="(a,3(3f8.4,tr2))") "     First:",h,hc,hd
             write(unit=*,fmt="(a,3(3f8.4,tr2))") "    Second:",k,kc,kd
             write(unit=*,fmt="(2(a,f10.5))")     "     Angle:",angle, " degrees  or  180-Angle: ",180.0-Angle
          else
             write(unit=*,fmt="(a)") " -> The unit cell must be given first !"
             call system("pause ")
             exit
          end if
          write(unit=*,fmt=*) " "
          call system("pause ")
       end do

    End Subroutine Menu_Geom_4

    !!----
    !!---- Subroutine Menu_Geom_5
    !!----
    !!
    Subroutine Menu_Geom_5()
       !---- Local Variables ----!
       character(len=80)      :: line
       real, dimension(3)     :: h,u,hc,uc,hd,ur
       real                   :: angle,mh,mu

       do
          call system('cls')
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "        GENERAL CRYSTALLOGRAPHY CALCULATOR "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "     Angle between two directions one in reciprocal space and the other in direct space"
          write(unit=*,fmt="(a)") " ==========================================================================================="
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " "

          if (Cell_Given) then
             write(unit=*,fmt="(a)",advance="no") " => Enter the indices of the reciprocal direction (<cr> for exit): "
             read(*,'(a)') line
             if (len_trim(line)==0) exit
             line=adjustl(line)
             read(unit=line,fmt=*,iostat=ierr) h
             if(ierr /= 0) then
                 write(unit=*,fmt="(a)") " -> Error in the indices ... retry!"
                 call system("pause ")
                 cycle
             end if
             write(unit=*,fmt="(a)",advance="no") " => Enter the indices of the direct direction (<cr> for exit): "
             read(*,'(a)') line
             if (len_trim(line)==0) exit
             line=adjustl(line)
             read(unit=line,fmt=*,iostat=ierr) u
             if(ierr /= 0) then
                 write(unit=*,fmt="(a)") " -> Error in the indices ... retry!"
                 call system("pause ")
                 cycle
             end if

             hc=Cart_Vector("R",h,Celda)
             hd=Matmul(Celda%GR,h)
             uc=Cart_Vector("D",u,Celda)
             ur=Matmul(Celda%GD,u)
             mh=sqrt(dot_Product(hc,hc))
             mu=sqrt(dot_Product(uc,uc))
             if(mh < eps .or. mu < eps) then
                 write(unit=*,fmt="(a)") " -> One of the directions is [0 0 0] ... retry!"
                 call system("pause ")
                 cycle
             end if
             angle=dot_Product(h,u)/mh/mu
             if(angle > 1.0) angle=1.0
             if(angle < -1.0) angle=-1.0
             angle=acosd(angle)
             write(unit=*,fmt="(a)")  "             Reciprocal space[A*]      Cartesian components[e]     Direct Components[A]"
             write(unit=*,fmt="(a,3(3f8.4,tr2))") "  h- First:",h,hc,hd
             write(unit=*,fmt="(a,3(3f8.4,tr2))") "  u-Second:",ur,uc,u
             write(unit=*,fmt="(2(a,f10.5))")     "     Angle:",angle, " degrees  or 180-Angle:",180.0-angle
          else
             write(unit=*,fmt="(a)") " -> The unit cell must be given first !"
             call system("pause ")
             exit
          end if
          write(unit=*,fmt=*) " "
          call system("pause ")
       end do

    End Subroutine Menu_Geom_5

    !!----
    !!---- Subroutine Menu_Geom_6
    !!----
    !!
    subroutine Menu_Geom_6()
       !---- Local Variables ----!
       character(len=80)      :: line
       real, dimension(3)     :: u,uc,vc,ur
       real                   :: mu,mv

       do
          call system('cls')
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "        GENERAL CRYSTALLOGRAPHY CALCULATOR "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "     Direction in direct space => unitary direction in reciprocal space"
          write(unit=*,fmt="(a)") " ==========================================================================="
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " "

          if (Cell_Given) then
             write(unit=*,fmt="(a)",advance="no") " => Enter the indices of direction (<cr> for exit): "
             read(*,'(a)') line
             if (len_trim(line)==0) exit
             line=adjustl(line)
             read(unit=line,fmt=*,iostat=ierr) u
             if(ierr /= 0) then
                 write(unit=*,fmt="(a)") " -> Error in the indices ... retry!"
                 call system("pause ")
                 cycle
             end if

             uc=Cart_Vector("D",u,Celda)
             mu=sqrt(dot_Product(uc,uc))
             ur=Matmul(Celda%GD,u)
             vc=Cart_Vector("R",ur,Celda)
             mv=sqrt(dot_Product(vc,vc))
             if(mu < eps) then
                 write(unit=*,fmt="(a)") " -> The directions is [0 0 0] ... retry!"
                 call system("pause ")
                 cycle
             end if
             uc=uc/mu
             ur=ur/mv
             write(unit=*,fmt="(2(a,3f9.4),a)") "                        [uvw] : [",u,"]  -> Unit vector: [",u/mu,"]"
             write(unit=*,fmt="(a,3f9.4,a)")    "   Cartesian Unit vector[uvw]c: [",uc,"]"
             write(unit=*,fmt="(a,3f9.4,a)")    "  Reciprocal Unit vector[hkl]*: [",ur,"]*"

          else
             write(unit=*,fmt="(a)") " -> The unit cell must be given first !"
             call system("pause ")
             exit
          end if
          write(unit=*,fmt=*) " "
          call system("pause ")
       end do


    End Subroutine Menu_Geom_6

    !!----
    !!---- Subroutine Menu_Geom_7
    !!----
    !!
    Subroutine Menu_Geom_7()
       !---- Local Variables ----!
       character(len=80)      :: line
       real, dimension(3)     :: u,v,uc,vc
       real                   :: mu,mv

       do
          call system('cls')
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "        GENERAL CRYSTALLOGRAPHY CALCULATOR "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "     Direction in reciprocal space => unitary direction in direct space"
          write(unit=*,fmt="(a)") " ==========================================================================="
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " "

          if (Cell_Given) then
             write(unit=*,fmt="(a)",advance="no") " => Enter the indices of direction (<cr> for exit): "
             read(*,'(a)') line
             if (len_trim(line)==0) exit
             line=adjustl(line)
             read(unit=line,fmt=*,iostat=ierr) u
             if(ierr /= 0) then
                 write(unit=*,fmt="(a)") " -> Error in the indices ... retry!"
                 call system("pause ")
                 cycle
             end if

             uc=Cart_Vector("R",u,Celda)
             mu=sqrt(dot_Product(uc,uc))
             v=Matmul(Celda%GR,u)
             vc=Cart_Vector("D",v,Celda)
             mv=sqrt(dot_Product(vc,vc))
             if(mu < eps) then
                 write(unit=*,fmt="(a)") " -> The directions is [0 0 0] ... retry!"
                 call system("pause ")
                 cycle
             end if
             uc=uc/mu
             v=v/mv
             write(unit=*,fmt="(2(a,3f9.4),a)") "                        [hkl]*: [",u,"]  -> Unit vector: [",u/mu,"]"
             write(unit=*,fmt="(a,3f9.4,a)")    "   Cartesian Unit vector[hkl]c: [",uc,"]"
             write(unit=*,fmt="(a,3f9.4,a)")    "      Direct Unit vector[uvw] : [",v,"]*"

          else
             write(unit=*,fmt="(a)") " -> The unit cell must be given first !"
             call system("pause ")
             exit
          end if
          write(unit=*,fmt=*) " "
          call system("pause ")
       end do
    End Subroutine Menu_Geom_7

    !!----
    !!---- Subroutine Menu_Geom_8
    !!----
    !!
    subroutine Menu_Geom_8()
       !---- Local Variables ----!
       character(len=80)        :: line,symb
       real, dimension(3,3)     :: mat
       type (Crystal_Cell_type) :: Celln
       integer                  :: i

       do
          call system('cls')
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "        GENERAL CRYSTALLOGRAPHY CALCULATOR "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "     Basis change matrix => get new unit cell parameters"
          write(unit=*,fmt="(a)") " =========================================================="
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " "

          if (Cell_Given) then
             write(unit=*,fmt="(a)",advance="no") " => Enter the transformation (<cr> for exit), e.g. 2a+b,c,-b: "
             read(*,'(a)') line
             if (len_trim(line)==0) exit
             line=L_Case(adjustl(line))
             Symb=line
             call Get_Mat_From_Symb(Symb,Mat,(/"a","b","c"/))
             if(Err_String) then
               write(unit=*,fmt="(a)") " => "//trim(Err_String_Mess)
               call system("pause ")
               cycle
             end if
             write(unit=*,fmt="(a)") " => Matrix corresponding to trasformation: "//trim(Symb)
             do i=1,3
                write(unit=*,fmt="(a,3f9.5)") "         ",Mat(i,:)
             end do
             call Change_Setting_Cell(Celda,Mat,Celln)
             write(unit=*,fmt="(/,a)") " ============================"
             write(unit=*,fmt="(a)")   " => New metric information <="
             write(unit=*,fmt="(a)")   " ============================"
             call Write_Crystal_Cell(Celln)
          else
             write(unit=*,fmt="(a)") " -> The unit cell must be given first !"
             call system("pause ")
             exit
          end if
          write(unit=*,fmt=*) " "
          call system("pause ")
       end do
    End Subroutine Menu_Geom_8


    subroutine Menu_Geom_9()
       !---- Local Variables ----!
       character(len=80)      :: line
       real,    dimension(3)  :: uc,vc
       integer, dimension(3)  :: u
       real,    dimension(:),allocatable  :: angs
       real,    dimension(:),allocatable  :: angles
       integer, dimension(:),allocatable  :: ii,jj,nang
       real                   :: dmin,angle,mu,mv,ca, tol
       integer                :: i,j,k,n,n_angles
       Type(Zone_Axis_Type)   :: Zone_Axis
       Type(Zone_Planes_Type) :: Zone_Planes
       Logical                :: ok

       do
          call system('cls')
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "                   GENERAL CRYSTALLOGRAPHY CALCULATOR "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "      Zone axis: list of zone planes and angles, special angles"
          write(unit=*,fmt="(a)") " ===================================================================="
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " "

          if (Cell_Given) then
             write(unit=*,fmt="(a)",advance="no") " => Enter the indices of direction (<cr> for exit): "
             read(*,'(a)') line
             if (len_trim(line)==0) exit
             line=adjustl(line)
             read(unit=line,fmt=*,iostat=ierr) u
             if(ierr /= 0) then
                 write(unit=*,fmt="(a)") " -> Error in the indices ... retry!"
                 call system("pause ")
                 cycle
             end if
             write(unit=*,fmt="(a)",advance="no") " => Enter the number of special inter-face angles to match (<cr> for exit): "
             read(*,*,iostat=ierr) n_angles
             if(ierr /= 0) then
                 write(unit=*,fmt="(a)") " -> Error in n_angles ... retry!"
                 call system("pause ")
                 cycle
             end if
             if(n_angles > 0) then
               if(allocated(angs)) deallocate(angs)
               allocate(angs(n_angles))
               do i=1,n_angles
                 write(unit=*,fmt="(a,i2,a)",advance="no") " => Enter the special angle : ",i,": "
                 read(*,*) angs(i)
               end do
               write(unit=*,fmt="(a)",advance="no") " => Angular tolerance in degrees (<cr> for exit): "
               read(unit=*,fmt=*,iostat=ierr) tol
               if(ierr /= 0) tol=1.0
               if(tol < 1.0) tol=1.0
             end if
             write(unit=*,fmt="(a)",advance="no") " => Enter the minimum d-spacing of planes (<cr> for exit): "
             read(*,'(a)') line
             if (len_trim(line)==0) exit
             line=adjustl(line)
             read(unit=line,fmt=*,iostat=ierr) dmin
             if(ierr /= 0) then
                 write(unit=*,fmt="(a)") " -> Error in dmin ... retry!"
                 call system("pause ")
                 cycle
             end if

             call Get_basis_from_uvw(dmin,u,celda,Zone_Axis,ok)

             if(.not. ok) then
                write(unit=*,fmt="(a,3i3,a)") " => Problem in getting the basis in the reciprocal plane perpendicular to [",u,"]"
                cycle
             end if

             call Get_Zone_Planes(Zone_Axis%rx,Zone_Axis%ry,dmin,Celda,Zone_Planes)
             write(unit=*,fmt="(a)")  "    Face indices    Face indices     Normal Angle     Face Angle"
             write(unit=*,fmt="(a)")  "    ============    ============     ============     =========="
             ! Calculation of all interplane angles
             if(n_angles > 0) then
               n=Zone_Planes%nplanes*(Zone_Planes%nplanes-1)/2
               if(allocated(angles)) deallocate(angles)
               if(allocated(ii)) deallocate(ii)
               if(allocated(jj)) deallocate(jj)
               if(allocated(nang)) deallocate(nang)
               allocate(angles(n),ii(n),jj(n),nang(n))
             end if
             n=0
             do i=1,Zone_Planes%nplanes-2
                uc=Zone_Planes%hc(:,i)
                mu=Zone_Planes%mh(i)
                do j=i+1,Zone_Planes%nplanes
                  vc=Zone_Planes%hc(:,j)
                  mv=Zone_Planes%mh(j)
                  angle=angle_val(uc,mu,vc,mv)
                  ca=180.0-angle
                  write(unit=*,fmt="(2(tr4,3i4),2(tr6,f10.4))") Zone_Planes%h(:,i),Zone_Planes%h(:,j),angle,ca
                  do k=1,n_angles
                    if(abs(ca-angs(k)) < tol) then
                      n=n+1
                      angles(n)=ca
                      ii(n)=i
                      jj(n)=j
                      nang(n)=k
                    end if
                  end do
                end do
             end do
             if(n_angles > 0) then
               write(unit=*,fmt="(/,a)")"    Special angles to compare with calculations"
               write(unit=*,fmt="(a)")  "    ============    ============     ==============     ================"
               write(unit=*,fmt="(a)")  "    Face indices    Face indices     Observed Angle     Calculated Angle"
               write(unit=*,fmt="(a)")  "    ============    ============     ==============     ================"
               do i=1,n
                 k=nang(i)
                 write(unit=*,fmt="(2(tr4,3i4),2(tr7,f10.4))") Zone_Planes%h(:,ii(i)),Zone_Planes%h(:,jj(i)),angs(k),angles(i)
               end do
             end if
          else
             write(unit=*,fmt="(a)") " -> The unit cell must be given first !"
             call system("pause ")
             exit
          end if
          write(unit=*,fmt=*) " "
          call system("pause ")
       end do
    End Subroutine Menu_Geom_9

    Subroutine Menu_Geom_10()
       !---- Local Variables ----!
       character(len=80)      :: line
       real,    dimension(3)  :: uc,vc,h1,h2,h3,h4
       real,    dimension(4)  :: angs,abest
       integer, dimension(3,8):: hb
       real,    dimension(:),allocatable  :: angles
       integer, dimension(:),allocatable  :: ii,jj,nang
       integer, dimension(3)  :: u,h,hi2,ho2,hi4,ho4
       real                   :: n_max,angle,mv1,mv2,mv3,mv4, tol, a13,a14,a23,a24, &
                                 a1,a2,a3,a4, mu,mv,rf, rmin
       integer                :: i,j,k,nedges,n,msol,nsol, nbest,neq
       Type(Zone_Axis_Type)   :: Zone_Axis
       Type(Zone_Planes_Type) :: Zone_Planes
       Logical                :: ok,v21,v41,v22,v42,vv

       do
          call system('cls')
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "        GENERAL CRYSTALLOGRAPHY CALCULATOR "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "          Indexing edges of trapezoids"
          write(unit=*,fmt="(a)") " ==============================================="
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " "

          if (Cell_Given) then
             write(unit=*,fmt="(a)",advance="no") " => Enter the indices of plane (<cr> for exit): "
             read(*,'(a)') line
             if (len_trim(line)==0) exit
             line=adjustl(line)
             read(unit=line,fmt=*,iostat=ierr) h
             if(ierr /= 0) then
                 write(unit=*,fmt="(a)") " -> Error in the indices ... retry!"
                 call system("pause ")
                 cycle
             end if
             nedges=4
             do i=1,nedges-1
               j=i+1
               write(unit=*,fmt="(2(a,i2),a)",advance="no") " => Enter the angle between edge: ",i," and edge ",j,": "
               read(*,*) angs(i)
             end do
             angs(4)=360.0-angs(1)-angs(2)-angs(3)

             write(unit=*,fmt="(a)",advance="no") " => Enter the maximun distance in angstrom for n_uvw (<cr> for exit): "
             read(*,'(a)') line
             if (len_trim(line)==0) exit
             line=adjustl(line)
             read(unit=line,fmt=*,iostat=ierr) n_max
             if(ierr /= 0) then
                 write(unit=*,fmt="(a)") " -> Error in n_max ... retry!"
                 call system("pause ")
                 cycle
             end if

             write(unit=*,fmt="(a)",advance="no") " => Angular tolerance in degrees (<cr> for exit): "
             read(unit=*,fmt=*,iostat=ierr) tol
             if(ierr /= 0) tol=1.0
             if(tol < 1.0) tol=1.0

             call Get_basis_from_uvw(n_max,h,celda,Zone_Axis,ok,"R")
             if(.not. ok) then
                write(unit=*,fmt="(a,3i3,a)") " => Problem in getting the basis in the direct plane perpendicular to [",h,"]*"
                cycle
             end if
             write(unit=*,fmt="(/,a,3i4,a)") "    Basis vector-a1: [",Zone_Axis%rx," ]"
             write(unit=*,fmt="(a,3i4,a)")   "    Basis vector-a2: [",Zone_Axis%ry," ]"


             call Get_Zone_Planes(Zone_Axis%rx,Zone_Axis%ry,n_max,Celda,Zone_Planes,"R")

           !  write(unit=*,fmt="(a)")  "    Edge indices    Edge indices     Normal Angle    Compl. Angle"
           !  write(unit=*,fmt="(a)")  "    ============    ============     ============     ==========="
             n=Zone_Planes%nplanes*(Zone_Planes%nplanes-1)
             if(allocated(ii)) deallocate(ii)
             if(allocated(jj)) deallocate(jj)
             if(allocated(angles)) deallocate(angles)
             if(allocated(nang)) deallocate(nang)
             allocate(ii(n),jj(n),angles(n),nang(n))
             ii=0; jj=0; angles=0.0;nang=0
             ! Calculation of all interplane angles
             msol=0
             do i=1,Zone_Planes%nplanes-2
                h1=Zone_Planes%h(:,i)
                uc=Zone_Planes%hc(:,i)
                mu=Zone_Planes%mh(i)
                do j=i+1,Zone_Planes%nplanes
                  h2=Zone_Planes%h(:,j)
                  vc=Zone_Planes%hc(:,j)
                  mv=Zone_Planes%mh(j)
                  angle=angle_val(uc,mu,vc,mv)
                  edges: do n=1,nedges
                    !if(abs(angle-angs(n)) < tol .or. abs(180.0-angle-angs(n)) < tol) then
                    if(abs(angle-angs(n)) < tol ) then
                      do k=1, msol
                        h3=Zone_Planes%h(:,ii(k))
                        h4=Zone_Planes%h(:,jj(k))
                        if( (equal_vector(h1,h3,3) .and. equal_vector(h2,h4,3)) .or.  &
                        (equal_vector(h1,h4,3) .and. equal_vector(h2,h3,3)) ) cycle edges
                      end do
                      msol=msol+1
                   !   write(unit=*,fmt="(2(tr4,3i4),2(tr6,f10.4),a,i5)") Zone_Planes%h(:,i),Zone_Planes%h(:,j),&
                   !                                                      angle,180-angle,"  #", msol
                      ii(msol)=i
                      jj(msol)=j
                      angles(msol)=angle
                      nang(msol) = n
                    end if
                  end do edges
                end do
             end do


              nsol=0
              rmin=999999.0
              do i=1,msol
                if(nang(i) /= 1) cycle
                h1= Zone_Planes%hc(:,ii(i))
                h2= Zone_Planes%hc(:,jj(i))
                mv1=Zone_Planes%mh(ii(i))
                mv2=Zone_Planes%mh(jj(i))
                a1=angle_val(h1,mv1,h2,mv2)
                do j=1,msol
                  if( i == j) cycle
                  if(nang(j) /= 3) cycle
                  !Test h3
                  h3=Zone_Planes%hc(:,ii(j))
                  mv3=Zone_Planes%mh(ii(j))

                  a13=angle_val(h1,mv1,h3,mv3)
                  v21=abs(a13-angs(2)) < tol
                  v41=abs(a13-angs(4)) < tol

                  a23=angle_val(h2,mv2,h3,mv3)
                  v22=abs(a23-angs(2)) < tol
                  v42=abs(a23-angs(4)) < tol

                  ok= v21 .or. v41 .or. v22 .or. v42
                  if(.not. ok) cycle

                  !Now knowing that h3 works test h4
                  !first determine with which vector we have to compare

                  h4=Zone_Planes%hc(:,jj(j))
                  mv4=Zone_Planes%mh(jj(j))
                  a3=angle_val(h3,mv3,h4,mv4)

                  if(v21) then
                    angle=angle_val(h4,mv4,h2,mv2)
                    vv=abs(angle-angs(4)) < tol
                    hi2=Zone_Planes%h(:,ii(i))
                    ho2=Zone_Planes%h(:,ii(j))
                    hi4=Zone_Planes%h(:,jj(i))
                    ho4=Zone_Planes%h(:,jj(j))
                    a2=a13
                    a4=angle
                  else if(v41) then
                    angle=angle_val(h4,mv4,h2,mv2)
                    vv=abs(angle-angs(2)) < tol
                    hi2=Zone_Planes%h(:,jj(i))
                    ho2=Zone_Planes%h(:,jj(j))
                    hi4=Zone_Planes%h(:,ii(i))
                    ho4=Zone_Planes%h(:,ii(j))
                    a2=angle
                    a4=a13
                  else if(v22) then
                    angle=angle_val(h4,mv4,h1,mv1)
                    vv=abs(angle-angs(4)) < tol
                    hi2=Zone_Planes%h(:,jj(i))
                    ho2=Zone_Planes%h(:,ii(j))
                    hi4=Zone_Planes%h(:,ii(i))
                    ho4=Zone_Planes%h(:,jj(j))
                    a2=a23
                    a4=angle
                  else if(v42) then
                    angle=angle_val(h4,mv4,h1,mv1)
                    vv=abs(angle-angs(2)) < tol
                    hi2=Zone_Planes%h(:,ii(i))
                    ho2=Zone_Planes%h(:,jj(j))
                    hi4=Zone_Planes%h(:,jj(i))
                    ho4=Zone_Planes%h(:,ii(j))
                    a2=angle
                    a4=a23
                 end if
                  if(.not. vv) cycle
                  !A solution has been found
                  nsol=nsol+1
                  rf=abs(a1-angs(1))+abs(a2-angs(2))+abs(a3-angs(3))+abs(a4-angs(4))
                  rf=100.0*rf/360.0
                  if(rf < rmin-0.001) then
                    rmin=rf
                    neq=0
                    nbest=nsol
                    hb(:,1)=Zone_Planes%h(:,ii(i))
                    hb(:,2)=Zone_Planes%h(:,jj(i))
                    hb(:,3)=hi2;  hb(:,4)=ho2
                    hb(:,5)=Zone_Planes%h(:,ii(j))
                    hb(:,6)=Zone_Planes%h(:,jj(j))
                    hb(:,7)=hi4;  hb(:,8)=ho4
                    abest=(/a1,a2,a3,a4/)
                  else if(abs(rf - rmin) < 0.001) then
                    neq=neq+1
                  end if
                  write(unit=*,fmt="(/,a,i5,tr4,3i3,a,3i3,f10.3,a,f7.3)") "   Sol# ",nsol, &
                          Zone_Planes%h(:,ii(i))," --- ",Zone_Planes%h(:,jj(i)),a1,"   R-factor(%)=",rf
                  write(unit=*,fmt="(a,tr9,3i3,a,3i3,f10.3)") "        ",  hi2," --- ",ho2,a2
                  write(unit=*,fmt="(a,tr9,3i3,a,3i3,f10.3)") "        ", &
                          Zone_Planes%h(:,ii(j))," --- ",Zone_Planes%h(:,jj(j)),a3
                  write(unit=*,fmt="(a,tr9,3i3,a,3i3,f10.3)") "        ",  hi4," --- ",ho4,a4
                end do
              end do
              if(nsol > 0) then
                write(unit=*,fmt="(/,a,i3,a,f7.3,/,a,i2,a)") " -> The best solution is that of number: ",nbest, &
                                    " of R-Factor(%)=",rmin, "    and ",neq, " additional equivalent solutions"
                write(unit=*,fmt="(a)")   "           u  v  w  ---  u  v  w   Angle(Obs)  Angle(Cal)    Obs-Cal"
                j=1
                do i=1,4
                  write(unit=*,fmt="(a,tr6,3i3,a,3i3,3f12.3)") "   ",  &
                      hb(:,j)," --- ",hb(:,j+1),angs(i),abest(i),angs(i)-abest(i)
                  j=j+2
                end do
              else
                write(unit=*,fmt="(a)") " => NO solution! ... increase tolerance?"
                call system("pause ")
                cycle
              end if
          else
             write(unit=*,fmt="(a)") " -> The unit cell must be given first !"
             call system("pause ")
             exit
          end if
          write(unit=*,fmt=*) " "
          call system("pause ")
       end do
    End Subroutine Menu_Geom_10

    !!----
    !!---- Subroutine Menu_Geom_11
    !!----
    !!
    Subroutine Menu_Geom_11()
       !---- Local Variables ----!
       character(len=80)      :: line
       integer, dimension(3)  :: h1,h2,u,uc
       real,    dimension(3)  :: c1,c2
       integer                :: i,j,k,n,hmax,kmax,lmax
       real                   :: ang,ang1,ang2,mc1,mc2,dmin

       do
          call system('cls')
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "                  GENERAL CRYSTALLOGRAPHY CALCULATOR "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "    List of planes intersecting a given one at a particular Zone Axis"
          write(unit=*,fmt="(a)") "   Forming an angle with the given plane within the interval [ang1,ang2]"
          write(unit=*,fmt="(a)") "    ================================================================="
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " "

          write(unit=*,fmt="(a)",advance="no") " => Enter the indices (hkl) of the Plane (<cr> for exit): "
          read(*,'(a)') line
          if (len_trim(line)==0) exit
          line=adjustl(line)
          read(unit=line,fmt=*,iostat=ierr) h1
          if(ierr /= 0) then
              write(unit=*,fmt="(a)") " -> Error in the indices ... retry!"
              call system("pause ")
              cycle
          end if
          write(unit=*,fmt="(a)",advance="no") " => Enter the indices [uvw] of the Zone Axis (<cr> for exit): "
          read(*,'(a)') line
          if (len_trim(line)==0) exit
          line=adjustl(line)
          read(unit=line,fmt=*,iostat=ierr) u
          if(ierr /= 0) then
              write(unit=*,fmt="(a)") " -> Error in the indices ... retry!"
              call system("pause ")
              cycle
          end if

          write(unit=*,fmt="(a)",advance="no") " => Enter the angular range to consider (<cr> for exit): "
          read(*,'(a)') line
          if (len_trim(line)==0) exit
          line=adjustl(line)
          read(unit=line,fmt=*,iostat=ierr) ang1,ang2
          if(ierr /= 0) then
              write(unit=*,fmt="(a)") " -> Error in ang1-ang2 ... retry!"
              call system("pause ")
              cycle
          end if

          write(unit=*,fmt="(a)",advance="no") " => Enter the minimum d-spacing of planes (<cr> for exit): "
          read(*,'(a)') line
          if (len_trim(line)==0) exit
          line=adjustl(line)
          read(unit=line,fmt=*,iostat=ierr) dmin
          if(ierr /= 0) then
              write(unit=*,fmt="(a)") " -> Error in dmin ... retry!"
              call system("pause ")
              cycle
          end if
          hmax=nint(Celda%cell(1)/dmin+1.0)
          kmax=nint(Celda%cell(2)/dmin+1.0)
          lmax=nint(Celda%cell(3)/dmin+1.0)
          c1=Cart_Vector("R",real(h1),Celda)
          mc1=sqrt(dot_Product(c1,c1))
          write(unit=*,fmt="(/,2(a,3i3),a)") " => List of planes intersecting (",h1,") at the zone axis: [",u,"]"
          n=0
          do i=0,kmax
            do j= kmax,-kmax,-1
              do k= lmax,-lmax,-1
                if(i == 0 .and. j == 0 .and. k == 0) cycle
                h2=(/i,j,k/)
                if(.not. co_prime(h2)) cycle
                uc=Cross_Product(h1,h2)
                if(co_linear(uc,u,3)) then
                  !Calculate the angle with the initial plane
                  c2=Cart_Vector("R",real(h2),Celda)
                  mc2=sqrt(dot_Product(c2,c2))
                  ang=180.0-Angle_val(c1,mc1,c2,mc2)
                  if(ang >= ang1 .and. ang <= ang2) then
                    n=n+1
                    write(unit=*,fmt="(a,i5,a,3i4,a,f9.3,a)") "    #",n,"  (hkl):",h2,"  Angle:",ang," degrees"
                  end if
                end if
              end do
            end do
          end do
          write(unit=*,fmt=*) " "
          call system("pause ")
       end do

    End Subroutine Menu_Geom_11

    Function Angle_val(h1,m1,h2,m2) result(angle)
     real, dimension(:), intent(in) :: h1,h2
     real,               intent(in) :: m1,m2
     real                           :: angle
     angle=dot_product(h1,h2)/m1/m2
     if(angle > 1.0) angle=1.0
     if(angle <-1.0) angle=-1.0
     angle=acosd(angle)
    End Function Angle_val

    !!----
    !!---- Subroutine Menu_Geom_12
    !!----
    !!
    Subroutine Menu_Geom_12()
       !---- Local Variables ----!
       character(len=80)      :: line
       real,    dimension(3)  :: u,v,w, h,c
       real                   :: angu,angv,angw,Chi,Phi,Theta,mc

       do
          call system('cls')
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "                  GENERAL CRYSTALLOGRAPHY CALCULATOR "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "           Enter Euler angles to get new Cartesian system u,v,w"
          write(unit=*,fmt="(a)") "    Calculate the angles of u,v,w with a direct and a reciprocal directions"
          write(unit=*,fmt="(a)") "    ======================================================================="
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " "

          write(unit=*,fmt="(a)",advance="no") " => Enter the Angles Chi, Phi, Theta (<cr> for exit): "
          read(*,'(a)') line
          if (len_trim(line)==0) exit
          line=adjustl(line)
          read(unit=line,fmt=*,iostat=ierr) Chi, Phi, Theta
          if(ierr /= 0) then
              write(unit=*,fmt="(a)") " -> Error in the angles ... retry!"
              call system("pause ")
              cycle
          end if
          ! Calculate the unitary vectors of the Cartesian system
           u(1)=cosd(phi)*cosd(Theta)*cosd(Chi)-sind(Phi)*sind(chi)
           u(2)=sind(phi)*cosd(Theta)*cosd(Chi)+cosd(Phi)*sind(chi)
           u(3)=-sind(Theta)*cosd(Chi)
           v(1)=-cosd(phi)*cosd(Theta)*sind(Chi)-sind(Phi)*cosd(chi)
           v(2)=-sind(phi)*cosd(Theta)*sind(Chi)+cosd(Phi)*cosd(chi)
           v(3)= sind(Theta)*sind(Chi)
           w=(/ cosd(Phi)*sind(Theta), sind(Phi)*sind(Theta),  cosd(Theta) /)
          Write(unit=*,fmt="(a)") " => The unit vectors of the reference system (u,v,w),"
          Write(unit=*,fmt="(a)") "    given in the Standard Cartesian system, are:"
          Write(unit=*,fmt="(a,3f10.5,a)") " => u =(",u,")"
          Write(unit=*,fmt="(a,3f10.5,a)") " => v =(",v,")"
          Write(unit=*,fmt="(a,3f10.5,a)") " => w =(",w,")"

          !Calculate the angles with a direct space direction

          write(unit=*,fmt="(a)",advance="no") " => Enter the components of a direct space vector (<cr> for exit): "
          read(*,'(a)') line
          if (len_trim(line)==0) exit
          line=adjustl(line)
          do
            read(unit=line,fmt=*,iostat=ierr) h
            if(ierr /= 0) then
                write(unit=*,fmt="(a)") " -> Error in the components ... retry!"
                call system("pause ")
                cycle
            else
                exit
            end if
          end do
          c=Cart_Vector("D",h,Celda)
          mc=sqrt(dot_Product(c,c))
          c=c/mc !Now I have a unitary vector
          angu=dot_product(c,u); if(angu > 1.0) angu=1.0; if(angu < -1.0) angu=-1.0
          angu=acosd(angu)
          angv=dot_product(c,v); if(angv > 1.0) angv=1.0; if(angv < -1.0) angv=-1.0
          angv=acosd(angv)
          angw=dot_product(c,w); if(angw > 1.0) angw=1.0; if(angw < -1.0) angw=-1.0
          angw=acosd(angw)
          write(unit=*,fmt="(a,3f10.5,a)") " => The angles of the direction [",h,"] with vectors u,v & w are:"
          Write(unit=*,fmt="(2(a,f10.5))") "    Angle(u,h) =",angu," Complementary:",180.0-angu
          Write(unit=*,fmt="(2(a,f10.5))") "    Angle(v,h) =",angv," Complementary:",180.0-angv
          Write(unit=*,fmt="(2(a,f10.5))") "    Angle(w,h) =",angw," Complementary:",180.0-angw

          !Calculate the angles with a reciprocal space direction
          write(unit=*,fmt="(a)",advance="no") " => Enter the components of a reciprocal space vector (<cr> for exit): "
          read(*,'(a)') line
          if (len_trim(line)==0) exit
          line=adjustl(line)
          do
            read(unit=line,fmt=*,iostat=ierr) h
            if(ierr /= 0) then
                write(unit=*,fmt="(a)") " -> Error in the components ... retry!"
                call system("pause ")
                cycle
            else
                exit
            end if
          end do

          c=Cart_Vector("R",h,Celda)
          mc=sqrt(dot_Product(c,c))
          c=c/mc !Now I have a unitary vector
          angu=dot_product(c,u); if(angu > 1.0) angu=1.0; if(angu < -1.0) angu=-1.0
          angu=acosd(angu)
          angv=dot_product(c,v); if(angv > 1.0) angv=1.0; if(angv < -1.0) angv=-1.0
          angv=acosd(angv)
          angw=dot_product(c,w); if(angw > 1.0) angw=1.0; if(angw < -1.0) angw=-1.0
          angw=acosd(angw)
          write(unit=*,fmt="(a,3f10.5,a)") " => The angles of the direction [",h,"]* with vectors u,v & w are:"
          Write(unit=*,fmt="(2(a,f10.5))") "    Angle(u,h) =",angu," Complementary:",180.0-angu
          Write(unit=*,fmt="(2(a,f10.5))") "    Angle(v,h) =",angv," Complementary:",180.0-angv
          Write(unit=*,fmt="(2(a,f10.5))") "    Angle(w,h) =",angw," Complementary:",180.0-angw
          write(unit=*,fmt=*) " "
          call system("pause ")
       end do

    End Subroutine Menu_Geom_12



    Subroutine Menu_Geom_12a()
       !---- Local Variables ----!
       character(len=80)      :: line
       real,    dimension(3)  :: uc,vc
       real,    dimension(:),allocatable  :: angs,angles
       integer, dimension(:),allocatable  :: ii,jj,nn,nang,nl
       integer, dimension(:,:),allocatable  :: p
       integer, dimension(3)  :: u,h,h1,h2,h3,h4,hi,ho
       real                   :: n_max,angle,mu,mv, tol
       integer                :: i,j,k,nedges,n,msol,nsol,k1,k2
       Type(Zone_Axis_Type)   :: Zone_Axis
       Type(Zone_Planes_Type) :: Zone_Planes
       Logical                :: ok,cond,v13,v14,v23,v24

       do
          call system('cls')
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "        GENERAL CRYSTALLOGRAPHY CALCULATOR "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "      Zone axis: list of zone planes and angles"
          write(unit=*,fmt="(a)") " ====================================================="
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " "

          if (Cell_Given) then
             write(unit=*,fmt="(a)",advance="no") " => Enter the indices of plane (<cr> for exit): "
             read(*,'(a)') line
             if (len_trim(line)==0) exit
             line=adjustl(line)
             read(unit=line,fmt=*,iostat=ierr) h
             if(ierr /= 0) then
                 write(unit=*,fmt="(a)") " -> Error in the indices ... retry!"
                 call system("pause ")
                 cycle
             end if
             write(unit=*,fmt="(a)",advance="no") " => Enter the number of edges (<cr> for exit): "
             read(*,*,iostat=ierr) nedges
             if(ierr /= 0) then
                 write(unit=*,fmt="(a)") " -> Error in nedges ... retry!"
                 call system("pause ")
                 cycle
             end if
             if(allocated(angs)) deallocate(angs)
             if(allocated(nl)) deallocate(nl)
             allocate(angs(nedges),nl(nedges))
             nl=0
             do i=1,nedges
               j=mod(i+1,nedges)
               if(j == 0) j=nedges
               write(unit=*,fmt="(2(a,i2),a)",advance="no") " => Enter the angle between edge: ",i," and edge ",j,": "
               read(*,*) angs(i)
             end do

             write(unit=*,fmt="(a)",advance="no") " => Enter the maximun distance in angstrom for n_uvw (<cr> for exit): "
             read(*,'(a)') line
             if (len_trim(line)==0) exit
             line=adjustl(line)
             read(unit=line,fmt=*,iostat=ierr) n_max
             if(ierr /= 0) then
                 write(unit=*,fmt="(a)") " -> Error in n_max ... retry!"
                 call system("pause ")
                 cycle
             end if

             write(unit=*,fmt="(a)",advance="no") " => Angular tolerance in degrees (<cr> for exit): "
             read(unit=*,fmt=*,iostat=ierr) tol
             if(ierr /= 0) tol=1.0
             if(tol < 1.0) tol=1.0

             call Get_basis_from_uvw(n_max,h,celda,Zone_Axis,ok,"R")
             write(unit=*,fmt="(a,3i4)") "  a1: ",Zone_Axis%rx
             write(unit=*,fmt="(a,3i4)") "  a2: ",Zone_Axis%ry

             if(.not. ok) then
                write(unit=*,fmt="(a,3i3,a)") " => Problem in getting the basis in the direct plane perpendicular to [",h,"]*"
                cycle
             end if

             call Get_Zone_Planes(Zone_Axis%rx,Zone_Axis%ry,n_max,Celda,Zone_Planes,"R")

             write(unit=*,fmt="(a)")  "    Edge indices    Edge indices     Normal Angle    Compl. Angle"
             write(unit=*,fmt="(a)")  "    ============    ============     ============     ==========="
             n=Zone_Planes%nplanes*(Zone_Planes%nplanes-1)
             if(allocated(ii)) deallocate(ii)
             if(allocated(jj)) deallocate(jj)
             if(allocated(angles)) deallocate(angles)
             if(allocated(nn)) deallocate(nn)
             if(allocated(nang)) deallocate(nang)
             if(allocated(p)) deallocate(p)
             allocate(ii(n),jj(n),angles(n),nn(n),nang(n),p(nedges,n))
             ii=0; jj=0; angles=0.0; nn=0 ; p=0; nang=0
             ! Calculation of all interplane angles
             msol=0
             do i=1,Zone_Planes%nplanes-2
                h1=Zone_Planes%h(:,i)
                uc=Zone_Planes%hc(:,i)
                mu=Zone_Planes%mh(i)
                do j=i+1,Zone_Planes%nplanes
                  h2=Zone_Planes%h(:,j)
                  vc=Zone_Planes%hc(:,j)
                  mv=Zone_Planes%mh(j)
                  angle=dot_product(uc,vc)/mu/mv
                  if(angle > 1.0) angle=1.0
                  if(angle <-1.0) angle=-1.0
                  angle=acosd(angle)
                  edges: do n=1,nedges
                    if(abs(angle-angs(n)) < tol .or. abs(180.0-angle-angs(n)) < tol) then
                      do k=1, msol
                        h3=Zone_Planes%h(:,ii(k))
                        h4=Zone_Planes%h(:,jj(k))
                        if( (equal_vector(h1,h3,3) .and. equal_vector(h2,h4,3)) .or.  &
                        (equal_vector(h1,h4,3) .and. equal_vector(h2,h3,3)) ) cycle edges
                      end do
                      msol=msol+1
                      write(unit=*,fmt="(2(tr4,3i4),2(tr6,f10.4),a,i5)") Zone_Planes%h(:,i),Zone_Planes%h(:,j),&
                                                                         angle,180-angle,"  #", msol
                      ii(msol)=i
                      jj(msol)=j
                      angles(msol)=angle
                      nang(msol) = n
                    end if
                  end do edges
                end do
             end do


              nsol=0
              do i=1,msol
                h1= Zone_Planes%h(:,ii(i))
                h2= Zone_Planes%h(:,jj(i))
                n=1
                nn(n)=i
                k=nang(i)
                !do j=i+1,msol
                do j=1,msol
                  if( i == j) cycle
                  h3=Zone_Planes%h(:,ii(j))
                  h4=Zone_Planes%h(:,jj(j))

                  v13=equal_vector(h3,h1,3)
                  v14=equal_vector(h4,h1,3)
                  v23=equal_vector(h3,h2,3)
                  v24=equal_vector(h4,h2,3)
                  cond= v13 .or. v14 .or. v23 .or. v24
                  if(.not. cond) cycle

                  !now compare if the present angle is consecutive to that of the
                  !previous pair of edges
                  k=k+1
                  if(k > nedges) k=1
                  ok= abs(angles(j)-angs(k)) < tol  .or.  abs(180-angles(j)-angs(k)) < tol
                  if(.not. ok)  then
                    k=k-1
                    cycle
                  end if

                  if(n == 1) then
                    if(v13 .or. v14) hi=h2
                    if(v23 .or. v24) hi=h1
                    if(v13 .or. v23) then
                      h1=h4
                      h2=h4
                    else
                      h1=h3
                      h2=h3
                    end if
                  else if(n == nedges-2) then
                    h1= hi
                    h2= ho
                  else if(n == nedges) then
                    if(.not.(v13 .and. v24)) then
                      if( .not. (v14 .and. v23) ) then
                        n=n-1
                        cycle
                      end if
                    end if
                  else
                    if(v13 .or. v23) ho=h4
                    if(v14 .or. v24) ho=h3
                    h1=ho
                    h2=ho
                  end if

                  n=n+1
                  if(n > nedges) then
                    n=n-1
                    exit
                  end if
                  nn(n)=j
                end do
                if(n /= nedges) cycle  !no solution
                nsol=nsol+1
                p(1:n,nsol)=nn(1:n)  !List of nedges pairs that are solutions
              end do

              write(unit=*,fmt="(/,a)")  "    List of solutions:"
              do i=1,nsol
                write(unit=*,fmt="(a,i5,a,40i5)")   "    Solutions #",i,":",p(:,i)
                do j=1,nedges
                  k=p(j,i)
                  k1=ii(k)
                  k2=jj(k)
                  write(unit=*,fmt="(a,3i3,a,3i3,f10.3)")"                ", &
                  Zone_Planes%h(:,k1)," --- ",Zone_Planes%h(:,k2),angles(k)
                end do
              end do


             !List of vectors
             write(unit=*,fmt="(/,a)")  "    Edge indices      n_uvw(Angstroms)"
             do i=1,Zone_Planes%nplanes
               write(unit=*,fmt="((tr4,3i4),tr4,f10.3)")  Zone_Planes%h(:,i),Zone_Planes%mh(i)
             end do
          else
             write(unit=*,fmt="(a)") " -> The unit cell must be given first !"
             call system("pause ")
             exit
          end if
          write(unit=*,fmt=*) " "
          call system("pause ")
       end do
    End Subroutine Menu_Geom_12a


   !!----  Get_Zone_Planes(u1,u2,dmin,Cell,Zone_Planes,Mode)
   !!----   !---- Arguments ----!
   !!----   integer, dimension(3),   intent(in) :: u1
   !!----   integer, dimension(3),   intent(in) :: u2
   !!----   real,                    intent(in) :: dmin
   !!----   Type(Crystal_Cell_Type), intent(in) :: Cell
   !!----   Type(Zone_Planes_Type),  intent(out):: Zone_Planes
   !!----   Character(len=*),optional,intent(in):: Mode
   !!----
   !!----   This subroutine provides the Zone planes up to dmin, with
   !!----   coprime indices. It needs as input two reciprocal lattice vectors
   !!----   u1 and u2 defining a primitive 2D-cell in the reciprocal plane
   !!----   perpendicular to the Zone axis. These two vectors can be obtained
   !!----   by a call to subroutine Get_basis_from_uvw from the CrysFML module:
   !!----   CFML_Crystal_Metrics.
   !!----   If Mode is provided, everything is interpreted in the direct space.
   !!----
   !!----   Created: February 2012 (JRC)
   !!----
   Subroutine Get_Zone_Planes(u1,u2,dmin,Cell,Zone_Planes,Mode)
      !---- Arguments ----!
      integer, dimension(3),   intent(in) :: u1
      integer, dimension(3),   intent(in) :: u2
      real,                    intent(in) :: dmin
      Type(Crystal_Cell_Type), intent(in) :: Cell
      Type(Zone_Planes_Type),  intent(out):: Zone_Planes
      Character(len=*),optional,intent(in):: Mode

      !---- Local variables ----!
      integer              :: i1,i2,n,i,j
      real                 :: smax,rm,rm1,rm2
      real                 :: cosgam,cos2,singam
      integer,dimension(3) :: h
      real,dimension(3,3)  :: mat,mg
      integer,dimension(:,:), allocatable :: cp_ind

      if(present(mode)) then
        call get_recip_to_Cart(real(u1),real(u2),dmin,cell,mat,"R")
        smax=dmin
        mg=cell%gd
      else
        call get_recip_to_Cart(real(u1),real(u2),dmin,cell,mat)
        smax=1.0/dmin
        mg=cell%gr
     end if


      !---- Determination of maximun linear combination indices ----!
      rm1=norm(u1,mg)
      rm2=norm(u2,mg)
      i1=nint(smax/rm1)+1
      i2=nint(smax/rm2)+1                         ! i1,i2: Maximum indices for lin. comb.

      n=(i1+1)*(2*i2+1)
      if (allocated(cp_ind)) deallocate(cp_ind)
      allocate(cp_ind(3,n))
      cp_ind=0

      n=0
      do i=-i1,i1
         do j=-i2,i2
            if (i == 0 .and. j == 0) cycle
            h = i*u1 + j*u2
            if(.not. co_prime(h)) cycle
            rm=norm(h,mg)
            if ( rm > smax) cycle
            n=n+1
            cp_ind(:,n)=h
         end do
      end do

      Zone_Planes%nplanes=n
      if(allocated(Zone_Planes%mh)) deallocate(Zone_Planes%mh)
      allocate(Zone_Planes%mh(n))
      if(allocated(Zone_Planes%h)) deallocate(Zone_Planes%h)
      allocate(Zone_Planes%h(3,n))
      if(allocated(Zone_Planes%hc)) deallocate(Zone_Planes%hc)
      allocate(Zone_Planes%hc(3,n))

      do i=1,n
         Zone_Planes%h(:,i)=cp_ind(:,i)
         Zone_Planes%hc(:,i)=matmul(mat,cp_ind(:,i))
         Zone_Planes%mh(i)=sqrt(dot_product(Zone_Planes%hc(:,i),Zone_Planes%hc(:,i)))
      end do
      return
   End Subroutine Get_Zone_Planes


   !!----Subroutine Get_Recip_to_Cart(u1,u2,dmin,Cell,Og,Mode)
   !!----   !---- Arguments ----!
   !!----   real, dimension(3),         intent(in) :: u1
   !!----   real, dimension(3),         intent(in) :: u2
   !!----   real,                       intent(in) :: Dmin
   !!----   type (Crystal_Cell_Type),   intent(in) :: Cell
   !!----   real, dimension(3,3),       intent(out):: Og
   !!----   Character(len=*), optional, intent(in) :: Mode
   !!----
   !!----   This subroutine defines a Cartesian frame in which the x-axis is
   !!----   along the reciprocal vector u1, the y-axis in within the plane
   !!----   defined by the vectors (u1,u2) and the z-axis completes the right-handed
   !!----   frame (z = x cross_product y). The output of the subroutine is the
   !!----   matrix OG, converting the indices of reciprocal lattice points to Cartesian
   !!----   coordinates with respect to the above frame: [hc]= OG [h]
   !!----   If Mode is provided we interpret everything above in the sense of direct
   !!----   space: u1 and u2 are lattice vectors, Dmin is n_max, etc.
   !!----
   !!----   Created: February 2006
   !!----   Updated: February 2012 (JRC)
   !!----
   Subroutine Get_Recip_to_Cart(u1,u2,dmin,Cell,Og,Mode)
      !---- Arguments ----!
      real, dimension(3),         intent(in) :: u1
      real, dimension(3),         intent(in) :: u2
      real,                       intent(in) :: Dmin
      type (Crystal_Cell_Type),   intent(in) :: cell
      real, dimension(3,3),       intent(out):: Og
      Character(len=*), optional, intent(in) :: Mode

      !---- Local variables ----!
      integer                          :: i1,i2,nref,n,ii,jj,j,k,iv
      real                             :: d,rm,rm1,rm2,rm3,ri,rj,s,smax
      real                             :: deter,cosgam,cos2,singam,cosbet,cost,vrec,anglp
      real, dimension(3)               :: h,u3,r,r1,r2
      real, dimension(3,3)             :: uc,b,tinv,mat
      real, parameter, dimension(3,16) :: rmore=reshape( (/1.,0.,0.,  0.,1.,0.,  0.,0.,1.,  1.,1.,0.,  1.,0.,1., &
                                                           0.,1.,1., -1.,1.,0., -1.,0.,1.,  0.,-1.,1., 1.,1.,1., &
                                                          -1.,1.,1.,  1.,-1.,1., 1.,1.,-1.,-1.,-1.,1.,-1.,1.,-1.,&
                                                           1.,-1.,-1./),(/3,16/))

      ! Initializing
      Og=0.0
      if(present(Mode)) then
        Smax=dmin
        mat=cell%gd
      else
        Smax=1.0/dmin
        mat=cell%gr
      end if
      !---- Determination of maximun linear combination indices ----!
      rm1=norm(u1,mat)
      rm2=norm(u2,mat)
      i1=nint(smax/rm1)+1
      i2=nint(smax/rm2)+1

      !---- Angle between the two reflections ----!
      cosgam=scalar(u1,u2,mat)/(rm1*rm2)
      cos2=cosgam*cosgam
      if (cos2 > 1.0) cos2=1.0
      singam= sqrt(1.0-cos2)

      !---- Cartesian coordinates in the plane: ----!
      !----    u1 is along x,  u2 is in the xy plane
      !----    r1,r2 cartesian components of u1 and u2
      !----
      r1=(/rm1,0.0,0.0/)
      r2=(/rm2*cosgam,rm2*singam,0.0/)
      uc(1,:)=r1
      uc(2,:)=r2

      !----  Determination of a third reciprocal lattice vector defining a primitive
      !----  basis in E* together with u1 and u2. The vector u3 is selected between
      !----  the vectors given in rmore(3,16) and calculate :
      !----                 (Rm1       0          0         )
      !----          [uc] = (Rm2*cosg  Rm2*sing   0         ) [e]
      !----                 (Rm3.cosb  Rm3.cosp   Rm3.cost  )
      do iv=1,16
         u3(:)=rmore(:,iv)
         deter=Determ_V(u1,u2,u3)
         if (deter > 0.1) then
            exit
         else if(-deter > 0.1) then
            u3(:)=-u3(:)
            exit
         end if
      end do
      deter= abs(deter)
      rm3=norm(u3,mat)
      cosbet=scalar(u1,u3,mat)/(rm1*rm3)
      uc(3,1)=rm3*cosbet
      vrec=sqrt(determ_a(mat))               !to obtain the volume of the reciprocal cell
      cost= deter*vrec/(rm1*rm2*rm3*singam)
      uc(3,3)=rm3*cost
      anglp=1.0-cost*cost-cosbet*cosbet
      if (anglp < 0) anglp=0.0
      uc(3,2)=rm3*sqrt(anglp)   !Matrix C  [u]=C[e]
      do j=1,3
         b(1,j)=u1(j)
         b(2,j)=u2(j)           !Matrix M      [u]=M[a*]
         b(3,j)=u3(j)           !Cartesian coord. of hkl: [X]=Ct Mt-1 [h] = OG [h]
      end do

      tinv=invert_a(b)     !M-1
      b=transpose(tinv)    !Mt-1  -->b
      tinv=transpose(uc)   !Ct  --> tinv

      og=matmul(tinv,b)    !Ct Mt-1 = OG

      return
   End Subroutine Get_Recip_to_Cart
    !!----
    !!---- Subroutine Get_basis_from_uvw(dmin,u,cell,ZoneB,ok,mode)
    !!----    real(kind=cp)             intent(in) :: dmin  !minimum d-spacing (smax=1/dmin)
    !!----    integer, dimension(3),    intent(in) :: u     !Zone axis indices
    !!----    type (Crystal_Cell_Type), intent(in) :: cell
    !!----    type (Zone_Axis_Type),    intent(out):: ZoneB !Object containing u and basis vector in the plane
    !!----    logical,                  intent(out):: ok
    !!----    character(len=*),optional,intent(in) :: mode
    !!----
    !!----  Subroutine to construct ZA of type Zone_Axis. This subroutine picks up two reciprocal
    !!----  lattice vectors satisfying the equation
    !!----                            hu+kv+lw=0
    !!----  The two reciprocal lattice vectors have no coprime factors and
    !!----  constitute the basis of a reciprocal lattice plane. They are
    !!----  obtained as the shortest two reciprocal lattice vectors satisfying
    !!----  the above equation. If mode is provided and mode="R", we interpret
    !!----  that the input zone axis is a reciprocal lattice vector and what we
    !!----  obtain is the basis of a direct plane in terms of lattice vectors.
    !!----  If mode="R", dmin corresponds n(uvw)max
    !!----  This subroutine has been imported from resvis_proc.f90.
    !!----
    !!----  Created: February 2006 (Imported from old programs for electron diffraction, Thesis JRC)
    !!----  Updated: February 2012 (JRC)
    !!----
!    Subroutine Get_basis_from_uvw(dmin,u,cell,ZoneB,ok,mode)
!       !--- Arguments ---!
!       real(kind=cp),            intent(in) :: dmin
!       integer, dimension(3),    intent(in) :: u
!       type (Crystal_Cell_Type), intent(in) :: cell
!       type (Zone_Axis_Type),    intent(out):: ZoneB
!       logical,                  intent(out):: ok
!       character(len=*),optional,intent(in) :: mode
!
!       !--- Local Variables ---!
!       integer                :: n,ik,il,um,iv,i1,i2,i,coun01,coun02,coun1,coun2
!       integer,dimension(1)   :: i0
!       integer                :: kmin,kmax,lmin,lmax
!       integer,dimension(3)   :: au,h,mu
!       real, dimension(2)     :: rm
!       real, dimension(3,3)   :: mat
!       integer,dimension(3,2) :: bas
!       real                   :: rv,s2max
!
!       ZoneB%nlayer=0
!       ZoneB%uvw=u
!       ok=.false.
!
!       au=abs(u)
!       um=3*maxval(au)
!       i0=maxloc(au)
!
!       i=i0(1)
!       iv=u(i)
!       mu=u
!       if (iv < 0) then
!         mu=-u
!         iv=-iv
!       end if
!
!       Select Case (i)
!         Case(1)
!           i1=2; i2=3
!         Case(2)
!           i1=1; i2=3
!         Case(3)
!           i1=1; i2=2
!       End Select
!
!       rm(1)=100000.0; rm(2)=rm(1)
!       bas(:,1) = (/ 71,121, 113/)
!       bas(:,2) = (/117, 91,-111/)
!
!       if(present(mode)) then
!         s2max=dmin*dmin   !here dmin is really n_max
!         kmax=nint(dmin/Cell%cell(i1)+1.0)
!         lmax=nint(dmin/Cell%cell(i2)+1.0)
!         kmax=min(um,kmax)
!         lmax=min(um,lmax)
!         mat=cell%gd
!       else
!         s2max=1.0/(dmin*dmin)
!         kmax=nint(Cell%cell(i1)/dmin+1.0)
!         lmax=nint(Cell%cell(i2)/dmin+1.0)
!         kmax=min(um,kmax)
!         lmax=min(um,lmax)
!         mat=cell%gr
!       end if
!
!       kmin=-kmax; lmin=-lmax
!       coun1=0; coun2=0
!       do ik=kmax,kmin,-1
!          do il=lmax,lmin,-1
!             if (ik == 0 .and. il == 0) cycle
!             n=-ik*mu(i1)-il*mu(i2)
!             if (mod(n,iv) == 0) then               !n is multiple of iv
!                h(i)= n/iv ; h(i1)=ik ; h(i2) = il  !h is solution of hu+kv+lw=0
!                rv=dot_product(real(h),matmul(mat,real(h)))
!                if (rv > s2max  .or. rv < 1.0e-20) cycle
!                if (rv < rm(1)) then
!                   if (.not. co_linear(bas(:,1),h,3) ) then
!                      bas(:,2)=bas(:,1)
!                      rm(2) = rm(1)
!                      if (coun1 >=1) coun2=coun2+1
!                   end if
!                   bas(:,1)=h
!                   rm(1) = rv
!                   coun1=coun1+1
!                else if (rv < rm(2) .and. .not. co_linear(bas(:,1),h,3) ) then
!                   bas(:,2)=h
!                   rm(2) = rv
!                   coun2=coun2+1
!                end if
!             end if
!          end do
!       end do
!       ZoneB%rx=bas(:,1)
!       ZoneB%ry=bas(:,2)
!       if (coun1 >= 1 .and. coun2 >=1) ok=.true.
!       coun01=0; coun02=0; coun1=0; coun2=0
!       do i=1,3
!          if (ZoneB%rx(i) < 0) coun1=coun1+1
!          if (ZoneB%ry(i) < 0) coun2=coun2+1
!          if (ZoneB%rx(i) == 0) coun01=coun01+1
!          if (ZoneB%ry(i) == 0) coun02=coun02+1
!       end do
!       if (coun1 >= 2 .or. (coun1 == 1 .and. coun01 == 2)) ZoneB%rx=-ZoneB%rx
!       if (coun2 >= 2 .or. (coun2 == 1 .and. coun02 == 2)) ZoneB%ry=-ZoneB%ry
!
!       return
!    End Subroutine Get_Basis_From_Uvw
 End Module Menu_5
