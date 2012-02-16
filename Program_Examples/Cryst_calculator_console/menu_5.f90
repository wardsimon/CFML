!!----
!!---- Menu: 2
!!---- Reflections
!!----
!!
 Module Menu_5
   !---- Use File ----!
   !use CFML_Crystallographic_Symmetry
   use CFML_Crystal_Metrics
   !use CFML_Math_3D
   use CFML_String_Utilities, only:  L_Case, pack_string, Get_Mat_From_Symb, Err_String, ERR_String_Mess


   !---- Variables ----!
   implicit none

   type (Crystal_Cell_type) :: Celda
   Logical                  :: cell_Given=.false.
   real, parameter          :: eps=1.0e-7
   integer                  :: ierr
   Character(len=1)         :: Cart_type="A"

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

          write(unit=*,fmt="(a)") "     GENERAL CRYSTRALLOGRAPHY CALCULATOR "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "     Geometrical Calculations Menu "
          write(unit=*,fmt="(a)") " ======================================"
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " [0] Back..."
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " [1] Enter unit cell parameters (or show current cell)"
          write(unit=*,fmt="(a)") " [2] Angle between two directions in direct space"
          write(unit=*,fmt="(a)") " [3] Angle between two directions in reciprocal space"
          write(unit=*,fmt="(a)") " [4] Direction in direct space => unitary direction in reciprocal space"
          write(unit=*,fmt="(a)") " [5] Direction in reciprocal space => unitary direction in direct space"
          write(unit=*,fmt="(a)") " [6] Basis change matrix => get new unit cell parameters"
          write(unit=*,fmt="(a)") " [7] Zone axis: list of zone planes and angles"
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

          end select
       end do

    End Subroutine Menu_Princ5

    !!----
    !!---- Subroutine Menu_Geom_1
    !!----
    !!
    Subroutine Menu_Geom_1()
       !---- Local Variables ----!
       character(len=20)        :: line
       character(len=1)         :: car
       real, dimension(3)       :: Cellv,Angl

       do
          call system('cls')
          write(unit=*,fmt="(a)") "     GENERAL CRYSTRALLOGRAPHY CALCULATOR "
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
                   Cart_Type="A"
               End Select

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
       character(len=40)      :: line
       real, dimension(3)     :: u,v,uc,vc,ur,vr
       real                   :: angle,mu,mv

       do
          call system('cls')
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "        GENERAL CRYSTRALLOGRAPHY CALCULATOR "
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
             write(unit=*,fmt="(a, f10.5,a)")     "     Angle:",angle, " degrees"
          else
             write(unit=*,fmt="(a)") " -> The unit cell must be given first !"
             call system("pause ")
             exit
          end if
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
       character(len=40)      :: line
       real, dimension(3)     :: h,k,hc,kc,hd,kd
       real                   :: angle,mh,mk

       do
          call system('cls')
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "        GENERAL CRYSTRALLOGRAPHY CALCULATOR "
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
             write(unit=*,fmt="(a)")              "             Reciprocal space[A*]      Cartesian components[e]     Direct Components[A]"
             write(unit=*,fmt="(a,3(3f8.4,tr2))") "     First:",h,hc,hd
             write(unit=*,fmt="(a,3(3f8.4,tr2))") "    Second:",k,kc,kd
             write(unit=*,fmt="(a, f10.5,a)")     "     Angle:",angle, " degrees"
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
    subroutine Menu_Geom_4()
       !---- Local Variables ----!
       character(len=40)      :: line
       real, dimension(3)     :: u,uc,vc,ur
       real                   :: mu,mv

       do
          call system('cls')
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "        GENERAL CRYSTRALLOGRAPHY CALCULATOR "
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


    End Subroutine Menu_Geom_4

    !!----
    !!---- Subroutine Menu_Geom_5
    !!----
    !!
    Subroutine Menu_Geom_5()
       !---- Local Variables ----!
       character(len=40)      :: line
       real, dimension(3)     :: u,v,uc,vc
       real                   :: mu,mv

       do
          call system('cls')
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "        GENERAL CRYSTRALLOGRAPHY CALCULATOR "
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
    End Subroutine Menu_Geom_5

    !!----
    !!---- Subroutine Menu_Geom_6
    !!----
    !!
    subroutine Menu_Geom_6()
       !---- Local Variables ----!
       character(len=40)        :: line,symb
       real, dimension(3,3)     :: mat
       type (Crystal_Cell_type) :: Celln
       integer                  :: i

       do
          call system('cls')
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "        GENERAL CRYSTRALLOGRAPHY CALCULATOR "
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
    End Subroutine Menu_Geom_6


    subroutine Menu_Geom_7()
       !---- Local Variables ----!


    End Subroutine Menu_Geom_7

!    Subroutine get_mat_from_symb(Symb,Mat,cod)
!      character(len=*),                intent(in)  :: Symb
!      real,dimension(3,3),             intent(out) :: Mat
!      character(len=1), dimension(3),  intent(in)  :: cod
!      !---- local variables ----!
!      integer :: i,j,k
!      character(len=len(Symb)), dimension(3) :: split
!
!      i=index(Symb,",")
!      j=index(Symb,",",back=.true.)
!      split(1)= pack_string(Symb(1:i-1))
!      split(2)= pack_string(Symb(i+1:j-1))
!      split(3)= pack_string(Symb(j+1:))
!      do i=1,3
!       call get_num_string(trim(split(i)), Mat(i,:),cod)
!      end do
!     contains
!
!        Subroutine get_num_string(string, v,cod)
!          character(len=*),                intent(in)  :: string
!          real,dimension(3),               intent(out) :: v
!          character(len=1), dimension(3),  intent(in)  :: cod
!          integer :: i
!          integer, dimension(3) :: j
!          character(len=len(String)), dimension(3) :: nsp
!
!          do i=1,3
!            j(i)=index(string,cod(i))
!          end do
!
!          if(j(1) == 0) then
!            v(1)=0.0
!            if(j(2) == 0) then
!              v(2)=0.0
!              !only the third component exist
!              call Read_Fract(string(1:j(3)-1), v(3))
!            else if(j(3) == 0) then
!              v(3)=0.0
!              !only the second component exist
!              call Read_Fract(string(1:j(2)-1), v(2))
!            else    !j(2) and j(3) are not zero
!              if(j(2) < j(3)) then
!                call Read_Fract(string(1:j(2)-1), v(2))
!                call Read_Fract(string(j(2)+1:j(3)-1), v(3))
!              else
!                call Read_Fract(string(1:j(3)-1), v(3))
!                call Read_Fract(string(j(3)+1:j(2)-1), v(2))
!              end if
!            end if
!
!          else if(j(2) == 0) then
!            v(2)=0.0
!            !j(1)/=0 now, just check j(3)
!            if(j(3) == 0) then
!              v(3)=0.0
!              !only the First component exist
!              call Read_Fract(string(1:j(1)-1), v(1))
!            else  ! only 1 and 3
!              if(j(1) < j(3)) then
!                call Read_Fract(string(1:j(1)-1), v(1))
!                call Read_Fract(string(j(1)+1:j(3)-1), v(3))
!              else
!                call Read_Fract(string(1:j(3)-1), v(3))
!                call Read_Fract(string(j(3)+1:j(1)-1), v(1))
!              end if
!
!            endif
!
!          else if(j(3) == 0) then
!            v(3)=0.0
!            if(j(1) < j(2)) then
!              call Read_Fract(string(1:j(1)-1), v(1))
!              call Read_Fract(string(j(1)+1:j(2)-1), v(2))
!            else
!              call Read_Fract(string(1:j(2)-1), v(2))
!              call Read_Fract(string(j(2)+1:j(1)-1), v(1))
!            end if
!
!          else !none of them is zero
!
!            if(j(1) < j(2) .and. j(2) < j(3)) Then !normal order a b c
!              call Read_Fract(string(1:j(1)-1),      v(1))
!              call Read_Fract(string(j(1)+1:j(2)-1), v(2))
!              call Read_Fract(string(j(2)+1:j(3)-1), v(3))
!            else if(j(2) < j(1) .and. j(1) < j(3)) then ! b a c
!              call Read_Fract(string(1:j(2)-1),      v(2))
!              call Read_Fract(string(j(2)+1:j(1)-1), v(1))
!              call Read_Fract(string(j(1)+1:j(3)-1), v(3))
!            else if(j(2) < j(3) .and. j(3) < j(1)) then ! b c a
!              call Read_Fract(string(1:j(2)-1),      v(2))
!              call Read_Fract(string(j(2)+1:j(3)-1), v(3))
!              call Read_Fract(string(j(3)+1:j(1)-1), v(1))
!            else if(j(3) < j(2) .and. j(2) < j(1)) then ! c b a
!              call Read_Fract(string(1:j(3)-1),      v(3))
!              call Read_Fract(string(j(3)+1:j(2)-1), v(2))
!              call Read_Fract(string(j(2)+1:j(1)-1), v(1))
!            else if(j(3) < j(1) .and. j(1) < j(2)) then ! c a b
!              call Read_Fract(string(1:j(3)-1),      v(3))
!              call Read_Fract(string(j(3)+1:j(1)-1), v(1))
!              call Read_Fract(string(j(1)+1:j(2)-1), v(2))
!            else if(j(1) < j(3) .and. j(3) < j(2)) then ! a c b
!              call Read_Fract(string(1:j(1)-1),      v(1))
!              call Read_Fract(string(j(1)+1:j(3)-1), v(3))
!              call Read_Fract(string(j(3)+1:j(2)-1), v(2))
!            else
!              !This is impossible
!              write(*,*) "  Something strange has occured!"
!            end if
!          end if
!
!        End Subroutine get_num_string
!
!        Subroutine Read_Fract(str,valu)
!         Character(len=*), intent(in) :: str
!         real,             intent(out):: valu
!         integer :: k
!         real    :: num,den
!         if(len_trim(str) == 0) then
!           valu=1.0
!           return
!         else if(len_trim(str) == 1) then
!           if(str == "+") then
!            valu=1.0
!            return
!           else if(str == "-") then
!            valu=-1.0
!            return
!           end if
!         end if
!         k=index(str,"/")
!         if(k == 0) then !a single number
!           read(unit=str,fmt=*,iostat=ierr) valu
!           if(ierr /= 0) then
!              valu=0.0
!              write(*,"(a)") " => Error, input string: "//str
!           end if
!         else !fraction
!           read(unit=str(1:k-1),fmt=*,iostat=ierr) num
!           if(ierr /= 0) then
!              valu=0.0
!              write(*,"(a)") " => Error, input string: "//str(1:k-1)
!              return
!           end if
!           read(unit=str(k+1:),fmt=*,iostat=ierr) den
!           if(ierr /= 0) then
!              valu=0.0
!              write(*,"(a)") " => Error, input string: "//str(k+1:)
!              return
!           end if
!           valu=num/den
!         end if
!        End Subroutine Read_Fract
!
!    End Subroutine get_mat_from_symb

 End Module Menu_5
