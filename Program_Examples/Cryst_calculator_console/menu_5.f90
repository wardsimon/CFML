!!----
!!---- Menu: 2
!!---- Reflections
!!----
!!
 Module Menu_5
   !---- Use File ----!
   use CFML_Crystallographic_Symmetry
   use CFML_Crystal_Metrics
   use CFML_Reflections_Utilities
   use CFML_Math_3D

   !---- Variables ----!
   implicit none

   type (Crystal_Cell_type) :: Celda
   Logical                  :: cell_Given=.false.
   real, parameter          :: eps=1.0e-7
   integer                  :: ierr

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
          write(unit=*,fmt="(a)") " [4] Direction in direct space => direction in reciprocal space"
          write(unit=*,fmt="(a)") " [5] Direction in reciprocal space => direction in direct space"
          write(unit=*,fmt="(a)") " [6] Basis change matrix => get new unit cell parameters"
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
               write(unit=*,fmt="(a)",advance="no") " => Enter a,b,c,alpha,beta,gamma: "
               read(*,'(a)') line
               if (len_trim(line)==0) exit
               read(unit=line,fmt=*,iostat=ierr) Cellv,Angl
               if(ierr /= 0) then
                  write(unit=*,fmt="(a)") " -> Error in the unit cell parameters ... retry!"
                  call system("pause ")
                  cycle
               end if
               call Set_Crystal_Cell(Cellv,Angl,Celda)
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
             uc=Matmul(u,Celda%Cr_Orth_cel)
             ur=Matmul(Celda%GD,u)
             vc=Matmul(v,Celda%Cr_Orth_cel)
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
             write(unit=*,fmt="(a)") "          Direction in direct space      Cartesian components     Reciprocal Components"
             write(unit=*,fmt="(a,3(3f8.4,tr2))") "     First:",u,uc,ur
             write(unit=*,fmt="(a,3(3f8.4,tr2))") "    Second:",v,vc,vr
             write(unit=*,fmt="(a, f10.5,a)")     "     Angle:",angle, " degrees"
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
             hc=Matmul(h,Celda%Cr_Orth_cel)
             hd=Matmul(Celda%GR,h)
             kc=Matmul(k,Celda%Cr_Orth_cel)
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
             write(unit=*,fmt="(a)") "          Direction in direct space      Cartesian components     Reciprocal Components"
             write(unit=*,fmt="(a,3(3f8.4,tr2))") "     First:",h,hc,hd
             write(unit=*,fmt="(a,3(3f8.4,tr2))") "    Second:",k,kc,kd
             write(unit=*,fmt="(a, f10.5,a)")     "     Angle:",angle, " degrees"
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


    End Subroutine Menu_Geom_4

    !!----
    !!---- Subroutine Menu_Geom_5
    !!----
    !!
    Subroutine Menu_Geom_5()


    End Subroutine Menu_Geom_5

    !!----
    !!---- Subroutine Menu_Geom_6
    !!----
    !!
    subroutine Menu_Geom_6()
       !---- Local Variables ----!


    End Subroutine Menu_Geom_6

 End Module Menu_5
