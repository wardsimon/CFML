!!----
!!---- CRYSTALLOGRAPHIC CALCULATOR
!!---- Version:1.00
!!----         Oct-2002
!!----
!!---- Authors: Juan Rodriguez-Carvajal
!!----          Javier Gonzalez-Platas
!!----
!!
 Program Calsym
    !---- Use files ----!
    use Menu_0
    use Menu_1
    use Menu_2
    use Menu_3
    use Menu_4
    use Menu_5
    use CFML_IO_Messages,  only: Wait_Message
    !---- Variables ----!
    implicit none

    character (len=2):: car

    if(OPS == 1) clear_string="cls" !Change for Windows

    !---- Menu Principal ----!
    open(unit=i_out,file=fileout,status="unknown",action="write",position="append")
    write(unit=i_out,fmt="(a)") " LOG-file of the program: CrysCalCon"
    write(unit=i_out,fmt="(a)") " ==================================="
    call Write_Date_Time(i_out)
    write(unit=i_out,fmt="(a)") "  "

    do
       call system(clear_string)

       write(unit=*,fmt="(a)") " "
       write(unit=*,fmt="(a)") "     GENERAL CRYSTALLOGRAPHY CALCULATOR "
       write(unit=*,fmt="(a)") " "
       write(unit=*,fmt="(a)") "     Principal Menu "
       write(unit=*,fmt="(a)") " ======================"
       write(unit=*,fmt="(a)") " "
       write(unit=*,fmt="(a)") " "
       write(unit=*,fmt="(a)") " [0] Exit"
       write(unit=*,fmt="(a)") " "
       write(unit=*,fmt="(a)") " [1] Space Groups"
       write(unit=*,fmt="(a)") " [2] Reflections"
       write(unit=*,fmt="(a)") " [3] Atomistics Calculations"
       write(unit=*,fmt="(a)") " [4] Chemical Information"
       write(unit=*,fmt="(a)") " [5] Geometry calculations"
       write(unit=*,fmt="(a)") " "
       write(unit=*,fmt="(a)",advance="no") " OPTION: "
       read(unit=*,fmt="(a)") car
       if (len_trim(car) == 0) exit
       car=adjustl(car)

       select case (car)
           case ("0 ")
              exit

           case ("1 ")
              call Menu_Princ1()

           case ("2 ")
              call Menu_Princ2()

           case ("3 ")
              call Menu_Princ3()

           case ("4 ")
              call Menu_Princ4()

           case ("5 ")
              call Menu_Princ5()

       end select
    end do
    close(unit=i_out)
    call Wait_Message(" => Press <enter> to close the program ...")
    call system(clear_string)
 End Program Calsym


