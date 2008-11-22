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
    use Menu_1
    use Menu_2
    use Menu_3
    use Menu_4

    !---- Variables ----!
    implicit none

    character (len=2):: car

    !---- Menu Principal ----!
    do
       call system("cls")

       print*,"     GENERAL CRYSTRALLOGRAPHY CALCULATOR "
       print*," "
       print*,"     Principal Menu "
       print*," ======================"
       print*," "
       print*," "
       print*," [0] Exit"
       print*," "
       print*," [1] Space Groups"
       print*," [2] Reflections"
       print*," [3] Atoms Calculations"
       print*," [4] Chemical Information"
       print*," "
       print*," OPTION: "
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

       end select
    end do

    call system("cls")

 End Program Calsym


