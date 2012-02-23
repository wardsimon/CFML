!!----
!!---- Menu: 4
!!---- Chemical Information
!!----
!!
 Module Menu_4
   !---- Use File ----!
   use CFML_scattering_chemical_tables
   use CFML_string_utilities

   !---- Variables ----!
   implicit none

 Contains

    !!----
    !!---- Subroutine Menu_Princ4
    !!----
    !!
    Subroutine Menu_Princ4()
       !---- Local Variables ----!
       character(len=2)  :: car

       do
          call system('cls')

          write(unit=*,fmt="(a)") "     GENERAL CRYSTALLOGRAPHY CALCULATOR "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "            Chemical Information "
          write(unit=*,fmt="(a)") "        ============================"
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " [0] Back..."
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " [1] International Table"
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)",advance="no") " OPTION: "
          read(*,'(a)') car
          if (len_trim(car) == 0) exit
          car=adjustl(car)

          select case (car)
             case ('0 ')
                exit

             case ('1 ')
                call Menu_CHEM_1()

             case ('2 ')

             case ('3 ')

             case ('4 ')

             case ('5 ')

          end select
       end do

    End Subroutine Menu_Princ4

    !!----
    !!---- Subroutine Menu_CHEM_1
    !!----
    !!
    Subroutine Menu_Chem_1()
       !---- Local Variables ----!
       character(len=20)     :: line
       integer               :: i, j, n, iv, ierr, npos
       integer, dimension(1) :: ivet
       real, dimension(1)    :: vet

       call Set_Chem_Info()

       do
          call system('cls')
          write(unit=*,fmt=*) "     GENERAL CRYSTALLOGRAPHY CALCULATOR "
          write(unit=*,fmt=*) " "
          write(unit=*,fmt=*) "     International Chemical Table "
          write(unit=*,fmt=*) " ===================================="
          write(unit=*,fmt=*) " "
          write(unit=*,fmt=*) " "
          write(unit=*,fmt="(a)",advance="no") " => Enter chemical symbol or atomic number of the Element: "
          read(unit=*,fmt='(a)') line
          if (len_trim(line)==0) exit
          line=adjustl(line)

          !---- Introduce Numero? ----!
          call getnum(line,vet,ivet,iv)
          if (iv == 1) then
             do i=1,num_chem_info
                if (ivet(1) /= chem_info(i)%z) cycle
                write(unit=*,fmt=*) " "
                write(unit=*,fmt=*) " "
                write(unit=*,fmt='(5x,a,a)')      "            Name: ",chem_info(i)%name
                write(unit=*,fmt='(5x,a,a)')      "          Symbol: ",chem_info(i)%symb
                write(unit=*,fmt='(5x,a,i3)')     "   Atomic Number: ",chem_info(i)%z
                write(unit=*,fmt='(5x,a,f10.5)')   "   Atomic Weight: ",chem_info(i)%atwe
                write(unit=*,fmt='(5x,a,f6.3)')" Covalent Radius: ",chem_info(i)%rcov
                write(unit=*,fmt='(5x,a,f6.3)')"    Waals Radius: ",chem_info(i)%rwaals

                n=count(chem_info(i)%oxid /= 0)
                line=" "
                npos=0
                do j=1,n
                   write(unit=line(npos+1:),fmt='(i3)') chem_info(i)%oxid(j)
                   npos=len_trim(line)
                end do
                write(unit=*,fmt='(5x,a,a)')     "Oxidation States: ",line
                write(unit=*,fmt=*) " "
                write(unit=*,fmt=*) " "
                exit
             end do
          else
             !---- Introduce notacion del Elemento ----!
             do i=1,num_chem_info
                if (u_case(line(1:2)) /= chem_info(i)%symb ) cycle
                write(unit=*,fmt=*) " "
                write(unit=*,fmt=*) " "
                write(unit=*,fmt='(5x,a,a)')      "            Name: ",chem_info(i)%name
                write(unit=*,fmt='(5x,a,a)')      "          Symbol: ",chem_info(i)%symb
                write(unit=*,fmt='(5x,a,i3)')     "   Atomic Number: ",chem_info(i)%z
                write(unit=*,fmt='(5x,a,f8.4)')   "   Atomic Weight: ",chem_info(i)%atwe
                write(unit=*,fmt='(5x,a,2x,f6.3)')" Covalent Radius: ",chem_info(i)%rcov
                write(unit=*,fmt='(5x,a,2x,f6.3)')"    Waals Radius: ",chem_info(i)%rwaals

                n=count(chem_info(i)%oxid /= 0)
                line=" "
                npos=0
                do j=1,n
                   write(unit=line(npos+1:),fmt='(i3)') chem_info(i)%oxid(j)
                   npos=len_trim(line)
                end do
                write(unit=*,fmt='(5x,a,a)')     "Oxidation States: ",line
                write(unit=*,fmt=*) " "
                write(unit=*,fmt=*) " "
                exit
             end do

          end if

          call system('pause')

       end do
    End Subroutine Menu_Chem_1

 End Module Menu_4
