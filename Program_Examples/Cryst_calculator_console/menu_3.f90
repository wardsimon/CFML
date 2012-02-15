!!----
!!---- Menu: 3
!!---- Atoms Calculations
!!----
!!
 Module Menu_3
    !---- Use File ----!
    use CFML_Crystallographic_Symmetry
    use CFML_Atom_TypeDef
    use CFML_string_utilities, only: getnum

    !---- Variables ----!
    implicit none

 Contains

    !!----
    !!---- Subroutine Menu_Princ3
    !!----
    !!
    Subroutine Menu_Princ3()
       !---- Local Variables ----!
       character(len=2)  :: car

       do
          call system('cls')

          write(unit=*,fmt="(a)") "     GENERAL CRYSTRALLOGRAPHY CALCULATOR "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "     Atoms Calculations "
          write(unit=*,fmt="(a)") " ============================"
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " [0] Back..."
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " [1] Multiplicity Position"
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)",advance="no") " OPTION: "
          read(*,'(a)') car
          if (len_trim(car) == 0) exit
          car=adjustl(car)

          select case (car)
             case ('0 ')
                exit

             case ('1 ')
                call Menu_ATOM_1()

             case ('2 ')

             case ('3 ')

             case ('4 ')

             case ('5 ')

          end select
       end do

    End Subroutine Menu_Princ3

    !!----
    !!---- Subroutine Menu_ATOM_1
    !!----
    !!
    Subroutine Menu_Atom_1()
       !---- Local Variables ----!
       character(len=20)     :: line, spgr
       integer               :: i, iv, ierr, npos
       integer, dimension(3) :: ivet
       real                  :: mlt
       real, dimension(3)    :: vet,xp
       type (Space_Group_type)    :: grp_espacial

       do
          call system('cls')
          write(unit=*,fmt="(a)") "     GENERAL CRYSTRALLOGRAPHY CALCULATOR "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "     Space Groups Information "
          write(unit=*,fmt="(a)") " ================================"
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)",advance="no") " Space Group (HM/Hall/Num): "

          read(*,'(a)') line
          if (len_trim(line)==0) exit
          line=adjustl(line)
          call set_spacegroup(line,grp_espacial)

          do
             call system('cls')
             write(unit=*,fmt="(a)") "     GENERAL CRYSTRALLOGRAPHY CALCULATOR "
             write(unit=*,fmt="(a)") " "
             write(unit=*,fmt="(a)") "     Multiplicity Position "
             write(unit=*,fmt="(a)") " ============================="
             write(unit=*,fmt="(a)") " "
             write(unit=*,fmt="(a)") " "
             write(unit=*,fmt="(a)",advance="no") " Position: "

             read(*,'(a)') line
             if (len_trim(line)==0) exit
             line=adjustl(line)

             call getnum(line,vet,ivet,iv)
             if (iv == 3) then
                npos=0
                do i=1,grp_espacial%multip
                   xp=ApplySO(grp_espacial%Symop(i), vet)
                   xp=mod(xp+10.0,1.0)

                   if (abs(xp(1)-vet(1)) .le. 0.001 .and. &
                       abs(xp(2)-vet(2)) .le. 0.001 .and. &
                       abs(xp(3)-vet(3)) .le. 0.001 ) npos =npos+1

                end do

                if (npos /= 0) then
                   mlt=1.0/real(npos)
                else
                   mlt=0.0
                end if

                write(unit=*,fmt=*) " "
                write(*,'(a,f7.4)') " Multiplicity: ",mlt
                call system('pause')
             end if
          end do
       end do

    End Subroutine Menu_Atom_1

end module Menu_3
