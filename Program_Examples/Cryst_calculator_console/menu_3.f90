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

          write(unit=*,fmt="(a)") "     GENERAL CRYSTALLOGRAPHY CALCULATOR "
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
       integer               :: i, iv, ierr, mlt
       integer, dimension(3) :: ivet
       real                  :: occ
       real, dimension(3)    :: vet,xp
       type (Space_Group_type)    :: grp_espacial

       do
          call system('cls')
          write(unit=*,fmt="(a)") "     GENERAL CRYSTALLOGRAPHY CALCULATOR "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "     Space Groups Information "
          write(unit=*,fmt="(a)") " ================================"
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)",advance="no") " Space Group (HM/Hall/Num): "

          read(unit=*,fmt='(a)') line
          if (len_trim(line)==0) exit
          line=adjustl(line)
          call set_spacegroup(line,grp_espacial)
          if(Err_Symm) then
            write(unit=*,fmt="(a)") trim(ERR_Symm_Mess)
            call system('pause')
            cycle
          end if
          do
             call system('cls')
             write(unit=*,fmt="(a)") "       GENERAL CRYSTALLOGRAPHY CALCULATOR "
             write(unit=*,fmt="(a)") " "
             write(unit=*,fmt="(a)") "     Multiplicity and Occupancy of Position "
             write(unit=*,fmt="(a)") "   ==========================================="
             write(unit=*,fmt="(a)") " "
             write(unit=*,fmt="(a)") " "
             write(unit=*,fmt="(a)",advance="no") " Position: "

             read(*,'(a)') line
             if (len_trim(line)==0) exit
             line=adjustl(line)
             vet=0.0
             call getnum(line,vet,ivet,iv)
             if (iv == 3) then
                mlt=Get_Multip_Pos(vet,grp_espacial)
                occ=real(mlt)/real(grp_espacial%Multip)
                write(unit=*,fmt=*) " "
                write(*,'(a,i4,a,f3.6)') " Multiplicity: ",mlt, "     Occupancy(SHELX/FullProf) proportional to: ",occ
                call system('pause')
             end if
          end do
       end do

    End Subroutine Menu_Atom_1

end module Menu_3
