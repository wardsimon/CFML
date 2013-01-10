!!----
!!---- Menu: 1
!!---- Space Groups
!!----
!!
 Module Menu_1
    !---- Use File ----!
    use CFML_Symmetry_Tables
    use CFML_Crystallographic_Symmetry
    use CFML_String_Utilities,          only: u_case, ucase, getnum
    use CFML_Math_General,              only: equal_vector
    use CFML_Math_3D,                   only: resolv_sist_3x3

    !---- Variables ----!
    implicit none

 Contains

    !!----
    !!---- Subroutine Menu_Princ1
    !!----
    !!
    Subroutine Menu_Princ1()
       !---- Local Variables ----!
       character(len=2)  :: car

       do
          call system('cls')

          write(unit=*,fmt="(a)") "     CRYSTRALLOGRAPHIC CALCULATOR "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "     Space Groups Menu "
          write(unit=*,fmt="(a)") " ========================="
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " [0] Back..."
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " [ 1] Space Group Information (Complete information on the Screen)"
          write(unit=*,fmt="(a)") " [ 2] List of the full list of Symmetry Operators (File/Screen)"
          write(unit=*,fmt="(a)") " [ 3] Construct a Space Group from a set of Generators"
          write(unit=*,fmt="(a)") " [ 4] Hall Symbol from a set of Generators"
          write(unit=*,fmt="(a)") " [ 5] Comparison of Two Space Groups"
          write(unit=*,fmt="(a)") " [ 6] Determination of the Laue class and Point Group"
          write(unit=*,fmt="(a)") " [ 7] Determination of the Symbol for Symmetry Operators"
          write(unit=*,fmt="(a)") " [ 8] Conversions: IT -> Kovalev, Miller&Love, Zack, etc..."
          write(unit=*,fmt="(a)") " [ 9] Wyckoff Positions (Testing...)"
          !write(unit=*,fmt="(a)") " [10] Wyckoff for all Groups (Testing...)"
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)",advance="no") " OPTION: "
          read(*,'(a)') car
          if (len_trim(car) == 0) exit
          car=adjustl(car)

          select case (car)
             case ('0 ')
                exit

             case ('1 ')
                call Menu_SPGR_1()

             case ('2 ')
                call Menu_SPGR_2()

             case ('3 ')
                call Menu_SPGR_3()

             case ('4 ')
                call Menu_SPGR_4()

             case ('5 ')
                call Menu_SPGR_5()

             case ('6 ')
                call Menu_SPGR_6()

             case ('7 ')
                call Menu_SPGR_7()

             case ('8 ')
                call Menu_SPGR_8()

             case ('9 ','10')
                call Menu_SPGR_9()

             !case ('10')
             !   call Menu_SPGR_10()
          end select
       end do

    End Subroutine Menu_Princ1

    !!----
    !!---- Subroutine Menu_SPGR_1
    !!----
    !!
    Subroutine Menu_Spgr_1()
       !---- Local Variables ----!
       character(len=20)      :: line, spgr
       integer                :: i, iv, ierr, npos
       integer, dimension(1)  :: ivet
       real, dimension(1)     :: vet
       type (Space_Group_Type):: grp_espacial

       do
          call system('cls')
          write(unit=*,fmt="(a)") "                   GENERAL CRYSTALLOGRAPHY CALCULATOR "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "     Space Groups Information (Complete information on the Screen)"
          write(unit=*,fmt="(a)") "   =================================================================="
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)",advance="no") " Space Group (HM/Hall/Num): "

          read(*,'(a)') line
          if (len_trim(line)==0) exit
          line=adjustl(line)

          call set_spacegroup(line,grp_espacial)
          if(Err_Symm) then
            write(unit=*,fmt="(a)") trim(ERR_Symm_Mess)
            call system('pause')
            cycle
          end if
          call write_spacegroup(grp_espacial,full=.true.)
          write(unit=*,fmt=*) " "
          call system("pause ")

       end do
    End Subroutine Menu_Spgr_1

    !!----
    !!---- Subroutine Menu_SPGR_2
    !!----
    !!
    Subroutine Menu_Spgr_2()
       !---- Local Variables ----!
       character(len=20)                :: line, spgr
       character(len=40)                :: str
       character(len=80), dimension(96) :: texto

       integer                          :: i, iv, ierr, lun, nlines, npos
       integer, dimension(1)            :: ivet

       real, dimension(1)               :: vet

       type (Space_Group_type)          :: grp_espacial

       do
          call system('cls')
          write(unit=*,fmt="(a)") "              GENERAL CRYSTALLOGRAPHY CALCULATOR "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "     List of the full list of Symmetry Operators (File/Screen) "
          write(unit=*,fmt="(a)") " =================================================================="
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)",advance="no") " Space Group (HM/Hall/Num): "

          read(*,'(a)') line
          if (len_trim(line)==0) exit
          line=adjustl(line)

          call set_spacegroup(line,grp_espacial)
          if(Err_Symm) then
            write(unit=*,fmt="(a)") trim(ERR_Symm_Mess)
            call system('pause')
            cycle
          end if

          !---- OutPut ----!
          write(unit=*,fmt=*) " "
          write(unit=*,fmt="(a)",advance="no") " Output on Screen or File ([S]/F)? "
          read(unit=*,fmt='(a)') line
          call ucase(line)
          if (len_trim(line) == 0) line='S'
          select case (line(1:1))
             case ('S')
                lun=6
             case ('F')
                lun=1
                open(unit=lun,file="symoper.dat",status="replace",action="write")
          end select

          write(unit=lun,fmt='(a)') " Full list of Symmetry Operators for space group: "//trim(grp_espacial%SPG_Symb)
          write(unit=lun,fmt='(a)') " ============================================================"
          write(unit=lun,fmt='(a,i4,/)') " General Multiplicity: ",grp_espacial%multip

          nlines= 1
          texto = " "
          do i=1,grp_espacial%multip
             call Get_SymSymb(grp_espacial%symop(i)%rot, &
                              grp_espacial%symop(i)%tr,str)
             if (mod(i,2) == 0) then
                write(texto(nlines)(40:80),'(a,i3,a,a)') &
                                           ' => SYMM(',i,'): ',trim(str)
                nlines=nlines+1
             else
                write(texto(nlines)( 1:39),'(a,i3,a,a)')  &
                                           ' => SYMM(',i,'): ',trim(str)
             end if
          end do

          do i=1,nlines
             write(unit=lun,fmt='(a)') texto(i)
          end do

          if (lun == 1) then
             close(1)
          else
             call system("pause")
          end if

          exit
       end do

    End Subroutine Menu_Spgr_2

    !!----
    !!---- Subroutine Menu_SPGR_3
    !!----
    !!
    Subroutine Menu_Spgr_3()
       !---- Local Variables ----!
       character (len=80)                :: line
       character (len=1)                 :: ans
       character (len=30)                :: spgr, spg, str
       character (len=30), dimension(10) :: gen
       character (len=140)               :: gener

       integer, dimension(3,3,24)        :: ss
       integer                           :: ng,i,istart,nlines,ier

       real, dimension(3,24)             :: ts

       type (Space_Group_type)                :: grp_espacial


       do
          call system('cls')
          write(unit=*,fmt="(a)") "     GENERAL CRYSTALLOGRAPHY CALCULATOR "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "     Space Group from Generators (Jones Faithful representation)"
          write(unit=*,fmt="(a)") " ==================================================================="
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)",advance="no") " Generators from International Tables ?(Y/[N]): "
          read(*,"(a)") ans
          call ucase(ans)
          if (ans == 'Y' ) then

             write(unit=*,fmt="(a)",advance="no") " Give the space group symbol: "
             read(*,'(a)') spgr
             call Get_Generators(Spgr,Gener)
             if(err_symtab) then
               write(unit=*,fmt="(a)") " => "//ERR_SymTab_Mess
               write(*,*) " "
               call system("pause")
               cycle
             end if
              ng=0
              do

                i=index(Gener,";")
                if(i /= 0 ) then
                  ng=ng+1
                  gen(ng)=Gener(1:i-1)
                  Gener=Gener(i+1:)
                else
                  if(len_trim(Gener) /= 0) then
                    ng=ng+1
                    gen(ng)= gener
                    exit
                  end if
                end if

              end do
                do i=1,ng
                  write(unit=*,fmt="(a,i3,a)")"  => Generator #", i, ": "//gen(i)
                end do
          else

             write(unit=*,fmt=*) " "
             write(unit=*,fmt="(a)",advance="no") " Give the number of generators (max 10): "
             read(*,*) ng
             if (ng == 0) exit
             if (ng > 10) ng=10
             istart=1
             do i=1,ng
                write(unit=*,fmt='(a,i1,a)',advance="no") " -> Give the generator number ",i,": "
                read(unit=*,fmt='(a)') gen(i)
                call Read_Xsym(gen(i),istart,ss(:,:,i),ts(:,i))
             end do
          end if

          call set_spacegroup(spgr,grp_espacial,gen,ng,mode='gen')
          if(Err_Symm) then
            write(unit=*,fmt="(a)") trim(ERR_Symm_Mess)
            call system('pause')
            cycle
          end if
          call write_spacegroup(grp_espacial,full=.true.)
          write(*,*)" "
          call system("pause")

          exit
       end do

    End Subroutine Menu_Spgr_3

    !!----
    !!---- Subroutine Menu_SPGR_4
    !!----
    !!
    Subroutine Menu_Spgr_4()
       !---- Local Variables ----!
       character (len=80)                :: line
       character (len=1)                 :: ans
       character (len=30)                :: spgr, spg, str
       character (len=16)                :: hall
       character (len=30), dimension(10) :: gen
       character (len=140)               :: gener
       integer                           :: ng,i,istart,nlines,ier
       type (Space_Group_type)           :: grp_espacial

       do
          call system('cls')
          write(unit=*,fmt="(a)") "     GENERAL CRYSTALLOGRAPHY CALCULATOR "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "     Hall Symbol from a set of Generators"
          write(unit=*,fmt="(a)") " ============================================"
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)",advance="no") " Generators from International Tables ?(Y/[N]): "
          read(*,"(a)") ans
          call ucase(ans)
          if (ans == 'Y' ) then
             write(unit=*,fmt="(a)",advance="no") " Give the number or the space group symbol : "
             read(unit=*,fmt='(a)') spgr
             call Get_Generators(Spgr,Gener)
             if(err_symtab) then
               write(unit=*,fmt="(a)") " => "//ERR_SymTab_Mess
               write(*,*) " "
               call system("pause")
               cycle
             end if

              ng=0
              do

                i=index(Gener,";")
                if(i /= 0 ) then
                  ng=ng+1
                  gen(ng)=Gener(1:i-1)
                  Gener=Gener(i+1:)
                else
                  if(len_trim(Gener) /= 0) then
                    ng=ng+1
                    gen(ng)= gener
                    exit
                  end if
                end if

              end do
                do i=1,ng
                  write(unit=*,fmt="(a,i3,a)")"  => Generator #", i, ": "//gen(i)
                end do
          else
             write(unit=*,fmt=*) " "
             write(unit=*,fmt="(a)",advance="no") " Give the number of generators (max 10): "
             read(*,*) ng
             if (ng == 0) exit
             if (ng > 10) ng=10
             istart=1
             do i=1,ng
                write(unit=*,fmt='(a,i1,a)',advance="no") " -> Give the generator number ",i,": "
                read(*,'(a)') gen(i)
             end do
          end if
          call set_spacegroup(spgr,grp_espacial,gen,ng,'gen')
          if(Err_Symm) then
            write(unit=*,fmt="(a)") trim(ERR_Symm_Mess)
            call system('pause')
            cycle
          end if
          call Get_HallSymb_from_Gener(grp_espacial,hall)
          write(unit=*,fmt=*)  " "
          write(unit=*,fmt=*)  "  Calculated Hall Symbol : "//hall
          write(unit=*,fmt=*) " "
          call system("pause")

          exit
       end do

    End Subroutine Menu_Spgr_4

    !!----
    !!---- Subroutine Menu_SPGR_5
    !!----
    !!
    Subroutine Menu_Spgr_5()
       !---- Local Variables ----!
       character (len=80)                :: line,spg
       character (len=20)                :: spgr
       character (len= 1)                :: car
       character (len=30), dimension(10) :: gen
       character (len=140)               :: gener

       integer                           :: i, ng, istart, npos, iv, ierr
       integer, dimension(1)             :: ivet
       integer, dimension(3,3,24)        :: ss
       integer, parameter                :: num_spgr_info=612

       real, dimension(1)                :: vet
       real, dimension(3,24)             :: ts

       type (Space_Group_type)           :: grp_espacial1, grp_espacial2

       do
          call system('cls')
          write(unit=*,fmt="(a)") "     GENERAL CRYSTALLOGRAPHY CALCULATOR "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "       Comparison of Two Space Groups"
          write(unit=*,fmt="(a)") "   ======================================"
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " [1] Generators <-> Generators"
          write(unit=*,fmt="(a)") " [2] Generators <-> International Tables"
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)",advance="no") " OPTION: "
          read(*,'(a)') car
          if (len_trim(car) == 0) exit
          car=adjustl(car)

          select case (car)
             case('1')
                write(unit=*,fmt=*) " "
                write(unit=*,fmt="(a)",advance="no") " Generators from International Tables SET1 ?(Y/[N]): "
                read(unit=*,fmt='(a)') car
                call ucase(car)
                if (car == 'Y' ) then
                  write(unit=*,fmt="(a)",advance="no") " Give the space group symbol (SET1)  : "
                  read(unit=*,fmt='(a)') spgr
                   call Get_Generators(Spgr,Gener)
                   if(err_symtab) then
                     write(unit=*,fmt="(a)") " => "//ERR_SymTab_Mess
                     write(*,*) " "
                     call system("pause")
                     cycle
                   end if
                   ng=0
                   do
                     i=index(Gener,";")
                     if(i /= 0 ) then
                       ng=ng+1
                       gen(ng)=Gener(1:i-1)
                       Gener=Gener(i+1:)
                     else
                       if(len_trim(Gener) /= 0) then
                         ng=ng+1
                         gen(ng)= gener
                         exit
                       end if
                     end if
                   end do
                   do i=1,ng
                     write(unit=*,fmt="(a,i3,a)")"  => Generator #", i, ": "//gen(i)
                   end do
                else
                   write(unit=*,fmt=*) " "
                   write(unit=*,fmt="(a)",advance="no") " Give the number of generators SET 1 (max 10): "
                   read(unit=*,fmt=*) ng
                   if (ng == 0) exit
                   if (ng > 10) ng=10
                   istart=1
                   do i=1,ng
                      write(unit=*,fmt='(a,i1,a)',advance="no") " -> Give the generator number ",i,": "
                      read(unit=*,fmt='(a)') gen(i)
                      call Read_Xsym(gen(i),istart,ss(:,:,i),ts(:,i))
                   end do
                end if

                call set_spacegroup(spgr,grp_espacial1,gen,ng,mode='gen ')
                if(Err_Symm) then
                  write(unit=*,fmt="(a)") trim(ERR_Symm_Mess)
                  call system('pause')
                  cycle
                end if

                write(unit=*,fmt=*) " "
                write(unit=*,fmt="(a)",advance="no") " Generators from International Tables SET2 ?(Y/[N]): "
                read(*,'(a)') car
                call ucase(car)
                if (car == 'Y' ) then
                   write(unit=*,fmt="(a)",advance="no") " Give the space group symbol (SET2)  : "
                   read(unit=*,fmt='(a)') spgr
                   call Get_Generators(Spgr,Gener)
                   if(err_symtab) then
                     write(unit=*,fmt="(a)") " => "//ERR_SymTab_Mess
                     write(*,*) " "
                     call system("pause")
                     cycle
                   end if
                   ng=0
                   do
                     i=index(Gener,";")
                     if(i /= 0 ) then
                       ng=ng+1
                       gen(ng)=Gener(1:i-1)
                       Gener=Gener(i+1:)
                     else
                       if(len_trim(Gener) /= 0) then
                         ng=ng+1
                         gen(ng)= gener
                         exit
                       end if
                     end if
                   end do
                   do i=1,ng
                     write(unit=*,fmt="(a,i3,a)")"  => Generator #", i, ": "//gen(i)
                   end do
                else
                   write(unit=*,fmt=*) " "
                   write(unit=*,fmt="(a)",advance="no") " Give the number of generators SET 2 (max 10): "
                   read(unit=*,fmt=*) ng
                   if (ng == 0) exit
                   if (ng > 10) ng=10
                   istart=1
                   do i=1,ng
                      write(unit=*,fmt='(a,i1,a)',advance="no") " -> Give the generator number ",i,": "
                      read(*,'(a)') gen(i)
                      call Read_Xsym(gen(i),istart,ss(:,:,i),ts(:,i))
                   end do
                end if
                call set_spacegroup(spgr,grp_espacial2,gen,ng,mode='gen ')

             case('2')
                write(unit=*,fmt=*) " "
                write(unit=*,fmt="(a)",advance="no") " Generators from International Tables SET1 ?(Y/[N]): "
                read(unit=*,fmt='(a)') car
                call ucase(car)
                if (car == 'Y' ) then
                  write(unit=*,fmt="(a)",advance="no") " Give the space group symbol (SET1)  : "
                  read(unit=*,fmt='(a)') spgr
                   call Get_Generators(Spgr,Gener)
                   if(err_symtab) then
                     write(unit=*,fmt="(a)") " => "//ERR_SymTab_Mess
                     write(*,*) " "
                     call system("pause")
                     cycle
                   end if
                   ng=0
                   do
                     i=index(Gener,";")
                     if(i /= 0 ) then
                       ng=ng+1
                       gen(ng)=Gener(1:i-1)
                       Gener=Gener(i+1:)
                     else
                       if(len_trim(Gener) /= 0) then
                         ng=ng+1
                         gen(ng)= gener
                         exit
                       end if
                     end if
                   end do
                   do i=1,ng
                     write(unit=*,fmt="(a,i3,a)")"  => Generator #", i, ": "//gen(i)
                   end do
                else
                   write(unit=*,fmt=*) " "
                   write(unit=*,fmt="(a)",advance="no") " Give the number of generators SET 1 (max 10): "
                   read(unit=*,fmt=*) ng
                   if (ng == 0) exit
                   if (ng > 10) ng=10
                   istart=1
                   do i=1,ng
                      write(unit=*,fmt='(a,i1,a)',advance="no") " -> Give the generator number ",i,": "
                      read(unit=*,fmt='(a)') gen(i)
                      call Read_Xsym(gen(i),istart,ss(:,:,i),ts(:,i))
                   end do
                end if
                call set_spacegroup(spgr,grp_espacial1,gen,ng,mode='gen ')

                write(unit=*,fmt=*) " "
                write(unit=*,fmt="(a)",advance="no") " Space Group SET2 (HM/Hall/Num): "
                read(*,'(a)') line
                if (len_trim(line)==0) exit
                line=adjustl(line)
                spgr=" "

                !---- Introduce Numero? ----!
                call getnum(line,vet,ivet,iv)
                if (iv == 1) then
                   do i=1,num_spgr_info
                      if (ivet(1) /= spgr_info(i)%n) cycle
                      spgr=spgr_info(i)%hm
                      call set_spacegroup(spgr,grp_espacial2,mode='hm')
                      exit
                   end do
                else
                   !---- Introduce notacion de HM ? ----!
                   do i=1,num_spgr_info
                      if (u_case(line(1:12)) /= spgr_info(i)%HM) cycle
                      spgr=spgr_info(i)%hm
                      npos=index(spgr,':')
                      if (npos == 0) then
                         call set_spacegroup(spgr,grp_espacial2,mode='hm')
                      else
                         spgr=spgr_info(i)%hall
                         call set_spacegroup(spgr,grp_espacial2,mode='hall')
                      end if
                      exit
                   end do

                   !---- Introduce notacion de Hall ? ----!
                   if (len_trim(spgr) == 0) then
                      do i=1,num_spgr_info
                         if (u_case(line(1:16)) /= u_case(spgr_info(i)%Hall) ) cycle
                         spgr=spgr_info(i)%hall
                         call set_spacegroup(spgr,grp_espacial2,mode='hall')
                         exit
                      end do
                   end if
                end if
                if (len_trim(spgr) == 0) cycle

             case default
                write(unit=*,fmt=*)  " Wrong Option"
                cycle
          end select

          if (SpGr_Equal(grp_espacial1, grp_espacial2)) then
             write(unit=*,fmt=*) " "
             write(unit=*,fmt=*) "  These Space Groups are EQUAL"
          else
             write(unit=*,fmt=*) " "
             write(unit=*,fmt=*) "  These Space Groups are NOT EQUAL"
          end if
          write(unit=*,fmt=*) " "
          call system("pause")

       end do

    End Subroutine Menu_Spgr_5

    !!----
    !!---- Subroutine Menu_SPGR_6
    !!----
    !!
    Subroutine Menu_Spgr_6()
       !---- Local Variables ----!
       character(len=20)     :: spgr
       character(len=5 )     :: laue_car, point_car
       type (Space_Group_type)    :: grp_espacial

       do
          call system('cls')
          write(unit=*,fmt="(a)") "     GENERAL CRYSTALLOGRAPHY CALCULATOR "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "     Determination of the Laue class and Point Group "
          write(unit=*,fmt="(a)") " ======================================================="
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)",advance="no") " Space Group (HM/Hall/Num): "

          read(*,'(a)') spgr
          if (len_trim(spgr)==0) exit
          spgr=adjustl(spgr)

          call set_spacegroup(spgr,grp_espacial)
          if(Err_Symm) then
            write(unit=*,fmt="(a)") trim(ERR_Symm_Mess)
            call system('pause')
            cycle
          end if

          call Get_Laue_PG(grp_espacial,laue_car,point_car)

          !---- Print Information ----!
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "             Table                       Calculated"
          write(unit=*,fmt="(a)") "    Laue Group   Point Group       Laue Group   Point Group"
          write(unit=*,fmt="(a)") "   ========================================================="
          write(*,'(6x,a5,10x,a5,13x,a5,10x,a5)') grp_espacial%laue,grp_espacial%pg, &
                                                 laue_car,point_car
          write(unit=*,fmt=*) " "
          call system("pause")

       end do

    End Subroutine Menu_Spgr_6

    !!----
    !!---- Subroutine Menu_SPGR_7
    !!----
    !!
    Subroutine Menu_Spgr_7()
       !---- Local Variables ----!
       character (len=80)                :: line,simbolo,symb
       character (len=12)                :: spgr,spg
       character (len= 1)                :: ans
       character (len=30), dimension(10) :: gen
       character (len=140)               :: gener

       integer, dimension(3,3,24)        :: ss
       integer                           :: ng,i,istart,ier

       real, dimension(3,24)             :: ts

       do
          call system('cls')
          write(unit=*,fmt="(a)") "     GENERAL CRYSTALLOGRAPHY CALCULATOR "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "     Derivation symbols for symmetry operations"
          write(unit=*,fmt="(a)") " =================================================="
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)",advance="no") " Generators from International Tables ?(Y/[N]): "
          read(*,"(a)") ans
          call ucase(ans)
          if (ans == 'Y' ) then
             write(unit=*,fmt="(a)",advance="no") " Give the space group symbol : "
             read(unit=*,fmt='(a)') spgr
             call Get_Generators(Spgr,Gener)
             if(err_symtab) then
               write(unit=*,fmt="(a)") " => "//ERR_SymTab_Mess
               write(*,*) " "
               call system("pause")
               cycle
             end if
             ng=0
             do
               i=index(Gener,";")
               if(i /= 0 ) then
                 ng=ng+1
                 gen(ng)=Gener(1:i-1)
                 Gener=Gener(i+1:)
               else
                 if(len_trim(Gener) /= 0) then
                   ng=ng+1
                   gen(ng)= gener
                   exit
                 end if
               end if
             end do
             do i=1,ng
               write(unit=*,fmt="(a,i3,a)")"  => Generator #", i, ": "//gen(i)
               call Read_Xsym(gen(i),1,ss(:,:,i),ts(:,i))
             end do

          else
             write(unit=*,fmt=*) " "
             write(unit=*,fmt="(a)",advance="no") " Give the number of generators (max 10): "
             read(unit=*,fmt=*) ng
             if (ng == 0) exit
             if (ng > 10) ng=10
             istart=1
             do i=1,ng
                write(unit=*,fmt='(a,i1,a)',advance="no") " -> Give the generator number ",i,": "
                read(unit=*,fmt='(a)') gen(i)
                call Read_Xsym(gen(i),istart,ss(:,:,i),ts(:,i))
             end do
          end if

          !---- Testeo de rutina ----!
          do i=1,ng
             call Get_SymSymb(ss(:,:,i),ts(:,i),simbolo)
             write(unit=*,fmt=*) " "
             write(unit=*,fmt='(a,i1,5x,a)') " Generator: ",i, ": "//gen(i)
             call Symmetry_Symbol(gen(i),Symb)
             write(unit=*,fmt='(a,6x,a)')    "    Symbol: ",trim(simbolo)//"   "//trim(symb)
          end do
          write(*,*) " "
          call system("pause")
       end do

    End Subroutine Menu_Spgr_7

    !!----
    !!---- Subroutine Menu_SPGR_8
    !!----
    !!
    Subroutine Menu_Spgr_8()
       !---- Local Variables ----!
       character(len=20)     :: spgr, cmd
       character(len=35)     :: Appl_trs, IAppl_trs
       type (Space_Group_type)    :: grp_espacial, grp_modified

       call Set_Spgr_Info()

       do
          call system('cls')
          write(unit=*,fmt="(a)") "     GENERAL CRYSTALLOGRAPHY CALCULATOR "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "     Conversion of Symmetry Operators of Space Groups"
          write(unit=*,fmt="(a)") "          ONLY CALCULATIONS ON SYMMETRY OPERATORS"
          write(unit=*,fmt="(a)") " ========================================================"
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)",advance="no") " Space Group (HM/Hall/Num): "

          read(*,'(a)') spgr
          if (len_trim(spgr)==0) exit
          spgr=adjustl(spgr)

          call set_spacegroup(spgr,grp_espacial)
          if(Err_Symm) then
            write(unit=*,fmt="(a)") trim(ERR_Symm_Mess)
            call system('pause')
            cycle
          end if

          !---- Systematic Conversions ----!
          call write_spacegroup(grp_espacial)
          write(unit=*,fmt=*) " "
          call system("pause")

          write(unit=*,fmt=*) " "
          write(unit=*,fmt=*) " IT -> ML"
          grp_modified=grp_espacial
          Call Setting_Change('IT','ML',grp_modified,Appl_trs, IAppl_trs)
          write(unit=*,fmt=*) "            Applied Transformation: ",Appl_trs
          write(unit=*,fmt=*) " Inverse of Applied Transformation: ",IAppl_trs
          call write_spacegroup(grp_modified)
          write(unit=*,fmt=*) " "
          call system("pause")

          write(unit=*,fmt=*) " "
          write(unit=*,fmt=*) " IT -> KO"
          grp_modified=grp_espacial
          Call Setting_Change('IT','KO',grp_modified,Appl_trs, IAppl_trs)
          write(unit=*,fmt=*) "            Applied Transformation: ",Appl_trs
          write(unit=*,fmt=*) " Inverse of Applied Transformation: ",IAppl_trs
          call write_spacegroup(grp_modified)
          write(unit=*,fmt=*) " "
          call system("pause")

          write(unit=*,fmt=*) " "
          write(unit=*,fmt=*) " IT -> ZA"
          grp_modified=grp_espacial
          Call Setting_Change('IT','ZA',grp_modified,Appl_trs, IAppl_trs)
          write(unit=*,fmt=*) "            Applied Transformation: ",Appl_trs
          write(unit=*,fmt=*) " Inverse of Applied Transformation: ",IAppl_trs
          call write_spacegroup(grp_modified)
          write(unit=*,fmt=*) " "
          call system("pause")

          write(unit=*,fmt=*) " "
          write(unit=*,fmt=*) " IT -> BC"
          grp_modified=grp_espacial
          Call Setting_Change('IT','BC',grp_modified,Appl_trs, IAppl_trs)
          write(unit=*,fmt=*) "            Applied Transformation: ",Appl_trs
          write(unit=*,fmt=*) " Inverse of Applied Transformation: ",IAppl_trs
          call write_spacegroup(grp_modified)
          write(unit=*,fmt=*) " "
          call system("pause")

          exit

       end do

    End Subroutine Menu_Spgr_8

    !!----
    !!---- Subroutine Menu_SPGR_9
    !!----
    !!
    Subroutine Menu_Spgr_9()
       !---- Local Variables ----!
       character(len=40)      :: line, spgr
       character(len=3)       :: car
       integer                :: i
       type (Space_Group_Type):: grp_espacial

       call Set_Spgr_Info()

       do
          call system('cls')
          write(unit=*,fmt="(a)") "     GENERAL CRYSTALLOGRAPHY CALCULATOR "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "             Wyckoff Information "
          write(unit=*,fmt="(a)") "        ==========================="
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
          call Get_Wyckoff(grp_espacial)
           write(unit=*,fmt=*) " "
          call system("pause ")
          exit
       end do

       return
    End Subroutine Menu_Spgr_9

    !!----
    !!---- Subroutine Menu_SPGR_10
    !!----
    !!
    Subroutine Menu_Spgr_10()
       !---- Local Variables ----!
       character(len=40)      :: line, spgr
       character(len=3)       :: car
       integer                :: i,ier
       type (Space_Group_Type):: grp_espacial

       call Set_Spgr_Info()

       do
          call system('cls')
          write(unit=*,fmt="(a)") "     GENERAL CRYSTALLOGRAPHY CALCULATOR "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "             Wyckoff Information "
          write(unit=*,fmt="(a)") "         ==========================="
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " "

          open(unit=1,file="itgen.spg",status="unknown",iostat=ier)
          if (ier /= 0) then
             write(unit=*,fmt="(a)") " "
             write(unit=*,fmt="(a)") " The file ITGEN.SPG is not present or some problems are present"
             write(unit=*,fmt="(a)") " "
             exit
          end if
          do
             read(1,'(a)',iostat=ier) line
             if (ier /= 0) exit
             if (line(1:4) /= "SPGR") cycle
             read(line(12:40),'(a)') spgr

             write(unit=*,fmt="(a)")  " Grupo Espacial: "//spgr

             call set_spacegroup(spgr,grp_espacial)
             call Get_Wyckoff(grp_espacial)

          end do

          close(unit=1)
           write(unit=*,fmt=*) " "
          call system("pause ")
          exit
       end do

       return
    End Subroutine Menu_Spgr_10


    !!----
    !!---- Subroutine Get_Wyckoff
    !!----
    !!
    Subroutine Get_Wyckoff(grp)
       !---- Arguments ----!
       type (Space_Group_Type)           :: grp

       !---- Local Variables ----!
       character (len=40)                :: carsymb, wyckoffcar, wyckoffcar1, wyckoffcar2
       integer                           :: i,j,k,n,m,ii,jj,kk
       integer                           :: iunit, num_orbita, nvar, npos, num_old, nzeros
       integer                           :: nx,ny,nz
       integer,dimension(3)              :: ix,ix1,ix2,ix11,ix22
       integer,dimension(3,3)            :: w,w1
       real,dimension(3)                 :: t,t1,ts,ts1,ts2,x,x1,x2
       type(wyckoff_type)                :: Wyck1
       type(wyck_pos_type)               :: wyckpos

       iunit=3
       open(unit=iunit,file="Wyckoff_Positions.inf",status="replace",action="write")
       num_orbita=0
       wyck1=wyckoff_Type(0,wyck_pos_type(0,' ',0,' ',' '))

       !---- Total Operators except Identity ----!
       do i=2,grp%multip
          call symmetry_symbol(grp%symop(i),carsymb)
          carsymb=adjustl(carsymb)

          !---- No traslation components ----!
          select case (carsymb(1:1))
             case ('t')
                cycle

             case ('a', 'b', 'c', 'd', 'n', 'g')
                cycle
          end select
          n=index(carsymb,'(')
          if (n /= 0) cycle

          !---- (R-1)X=0 ----!
          w=grp%symop(i)%rot
          do j=1,3
            w(j,j)=w(j,j)-1
          end do

          !---- Cell traslations ----!
          do ii=0,1
             do jj=0,1
                uno:do kk=0,1

                   t(1)=grp%symop(i)%tr(1) + real(ii)
                   t(2)=grp%symop(i)%tr(2) + real(jj)
                   t(3)=grp%symop(i)%tr(3) + real(kk)
                   call resolv_sist_3x3(w,-t,ts,x,ix)
                   ts=mod(ts+10.0,1.0)

                   !---- write Wyckoff String ----!
                   if (any(abs(x) > 0.0 .and. abs(x) < 1.0)) x=2.0*x
                   call get_string_resolv(ts,x,ix,wyckoffcar)

                   !---- Reduce Wyckoff String ----!
                   call Wyckoff_represen(ts,x,ix)
                   call get_string_resolv(ts,x,ix,wyckoffcar1)

                   if (num_orbita /= 0) then
                      !---- Searching in the actual List ----!
                      do n=1,num_orbita
                         do j=1,wyck1%orbit(n)%norb
                            if (wyckoffcar1 == wyck1%orbit(n)%str_orbit(j)) cycle uno
                         end do
                      end do

                      !---- Lattice traslations ----!
                      do j=1,grp%numlat

                         !---- Control 1 ----!
                         ts1=ts+grp%latt_trans(:,j)
                         ts1=mod(ts1+10.0,1.0)
                         x1=x
                         ix1=ix
                         call wyckoff_represen(ts1,x1,ix1)
                         call get_string_resolv(ts1,x1,ix1,wyckoffcar2)
                         do n=1,num_orbita
                            do k=1,wyck1%orbit(n)%norb
                               if (wyckoffcar2 == wyck1%orbit(n)%str_orbit(k)) cycle uno
                            end do
                         end do

                         !---- Control 1A ----!
                         nzeros=count(ix==0)
                         select case (nzeros)
                            case (2)
                               do n=1,3
                                  if (ix1(n) ==0) cycle
                                  x1(n)=-x1(n)
                               end do

                            case (1)
                            case (0)
                         end select
                         call wyckoff_represen(ts1,x1,ix1)
                         call get_string_resolv(ts1,x1,ix1,wyckoffcar2)
                         do n=1,num_orbita
                            do k=1,wyck1%orbit(n)%norb
                               if (wyckoffcar2 == wyck1%orbit(n)%str_orbit(k)) cycle uno
                            end do
                         end do

                         !---- Control 2 ----!
                         nzeros=count(ix==0)
                         if (nzeros /=3) then
                            ts1=-ts+grp%latt_trans(:,j)
                            ts1=mod(ts1+10.0,1.0)
                            x1=-x
                            ix1=ix
                            call get_string_resolv(ts1,x1,ix1,wyckoffcar2)
                            do n=1,num_orbita
                               do k=1,wyck1%orbit(n)%norb
                                  if (wyckoffcar2 == wyck1%orbit(n)%str_orbit(k)) cycle uno
                               end do
                            end do

                            !---- Control 2A ----!
                            select case (nzeros)
                               case (2)
                                  do n=1,3
                                     if (ix1(n) ==0) cycle
                                     x1(n)=-x1(n)
                                  end do

                               case (1)
                               case (0)
                            end select
                            call wyckoff_represen(ts1,x1,ix1)
                            call get_string_resolv(ts1,x1,ix1,wyckoffcar2)
                            do n=1,num_orbita
                               do k=1,wyck1%orbit(n)%norb
                                  if (wyckoffcar2 == wyck1%orbit(n)%str_orbit(k)) cycle uno
                               end do
                            end do
                         end if

                         !---- Control 3----!
                         ts1=ts+grp%latt_trans(:,j)
                         ts1=mod(ts1+10.0,1.0)
                         x1=x
                         ix1=ix
                         call wyckoff_represen(ts1,x1,ix1)
                         nvar=count(ix1 /=0)
                         nx=count(ix1 == 1)
                         ny=count(ix1 == 2)
                         nz=count(ix1 == 3)

                         wyckpos=wyck_pos_type(0,' ',0,' ',' ')
                         call get_string_resolv(ts1,x1,ix1,wyckoffcar2)
                         call wyckoff_Orbit(Grp,WyckoffCar2,wyckpos%norb,wyckpos%str_orbit)

                         do n=1,num_orbita
                            !if (wyckpos%norb /= wyck1%orbit(n)%norb) cycle

                            do k=1,wyck1%orbit(n)%norb
                               call Read_RSymm(wyck1%orbit(n)%str_orbit(k),ts2,x2,ix2)

                               select case (nvar)
                                  case (0) ! Punto fijo

                                  case (1) ! (v,f,f), (f,v,f), (f,f,v)
                                     if (equal_vector(ix1,ix2,3) ) then
                                        do m=1,3
                                           if (ix1(m) ==0) cycle
                                           ts2(m)=0.0
                                        end do
                                        if (equal_vector(ts1,ts2,3) ) cycle uno
                                     end if

                                     if (equal_vector(ts1,ts2,3) ) then
                                        ix11=0
                                        ix22=0
                                        do m=1,3
                                           if (ix1(m) /=0) ix11(m)=1
                                           if (ix2(m) /=0) ix22(m)=1
                                        end do
                                        if (equal_vector(ix11,ix22,3) ) cycle uno
                                     end if

                                  case (2) ! (v,v,f), (v,f,v), (f,v,v)
                                     if (nx==2 .or. ny==2 .or. nz==2) then
                                        if (equal_vector(ts1,ts2,3) ) then
                                           ix11=0
                                           ix22=0
                                           do m=1,3
                                              if (ix1(m) /=0) ix11(m)=1
                                              if (ix2(m) /=0) ix22(m)=1
                                           end do
                                           if (equal_vector(ix11,ix22,3) ) cycle uno
                                        end if

                                        if (equal_vector(ix1,ix2,3) ) then
                                           do m=1,3
                                              if (ix1(m) ==0) cycle
                                              ts2(m)=0.0
                                           end do
                                           if (equal_vector(ts1,ts2,3) ) cycle uno
                                        end if
                                     end if

                                     if (nx==1 .or. ny==1 .or. nz==1) then

                                        if (equal_vector(ix1,ix2,3) ) then
                                           do m=1,3
                                              if (ix1(m) ==0) cycle
                                              ts2(m)=0.0
                                           end do
                                           if (equal_vector(ts1,ts2,3) ) cycle uno
                                        end if

                                        if (equal_vector(ts1,ts2,3) ) then
                                           ix11=0
                                           ix22=0
                                           do m=1,3
                                              if (ix1(m) /=0) ix11(m)=1
                                              if (ix2(m) /=0) ix22(m)=1
                                           end do
                                           if (equal_vector(ix11,ix22,3) ) cycle uno
                                        end if
                                     end if

                                  case (3) ! (v,v,v)
                               end select

                            end do
                         end do

                      end do

                   end if

                   !---- Add Orbite ----!
                   num_orbita=num_orbita+1
                   wyck1%num_orbit=num_orbita
                   call wyckoff_Orbit(Grp,WyckoffCar1,wyck1%orbit(num_orbita)%norb,wyck1%orbit(num_orbita)%str_orbit)
                   wyck1%Orbit(num_orbita)%multp=wyck1%Orbit(num_orbita)%norb*grp%numlat
                end do uno
             end do
          end do

       end do

       !----
       !---- Special points could be added to Wyckoff Positions ----!
       !----
       num_old=num_orbita

       do i=1, num_old
          if (wyck1%orbit(i)%multp ==1) cycle
          call Read_Rsymm(wyck1%orbit(i)%str_orbit(1),t,x,ix)

          !---- Cuantas columnas son diferentes de zero ----!
          nvar=count(ix /= 0)
          select case (nvar)
            case (0)
               !---- Fixed Point ----!
               ! No hay que hacer nada

            case (1)
               !---- One Free Variable ----!
               do j=1,3
                  if (ix(j) /= 0) exit
               end do

               t1=t
               dos: do k=0,7
                  t1(j)=real(k)/8.0

                  npos=0
                  do n=1, grp%multip
                     ts1=ApplySO(grp%symop(n),t1)
                     ts1=mod(ts1+10.0,1.0)
                     if (equal_vector(t1,ts1,3)) npos=npos+1
                  end do
                  npos=grp%multip/npos
                  if (npos >= wyck1%orbit(i)%multp) cycle

                  x1=0.0
                  ix1=0
                  call get_string_resolv(t1,x1,ix1,wyckoffcar)

                  !---- Searching in the actual List ----!
                  do n=1,num_orbita
                     do jj=1,wyck1%orbit(n)%norb
                        if (wyckoffcar == wyck1%orbit(n)%str_orbit(jj)) cycle dos
                     end do
                  end do

                  !---- Lattice traslations ----!
                  do jj=1,grp%numlat
                     !---- Control 1 ----!
                     ts1=t1+grp%latt_trans(:,jj)
                     ts1=mod(ts1+10.0,1.0)
                     x1=0.0
                     ix1=0
                     call get_string_resolv(ts1,x1,ix1,wyckoffcar1)
                     do n=1,num_orbita
                        do kk=1,wyck1%orbit(n)%norb
                           if (wyckoffcar1 == wyck1%orbit(n)%str_orbit(kk)) cycle dos
                        end do
                     end do
                  end do

                  num_orbita=num_orbita+1
                  wyck1%num_orbit=num_orbita
                  call wyckoff_Orbit(Grp,wyckoffCar,wyck1%orbit(num_orbita)%norb,wyck1%orbit(num_orbita)%str_orbit)
                  wyck1%Orbit(num_orbita)%multp=wyck1%Orbit(num_orbita)%norb*grp%numlat
               end do dos

            case (2)
               !---- Two Free Variables ----!

            case (3)
               !---- Three Free Variables ----!

         end select
      end do

      !---- Print ----!
      !call Write_Wyckoff(Wyck1,Grp%SPG_Symb,iunit)
      if (grp%wyckoff%num_orbit /= wyck1%num_orbit) then
         call Write_Wyckoff(grp%wyckoff,Grp%SPG_Symb,iunit)
         call Write_Wyckoff(Wyck1,Grp%SPG_Symb,iunit)
      end if
      close(unit=iunit)

      return
   end subroutine Get_Wyckoff


   !!----
   !!---- Subroutine
   !!----
   Subroutine Wyckoff_Represen(ts,x,ix)
      !---- Arguments ----!
      real, dimension(3),    intent(in out) :: ts
      real, dimension(3),    intent(in out) :: x
      integer, dimension(3), intent(in out) :: ix

      !---- Local variables ----!
      integer :: i, nzeros

      !---- (x+1/2,0,0)  -> (x,0,0) ----!
      !---- (x,-y+1/2,0) -> (x,y,0) ----!
      do i=1,3
         if (ix(i) /=0) then
            x(i)=1.0
            ts(i)=0.0
         end if
      end do

      nzeros=count(ix==0)
      select case (nzeros)
         case (0)
         case (1)
         case (2)
            do i=1,3
               if (ix(i) == 0) cycle
               ix(i)=i
               x(i)=1.0
            end do

         case (3)
      end select

      return
   End Subroutine Wyckoff_Represen


   !!----
   !!----
   !!----
   Subroutine Read_RSymm(Symb,T,X,Ix)
      !---- Arguments ----!
      character(len=*),      intent(in ) :: Symb
      real, dimension(3),    intent(out) :: T
      real, dimension(3),    intent(out) :: X
      integer, dimension(3), intent(out) :: Ix

      !---- Local Variables ----!
      integer                   :: i,j,k
      integer,dimension(3,3)    :: w

      t=0.0
      x=0.0
      ix=0

      call Read_Xsym(Symb,1,w,t)
      do i=1,3
         do j=1,3
            if (w(i,j) == 0) cycle
            ix(i)=j
            x(i)=real(w(i,j))
         end do
      end do

      return
   End Subroutine Read_RSymm

 End module Menu_1
