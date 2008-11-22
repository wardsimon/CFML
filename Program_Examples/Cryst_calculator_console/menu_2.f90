!!----
!!---- Menu: 2
!!---- Reflections
!!----
!!
 Module Menu_2
   !---- Use File ----!
   use Crystallographic_Symmetry
   use Reflections_Utilities
   use Crystal_Types, only: Crystal_Cell_Type, Write_Crystal_Cell

   !---- Variables ----!
   implicit none

 Contains

    !!----
    !!---- Subroutine Menu_Princ2
    !!----
    !!
    Subroutine Menu_Princ2()
       !---- Local Variables ----!
       character(len=2)  :: car

       do
          call system('cls')

          print*,"     GENERAL CRYSTRALLOGRAPHY CALCULATOR "
          print*," "
          print*,"     Reflections Menu "
          print*," ========================"
          print*," "
          print*," "
          print*," [0] Back..."
          print*," "
          print*," [1] Asymmetric Unit in Reciprocal Space"
          print*," [2] Change an arbitrary reflection to the Asymmetric Unit"
          print*," [3] Systematic Absences of Reflections"
          print*," [4] List of Equivalent Reflections"
          print*," [5] List of unique reflections within a range of S or D"
          print*," [6] Multiplicity of reflections"
          print*," "
          print*," OPTION: "
          read(*,'(a)') car
          if (len_trim(car) == 0) exit
          car=adjustl(car)

          select case (car)
             case ('0 ')
                exit

             case ('1 ')
                call Menu_REFL_1()

             case ('2 ')
                call Menu_REFL_2()

             case ('3 ')
                call Menu_REFL_3()

             case ('4 ')
                call Menu_REFL_4()

             case ('5 ')
                call Menu_REFL_5()

             case ('6 ')
                call Menu_REFL_6()

          end select
       end do

    End Subroutine Menu_Princ2

    !!----
    !!---- Subroutine Menu_REFL_1
    !!----
    !!
    Subroutine Menu_Refl_1()
       !---- Local Variables ----!
       character(len=20)       :: line, spgr
       integer                 :: i, iv, ierr, npos
       integer, dimension(1)   :: ivet
       integer, dimension(3)   :: h,k
       real, dimension(1)      :: vet
       type (Space_Group_type) :: grp_espacial

       do
          call system('cls')
          print*,"     GENERAL CRYSTRALLOGRAPHY CALCULATOR "
          print*," "
          print*,"     Asymmetric Unit in Reciprocal Space"
          print*," ==========================================="
          print*," "
          print*," "
          print*," Space Group (HM/Hall/Num): "

          read(*,'(a)') line
          if (len_trim(line)==0) exit
          spgr=adjustl(line)
          call set_spacegroup(spgr,grp_espacial)
          call write_asu(grp_espacial)
          print*," "
          call system("pause ")
          do
             print*," Enter a reflection [3 integers] <cr> -> stops:"
             read(*,'(a)') line
             if (len_trim(line)==0) exit
             read(unit=line,fmt=*,iostat=ierr) h
             if(ierr /= 0) cycle
             k=Get_Hequiv_Asu(h,grp_espacial)
             write(unit=*,fmt="(2(a,3i4))") " Input reflection: ",h,"   Equivalent reflection in A.U.: ",k
          end do
       end do
    End Subroutine Menu_Refl_1

    !!----
    !!---- Subroutine Menu_REFL_2
    !!----
    !!
    Subroutine Menu_Refl_2()
       !---- Local Variables ----!
       logical                                          :: info
       character(len=20)                                :: spgr,line
       integer, dimension(3)                            :: h,k,nulo
       integer, dimension(1)                            :: ivet
       integer                                          :: i,j,num,iv,ierr,npos
       real                                             :: fase,fase1,fase2
       real, dimension(1)                               :: vet
       type (Space_Group_type)                          :: grp_espacial
       type (Reflection_type),allocatable, dimension(:) :: reflexiones

       nulo=0
       info=.true.
       num=0

       do
          call system('cls')
          print*," "
          print*,"     GENERAL CRYSTRALLOGRAPHY CALCULATOR "
          print*," "
          print*,"     Put Reflections in the Asymmetric Unit"
          print*," =============================================="
          print*," "
          print*," "

          if (info) then
             print*," Space Group (HM/Hall/Num): "

             read(*,'(a)') line
             if (len_trim(line)==0) exit
             line=adjustl(line)
             call set_spacegroup(line,grp_espacial)
             info=.false.
          end if

          print*," "
          print*,"  Reflection hkl and phase: "
          read(*,*) h(1),h(2),h(3),fase
          print*," "
          print*," "
          if (hkl_equal(h,nulo)) exit

          if (.not. allocated(reflexiones)) allocate (reflexiones(grp_espacial%multip))

          num=1
          reflexiones(num)%h=h
          reflexiones(num)%phase=fase
          uno:do i=2,grp_espacial%multip
             call hkl_RP(h,fase,grp_espacial%Symop(i),k,fase1)
             do j=1,num
                if (hkl_equal(k,reflexiones(j)%h)) cycle uno
             end do
             num=num+1
             reflexiones(num)%h=k
             reflexiones(num)%phase=fase1
          end do uno

          do i=1,num
             k=reflexiones(i)%h
             fase1=reflexiones(i)%phase
             h=asu_hkl(k,grp_espacial)
             if (hkl_equal(h,nulo)) cycle
             if (hkl_equal(h,-k)) cycle
             write(*,'(a,3i4,f8.1)') " Reflection in asymmetric unit: ", h,fase1
          end do

          print*," "
          call system("pause ")
       end do

    End Subroutine Menu_Refl_2

    !!----
    !!---- Subroutine Menu_REFL_3
    !!----
    !!
    Subroutine Menu_Refl_3()
       !---- Local Variables ----!
       logical               :: info
       character(len=20)     :: line, spgr
       integer               :: i, iv, ierr, npos
       integer, dimension(1) :: ivet
       integer, dimension(3) :: h,k,nulo
       real, dimension(1)    :: vet
       type (Space_Group_type)    :: grp_espacial

       nulo=0
       info=.true.

       do
          call system('cls')
          print*," "
          print*,"     GENERAL CRYSTRALLOGRAPHY CALCULATOR "
          print*," "
          print*,"     Systematic Absences"
          print*," ==========================="
          print*," "
          print*," "

          if (info) then
             print*," Space Group (HM/Hall/Num): "

             read(*,'(a)') line
             if (len_trim(line)==0) exit
             line=adjustl(line)
             call set_spacegroup(line,grp_espacial)
             info=.false.
          end if

          print*," "
          print*," "
          print*,"  Reflection hkl: "
          read(*,*) h(1),h(2),h(3)
          print*," "
          if (hkl_equal(h,nulo)) exit
          if (hkl_absent(h,grp_espacial)) then
             print*," This reflexion IS a Systematic Absence"
          else
             print*," This reflexion IS NOT a Systematic Absence"
          end if
          print*," "
          call system("pause ")
       end do

    End Subroutine Menu_Refl_3

    !!----
    !!---- Subroutine Menu_REFL_4
    !!----
    !!
    subroutine Menu_REFL_4()
       !---- Local Variables ----!
       logical                                     :: info
       character(len=20)                           :: spgr,line
       integer, dimension(3)                       :: h,k,l,nulo
       integer, dimension(1)                       :: ivet
       integer                                     :: i,j,num,iv,ierr,npos
       real                                        :: fase,fase1,fase2
       real, dimension(1)                          :: vet
       type (Space_Group_type)                          :: grp_espacial
       type (Reflection_type),allocatable, dimension(:) :: reflexiones

       nulo=0
       info=.true.
       num=0

       do
          call system('cls')
          print*," "
          print*,"     GENERAL CRYSTRALLOGRAPHY CALCULATOR "
          print*," "
          print*,"     List of equivalent reflections "
          print*," ======================================"
          print*," "
          print*," "

          if (info) then
             print*," Space Group (HM/Hall/Num): "

             read(*,'(a)') line
             if (len_trim(line)==0) exit
             line=adjustl(line)
             call set_spacegroup(line,grp_espacial)
             info=.false.
          end if

          print*," "
          print*," Reflection hkl and phase: "
          read(*,*) h(1),h(2),h(3),fase
          print*," "
          if (hkl_equal(h,nulo)) exit
          if (.not. allocated(reflexiones)) allocate (reflexiones(grp_espacial%multip))

          num=1
          reflexiones(num)%h=h
          reflexiones(num)%phase=fase
          uno:do i=2,grp_espacial%multip
             call hkl_RP(h,fase,grp_espacial%Symop(i),k,fase1)
             do j=1,num
                if (hkl_equal(k,reflexiones(j)%h)) cycle uno
             end do
             num=num+1
             reflexiones(num)%h=k
             reflexiones(num)%phase=fase1
          end do uno

          j=0
          do i=1,num,2
             k=reflexiones(i)%h
             fase1=reflexiones(i)%phase
             if ((i+1) <= num) then
                l=reflexiones(i+1)%h
                fase2=reflexiones(i+1)%phase
                write(*,'(5x,3i4,f8.1,10x,3i4,f8.1)')k,fase1,l,fase2
             else
                write(*,'(5x,3i4,f8.1)')k,fase1
             end if
             j=j+1
             j=mod(j,12)
             if (j==0) then
                print*," "
                call system("pause ")
                print*," "
             end if
          end do

          print*," "
          call system("pause ")
       end do

    End Subroutine Menu_Refl_4

    !!----
    !!---- Subroutine Menu_REFL_5
    !!----
    !!
    Subroutine Menu_Refl_5()
       !---- Local Variables ----!
       logical                                     :: info
       character(len=20)                           :: spgr,line
       character(len=80)                           :: name_file
       integer, dimension(3)                       :: h,k,l,nulo
       integer, dimension(1)                       :: ivet
       integer                                     :: i,j,num,iv,ierr,npos
       real                                        :: fase,fase1,fase2,val1,val2
       real, dimension(1)                          :: vet
       real, dimension(3)                          :: celda, angulo
       type (Space_Group_type)                          :: grp_espacial
       type (Reflect_type),allocatable, dimension(:) :: reflexiones
       type (Crystal_Cell_type)                         :: celdilla

       nulo=0
       info=.true.
       num=0

       do
          call system('cls')
          print*," "
          print*,"     GENERAL CRYSTRALLOGRAPHY CALCULATOR "
          print*," "
          print*,"     List of Unique Reflections  "
          print*," ======================================"
          print*," "
          print*," "

          if (info) then
             print*," Space Group (HM/Hall/Num): "

             read(*,'(a)') line
             if (len_trim(line)==0) exit
             line=adjustl(line)
             call set_spacegroup(line,grp_espacial)
             info=.false.
          end if
          print*," "
          print*," INTERVAL in Sin_Theta/Lambda: "
          read(*,*) val1,val2
          if (val1*val2 < 0.0001) exit
          print*," "
          if (.not. allocated(reflexiones)) allocate (reflexiones(10000))
          call hkl_uni(celdilla,grp_espacial,.true.,val1,val2,'s',num,reflexiones)
          print*," "
          print*," Name of the output file: "
          read(*,'(a)') name_file
          open(1,file=name_file)
          write(1,'(a)') "    LIST OF UNIQUE REFLECTIONS"
          write(1,'(a)') " ================================"
          call write_spacegroup(grp_espacial,1)
          call Write_Crystal_Cell(Celdilla,1)
          write(1,'(/a,2f8.4,a/)') " => List of reflections within: ",val1,val2," 1/Ang"
          do i=1,num
             write(1,'(3i4,i5,f10.5,i8)') reflexiones(i)%h, reflexiones(i)%mult, &
                                          reflexiones(i)%S, i
          end do
          close(1)

       end do

    End Subroutine Menu_Refl_5

    !!----
    !!---- Subroutine Menu_REFL_6
    !!----
    !!
    subroutine Menu_REFL_6()
       !---- Local Variables ----!
       logical                                     :: info
       character(len=20)                           :: spgr,line
       integer, dimension(3)                       :: h,k,l,nulo
       integer, dimension(1)                       :: ivet
       integer                                     :: mul
       integer                                     :: i,j,iv,ierr,npos
       real, dimension(1)                          :: vet
       type (Space_Group_type)                     :: grp_espacial

       nulo=0
       info=.true.

       do
          call system('cls')
          print*," "
          print*,"     GENERAL CRYSTRALLOGRAPHY CALCULATOR "
          print*," "
          print*,"     Multiplicity of the Reflection"
          print*," ======================================"
          print*," "
          print*," "

          if (info) then
             print*," Space Group (HM/Hall/Num): "

             read(*,'(a)') line
             if (len_trim(line)==0) exit
             line=adjustl(line)
             call set_spacegroup(line,grp_espacial)
             info=.false.
          end if

          print*," "
          print*," Reflection hkl: "
          read(*,*) h(1),h(2),h(3)
          print*," "
          if (hkl_equal(h,nulo)) exit
          mul=hkl_mult(h,grp_espacial,.true.)
          write(*,'(6x,a,i6)') "Multiplicity of reflection: ",mul
          print*," "
          call system("pause ")
       end do

    End Subroutine Menu_Refl_6

 End Module Menu_2
