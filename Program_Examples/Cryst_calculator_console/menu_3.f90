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
    use CFML_Crystal_Metrics,  only: Crystal_Cell_Type, Write_Crystal_Cell, Cart_Vector
    use CFML_IO_Formats,       only: Readn_set_Xtal_Structure,err_form_mess,err_form,file_list_type
    use Menu_1,                only: Get_Wyckoff
    !---- Variables ----!
    implicit none

    logical :: structure_read=.false.
    type (space_group_type)     :: SpG
    type (Atom_list_Type)       :: A
    type (Crystal_Cell_Type)    :: Cell
    real, parameter :: eps=0.00001

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
          write(unit=*,fmt="(a)") "     Atomistic Calculations "
          write(unit=*,fmt="(a)") " =============================="
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " [0] Back..."
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " [1] Multiplicity of an atom position"
          write(unit=*,fmt="(a)") " [2] Read a CFL or CIF file to load a crystal structure"
          write(unit=*,fmt="(a)") " [3] Show current structure information"
          write(unit=*,fmt="(a)") " [4] Calculate polarisation in e.A per unit cell"
          write(unit=*,fmt="(a)",advance="no") " OPTION: "
          read(*,'(a)') car
          if (len_trim(car) == 0) exit
          car=adjustl(car)

          select case (car)
             case ('0 ')
                exit

             case ('1 ')
                call Menu_ATOM_1()   !Multiplicity

             case ('2 ')
                call Menu_ATOM_2()   !Reading a CFL or CIF file

             case ('3 ')
                call Menu_ATOM_3()   !Showing the crystal structure

             case ('4 ')
                call Menu_ATOM_4()   !Calculating polarisation if charges are given

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
    !!----
    !!---- Subroutine Menu_ATOM_2  Reading a CFL or CIF file
    !!----
    !!
    Subroutine Menu_Atom_2()
       !---- Local Variables ----!
       character(len=20)     :: line
       integer               :: i
       logical               :: esta


       do
          call system('cls')
          write(unit=*,fmt="(a)") "     GENERAL CRYSTALLOGRAPHY CALCULATOR "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "     Atomistic Calculations "
          write(unit=*,fmt="(a)") " ================================"
          write(unit=*,fmt="(a)") " Reading a CFL or CIF file ..."
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)",advance="no") " => Enter the full name of the file: "

          read(unit=*,fmt='(a)') line
          if (len_trim(line)==0) exit

          line=adjustl(line)
          i=index(line,".cif")
          if(i /= 0) then
             inquire(file=trim(line),exist=esta)
             if (esta) then
                call Readn_set_Xtal_Structure(trim(line),Cell,SpG,A,Mode="CIF")
             else
                write(unit=*,fmt="(a)") " => File: "//trim(line)//"  does not exist!"
                call system('pause')
                cycle
             end if
          else
            i=index(line,".cfl")
            if(i /= 0) then
               inquire(file=trim(line),exist=esta)
               if (esta) then
                  call Readn_set_Xtal_Structure(trim(line),Cell,SpG,A,Mode="CFL")
               else
                  write(unit=*,fmt="(a)") " => File: "//trim(line)//"  does not exist!"
                  call system('pause')
                  cycle
               end if
            else
               write(unit=*,fmt="(a)") " => Illegal File name: "//trim(line)
               call system('pause')
               cycle
            end if
          end if
          if (err_form) then
             write(unit=*,fmt="(a)") trim(err_form_mess)
          else
             write(unit=*,fmt="(a)") " => File: "//trim(line)//" successfully read!"
             structure_read=.true.
             call system('pause')
             exit
          end if

       end do

    End Subroutine Menu_Atom_2

    !!----
    !!---- Subroutine Menu_ATOM_3  Showing Cell, Space Group and Atom positions
    !!----
    !!
    Subroutine Menu_Atom_3()
       !---- Local Variables ----!
       character(len=20)     :: line
       integer               :: i

       call system('cls')
       write(unit=*,fmt="(a)") "     GENERAL CRYSTALLOGRAPHY CALCULATOR "
       write(unit=*,fmt="(a)") " "
       write(unit=*,fmt="(a)") "     Atomistic Calculations "
       write(unit=*,fmt="(a)") " ================================"
       write(unit=*,fmt="(a)") " Showing Cell, Space Group and Atom positions ..."
       write(unit=*,fmt="(a)") " "
       if(.not. structure_read) then
         write(unit=*,fmt="(/a/)") " => A CIF or CFL file must read before using this option! "
         call system('pause')
         return
       end if
       write(unit=*,fmt="(/a/)") " => UNIT CELL information: "

       call Write_Crystal_Cell(Cell)
       call system('pause')

       call system('cls')
       write(unit=*,fmt="(/a/)") " => SPACE GROUP information: "
       call Write_SpaceGroup(SpG,Full=.true.)
       call system('pause')

       call system('cls')
       write(unit=*,fmt="(/a/)") " => ATOMS information: "
       call Write_Atom_List(A,level=1)
       call system('pause')

    End Subroutine Menu_Atom_3

    Subroutine Menu_Atom_4()
       !---- Local Variables ----!
       character(len=20)     :: line
       integer               :: i,j,Mult
       real                  :: q, qp,qm,pol,ang
       real, dimension(3)    :: pos,cpos,r_plus, r_minus, r_pol,u_pol, fr_pol
       real, dimension(3,192):: orb
       logical               :: calc_possible=.true.

       call system('cls')
       write(unit=*,fmt="(a)") "     GENERAL CRYSTALLOGRAPHY CALCULATOR "
       write(unit=*,fmt="(a)") " "
       write(unit=*,fmt="(a)") "     Atomistic Calculations "
       write(unit=*,fmt="(a)") " ================================"
       write(unit=*,fmt="(a)") " Calculating the polarisation per unit cell ..."
       write(unit=*,fmt="(a)") " "
       do i=1,A%natoms
         if(abs(A%atom(i)%charge) <= 0.001)  then
           calc_possible=.false.
           exit
         end if
       end do
       if( .not. calc_possible) then
         write(unit=*,fmt="(a)") " => Calculation of P impossible. No charges have been provided! "
         call system('pause')
         return
       end if
       r_plus=0.0; r_minus=0.0
       qp=0.0; qm=0.0
       do i=1,A%natoms
          pos=A%atom(i)%x
          q=A%atom(i)%charge
          call Orbit(pos,SpG,Mult,orb)
          do j=1,Mult
            pos=orb(:,j)
            cpos=Cart_Vector("D",pos,Cell)
            if( q > 0.0) then
               qp=qp+q
               r_plus=r_plus+cpos*q
            else
               qm=qm+abs(q)
               r_minus=r_minus+cpos*abs(q)
            end if
          end do
       end do
       r_plus=r_plus/qp     !Cartesian vectors giving the c.o.g of charges
       r_minus=r_minus/qm
       r_pol=r_plus-r_minus  !Polarisation vector
       pol=sqrt(dot_product(r_pol,r_pol))
       !fr_pol = matmul(cell%Orth_Cr_cel,r_pol)

       write(unit=*,fmt="(a,f14.5)")     " => Polarisation (electron.Angstrom): ",pol
       write(unit=*,fmt="(a,3f14.5,a)")  " => Cartesian Polarisation vector: (",r_pol," )"
       if(pol > eps) then
         cpos=Cart_Vector("D",[1.0,0.0,0.0],Cell)
         ang=Angle_Vect(cpos, r_pol)
         write(unit=*,fmt="(a,f14.5,a)")  " => Angle of Polarisation vector with a-axis: ",ang," degrees"
         cpos=Cart_Vector("D",[0.0,1.0,0.0],Cell)
         ang=Angle_Vect(cpos, r_pol)
         write(unit=*,fmt="(a,f14.5,a)")  " => Angle of Polarisation vector with b-axis: ",ang," degrees"
         cpos=Cart_Vector("D",[0.0,0.0,1.0],Cell)
         ang=Angle_Vect(cpos, r_pol)
         write(unit=*,fmt="(a,f14.5,a)")  " => Angle of Polarisation vector with c-axis: ",ang," degrees"
       end if
       call system('pause')

    End Subroutine Menu_Atom_4

    Subroutine Orbit(x,Spg,Mult,orb)
       !---- Arguments ----!
       real, dimension(3),    intent (in) :: x
       type(Space_Group_type),intent (in) :: spg
       integer,               intent(out) :: mult
       real,dimension(:,:),   intent(out) :: orb

       !---- Local variables ----!
       integer            :: j, nt
       real, dimension(3) :: xx,v
       character(len=1)   :: laty

       laty="P"
       mult=1
       orb(:,1)=x(:)
       ext: do j=2,Spg%multip
          xx=ApplySO(Spg%SymOp(j),x)
          do nt=1,mult
             v=orb(:,nt)-xx(:)
             if (Lattice_trans(v,Spg%spg_lat)) then
               if (.not. Lattice_trans(v,laty)) cycle  !Count in orbit the centred related atoms
               cycle ext
             end if
          end do
          mult=mult+1
          orb(:,mult)=xx(:)
       end do ext

       return
    End Subroutine Orbit



    Function Angle_vect(u,v) Result(angle)
      real, dimension(:), intent(in) :: u,v
      real :: angle
      real :: mu, mv
      mu=sqrt(dot_Product(u,u))
      mv=sqrt(dot_Product(v,v))
      if(mu < eps .or. mv < eps) then
          write(unit=*,fmt="(a)") " -> One of the directions is [0 0 0] ... retry!"
          call system("pause ")
          return
      end if
      angle=dot_Product(u,v)/mu/mv
      if(angle > 1.0) angle=1.0
      if(angle < -1.0) angle=-1.0
      angle=acosd(angle)

    End Function Angle_vect
end module Menu_3
