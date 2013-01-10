!!----
!!---- Menu: 3
!!---- Atoms Calculations
!!----
!!
 Module Menu_3
    !---- Use File ----!
    use Menu_0
    use CFML_Crystallographic_Symmetry
    use CFML_Atom_TypeDef
    use CFML_IO_Messages,      only: Wait_Message
    use CFML_string_utilities, only: getnum
    use CFML_Math_General,     only: cosd, acosd, sind, asind
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
          call system(clear_string)

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
          write(unit=*,fmt="(a)") " [4] Calculate Ionic Dipolar moment & polarisation"
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
          call system(clear_string)
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
            call Wait_Message(" => Press <enter> to continue ...")
            cycle
          end if
          do
             call system(clear_string)
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
                call Wait_Message(" => Press <enter> to continue ...")
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
          call system(clear_string)
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
                call Wait_Message(" => Press <enter> to continue ...")
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
                  call Wait_Message(" => Press <enter> to continue ...")
                  cycle
               end if
            else
               write(unit=*,fmt="(a)") " => Illegal File name: "//trim(line)
               call Wait_Message(" => Press <enter> to continue ...")
               cycle
            end if
          end if
          if (err_form) then
             write(unit=*,fmt="(a)") trim(err_form_mess)
          else
             write(unit=*,fmt="(a)") " => File: "//trim(line)//" successfully read!"
             structure_read=.true.
             call Wait_Message(" => Press <enter> to continue ...")
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

       call system(clear_string)
       write(unit=*,fmt="(a)") "     GENERAL CRYSTALLOGRAPHY CALCULATOR "
       write(unit=*,fmt="(a)") " "
       write(unit=*,fmt="(a)") "     Atomistic Calculations "
       write(unit=*,fmt="(a)") " ================================"
       write(unit=*,fmt="(a)") " Showing Cell, Space Group and Atom positions ..."
       write(unit=*,fmt="(a)") " "
       if(.not. structure_read) then
         write(unit=*,fmt="(/a/)") " => A CIF or CFL file must read before using this option! "
         call Wait_Message(" => Press <enter> to continue ...")
         return
       end if
       write(unit=*,fmt="(/a/)") " => UNIT CELL information: "

       call Write_Crystal_Cell(Cell)
       call Wait_Message(" => Press <enter> to continue ...")

       call system(clear_string)
       write(unit=*,fmt="(/a/)") " => SPACE GROUP information: "
       call Write_SpaceGroup(SpG,Full=.true.)
       call Wait_Message(" => Press <enter> to continue ...")

       call system(clear_string)
       write(unit=*,fmt="(/a/)") " => ATOMS information: "
       call Write_Atom_List(A,level=1)
       call Wait_Message(" => Press <enter> to continue ...")

    End Subroutine Menu_Atom_3

    Subroutine Menu_Atom_4()
       !---- Local Variables ----!
       character(len=20)     :: line
       integer               :: i,j,Mult
       real                  :: q, qp,qm,pol,ang
       real, dimension(3)    :: pos,cpos,r_plus, r_minus, r_pol,u_pol, rt,vc
       real, dimension(3,192):: orb
       logical               :: calc_possible=.true.

       call system(clear_string)
       write(unit=*,fmt="(a)") "     GENERAL CRYSTALLOGRAPHY CALCULATOR "
       write(unit=*,fmt="(a)") " "
       write(unit=*,fmt="(a)") "     Atomistic Calculations "
       write(unit=*,fmt="(a)") " ================================"
       write(unit=*,fmt="(a)") " Calculating the ionic polarisation ..."
       write(unit=*,fmt="(a)") " "
       do i=1,A%natoms
         if(abs(A%atom(i)%charge) <= 0.001)  then
           calc_possible=.false.
           exit
         end if
       end do
       if( .not. calc_possible) then
         write(unit=*,fmt="(a)") " => Calculation of P impossible. No charges have been provided! "
         call Wait_Message(" => Press <enter> to continue ...")
         return
       end if
       !if(SpG%Centred /= 1) then
       !  write(unit=*,fmt="(a)") " => Polarisation = 0.0 Centrosymmetric Crystal! "
       !  call Wait_Message(" => Press <enter> to continue ...")
       !  return
       !end if
       r_plus=0.0; r_minus=0.0
       qp=0.0; qm=0.0
       vc=[0.5,0.5,0.5] !Centre of the unit cell with respet to which we take the vector positions
       if(SpG%Centred == 1) vc=0.0
       write(unit=*,fmt="(a,3f7.2,a)") " => Origin of vector positions taken at: (",vc," )"

       do i=1,A%natoms
          pos=A%atom(i)%x
          q=A%atom(i)%charge
          call Get_Orbit(pos,SpG,Mult,orb) !here the orbit has always positive coordinates
          rt=0.0
          do j=1,Mult
            rt=rt + orb(:,j)-vc
          end do
          if (.not. Lattice_trans(rt,Spg%spg_lat)) rt=0.0 !If non-centrosymmetric no problem
          orb(:,mult)= orb(:,mult)-rt  !Apply a lattice translation to have vector sum = 0
          write(unit=*,fmt="(a,f8.2,a,3f10.5)") " => Atom: "//trim(A%Atom(i)%Lab)//",  Charge: ", A%Atom(i)%Charge
          do j=1,Mult
            pos=orb(:,j)-vc  !calculate positions w.r.t. centre of the cell
            cpos=Cart_Vector("D",pos,Cell)
            if( q > 0.0) then
               qp=qp+q
               r_plus=r_plus+cpos*q
            else
               qm=qm+abs(q)
               r_minus=r_minus+cpos*abs(q)
            end if
            write(unit=*,fmt="(a,i3,a,3f9.5,3f11.4)") "   ",j,": ",pos,cpos
          end do
       end do
       if(abs(qp-qm) > eps) then
         write(unit=*,fmt="(a,2f10.2)") " => Warning! The crystal is NOT NEUTRAL! -> Q+, Q-: ",qp,qm
       end if
       r_plus=r_plus/qp     !Cartesian vectors giving the c.o.g of charges
       r_minus=r_minus/qm
       r_pol=qp*(r_plus-r_minus)  !Dipole moment vector
       pol=sqrt(dot_product(r_pol,r_pol))

       write(unit=*,fmt="(a,3f12.5,a)")  " => Q+ charge centre : (",r_plus," )"
       write(unit=*,fmt="(a,3f12.5,a)")  " => Q- charge centre : (",r_minus," )"

       write(unit=*,fmt="(a,f12.5)")     " => Ionic Dipolar Moment (electron.Angstrom): ",pol
       write(unit=*,fmt="(a,g12.5)")     " => Ionic Dipolar Moment (Coulomb.Metre): ",1.60217646e-29*pol
       write(unit=*,fmt="(a,3f12.5,a)")  " => Cartesian Dipolar Moment vector :(",r_pol," )"
       write(unit=*,fmt="(a,f12.5)")     " => Ionic Polarisation (Coulomb/m^2): ",16.0217646*pol/Cell%CellVol  !/real(SpG%Numlat)
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
       call Wait_Message(" => Press <enter> to continue ...")

    End Subroutine Menu_Atom_4


    Function Angle_vect(u,v) Result(angle)
      real, dimension(:), intent(in) :: u,v
      real :: angle
      real :: mu, mv
      mu=sqrt(dot_Product(u,u))
      mv=sqrt(dot_Product(v,v))
      if(mu < eps .or. mv < eps) then
          write(unit=*,fmt="(a)") " -> One of the directions is [0 0 0] ... retry!"
          call Wait_Message(" => Press <enter> to continue ...")
          return
      end if
      angle=dot_Product(u,v)/mu/mv
      if(angle > 1.0) angle=1.0
      if(angle < -1.0) angle=-1.0
      angle=acosd(angle)
      return
    End Function Angle_vect

end module Menu_3
