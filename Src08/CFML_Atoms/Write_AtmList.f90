!!----
!!----
!!----
SubModule (CFML_IOForm)  Write_Atoms
   Contains
   !!----
   !!---- WRITE_INFO_ATOMS_LIST
   !!----    Write the atoms in the asymmetric unit
   !!----
   !!---- 12/06/2019
   !!
   Module Subroutine Write_Atom_List(A, Iunit)
      !---- Arguments ----!
      type(atlist_type), intent(in) :: A        ! Atom list object
      integer, optional, intent(in) :: IUnit    ! Logical unit

      !---- Local Variables ----!
      character(len=1)             :: car
      character(len=:),allocatable :: car2
      character(len=:),allocatable :: fmt1,fmt2
      character(len=:),allocatable :: line
      integer                      :: n, lun, k, j

      !> Init
      lun=6
      if (present(iunit)) lun=iunit

      !> Header
      write(unit=lun, fmt="(/,a)")    "  ===="
      write(unit=lun, fmt="(a)")      "  ====  Atoms information"
      write(unit=lun, fmt="(a,/)")    "  ===="

      !> Check number of atoms
      if (a%natoms == 0) then
         write(unit=lun,fmt="(/,a,/)") "  => No atoms provided!"
         return
      end if

      car2="ISO"
      if(any(A%Atom(:)%Thtype == "ANI")) car2="ANI"

      if(car2 == "ISO") then
          line="    Atom  Scatt/Chem    Mult      x/a       y/b       z/c      B[iso]     Occ"
      else
          car2=trim(u_case(A%Atom(1)%Utype))
          select case (car2)
             case ("BETA")
                line="    Atom  Scatt/Chem    Mult      x/a       y/b       z/c     B[iso]      Occ  "
                line=line//"  beta_11   beta_22   beta_33   beta_12   beta_13   beta_23"
             case ("U")
                line="    Atom  Scatt/Chem    Mult      x/a       y/b       z/c    U[iso/eq]   Occ  "
                line=line//"     U_11      U_22      U_33      U_12      U_13      U_23"
             case ("B")
                line="    Atom  Scatt/Chem    Mult      x/a       y/b       z/c     B[iso]      Occ"
                line=line//"      B_11      B_22      B_33      B_12      B_13      B_23"
             case default
                line="    Atom  Scatt/Chem    Mult      x/a       y/b       z/c     B[iso]      Occ"
          end select
      end if
           !1234567890123456781234567
      write(unit=lun,fmt="(T3,a)") trim(line)

      line=repeat("=", len_trim(line))
      write(unit=lun,fmt="(T3,a)") trim(line)

      fmt1="(T3,a,T6,a,T18,a,T28,i3,5f10.5)"     ! Iso
      fmt2="(T3,a,T6,a,T18,a,T28,i3,11f10.5)"    ! Aniso

      select type (Atm => A%atom)
         type is (atm_type)
            do n=1,A%natoms
               car2=Atm(n)%ThType
               car2=u_case(car2)
               car=car2(1:1)
               if (.not. A%active(n)) car='-'
               select case (trim(car2))
                  case ('ISO')
                     write(unit=lun,fmt=fmt1)  car, trim(Atm(n)%Lab), trim(Atm(n)%SfacSymb), &
                          Atm(n)%mult,  Atm(n)%X, Atm(n)%U_iso, Atm(n)%Occ
                  case ('ANI')
                     write(unit=lun,fmt=fmt2)  car, trim(Atm(n)%Lab), trim(Atm(n)%SfacSymb), &
                          Atm(n)%mult,  Atm(n)%X, Atm(n)%U_iso, Atm(n)%Occ, Atm(n)%U
               end select
               if(Atm(n)%magnetic)  then
                    Select Case (trim(u_case(a%mcomp)))
                      Case ("CRYSTAL")
                         write(unit=lun,fmt="(T20,a,3f10.5,a)")  "Moment(uB):", Atm(n)%moment, "  Crystal Components:   Mx(a)  My(b)  Mz(b)"
                      Case ("SPHERICAL")
                         write(unit=lun,fmt="(T20,a,3f10.5,a)")  "Moment(uB):", Atm(n)%moment, "  Spherical Components:  Moment  Phi  Theta"
                      Case ("CARTESIAN")
                         write(unit=lun,fmt="(T20,a,3f10.5,a)")  "Moment(uB):", Atm(n)%moment, "  Cartesian Components:  MxC  MyC  MzC"
                    End Select
               end if
            end do

         type is (atm_std_type)
            do n=1,A%natoms
               car=" "
               if (.not. A%active(n)) car='-'
               car2=Atm(n)%ThType
               car2=u_case(car2)
               select case (trim(car2))
                  case ('ISO')
                     write(unit=lun,fmt=fmt1)  car, trim(Atm(n)%Lab), trim(Atm(n)%chemSymb), &
                          Atm(n)%mult,  Atm(n)%X, Atm(n)%U_iso, Atm(n)%Occ
                  case ('ANI')
                     write(unit=lun,fmt=fmt2)  car, trim(Atm(n)%Lab), trim(Atm(n)%chemSymb), &
                          Atm(n)%mult,  Atm(n)%X, Atm(n)%U_iso, Atm(n)%Occ, Atm(n)%U
               end select
               if(Atm(n)%magnetic)  then
                    Select Case (trim(u_case(a%mcomp)))
                      Case ("CRYSTAL")
                         write(unit=lun,fmt="(T20,a,3f10.5,a)")  "Moment(uB):", Atm(n)%moment, "  Crystal Components:   Mx(a)  My(b)  Mz(b)"
                      Case ("SPHERICAL")
                         write(unit=lun,fmt="(T20,a,3f10.5,a)")  "Moment(uB):", Atm(n)%moment, "  Spherical Components:  Moment  Phi  Theta"
                      Case ("CARTESIAN")
                         write(unit=lun,fmt="(T20,a,3f10.5,a)")  "Moment(uB):", Atm(n)%moment, "  Cartesian Components:  MxC  MyC  MzC"
                    End Select
               end if
            end do

         class is (MAtm_std_Type)
            do n=1,A%natoms
               car=" "
               if (.not. A%active(n)) car='-'
               car2=atm(n)%ThType
               car2=u_case(car2)
               select case (trim(car2))
                  case ('ISO')
                     write(unit=lun,fmt=fmt1)  car, trim(Atm(n)%Lab), trim(Atm(n)%chemSymb), &
                          Atm(n)%mult,  Atm(n)%X, Atm(n)%U_iso, Atm(n)%Occ
                  case ('ANI')
                     write(unit=lun,fmt=fmt2)  car, trim(Atm(n)%Lab), trim(Atm(n)%chemSymb), &
                          Atm(n)%mult,  Atm(n)%X, Atm(n)%U_iso, Atm(n)%Occ, Atm(n)%U
               end select
               if(Atm(n)%magnetic)  then
                    Select Case (trim(u_case(a%mcomp)))
                      Case ("CRYSTAL")
                         write(unit=lun,fmt="(T20,a,3f10.5,a)")  "Moment(uB):", Atm(n)%moment, "  Crystal Components:   Mx(a)  My(b)  Mz(b)"
                      Case ("SPHERICAL")
                         write(unit=lun,fmt="(T20,a,3f10.5,a)")  "Moment(uB):", Atm(n)%moment, "  Spherical Components:  Moment  Phi  Theta"
                      Case ("CARTESIAN")
                         write(unit=lun,fmt="(T20,a,3f10.5,a)")  "Moment(uB):", Atm(n)%moment, "  Cartesian Components:  MxC  MyC  MzC"
                    End Select
               end if
            end do

            Select Case (trim(u_case(a%mcomp)))
              Case ("CRYSTAL")
                 line="    Atom        Modulation  Harmonic             Cos_x     Cos_y     Cos_z     Sin_x     Sin_y     Sin_z"
              Case ("SPHERICAL")
                 line="    Atom        Modulation  Harmonic           Cos_Mod   Cos_Phi Cos_Theta   Sin_Mod   Sin_Phi Sin_Theta"
              Case ("CARTESIAN")
                 line="    Atom        Modulation  Harmonic            Cos_xC    Cos_yC    Cos_zC    Sin_xC    Sin_yC    Sin_zC"
            End Select
            write(unit=lun,fmt="(/,/,T3,a)") line
            line=repeat("=", len_trim(line))
            write(unit=lun,fmt="(T3,a)") line
            do n=1,A%natoms
               do j=1,atm(n)%n_mc
                 k=atm(n)%pmc_q(j)
                 write(unit=lun,fmt="(T7,a,t15,a,i3,a,t47,6f10.5)")   trim(Atm(n)%Lab),           "        Moment    [",k,"]", Atm(n)%Mcs(:,j)
               end do
               do j=1,atm(n)%n_dc
                 k=atm(n)%pdc_q(j)
                 write(unit=lun,fmt="(T7,a,t15,a,i3,a,t47,6f10.5)")   trim(Atm(n)%Lab),           "  Displacement    [",k,"]", Atm(n)%Dcs(:,j)
               end do
               do j=1,atm(n)%n_oc
                 k=atm(n)%poc_q(j)
                 write(unit=lun,fmt="(T7,a,t15,a,i3,a,t47,f10.5,tr20,f10.5)")   trim(Atm(n)%Lab), "     Occupancy    [",k,"]", Atm(n)%Ocs(:,j)
               end do
            end do

            if(any(atm(:)%n_uc > 0)) then
               line="    Atom        Modulation       Harmonic   U_11_Cos  U_22_Cos  U_33_Cos  U_12_Cos  U_13_Cos  U_23_Cos  U_11_Sin  U_22_Sin  U_33_Sin  U_12_Sin  U_13_Sin  U_23_Sin"
               write(unit=lun,fmt="(/,/,T3,a)") line
               line=repeat("=", len_trim(line))
               write(unit=lun,fmt="(T3,a)") line
               do n=1,A%natoms
                  do j=1,atm(n)%n_uc
                    k=atm(n)%puc_q(j)
                    write(unit=lun,fmt="(T7,a,t15,a,i3,a,t45,12f10.5)")   trim(Atm(n)%Lab), "Thermal Displacement [",k,"]", Atm(n)%Ucs(:,j)
                  end do
               end do
            end if

         type is (Atm_ref_Type)
            do n=1,A%natoms
               car=" "
               if (.not. A%active(n)) car='-'
               car2=atm(n)%ThType
               car2=u_case(car2)
               select case (trim(car2))
                  case ('ISO')
                     write(unit=lun,fmt=fmt1)  car, trim(Atm(n)%Lab), trim(Atm(n)%chemSymb), &
                          Atm(n)%mult,  Atm(n)%X, Atm(n)%U_iso, Atm(n)%Occ
                  case ('ANI')
                     write(unit=lun,fmt=fmt2)  car, trim(Atm(n)%Lab), trim(Atm(n)%chemSymb), &
                          Atm(n)%mult,  Atm(n)%X, Atm(n)%U_iso, Atm(n)%Occ, Atm(n)%U
               end select
               if(Atm(n)%magnetic)  then
                    Select Case (trim(u_case(a%mcomp)))
                      Case ("CRYSTAL")
                         write(unit=lun,fmt="(T20,a,3f10.5,a)")  "Moment(uB):", Atm(n)%moment, "  Crystal Components:   Mx(a)  My(b)  Mz(b)"
                      Case ("SPHERICAL")
                         write(unit=lun,fmt="(T20,a,3f10.5,a)")  "Moment(uB):", Atm(n)%moment, "  Spherical Components:  Moment  Phi  Theta"
                      Case ("CARTESIAN")
                         write(unit=lun,fmt="(T20,a,3f10.5,a)")  "Moment(uB):", Atm(n)%moment, "  Cartesian Components:  MxC  MyC  MzC"
                    End Select
               end if
            end do
      end select

   End Subroutine Write_Atom_List

End SubModule Write_Atoms