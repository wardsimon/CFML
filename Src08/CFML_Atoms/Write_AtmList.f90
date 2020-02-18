!!----
!!----
!!----
SubModule (CFML_Atoms) Write_Atoms
   Contains
   !!----
   !!---- WRITE_INFO_ATOMS_LIST
   !!----    Write the atoms in the asymmetric unit
   !!----
   !!---- 12/06/2019
   !!
   Module Subroutine Write_Info_Atom_List(A, Iunit)
      !---- Arguments ----!
      type(atlist_type), intent(in) :: A        ! Atom list object
      integer, optional, intent(in) :: IUnit    ! Logical unit

      !---- Local Variables ----!
      character(len=1)             :: car
      character(len=4)             :: car2
      character(len=:),allocatable :: fmt1,fmt2
      character(len=:),allocatable :: line
      integer                      :: n, lun

      !> Local Types
      type (atm_type)      :: atm
      type (atm_std_type)  :: atms
      type (matm_std_type) :: matm
      type (atm_ref_type)  :: atr

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

      associate (Aat => A%Atom)
         car2=aat(1)%UType
      end associate

      line=" Atom    Chem     Mult      x/a       y/b       z/c      B[iso/eq]    Occ  "
      select case (trim(car2))
         case ("BETA")
            line=" Atom    Chem     Mult      x/a       y/b       z/c     beta[iso/eq]    Occ  "
            line=line//"   beta_11      beta_22      beta_33      beta_12      beta_13      beta_23"
         case ("U")
            line=" Atom    Chem     Mult      x/a       y/b       z/c      U[iso/eq]    Occ  "
            line=line//"      U_11         U_22         U_33         U_12         U_13         U_23"
         case ("B")
            line=" Atom    Chem     Mult      x/a       y/b       z/c      B[iso/eq]    Occ  "
            line=line//"      B_11         B_22         B_33         B_12         B_13         B_23"
         case default
            line=" Atom    Chem     Mult      x/a       y/b       z/c      B[iso/eq]    Occ  "
      end select
           !1234567890123456781234567
      write(unit=lun,fmt="(T3,a)") trim(line)

      line=repeat("=", len_trim(line))
      write(unit=lun,fmt="(T3,a)") trim(line)

      fmt1="(T3,a,T6,a,T18,a,T28,i3,5f10.5)"     ! Iso
      fmt2="(T3,a,T6,a,T18,a,T28,i3,11f10.5)"    ! Aniso

      select type (aat => A%atom)
         type is (atm_type)
            do n=1,A%natoms
               car=" "
               if (.not. A%active(n)) car='-'
               atm=aat(n)
               car2=Atm%ThType
               car2=u_case(car2)
               select case (trim(car2))
                  case ('ISO')
                     write(unit=lun,fmt=fmt1)  car, trim(Atm%Lab), trim(Atm%chemSymb), &
                          Atm%mult,  Atm%X, Atm%U_iso, Atm%Occ
                  case ('ANI')
                     write(unit=lun,fmt=fmt2)  car, trim(Atm%Lab), trim(Atm%chemSymb), &
                          Atm%mult,  Atm%X, Atm%U_iso, Atm%Occ, Atm%U
               end select
            end do

         type is (atm_std_type)
            do n=1,A%natoms
               car=" "
               if (.not. A%active(n)) car='-'

               atms=aat(n)
               car2=Atms%ThType
               car2=u_case(car2)
               select case (trim(car2))
                  case ('ISO')
                     write(unit=lun,fmt=fmt1)  car, trim(Atms%Lab), trim(Atms%chemSymb), &
                          Atms%mult,  Atms%X, Atms%U_iso, Atms%Occ
                  case ('ANI')
                     write(unit=lun,fmt=fmt2)  car, trim(Atms%Lab), trim(Atms%chemSymb), &
                          Atms%mult,  Atms%X, Atms%U_iso, Atms%Occ, Atms%U
               end select
            end do

         type is (MAtm_std_Type)
            do n=1,A%natoms
               car=" "
               if (.not. A%active(n)) car='-'

               matm=aat(n)
               car2=matm%ThType
               car2=u_case(car2)
               select case (trim(car2))
                  case ('ISO')
                     write(unit=lun,fmt=fmt1)  car, trim(matm%Lab), trim(matm%chemSymb), &
                          matm%mult,  matm%X, matm%U_iso, matm%Occ
                  case ('ANI')
                     write(unit=lun,fmt=fmt2)  car, trim(matm%Lab), trim(matm%chemSymb), &
                          matm%mult,  matm%X, matm%U_iso, matm%Occ, matm%U
               end select
            end do

         type is (Atm_ref_Type)
            do n=1,A%natoms
               car=" "
               if (.not. A%active(n)) car='-'

               atr=aat(n)
               car2=atr%ThType
               car2=u_case(car2)
               select case (trim(car2))
                  case ('ISO')
                     write(unit=lun,fmt=fmt1)  car, trim(atr%Lab), trim(atr%chemSymb), &
                          atr%mult,  atr%X, atr%U_iso, atr%Occ
                  case ('ANI')
                     write(unit=lun,fmt=fmt2)  car, trim(atr%Lab), trim(atr%chemSymb), &
                          atr%mult,  atr%X, atr%U_iso, atr%Occ,  atr%U
               end select
            end do
      end select

   End Subroutine Write_Info_Atom_List

End SubModule Write_Atoms