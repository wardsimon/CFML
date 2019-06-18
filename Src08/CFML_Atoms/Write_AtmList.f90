!!----
!!----
!!----
SubModule (CFML_Atoms) Atm_003
   Contains
   !!----
   !!---- WRITE_INFO_ATOMS_LIST
   !!----    Write the atoms in the asymmetric unit
   !!----
   !!---- 12/06/2019 
   !!
   Module Subroutine Write_Info_Atom_List(A, Iunit)
      !---- Arguments ----!
      class(alist_type),              intent(in) :: A        ! Atom list object
      integer,              optional, intent(in) :: IUnit    ! Logical unit

      !---- Local Variables ----!
      character(len=1)   :: car
      character(len=4)   :: car2
      character(len=80)  :: fmt1,fmt2
      character(len=132) :: line
      integer            :: n, lun

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
      
      select type(A)
         type is (alist_type)
            err_CFML%Ierr=1
            err_CFML%Msg=" Write_Info_Atom_List@CFML_ATOMS: Incompatible argument for this peocedure!"
            return
            
         type is (atm_list_type)
            car2=u_case(A%Atom(1)%Utype)
            
         type is (atm_std_list_type)
            car2=u_case(A%Atom(1)%Utype)  
            
         type is (MAtm_std_List_Type)
            car2=u_case(A%Atom(1)%Utype)  
            
         type is (Atm_ref_List_Type)
            car2=u_case(A%Atom(1)%Utype)         
      end select
      
      line=" "
      select case (trim(car2))
         case ("BETA")
            line="beta[iso/eq]    beta_11      beta_22      beta_33      beta_12      beta_13      beta_23"
         case ("U")
            line="U[iso/eq]        U_11         U_22         U_33         U_12         U_13         U_23"
         case ("B") 
            line="U[iso/eq]        B_11         B_22         B_33         B_12         B_13         B_23"
      end select
            
      line="A   Atom      Chem     Mult     x/a       y/b       z/c       Occ     "//trim(line)      
      write(unit=lun,fmt="(T3,a)") trim(line)
      
      line=repeat("=", len_trim(line))
      write(unit=lun,fmt="(T3,a)") trim(line)
      
      fmt1="(T3,a,T11,a,T20,a,T28,i3,5f10.4)"     ! Iso
      fmt2="(T3,a,T11,a,T20,a,T28,i3,11f10.4)"    ! Aniso 
      
      select type (A)
         type is (atm_list_type)
            do n=1,A%natoms
               car=" "
               if (.not. A%active(n)) car='-'
               car2=A%Atom(n)%ThType
               car2=u_case(car2)      
               select case (trim(car2))
                  case ('ISO')
                     write(unit=lun,fmt=fmt1)  car, trim(A%Atom(n)%Lab), trim(A%Atom(n)%chemSymb), &
                          A%Atom(n)%mult,  A%Atom(n)%X, A%Atom(n)%Occ, A%Atom(n)%U_iso
                  case ('ANI')    
                     write(unit=lun,fmt=fmt2)  car, trim(A%Atom(n)%Lab), trim(A%Atom(n)%chemSymb), &
                          A%Atom(n)%mult,  A%Atom(n)%X, A%Atom(n)%Occ, A%Atom(n)%U_iso, A%Atom(n)%U
               end select
            end do   
               
         type is (atm_std_list_type)
            do n=1,A%natoms
               car=" "
               if (.not. A%active(n)) car='-'
               car2=A%Atom(n)%ThType
               car2=u_case(car2)      
               select case (trim(car2))
                  case ('ISO')
                     write(unit=lun,fmt=fmt1)  car, trim(A%Atom(n)%Lab), trim(A%Atom(n)%chemSymb), &
                          A%Atom(n)%mult,  A%Atom(n)%X, A%Atom(n)%Occ, A%Atom(n)%U_iso
                  case ('ANI')    
                     write(unit=lun,fmt=fmt2)  car, trim(A%Atom(n)%Lab), trim(A%Atom(n)%chemSymb), &
                          A%Atom(n)%mult,  A%Atom(n)%X, A%Atom(n)%Occ, A%Atom(n)%U_iso, A%Atom(n)%U
               end select
            end do 
            
         type is (MAtm_std_List_Type)
            do n=1,A%natoms
               car=" "
               if (.not. A%active(n)) car='-'
               car2=A%Atom(n)%ThType
               car2=u_case(car2)      
               select case (trim(car2))
                  case ('ISO')
                     write(unit=lun,fmt=fmt1)  car, trim(A%Atom(n)%Lab), trim(A%Atom(n)%chemSymb), &
                          A%Atom(n)%mult,  A%Atom(n)%X, A%Atom(n)%Occ, A%Atom(n)%U_iso
                  case ('ANI')    
                     write(unit=lun,fmt=fmt2)  car, trim(A%Atom(n)%Lab), trim(A%Atom(n)%chemSymb), &
                          A%Atom(n)%mult,  A%Atom(n)%X, A%Atom(n)%Occ, A%Atom(n)%U_iso, A%Atom(n)%U
               end select
            end do 
            
         type is (Atm_ref_List_Type)
            do n=1,A%natoms
               car=" "
               if (.not. A%active(n)) car='-'
               car2=A%Atom(n)%ThType
               car2=u_case(car2)      
               select case (trim(car2))
                  case ('ISO')
                     write(unit=lun,fmt=fmt1)  car, trim(A%Atom(n)%Lab), trim(A%Atom(n)%chemSymb), &
                          A%Atom(n)%mult,  A%Atom(n)%X, A%Atom(n)%Occ, A%Atom(n)%U_iso
                  case ('ANI')    
                     write(unit=lun,fmt=fmt2)  car, trim(A%Atom(n)%Lab), trim(A%Atom(n)%chemSymb), &
                          A%Atom(n)%mult,  A%Atom(n)%X, A%Atom(n)%Occ, A%Atom(n)%U_iso, A%Atom(n)%U
               end select
            end do 
      end select          
               
   End Subroutine Write_Info_Atom_List
   
End SubModule Atm_003   