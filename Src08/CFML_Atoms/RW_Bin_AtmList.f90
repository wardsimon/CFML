!!----
!!----
!!----
SubModule (CFML_Atoms) Atm_004
   Contains
   !!----
   !!---- READ_BIN_ATOM_LIST
   !!----
   !!----    Reads the atoms in the asymmetric unit in a binary file.
   !!---- 
   !!----    The procedure reads in the given order a series of bytes corresponding to the
   !!----    components of the type Ats. The full structure is re-allocated inside the procedure
   !!----    before reading the components.
   !!----
   !!----    The number of atoms is the first element read in the file.
   !!----
   !!---- 12/06/2019 
   !!
   Module Subroutine Read_Bin_Atom_List(filename, A)
      !---- Arguments ----!
      character(len=*),  intent(in)  :: filename
      class(alist_type), intent(out) :: A
      
      !---- Local Variables ----!
      integer                        :: i,n,ierr,lun

      !> Init
      call clear_error()
      
      !> Exits file?
      open(newunit=lun,file=trim(filename), access="stream", status="old", iostat=ierr)
      if (ierr /=0) then
         err_CFML%IErr=1
         err_CFML%Msg="Read_Bin_Atoms_List@CFML_ATOMS: Error opening the binary file "//trim(filename)
         return
      end if   

      !> First: read number of atoms
      n=0       
      call Allocate_Atom_List(N, A)
      
      read(unit=lun,iostat=ierr) n
      if (ierr /= 0) then
         err_CFML%IErr=1
         err_CFML%Msg="Read_Bin_Atoms_List@CFML_ATOMS: Error reading number of atoms in the binary file "//trim(filename)
         close(unit=lun)
         return
      end if
      if (n <= 0) return
      
      !> Allocating atom_list
      call Allocate_Atom_List(N, A)
      
      !> Read active
      read(unit=lun,iostat=ierr)  A%Active 
      if (ierr /= 0) then
         err_CFML%IErr=1
         err_CFML%Msg="Read_Bin_Atoms_List@CFML_ATOMS: Error reading active atoms!"
         close(unit=lun)
         return
      end if

      !> Load information 
      select type (A)
         type is (Atm_List_Type)
            do i=1,n
               read(unit=lun,iostat=ierr)  A%Atom(i)%Lab,        &
                                           A%Atom(i)%ChemSymb,   &     
                                           A%Atom(i)%Z,          &
                                           A%Atom(i)%Mult,       &
                                           A%Atom(i)%X,          &
                                           A%Atom(i)%Occ,        &
                                           A%Atom(i)%UType,      &
                                           A%Atom(i)%ThType,     &
                                           A%Atom(i)%U_iso,      &
                                           A%Atom(i)%U,          &
                                           A%Atom(i)%Magnetic,   &
                                           A%Atom(i)%Moment,     &
                                           A%Atom(i)%SfacSymb,   &
                                           A%Atom(i)%Ind_ff

                if (ierr /=0) then
                   err_CFML%IErr=1
                   err_CFML%Msg="Read_Bin_Atoms_List@CFML_ATOMS: Error reading atoms information!"
                   exit
                end if                            
            end do
               
         type is (Atm_std_List_Type)   
            do i=1,n
               read(unit=lun,iostat=ierr)  A%Atom(i)%Lab,        &
                                           A%Atom(i)%ChemSymb,   &     
                                           A%Atom(i)%Z,          &
                                           A%Atom(i)%Mult,       &
                                           A%Atom(i)%X,          &
                                           A%Atom(i)%Occ,        &
                                           A%Atom(i)%UType,      &
                                           A%Atom(i)%ThType,     &
                                           A%Atom(i)%U_iso,      &
                                           A%Atom(i)%U,          &
                                           A%Atom(i)%Magnetic,   &
                                           A%Atom(i)%Moment,     &
                                           A%Atom(i)%SfacSymb,   &
                                           A%Atom(i)%Ind_ff,     &
                                           A%Atom(i)%X_Std,      &      
                                           A%Atom(i)%Occ_Std,    &    
                                           A%Atom(i)%U_iso_Std,  & 
                                           A%Atom(i)%U_Std,      &
                                           A%Atom(i)%Moment_std     
                                           
                if (ierr /=0) then
                   err_CFML%IErr=1
                   err_CFML%Msg="Read_Bin_Atoms_List@CFML_ATOMS: Error reading atoms information!"
                   exit
                end if                            
            end do
            
         type is (MAtm_std_List_Type)   
            do i=1,n
               read(unit=lun,iostat=ierr)  A%Atom(i)%Lab,        &
                                           A%Atom(i)%ChemSymb,   &     
                                           A%Atom(i)%Z,          &
                                           A%Atom(i)%Mult,       &
                                           A%Atom(i)%X,          &
                                           A%Atom(i)%Occ,        &
                                           A%Atom(i)%UType,      &
                                           A%Atom(i)%ThType,     &
                                           A%Atom(i)%U_iso,      &
                                           A%Atom(i)%U,          &
                                           A%Atom(i)%Magnetic,   &
                                           A%Atom(i)%Moment,     &
                                           A%Atom(i)%SfacSymb,   &
                                           A%Atom(i)%Ind_ff,     &
                                           A%Atom(i)%X_Std,      &      
                                           A%Atom(i)%Occ_Std,    &    
                                           A%Atom(i)%U_iso_Std,  & 
                                           A%Atom(i)%U_Std,      &
                                           A%Atom(i)%Moment_std, &
                                           A%Atom(i)%wyck,       &
                                           A%Atom(i)%n_mc,       &
                                           A%Atom(i)%n_dc,       &
                                           A%Atom(i)%Mcs,        &
                                           A%Atom(i)%Mcs_std,    &
                                           A%Atom(i)%Dcs,        &
                                           A%Atom(i)%Dcs_std
      
                                           
                if (ierr /=0) then
                   err_CFML%IErr=1
                   err_CFML%Msg="Read_Bin_Atoms_List@CFML_ATOMS: Error reading atoms information!"
                   exit
                end if                            
            end do
            
         type is (Atm_Ref_List_Type)
            do i=1,n
               read(unit=lun,iostat=ierr)  A%Atom(i)%Lab,        &
                                           A%Atom(i)%ChemSymb,   &     
                                           A%Atom(i)%Z,          &
                                           A%Atom(i)%Mult,       &
                                           A%Atom(i)%X,          &
                                           A%Atom(i)%Occ,        &
                                           A%Atom(i)%UType,      &
                                           A%Atom(i)%ThType,     &
                                           A%Atom(i)%U_iso,      &
                                           A%Atom(i)%U,          &
                                           A%Atom(i)%Magnetic,   &
                                           A%Atom(i)%Moment,     &
                                           A%Atom(i)%SfacSymb,   &
                                           A%Atom(i)%Ind_ff,     &
                                           A%Atom(i)%X_Std,      &      
                                           A%Atom(i)%Occ_Std,    &    
                                           A%Atom(i)%U_iso_Std,  & 
                                           A%Atom(i)%U_Std,      &
                                           A%Atom(i)%Moment_std, &
                                           A%Atom(i)%LX,         &
                                           A%Atom(i)%LOcc,       &
                                           A%Atom(i)%LU_iso,     &
                                           A%Atom(i)%LU,         &
                                           A%Atom(i)%MX,         &
                                           A%Atom(i)%MOcc,       &
                                           A%Atom(i)%MU_iso,     &
                                           A%Atom(i)%MU      
                                           
                if (ierr /=0) then
                   err_CFML%IErr=1
                   err_CFML%Msg="Read_Bin_Atoms_List@CFML_ATOMS: Error reading atoms information!"
                   exit
                end if                            
            end do
              
      end select
      
      close(unit=lun)
   End Subroutine Read_Bin_atom_list
   
   !!----
   !!---- WRITE_BIN_ATOM_LIST
   !!----
   !!----    Write the atoms in the asymmetric unit in a binary file.
   !!----    The file should have been opened with the access="stream" attribute. 
   !!----
   !!---- 12/06/2019 
   !!
   Module Subroutine Write_Bin_Atom_List(filename, A)
      !---- Arguments ----!
      character(len=*),  intent(in) :: filename
      class(alist_type), intent(in) :: A

      !---- Local Variables ----!
      integer                        :: i,lun

      !> Init
      call clear_error()
      if (A%natoms ==0) return
      
      !> Exits file?
      open(newunit=lun,file=trim(filename), access="stream", status="replace", iostat=ierr)
      if (ierr /=0) then
         err_CFML%IErr=1
         err_CFML%Msg="Write_Bin_Atoms_List@CFML_ATOMS: Error opening the binary file "//trim(filename)
         return
      end if  
      
      !> First: Write number of atoms
      write(unit=lun) a%natoms  
      
      !> Second: Write active list
      write(unit=lun) a%active  
      
      select type(A)
         type is (Atm_List_Type)
            do i=1,a%natoms
               write(unit=lun,iostat=ierr)  A%Atom(i)%Lab,        &
                                            A%Atom(i)%ChemSymb,   &     
                                            A%Atom(i)%Z,          &
                                            A%Atom(i)%Mult,       &
                                            A%Atom(i)%X,          &
                                            A%Atom(i)%Occ,        &
                                            A%Atom(i)%UType,      &
                                            A%Atom(i)%ThType,     &
                                            A%Atom(i)%U_iso,      &
                                            A%Atom(i)%U,          &
                                            A%Atom(i)%Magnetic,   &
                                            A%Atom(i)%Moment,     &
                                            A%Atom(i)%SfacSymb,   &
                                            A%Atom(i)%Ind_ff
                if (ierr /=0) then
                   err_CFML%IErr=1
                   err_CFML%Msg="Write_Bin_Atoms_List@CFML_ATOMS: Error writting atoms information!"
                   exit
                end if                            
            end do
               
         type is (Atm_std_List_Type)   
            do i=1,a%natoms
               write(unit=lun,iostat=ierr)  A%Atom(i)%Lab,        &
                                            A%Atom(i)%ChemSymb,   &     
                                            A%Atom(i)%Z,          &
                                            A%Atom(i)%Mult,       &
                                            A%Atom(i)%X,          &
                                            A%Atom(i)%Occ,        &
                                            A%Atom(i)%UType,      &
                                            A%Atom(i)%ThType,     &
                                            A%Atom(i)%U_iso,      &
                                            A%Atom(i)%U,          &
                                            A%Atom(i)%Magnetic,   &
                                            A%Atom(i)%Moment,     &
                                            A%Atom(i)%SfacSymb,   &
                                            A%Atom(i)%Ind_ff,     &
                                            A%Atom(i)%X_Std,      &      
                                            A%Atom(i)%Occ_Std,    &    
                                            A%Atom(i)%U_iso_Std,  & 
                                            A%Atom(i)%U_Std,      &
                                            A%Atom(i)%Moment_std   
               
                if (ierr /=0) then
                   err_CFML%IErr=1
                   err_CFML%Msg="Write_Bin_Atoms_List@CFML_ATOMS: Error writting atoms information!"
                   exit
                end if                            
            end do
            
         type is (MAtm_std_List_Type) 
            do i=1,a%natoms
               write(unit=lun,iostat=ierr)  A%Atom(i)%Lab,        &
                                            A%Atom(i)%ChemSymb,   &     
                                            A%Atom(i)%Z,          &
                                            A%Atom(i)%Mult,       &
                                            A%Atom(i)%X,          &
                                            A%Atom(i)%Occ,        &
                                            A%Atom(i)%UType,      &
                                            A%Atom(i)%ThType,     &
                                            A%Atom(i)%U_iso,      &
                                            A%Atom(i)%U,          &
                                            A%Atom(i)%Magnetic,   &
                                            A%Atom(i)%Moment,     &
                                            A%Atom(i)%SfacSymb,   &
                                            A%Atom(i)%Ind_ff,     &
                                            A%Atom(i)%X_Std,      &      
                                            A%Atom(i)%Occ_Std,    &    
                                            A%Atom(i)%U_iso_Std,  & 
                                            A%Atom(i)%U_Std,      &
                                            A%Atom(i)%Moment_std, &
                                            A%Atom(i)%wyck,       &
                                            A%Atom(i)%n_mc,       &
                                            A%Atom(i)%n_dc,       &
                                            A%Atom(i)%Mcs,        &
                                            A%Atom(i)%Mcs_std,    &
                                            A%Atom(i)%Dcs,        &
                                            A%Atom(i)%Dcs_std
                                           
                if (ierr /=0) then
                   err_CFML%IErr=1
                   err_CFML%Msg="Write_Bin_Atoms_List@CFML_ATOMS: Error reading atoms information!"
                   exit
                end if                            
            end do 
             
         type is (Atm_Ref_List_Type)
            do i=1,a%natoms
               write(unit=lun,iostat=ierr)  A%Atom(i)%Lab,        &
                                            A%Atom(i)%ChemSymb,   &     
                                            A%Atom(i)%Z,          &
                                            A%Atom(i)%Mult,       &
                                            A%Atom(i)%X,          &
                                            A%Atom(i)%Occ,        &
                                            A%Atom(i)%UType,      &
                                            A%Atom(i)%ThType,     &
                                            A%Atom(i)%U_iso,      &
                                            A%Atom(i)%U,          &
                                            A%Atom(i)%Magnetic,   &
                                            A%Atom(i)%Moment,     &
                                            A%Atom(i)%SfacSymb,   &
                                            A%Atom(i)%Ind_ff,     &
                                            A%Atom(i)%X_Std,      &      
                                            A%Atom(i)%Occ_Std,    &    
                                            A%Atom(i)%U_iso_Std,  & 
                                            A%Atom(i)%U_Std,      &
                                            A%Atom(i)%Moment_std, &
                                            A%Atom(i)%LX,         &
                                            A%Atom(i)%LOcc,       &
                                            A%Atom(i)%LU_iso,     &
                                            A%Atom(i)%LU,         &
                                            A%Atom(i)%MX,         &
                                            A%Atom(i)%MOcc,       &
                                            A%Atom(i)%MU_iso,     &
                                            A%Atom(i)%MU 
               
                if (ierr /=0) then
                   err_CFML%IErr=1
                   err_CFML%Msg="Write_Bin_Atoms_List@CFML_ATOMS: Error writting atoms information!"
                   exit
                end if                            
            end do  
      end select
      
      !> Close file
      close(unit=lun)
   End Subroutine Write_Bin_atom_list
   
End SubModule Atm_004   