!!----
!!----
!!----
SubModule (CFML_Atoms) ReadnWrite_Atoms_Binary
  implicit none
   Contains
   !!----
   !!---- READ_BIN_ATOM_LIST
   !!----    Reads the atoms in the asymmetric unit in a binary file.
   !!----
   !!----    The procedure reads in the given order a series of bytes corresponding to the
   !!----    components of the type Ats. The full structure is re-allocated inside the procedure
   !!----    before reading the components.
   !!----
   !!----    The number of atoms is the first element read in the file.
   !!----
   !!----    Note: the procedure is general but the user need to know which type of atom will be
   !!----          read.
   !!----
   !!---- 12/06/2019
   !!
   Module Subroutine Read_Bin_Atom_List(filename, A, Type_Atm) !, d)
      !---- Arguments ----!
      character(len=*),   intent(in)     :: filename
      type(atlist_type),  intent(in out) :: A
      character(len=*),   intent(in)     :: Type_Atm
      !integer,            intent(in)     :: d !Number of k-vectors
      !---- Local Variables ----!
      integer                            :: i,n,ierr,lun
      type (atm_type)      :: atm
      type (atm_std_type)  :: atms
      !type (matm_std_type) :: matm
      type (atm_ref_type)  :: atr

      !> Init
      call clear_error()

      !> Exits file?
      open(newunit=lun,file=trim(filename), access="stream", status="old", iostat=ierr)
      if (ierr /=0) then
         err_CFML%IErr=1
         err_CFML%Msg="Read_Bin_Atom_List@CFML_ATOMS: Error opening the binary file "//trim(filename)
         return
      end if

      !> First: read number of atoms
      n=0
      call Allocate_Atom_List(N, A,Type_Atm,3)

      read(unit=lun,iostat=ierr) n
      if (ierr /= 0) then
         err_CFML%IErr=1
         err_CFML%Msg="Read_Bin_Atom_List@CFML_ATOMS: Error reading number of atoms in the binary file "//trim(filename)
         close(unit=lun)
         return
      end if
      if (n <= 0) return

      !> Allocating
      call Allocate_Atom_List(N, A,Type_Atm,3)

      !> Read active
      read(unit=lun,iostat=ierr)  A%active
      if (ierr /= 0) then
         err_CFML%IErr=1
         err_CFML%Msg="Read_Bin_Atoms_List@CFML_ATOMS: Error reading active atoms!"
         close(unit=lun)
         return
      end if

      !> Load information
      select type (aat => A%Atom)
         type is (atm_type)
            do i=1,n
               read(unit=lun,iostat=ierr) Atm
               if (ierr /=0) then
                  err_CFML%IErr=1
                  err_CFML%Msg="Write_Bin_Atom_List@CFML_ATOMS: Error writting atoms information!"
                  exit
               end if
               aat(i)=atm
            end do

         type is (atm_std_type)
            do i=1,n
               read(unit=lun,iostat=ierr) Atms
               if (ierr /=0) then
                  err_CFML%IErr=1
                  err_CFML%Msg="Write_Bin_Atom_List@CFML_ATOMS: Error writting atoms information!"
                  exit
               end if
               aat(i)=atms
            end do

         !type is (matm_std_type)
         !   do i=1,n
         !      read(unit=lun,iostat=ierr) matm
         !      if (ierr /=0) then
         !         err_CFML%IErr=1
         !         err_CFML%Msg="Write_Bin_Atom_List@CFML_ATOMS: Error writting atoms information!"
         !         exit
         !      end if
         !      aat(i)=matm
         !   end do

         type is (atm_ref_type)
            do i=1,n
               write(unit=lun,iostat=ierr) atr
               if (ierr /=0) then
                  err_CFML%IErr=1
                  err_CFML%Msg="Write_Bin_Atom_List@CFML_ATOMS: Error writting atoms information!"
                  exit
               end if
               aat(i)=atr
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
      character(len=*),   intent(in) :: filename
      type(atlist_type), intent(in) :: A

      !---- Local Variables ----!
      integer                        :: i,n,lun,ierr
      type (atm_type)      :: atm
      type (atm_std_type)  :: atms
      !type (matm_std_type) :: matm
      type (atm_ref_type)  :: atr

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
      n=a%natoms
      write(unit=lun) n

      !> Second: Write active list
      write(unit=lun) a%active

      select type (aat => A%Atom)
         type is (atm_type)
            do i=1,n
               atm=aat(i)
               write(unit=lun,iostat=ierr) Atm
               if (ierr /=0) then
                  err_CFML%IErr=1
                  err_CFML%Msg="Write_Bin_Atom_List@CFML_ATOMS: Error writting atoms information!"
                  exit
               end if
            end do

         type is (atm_std_type)
            do i=1,n
               atms=aat(i)
               write(unit=lun,iostat=ierr) Atms
               if (ierr /=0) then
                  err_CFML%IErr=1
                  err_CFML%Msg="Write_Bin_Atom_List@CFML_ATOMS: Error writting atoms information!"
                  exit
               end if
            end do

         !type is (matm_std_type)
         !   do i=1,n
         !      matm=aat(i)
         !      write(unit=lun,iostat=ierr) matm
         !      if (ierr /=0) then
         !         err_CFML%IErr=1
         !         err_CFML%Msg="Write_Bin_Atom_List@CFML_ATOMS: Error writting atoms information!"
         !         exit
         !      end if
         !   end do

         type is (atm_ref_type)
            do i=1,n
               atr=aat(i)
               write(unit=lun,iostat=ierr) atr
               if (ierr /=0) then
                  err_CFML%IErr=1
                  err_CFML%Msg="Write_Bin_Atom_List@CFML_ATOMS: Error writting atoms information!"
                  exit
               end if
            end do
      end select

      !> Close file
      close(unit=lun)
   End Subroutine Write_Bin_atom_list

End SubModule ReadnWrite_Atoms_Binary