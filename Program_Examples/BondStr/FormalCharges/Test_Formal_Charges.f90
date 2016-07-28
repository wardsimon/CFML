!!----
!!---- This program uses the Set_Charges subroutine to
!!---- set atomic charges in a crystal structure read
!!---- from a .cif .or. cfl file. It generates a .cfl file
!!---- called xxx_auto.cfl
!!----
Program Test_Formal_Charges

  !---- Use Modules ----!
  Use CFML_Crystal_Metrics,             Only: Crystal_Cell_Type
  Use CFML_Crystallographic_Symmetry,   Only: Space_Group_Type
  Use CFML_Atom_TypeDef,                Only: Atom_List_Type
  Use CFML_IO_Formats,                  Only: Readn_set_Xtal_Structure,ERR_Form_Mess,err_form,File_List_Type,Write_CFL
  Use CFML_BVS_Energy_Calc,             Only: Err_Conf, Err_Conf_Mess
  Use CFML_Formal_Charges

  !---- Variables ----!
  Implicit None

  Type (Space_Group_Type)                    :: SpGr
  Type (Crystal_Cell_Type)                   :: Cell
  Type (Atom_List_Type)                      :: A
  Type (File_List_Type)                      :: Fich_cfl,Fich_cif

  Character(len=132)                         :: line
  Character(len=256)                         :: filcod

  Integer                                    :: i,j,k
  Integer                                    :: narg,lun=19
  
  Real(kind=cp), Dimension(:), Allocatable   :: charges
  Logical                                    :: arggiven=.False.,esta,cif=.False.,out_cif=.False.,input_charges=.False.
  Logical                                    :: identical=.True.

  ! Arguments on the command line

  narg=Command_Argument_Count()

  If (narg > 0) Then
     Call GET_COMMAND_ARGUMENT(1,filcod)
     i=Index(filcod,".cfl")
     If(i /= 0) filcod=filcod(1:i-1)
     i=Index(filcod,".cif")
     If(i /= 0) filcod=filcod(1:i-1)
     arggiven=.True.
  End If

  If(.Not. arggiven) Then
     Write(unit=*,fmt="(a)",advance='no') " => Code of the file xx.cfl(cif) (give xx): "
     Read(unit=*,fmt="(a)") filcod
     If(Len_trim(filcod) == 0) Stop
  End If

  Inquire(file=Trim(filcod)//".cfl",exist=esta)

  If(esta) Then
     Call Readn_set_Xtal_Structure(Trim(filcod)//".cfl",Cell,SpGr,A,Mode="CFL",file_list=fich_cfl)
     cif=.False.
     out_cif=.True.
  Else
     Inquire(file=Trim(filcod)//".cif",exist=esta)
     If(.Not. esta) Then
        Write(unit=*,fmt="(a)") " File: "//Trim(filcod)//".cfl (or .cif) doesn't exist!"
        Stop
     End If
     Call Readn_set_Xtal_Structure(Trim(filcod)//".cif",Cell,SpGr,A,Mode="CIF",file_list=fich_cif)
     ! Store oxidation states, if any
     i=1
     Do 
        line=adjustl(fich_cif%line(i))
        If (line(1:27) == "_atom_type_oxidation_number") Then 
           Do j = 1 , A%Natoms
              line=Adjustl(fich_cif%line(i+j))
              k = Index(line," ")
              Read(line(k:),*) A%Atom(j)%Charge
           End Do
           Exit
        End If
        i=i+1
        If (i > fich_cif%nlines) Exit
     End Do
     cif=.True.
  End If

  If (Err_Form) Then
     Write(unit=*,fmt="(a)") Trim(ERR_Form_Mess)
     Write(unit=*,fmt="(/,a)") " => PROGRAM TEST_FORMAL_CHARGES finished in error!"
     Stop
  End If

  ! Objects Cell, SpGr and A must be built before calling this subroutine

  Allocate(charges(A%Natoms))
  Call Set_Charges(SpGr,Cell,A,charges,filcod)

  If (Err_Conf) Then
     Write(unit=*,fmt="(a)") Trim(Err_Conf_Mess)
     Write(unit=*,fmt="(/,a)") " => PROGRAM TEST_FORMAL_CHARGES finished in error!"
     Stop
  End If

  If (Err_Char) Then
     Write(unit=*,fmt="(a)") Trim(Err_Char_Mess)
     Write(unit=*,fmt="(/,a)") " => PROGRAM TEST_FORMAL_CHARGES finished in error!"
     Stop
  End If

  If (WARN_Char) Write(unit=*,fmt="(a)") Trim(WARN_Char_Mess)

  ! Check if charges have given in the input file
  
  i = 1
  Do 
     If (A%Atom(i)%Charge /= 0) input_charges = .True.
     If (input_charges) Exit
     i = i + 1
     If (i > A%Natoms) Exit
  End Do

  ! Compare input charges with the automatic setting

  If (input_charges) Then
     i = 1
     Do
        If (A%Atom(i)%Charge /= charges(i)) identical = .False.
        If (.Not. identical) Exit
        i = i + 1
        If (i > A%Natoms) Exit
     End Do
  End If

  ! Write CFL files

  If (input_charges .And. cif) Then
     Open(unit=lun,file=Trim(filcod)//".cfl",status="replace",action="write")
     Call Write_CFL(lun,Cell,SpGr,A)
     Close(unit=lun)
  End If

  Open(unit=lun,file=Trim(filcod)//"_auto.cfl",status="replace",action="write")
  Do i = 1 , A%Natoms
     A%Atom(i)%Charge = charges(i)
  End Do
  Call Write_CFL(lun,Cell,SpGr,A)
  Close(unit=lun)

End Program Test_Formal_Charges
