!!----
!!---- This program uses the Set_Formal_Charges subroutine to
!!---- set atomic charges in a crystal structure read
!!---- from a .cif .or. cfl file. It generates a .cfl file
!!---- called xxx_auto.cfl
!!----
Program Test_Formal_Charges

  !---- Use Modules ----!
  Use CFML_GlobalDeps
  Use CFML_Crystal_Metrics,             Only: Crystal_Cell_Type
  Use CFML_Crystallographic_Symmetry,   Only: Space_Group_Type
  Use CFML_Atom_TypeDef,                Only: Atom_List_Type, write_atom_list
  Use CFML_IO_Formats,                  Only: Readn_set_Xtal_Structure,ERR_Form_Mess,err_form,File_List_Type,Write_CFL
  Use CFML_BVS_Energy_Calc,             Only: Err_Conf, Err_Conf_Mess, Set_Formal_Charges, Err_char,Err_Char_Mess, &
                                              Warn_Char,Warn_Char_Mess
!  Use Formal_charges_mod

  !---- Variables ----!
  Implicit None

  Type (Space_Group_Type)                    :: SpGr
  Type (Crystal_Cell_Type)                   :: Cell
  Type (Atom_List_Type)                      :: A
  Type (File_List_Type)                      :: Fich_cfl,Fich_cif

  Character(len=132)                         :: line
  Character(len=256)                         :: filcod
  character(len=6)                           :: exten,speciess=" "
  character(len=2)                           :: aux

  Integer                                    :: i,j,k,ier,fcount
  Integer                                    :: narg,lun=19, i_buf=4,i_cfl=2

  Logical                                    :: arggiven=.false.,esta,cif=.False.,out_cif=.False.
  Logical                                    :: buffer_file=.false.

  Real(kind=cp) :: tini, tfin
  ! Arguments on the command line

  narg=Command_Argument_Count()

  If (narg > 0) Then
     Call GET_COMMAND_ARGUMENT(1,filcod)
     i=Index(filcod,".cfl")
     If(i /= 0) then
       filcod=filcod(1:i-1)
       exten=".cfl"
     end if
     i=Index(filcod,".cif")
     If(i /= 0) then
       filcod=filcod(1:i-1)
       exten=".cif"
     end if
     i=Index(filcod,".buf")
     If(i /= 0) buffer_file=.true.
     arggiven=.True.
     if(narg == 2) Call GET_COMMAND_ARGUMENT(2,speciess)
  End If
  if(buffer_file) then
    open(unit=i_buf,file=trim(filcod),status="old",action="read",position="rewind")
  end if
  call cpu_time(tini)
  fcount=0
  do !external loop for treating a buffer file

     if(buffer_file) then
       read(unit=i_buf,fmt="(a)",iostat=ier) filcod
       if(ier /= 0) exit
       i=index(filcod,".",back=.true.)
       exten=filcod(i:)
       filcod=filcod(1:i-1)
     end if


     If(.Not. arggiven) Then
        Write(unit=*,fmt="(a)",advance='no') " => Code of the file xx.cfl(cif) (give xx): "
        Read(unit=*,fmt="(a)") filcod
        If(Len_trim(filcod) == 0) Stop
     End If

     Inquire(file=Trim(filcod)//".cfl",exist=esta)
     fich_cif%nlines=0 !necessary to initialize the list of lines inside Readn_set_Xtal_Structure
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

     ! Opening the file for writing
     fcount=fcount+1
     Open(unit=lun,file=Trim(filcod)//".fch",status="replace",action="write")
     write(*,"(a,t50,a,i6)") "  => Treating the file: "//trim(filcod)//trim(exten),"#",fcount

     If (Err_Form) Then
        Write(unit=*,fmt="(a)") Trim(ERR_Form_Mess)
        Write(unit=*,fmt="(/,a)") " => PROGRAM FORMAL_CHARGES finished in error!"
        Stop
     End If
     ! Objects Cell, SpGr and A must be built before calling this subroutine

     Call Set_Formal_Charges(SpGr,Cell,A,eps_val=0.002,iwrt=lun)

     If (Err_Conf) Then
        Write(unit=*,fmt="(a)") Trim(Err_Conf_Mess)
        Write(unit=*,fmt="(/,a)") " => PROGRAM FORMAL_CHARGES finished in error!"
        Stop
     End If

     If (Err_Char) Then
        Write(unit=*,fmt="(a)") Trim(Err_Char_Mess)
        Write(unit=*,fmt="(/,a)") " => PROGRAM FORMAL_CHARGES finished in error!"
        Stop
     End If

     If (WARN_Char) Write(unit=*,fmt="(a)") Trim(WARN_Char_Mess)
     Write(unit=lun,fmt="(a,/)") " => Formal Charges: "
     Do i = 1 , A%Natoms
        Write(unit=lun,fmt="(4x,a4,a1,1x,f5.2)") A%Atom(i)%Lab, ":", A%Atom(i)%Charge
     End Do

     if(len_trim(speciess) == 0) then !look for the cation probable species to be searched
       j=1
       Do i = 1 , A%Natoms
         if(A%Atom(i)%charge < 0.0) cycle
         aux=A%Atom(i)%ChemSymb
         if( aux == "Li") then
          j=i
          exit
         else if( aux == "Na") then
          j=i
          exit
         else if( aux == "H ") then
          j=i
          exit
         end if
       End Do
       i=nint(A%Atom(j)%charge)
       if(i < 0)  then
         write(unit=speciess,fmt="(a,i2)") trim(A%Atom(j)%ChemSymb),i
       else
         write(unit=speciess,fmt="(a,i1)") trim(A%Atom(j)%ChemSymb)//"+",i
       end if
     end if

     Open(unit=i_cfl,file=Trim(filcod)//"_auto.cfl",status="replace",action="write")
     Write(unit=i_cfl,fmt="(a)") "Title  CFL-file generated from by Formal_Charges.f90"
     Call Write_CFL(i_cfl,Cell,SpGr,A)
     Write(unit=i_cfl,fmt="(/,a)") "!   Bond_STR instructions"
     Write(unit=i_cfl,fmt="(a)")   "!   Nx, Ny, Nz, Species, Dmax, Delta(eV): Values within Emin+delta are counted for fractional volume estimation"
     Write(unit=i_cfl,fmt="(a,3i5,a,2f8.2)") "BVEL ",nint(Cell%cell(1:3)*10.0), " "//speciess, 10.0, 3.0
     Write(unit=i_cfl,fmt="(a)") "PERCOLATION   3.5"
     Close(unit=i_cfl)
     Close(unit=lun)
     if(.not. buffer_file) exit
  end do
  call cpu_time(tfin)
  write(unit=*,fmt="(/a,i5)")      "  Total number of treated files: ",fcount
  write(unit=*,fmt="(a,f10.2,a)")  "  CPU-Time: ", tfin-tini," seconds"
  write(unit=*,fmt="(a,f10.2,a)")  "  CPU-Time: ", (tfin-tini)/60.0," minutes"

End Program Test_Formal_Charges
