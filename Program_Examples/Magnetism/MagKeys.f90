Program Get_Read_Write_Magkeys
   !---- Auxiliary program to check that magnetic keys and codes work
   !  Key_Code=(/ "Rxyz","Ixyz","Mxyz", C1-12/)
   !  Code_Nam=(/ "Rx_", "Ry_", "Rz_", "Ix_", "Iy_", "Iz_", also spherical, "MagPh_"/)

   !---- Use Modules ----!
   !use f2kcli  !comment for non-Lahey
   use CFML_crystallographic_symmetry, only: space_group_type, Write_SpaceGroup
   use CFML_string_utilities,          only: u_case, Get_LogUnit
   use CFML_Atom_TypeDef,              only: Atom_List_Type, Write_Atom_List, MAtom_list_Type
   use CFML_crystal_Metrics,           only: Crystal_Cell_Type, Write_Crystal_Cell
   use CFML_IO_Formats,                only: Readn_set_Xtal_Structure,err_form_mess,err_form,     &
                                             file_list_type
   use CFML_Keywords_Code_Parser,      only: Init_RefCodes,Read_RefCodes_File,Write_Info_RefCodes,&
                                             Allocate_VParam,Err_RefCodes, err_refcodes_mess,     &
                                             NP_Max,NP_Refi

   use CFML_Magnetic_Symmetry

   !---- Local Variables ----!
   implicit none
   type (file_list_type)       :: fich_cfl
   type (Space_Group_Type)     :: SpG
   type (MagSymm_k_Type)       :: MGp
   type ( Atom_list_Type)      :: A
   type (mAtom_list_Type)      :: mA
   type (Crystal_Cell_Type)    :: Cell

   character(len=256)                 :: filcod     !Name of the input file
   integer                            :: narg, lun=1, n_ini,n_end, i
   logical                            :: esta, arggiven=.false., rest_file=.false.

    narg=COMMAND_ARGUMENT_COUNT()

    if(narg > 0) then
            call GET_COMMAND_ARGUMENT(1,filcod)
            arggiven=.true.
    end if

    write (unit=*,fmt="(/,/,6(a,/))")                                                  &
            "                 ------ Program for Reading MagKeyS ------"            , &
            "                    ---- Version 0.0 December-2011 ----"                          , &
            "    *******************************************************************************"  , &
            "    * Reads magnetic keys from *.CFL file and writes them into *.OUT file  "  , &
            "    *******************************************************************************"  , &
            "                             (OZ- December-2011 )"
   write (unit=*,fmt=*) " "

   if(.not. arggiven) then
     write(unit=*,fmt="(a)",advance='no') " => Code of the file xx.cfl (give xx): "
     read(unit=*,fmt="(a)") filcod
     if(len_trim(filcod) == 0) stop
   end if

   open(unit=lun,file=trim(filcod)//".out", status="replace",action="write")
    write(unit=lun,fmt="(/,/,6(a,/))")                                                  &
            "                 ------ Program for Reading MagKeyS ------"            , &
            "                    ---- Version 0.0 December-2011 ----"                          , &
            "    *******************************************************************************"  , &
            "    * Reads magnetic keys from *.CFL file and writes them into *.OUT file  "  , &
            "    *******************************************************************************"  , &
            "                             (OZ- December-2011 )"

   inquire(file=trim(filcod)//".cfl",exist=esta)
   if( .not. esta) then
     write(unit=*,fmt="(a)") " File: "//trim(filcod)//".cfl doesn't exist!"
     stop
   end if

   call Readn_Set_Xtal_Structure(trim(filcod)//".cfl",Cell,SpG,A,Mode="CFL",file_list=fich_cfl)

   If(err_form) then
     write(unit=*,fmt="(a)") trim(err_form_mess)
   Else
     write(unit=lun,fmt="(/,a)") "    =========================="
     write(unit=lun,fmt="( a )") "    Text of the input CFL file"
     write(unit=lun,fmt="(a,/)") "    =========================="

     do i=1,fich_cfl%nlines
       write(unit=lun,fmt="(a,i5,a)") " Line:",i,"  "//fich_cfl%line(i)
     end do

     call Write_Crystal_Cell(Cell,lun)
     call Write_SpaceGroup(SpG,lun)
     call Write_Atom_List(A,lun=lun)

     n_ini=1
     n_end=fich_cfl%nlines
     call Readn_Set_Magnetic_Structure(fich_cfl,n_ini,n_end,MGp,mA,Cell=Cell)

       if(err_MagSym) then
         write(unit=*,fmt="(a)") "   "//err_MagSym_Mess
         stop
       end if
     call Write_Magnetic_Structure(lun,MGp,mA)

     NP_Max= mA%natoms*17 ! Rxyz,Ixyz (also Rm, Rphi, Rth), mPhase, C1-12
     call Allocate_Vparam(NP_Max)

     !Determine the mag refinement codes from vary/fix instructions
     call Init_RefCodes(mA)
     call Read_RefCodes_File(fich_cfl,1,fich_cfl%nlines,mA,Spg)

     if(Err_RefCodes) then
       write(unit=*,fmt="(a)") trim(err_refcodes_mess)
       write(*,*) 'Err_RefCodes'
     end if

     call Write_Info_RefCodes(mA,mGp,lun)

     write(unit=*,fmt="(a)") " Normal End of: PROGRAM FOR READING MagKeyS "
     write(unit=*,fmt="(a)") " Results in File: "//trim(filcod)//".out"
   End If

   close(unit=lun)
   stop
End Program Get_Read_Write_Magkeys
