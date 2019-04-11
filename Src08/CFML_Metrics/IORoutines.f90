Submodule (CFML_Metrics) IORoutines
   !---- Variables ----!
   implicit none
   
 Contains
    !!----
    !!---- WRITE_CRYSTAL_CELL
    !!----    Writes the information of the Cell
    !!----
    !!---- 10/04/2019 
    !!
    Module Subroutine Write_Crystal_Cell(Cell, Iunit)
       !---- Arguments ----!
       class(Cell_Type),  intent(in) :: Cell         ! Cell object
       integer, optional, intent(in) :: Iunit

       !---- Local variables ----!
       integer :: lun, i, j

       !> Init
       lun=6
       if (present(iunit)) lun=iunit
       
       !> Print Zone 
       Write(unit=lun,fmt="(/,a)")    "        Metric information:"
       Write(unit=lun,fmt="(a,/)")    "        -------------------"
       Write(unit=lun,fmt="(a,/)")    " => Direct cell parameters:"
       Write(unit=lun,fmt="(3(a,f12.4))")"         a = ", Cell%cell(1),"      b = ", Cell%cell(2), "      c = ", Cell%cell(3)
       Write(unit=lun,fmt="(3(a,f12.3))")"     alpha = ", Cell%ang(1) ,"   beta = ", Cell%ang(2) , "  gamma = ", Cell%ang(3)
       Write(unit=lun,fmt="(a,f12.4)")   "                        Direct Cell Volume = ",Cell%Vol
       
       select type (cell)
          class is (cell_g_type)
             Write(unit=lun,fmt="(/,a,/)")     " => Reciprocal cell parameters:"
             Write(unit=lun,fmt="(3(a,f12.6))")"         a*= ", Cell%rcell(1),"      b*= ",Cell%rcell(2),"      c*= ", Cell%rcell(3)
             Write(unit=lun,fmt="(3(a,f12.3))")"     alpha*= ", Cell%rang(1) ,"   beta*= ",Cell%rang(2) ,"  gamma*= ", Cell%rang(3)
             Write(unit=lun,fmt="(a,f12.8)")   "                    Reciprocal Cell Volume = ",1.0_cp/cell%vol
             Write(unit=lun,fmt="(/,a,/)")     " => Direct and Reciprocal Metric Tensors:"
             Write(unit=lun,fmt="(a)")         "                   GD                                       GR"

             do i=1,3
                Write(unit=lun,fmt="(3f12.4,a,3f12.6)") (Cell%GD(i,j),j=1,3),"      ", (Cell%GR(i,j),j=1,3)
             end do

             if (Cell%CartType == "A") then
                Write(unit=lun,fmt="(/,a,/)") " =>  Cartesian frame: x // a; y is in the ab-plane; z is x ^ y   "
             else
                Write(unit=lun,fmt="(/,a,/)") " =>  Cartesian frame: z // c; y is in the bc-plane; x is y ^ z   "
             end if

             Write(unit=lun,fmt="(a)")       "     Crystal_to_Orthonormal_Matrix              Orthonormal_to_Crystal Matrix"
             Write(unit=lun,fmt="(a)")       "              Cr_Orth_cel                               Orth_Cr_cel  "
       
             do i=1,3
                Write(unit=lun,fmt="(3f12.4,a,3f12.6)") (Cell%Cr_Orth_cel(i,j),j=1,3),"      ", (Cell%Orth_Cr_cel(i,j),j=1,3)
             end do

             Write(unit=lun,fmt="(/,a)")     "     Busing-Levy B-matrix: Hc=B.H            Inverse of the Busing-Levy B-matrix"
             Write(unit=lun,fmt="(a)")       "                BL_M                                      BL_Minv  "
       
             do i=1,3
                Write(unit=lun,fmt="(3f12.6,a,3f12.4)") (Cell%BL_M(i,j),j=1,3),"      ", (Cell%Inv_BL_M(i,j),j=1,3)
             end do
            
          class is (cell_ls_type)
             Write(unit=lun,fmt="(/,a,/)")     " => Refinement codes for cell parameters:"
             Write(unit=lun,fmt="(3(a,i12))")"         a = ", Cell%lcell(1),"      b = ", Cell%lcell(2), "      c = ", Cell%lcell(3)
             Write(unit=lun,fmt="(3(a,i12))")"     alpha = ", Cell%lang(1) ,"   beta = ", Cell%lang(2) , "  gamma = ", Cell%lang(3)
       
       end select
       
       return
    End Subroutine Write_Crystal_Cell
    
    
    !!----
    !!---- READ_CRYSTAL_CELL
    !!----    Read the cell characteristics from a binary file associated to the
    !!----    logical unit lun. 
    !!----    The file is supposed to be opened with form="unformatted",
    !!----    access="stream" or equivalent
    !!----
    !!---- 10/04/2019 
    !!
    Module Subroutine Read_Bin_Crystal_Cell(Cell,Iunit)
       !---- Arguments ----!
       class(Cell_Type),  intent(out) :: Cell       ! Cell object
       Integer,           intent(in)  :: Iunit
       
       !---- Local Variables ----!
       logical :: info
       integer :: ier
        
       !> Control
       inquire(unit=Iunit, opened=info)
       if (.not. info) then
          err_CFML%IErr=1
          err_CFML%Msg="READ_CRYSTAL_CELL@METRICS: Error reading Cell parameters from a closed binary file!"
          return 
       end if 
       
       select type (cell)
          type is (cell_type)
             read(unit=iunit,iostat=ier) Cell%cell, Cell%ang, Cell%scell, Cell%sang, Cell%vol, Cell%svol
             
          type is (cell_G_type)
             read(unit=iunit,iostat=ier) Cell%cell,  Cell%ang,  Cell%scell, Cell%sang, Cell%vol, Cell%svol, &
                                         Cell%rcell, Cell%rang, Cell%rvol,  Cell%GD,   Cell%GR,             &
                                         Cell%Cr_Orth_cel, Cell%Orth_Cr_cel, Cell%BL_M, Cell%Inv_BL_M,       &
                                         Cell%CartType          

          type is (cell_LS_type)
             read(unit=iunit,iostat=ier) Cell%cell, Cell%ang, Cell%scell, Cell%sang, Cell%vol, Cell%svol,   &
                                         Cell%lcell,Cell%lang
                               
          type is (cell_GLS_type)
             read(unit=iunit,iostat=ier) Cell%cell,  Cell%ang,  Cell%scell, Cell%sang, Cell%vol, Cell%svol, &
                                         Cell%rcell, Cell%rang, Cell%rvol,  Cell%GD,   Cell%GR,             &
                                         Cell%Cr_Orth_cel, Cell%Orth_Cr_cel, Cell%BL_M, Cell%Inv_BL_M,       &
                                         Cell%CartType, cell%lcell, cell%lang 
       end select
       if (ier /= 0) then
          Err_CFML%IErr=1
          Err_CFML%Msg="READ_CRYSTAL_CELL@METRICS:  Error reading Cell parameters from a binary file!"
       end if 
       
       return
    End Subroutine Read_Bin_Crystal_Cell
    
    !!----
    !!---- WRITE_CRYSTAL_CELL
    !!----    Writes the cell characteristics in a binary file associated to the
    !!----    logical unit lun. The file is supposed to be opened with form="unformatted",
    !!----    access="stream" or equivalent
    !!----
    !!---- 10/04/2019 
    !!
    Module Subroutine Write_Bin_Crystal_Cell(Cell,Iunit)
       !---- Arguments ----!
       class(Cell_Type),  intent(in) :: Cell       ! Cell object
       Integer,           intent(in) :: Iunit
       
       !---- Local Variables ----!
       logical :: info
        
       !> Control
       inquire(unit=Iunit, opened=info)
       if (.not. info) return
       
       select type (cell)
          type is (cell_type)
             write(unit=iunit) Cell%cell, Cell%ang, Cell%scell, Cell%sang, Cell%vol, Cell%svol
             
          type is (cell_G_type)
             write(unit=iunit) Cell%cell,  Cell%ang,  Cell%scell, Cell%sang, Cell%vol, Cell%svol, &
                               Cell%rcell, Cell%rang, Cell%rvol,  Cell%GD,   Cell%GR,             &
                               Cell%Cr_Orth_cel, Cell%Orth_Cr_cel, Cell%BL_M, Cell%Inv_BL_M,      &
                               Cell%CartType          

          type is (cell_LS_type)
             write(unit=iunit) Cell%cell, Cell%ang, Cell%scell, Cell%sang, Cell%vol, Cell%svol,   &
                               Cell%lcell,Cell%lang
                               
          type is (cell_GLS_type)
             write(unit=iunit) Cell%cell,  Cell%ang,  Cell%scell, Cell%sang, Cell%vol, Cell%svol, &
                               Cell%rcell, Cell%rang, Cell%rvol,  Cell%GD,   Cell%GR,             &
                               Cell%Cr_Orth_cel, Cell%Orth_Cr_cel, Cell%BL_M, Cell%Inv_BL_M,      &
                               Cell%CartType, cell%lcell, cell%lang 
       end select
       
       return
    End Subroutine Write_Bin_Crystal_Cell
    
End Submodule IORoutines  