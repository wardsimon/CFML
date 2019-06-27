!!----
!!----
!!----
SubModule (CFML_IOForm) IOF_002
   Contains
   
   !!----
   !!---- READ_CFL_ATOM
   !!----    Subroutine to read atoms parameters 
   !!----   
   !!----    Format: ATOM   Label  ChemSymb   x y z B Occ Us or Moment 
   !!----      
   !!----     For charge obtained from Label: Label[+/-][number]
   !!----       
   !!----
   !!---- 25/06/2019 
   !!
   Module Subroutine Read_CFL_Atom(lines,n_ini, n_end, At_List)
      !---- Arguments ----!
      character(len=*), dimension(:), intent(in)     :: lines
      integer,                        intent(in out) :: n_ini
      integer,                        intent(in)     :: n_end
      Type (AtList_Type),             intent(out)    :: At_List

      !---- Local variables -----!
      character(len=12), parameter      :: DIGPM="0123456789+-"
      character(len=:), allocatable     :: line
      character(len=180), dimension(1)  :: extline
      character(len=10)                 :: dire
      character(len=20)                 :: label
      character(len=2)                  :: Chemsymb
      integer                           :: i, n, na, npos, nlong, iv
      real(kind=cp), dimension (20)     :: vet1, vet2
      
       
      !> Init
      call clear_error()
      call Allocate_Atom_List(0, At_list)
      
      !> Calculate number of Atoms
      na=0
      do i=n_ini, n_end
         !> Checks
         if (len_trim(lines(i)) <=0) cycle                 ! Empty line
         line=adjustl(lines(i))
         if (line(1:1) =='!' .or. line(1:1) =='#') cycle   ! Comment line
         npos=index(u_case(line),'ATOM')
         if (npos > 0) na=na+1
      end do
      if (na == 0) return             ! No atoms in the lines 
      
      !> Allocate List
      call Allocate_Atom_List(na, At_list)                             
      
      !> Read Information
      na=0
      do i=n_ini, n_end
         !> Checks
         if (len_trim(lines(i)) <=0) cycle                 ! Empty line
         line=adjustl(lines(i))
         if (line(1:1) =='!' .or. line(1:1) =='#') cycle   ! Comment line

         !> Truncate line from symbols: # and !         
         npos=index(line,'!')                             ! Cut line according to ! 
         if (npos > 0) line=line(:npos-1)
         npos=index(line,'#')
         if (npos > 0) line=line(:npos-1)
         
         !> Tabs
         do
            npos=index(line,TAB)
            if (npos == 0) exit
            line(npos:npos)=' ' 
         end do
         
         !> ATOM Directive
         call cut_string(line,nlong,dire)
         if (u_case(trim(dire)) /= 'ATOM') cycle
         
         na=na+1
         
         !> Label
         call cut_string(line,nlong,label)
         At_list%atom(na)%Lab=trim(label)
         
         !> Charge / Oxidation state
         !> Anions
         npos=index(label,'-')
         if (npos > 0) then
            iv=index(DIGPM,label(npos+1:npos+1))
            if (iv ==0 .or. iv > 10) then
               err_CFML%Ierr=1
               err_CFML%Msg='Read_CFL_Atom@CFML_IOFORM: Error reading Charge information: '//trim(line)
               At_list%Natoms=Max(0,na-1)
               return
            end if 
            At_list%atom(na)%charge=-(iv-1)  
         end if   
         
         !> Cations
         npos=index(label,'+')
         if (npos > 0) then
            iv=index(DIGPM,label(npos+1:npos+1))
            if (iv ==0 .or. iv > 10) then
               err_CFML%Ierr=1
               err_CFML%Msg='Read_CFL_Atom@CFML_IOFORM: Error reading Charge information: '//trim(line)
               At_list%Natoms=Max(0,na-1)
               return
            end if
            At_list%atom(na)%charge=iv-1 
         end if   
         
         
         !> Chemical symbol / Scattering Factor
         call cut_string(line,nlong,dire)
         npos=index(DIGPM,dire(2:2))
         if (npos > 0) then
            Chemsymb=u_case(dire(1:1))
         else
            Chemsymb=u_case(dire(1:1))//l_case(dire(2:2))
         end if
         At_list%atom(na)%ChemSymb=ChemSymb
         At_list%atom(na)%SFacSymb=dire
         
         !> Parameters 
         extline='atm '//trim(line) 
         n=1  
         call Read_Key_ValueSTD(extline,n,n,"atm",vet1,vet2,iv) 
         if (iv <= 0) then
            err_CFML%Ierr=1
            err_CFML%Msg='Read_CFL_Atom@CFML_IOFORM: Error reading Atom information: '//trim(line)
            At_list%Natoms=Max(0,na-1)
            return
         end if 
         
         if (iv < 3) then
            err_CFML%Ierr=1
            err_CFML%Msg='Read_CFL_Atom@CFML_IOFORM: Error reading Atom coordinates: '//trim(line)
            At_list%Natoms=Max(0,na-1)
            return
         end if   
         
         !> Coordinates
         At_List%atom(na)%x=vet1(1:3) 
         select type (atm => at_list%atom)
            class is (Atm_Std_Type)
               atm(na)%x_std=vet2(1:3)
         end select 
         
         !> Biso
         if (iv > 3) then
            At_List%atom(na)%U_iso=vet1(4)
            select type (atm => at_list%atom)
               class is (Atm_Std_Type)
                  atm(na)%u_iso_std=vet2(4)
            end select  
         end if 
         
         !> Occ
         if (iv > 4) then
            At_List%atom(na)%Occ=vet1(5)
            select type (atm => at_list%atom)
               class is (Atm_Std_Type)
                  atm(na)%Occ_std=vet2(5)
            end select  
         end if  
         
         select case (iv)
            case (8)  ! Only moments
               At_List%atom(na)%Moment=vet1(6:8)
               select type (atm => at_list%atom)
                  class is (Atm_Std_Type)
                     atm(na)%Moment_std=vet2(6:8)
               end select 
               At_List%atom(na)%Magnetic=.true.
               
            case (11) ! Only anisotropic parameters
               At_List%atom(na)%U=vet1(6:11)
               select type (atm => at_list%atom)
                  class is (Atm_Std_Type)
                     atm(na)%U_std=vet2(6:11)
               end select
               
            case (14) ! Anisotropic + Moments
               At_List%atom(na)%U=vet1(6:11)
               At_List%atom(na)%Moment=vet1(12:14)
               select type (atm => at_list%atom)
                  class is (Atm_Std_Type)
                     atm(na)%U_std=vet2(6:11)
                     atm(na)%Moment_std=vet2(12:14)
               end select
               
            case default
               err_CFML%Ierr=1
               err_CFML%Msg='Read_CFL_Atom@CFML_IOFORM: Error reading Atom information: '//trim(line)
               At_list%Natoms=Max(0,na-1)
               return       
         end select   
         
      end do   
   End Subroutine Read_CFL_Atom
   
   !!----
   !!---- READ_CFL_CELL
   !!----    Obtaining Cell Parameter from CFL Format
   !!----
   !!---- 27/06/2019 
   !!
   Module Subroutine Read_CFL_Cell(lines, n_ini, n_end, Cell)
      !---- Arguments ----!
      character(len=*), dimension(:),  intent(in)     :: lines   ! Containing information
      integer,                         intent(in out) :: n_ini   ! Index to start
      integer,                         intent(in)     :: n_end   ! Index to Finish
      class(Cell_Type),                intent(out)    :: Cell    ! Cell object

      !---- Local variables -----!
      integer                              :: nlong, iv
      real(kind=cp), dimension (6)         :: vcell, std
      
      character(len=132)                   :: line
      character(len=4)                     :: dire

      !> Init
      call clear_error()

      call Read_Key_ValueSTD(lines,n_ini,n_end,"CELL",vcell,std,iv)
      if (iv /= 6) then
         err_CFML%Ierr=1
         err_CFML%Msg="Read_CFL_Cell@CFML_IOForm: Problems reading cell parameters!"  
         return
      end if
      
      call Set_Crystal_Cell(vcell(1:3),vcell(4:6), Cell, Vscell=std(1:3), Vsang=std(4:6))
   End Subroutine Read_CFL_Cell
   
   
    
End SubModule IOF_002   