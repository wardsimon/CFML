!!-------------------------------------------------------
!!---- Crystallographic Fortran Modules Library (CrysFML)
!!-------------------------------------------------------

!!----
!!
 Module CFML_IO_Formats

    !---- Use modules ----!
    Use CFML_GlobalDeps,                only: cp,sp,dp,pi,eps,Write_Date_Time
    Use CFML_Math_General,              only: sind,equal_matrix,Sort
    Use CFML_Math_3D,                   only: determ_a
    Use CFML_String_Utilities
    Use CFML_Crystal_Metrics,           only: Crystal_Cell_Type, Set_Crystal_Cell, Convert_U_Betas, &
                                              Convert_B_Betas, U_Equiv, Convert_Betas_U
    Use CFML_Crystallographic_Symmetry, only: Space_Group_Type, Magnetic_Space_Group_Type,Set_SpaceGroup,     &
                                              Init_Magnetic_Space_Group_Type,Get_Multip_Pos,Get_MagMatSymb,   &
                                              Read_Xsym,Read_Msymm, Setting_Change, get_symsymb,Sym_Oper_type,&
                                              Msym_Oper_Type,is_Lattice_vec,Get_Stabilizer,Get_mOrbit
    Use CFML_Atom_TypeDef,              only: Atom_Type, Init_Atom_Type,atom_list_type,         &
                                              Allocate_atom_list, Deallocate_atom_list
    Use CFML_Molecular_Crystals,        only: Err_Molec, Err_Molec_Mess,Molecular_Crystal_Type, &
                                              Read_Molecule, Set_Euler_Matrix, Write_Molecule
    Use CFML_Geometry_Calc,             only: Point_List_Type, Get_Euler_from_Fract
    Use CFML_Diffraction_Patterns,      only: Diffraction_Pattern_type
    Use CFML_Scattering_Chemical_Tables,only: Set_Magnetic_Form, Remove_Magnetic_Form, num_mag_form, &
                                              Magnetic_Form, get_magnetic_form_factor
    Use CFML_Magnetic_Groups
    Use CFML_EisPack,                   only: rg_ort

    !---- Variables ----!
    implicit none

    private


    !---- List of public functions ----!

    !---- List of public subroutines ----!
    public :: Init_Err_Form, Read_Atom, Read_Cell, Read_Cif_Atom, Read_Cif_Cell,                 &
              Read_Cif_Cont, Read_Cif_Hall, Read_Cif_Hm, Read_Cif_Lambda, Read_Cif_Symm,         &
              Read_Cif_Title, Read_Cif_Z, Read_File_Atom, Read_File_Spg, Read_Cif_ChemicalName,  &
              Read_File_Transf, Read_Shx_Atom, Read_Shx_Cell, Read_Shx_Cont, Read_Shx_Fvar,      &
              Read_Shx_Latt, Read_Shx_Symm, Read_Shx_Titl, Read_Uvals, Write_Cif_Powder_Profile, &
              Write_Cif_Template, Write_Shx_Template, Read_File_rngSINTL, Read_File_Lambda,      &
              Get_job_info, File_To_FileList, Get_Phases_File, Read_Cif_Pressure,                &
              Read_Cif_Temp,Readn_Set_Magnetic_Space_Group, Set_Magnetic_Space_Group,            &
              Cleanup_Symmetry_Operators,Write_mCIF, Get_Refinement_Codes, Get_moment_ctr,       &
              Readn_Set_Magnetic_Structure_MCIF

    !---- List of public overloaded procedures: subroutines ----!
    public :: Read_File_Cell, Readn_Set_Xtal_Structure, Write_Atoms_CFL, Write_CFL

    !---- List of private functions ----!

    !---- List of private subroutines ----!
    private:: Read_File_Cellc, Read_File_Cellt, Read_File_Atomlist,Read_File_Pointlist,               &
              Readn_Set_Xtal_CFL, Readn_Set_Xtal_CIF, Readn_Set_Xtal_PCR,Readn_Set_Xtal_SHX,          &
              Readn_Set_Xtal_CFL_Molec, Readn_Set_Xtal_Structure_Split,                               &
              Readn_Set_Xtal_Structure_Molcr, Get_NPhases_CIFFile,Get_NPHases_PCRFile,                &
              Write_CFL_Molcrys, Write_CFL_Atom_List_Type, Write_Atoms_CFL_ATM, Write_Atoms_CFL_MOLX, &
              Write_Atoms_CFL_MOLX_orig, Readn_Set_XTal_CFL_Shub

    !---- Definitions ----!


    character (len=1), dimension(26),parameter, private   :: &
    cdd=(/'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r', &
        's','t','u','v','w','x','y','z'/)
    real(kind=dp), parameter, private :: epps=0.000001_dp


    !---- Interfaces - Overloaded procedures--!
    Interface  Read_File_Cell
       Module Procedure Read_File_Cellc  !Last Output Argument Vector Of Six Component With The Cell Parameters
       Module Procedure Read_File_Cellt  !Last output argument object of type Crystal_cell_type
    End interface

    Interface Read_File_Atom
       Module Procedure Read_File_Atomlist   !Last Output Argument of type Atom_list_type
       Module Procedure Read_File_Pointlist  !Last output argument of type Point_list_type
    End Interface

    Interface Readn_Set_Xtal_Structure
       Module Procedure Readn_Set_Xtal_Structure_Molcr ! For Molecular Crystal Type
       Module Procedure Readn_Set_Xtal_Structure_Split ! For Cell, Spg, A types
       Module Procedure Readn_Set_Xtal_Structure_Magn  ! Use Shubnikov groups
    End Interface

    Interface Write_CFL
       Module Procedure Write_CFL_Molcrys        ! For Molecular Crystal Type
       Module Procedure Write_CFL_Atom_List_Type ! For Cell, Spg, A Types
    End Interface

    Interface Write_Atoms_CFL
       Module Procedure Write_Atoms_CFL_MOLX ! For Molecular Crystal Type
       Module Procedure Write_Atoms_CFL_ATM  ! For Cell, Spg, A Types
    End Interface

 Contains

    !---- Functions ----!

    !---- Subroutines ----!





    !!----
    !!---- Subroutine Get_Job_Info(file_dat,i_ini,i_end,Job_info)
    !!----   character(len=*), dimension(:), intent( in) :: file_dat     !Lines of text (content of a file)
    !!----   integer,                        intent( in) :: i_ini,i_end  !Lines to explore
    !!----   type(job_info_type),            intent(out) :: Job_info     !Object to be constructed here
    !!----
    !!----
    !!----    Constructor of the object Job_info. The arrary of strings file_dat
    !!----    have to be provided as input. It contains lines corresponding to the
    !!----    input control file. The analysis of the command lines is not given here.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Get_Job_Info(file_dat,i_ini,i_end,Job_info)
       !---- Arguments ----!
       character(len=*), dimension(:), intent( in) :: file_dat
       integer,                        intent( in) :: i_ini,i_end
       type(job_info_type),            intent(out) :: Job_info

       !---- Local Variables ----!
       integer                           :: i,nphas, ncmd,n_pat,ier, j
       integer, dimension(i_end-i_ini+1) :: ip,ic,ipt
       real(kind=sp)                     :: a1,a2,a3,a4,a5
       character(len=120)                :: line, fmtfields, fmtformat

       !--- Initialize FindFMT
       call Init_FindFMT(i_ini)
       nphas=0
       ncmd=0
       n_pat=0
       ip=i_end
       ic=0
       ipt=0
       Job_info%title=" General Job: CrysFML"
       Job_info%Num_Patterns=1

       do i=i_ini,i_end
          line=u_case(adjustl(file_dat(i)))
          if (line(1:5) == "TITLE") Job_info%title=line(7:)
          if (line(1:5) == "NPATT") then
             read(unit=line(7:), fmt=*,iostat=ier) Job_info%Num_Patterns
             if (ier /= 0) Job_info%Num_Patterns=1
          end if
          if (line(1:6) == "PHASE_") then
             nphas=nphas+1
             ip(nphas)=i
          end if
          if (line(1:4) == "CMDL") then
             ncmd=ncmd+1
             ic(ncmd)=i
          end if
          if (line(1:5) == "PATT_") then
             n_pat=n_pat+1
             ipt(n_pat)=i
          end if
       end do

       if (nphas == 0) then
          nphas=1
          ip(nphas)=0
       end if
       if (n_pat == 0) then
          n_pat=1
          ipt(n_pat) = 0
       end if

       if (Job_info%Num_Patterns /= n_pat) Job_info%Num_Patterns = n_pat
       Job_info%Num_Phases=nphas
       Job_info%Num_Cmd=ncmd

       if (allocated(Job_Info%Patt_typ)) deallocate(Job_Info%Patt_typ)
       allocate(Job_Info%Patt_typ(n_pat))

       if (allocated(Job_Info%Phas_nam)) deallocate(Job_Info%Phas_nam)
       allocate(Job_Info%Phas_nam(nphas))

       if (allocated(Job_Info%range_stl)) deallocate(Job_Info%range_stl)
       allocate(Job_Info%range_stl(n_pat))

       if (allocated(Job_Info%range_q)) deallocate(Job_Info%range_q)
       allocate(Job_Info%range_q(n_pat))

       if (allocated(Job_Info%range_d)) deallocate(Job_Info%range_d)
       allocate(Job_Info%range_d(n_pat))

       if (allocated(Job_Info%range_2theta)) deallocate(Job_Info%range_2theta)
       allocate(Job_Info%range_2theta(n_pat))

       if (allocated(Job_Info%range_energy)) deallocate(Job_Info%range_energy)
       allocate(Job_Info%range_energy(n_pat))

       if (allocated(Job_Info%range_tof)) deallocate(Job_Info%range_tof)
       allocate(Job_Info%range_tof(n_pat))

       if (allocated(Job_Info%lambda)) deallocate(Job_Info%lambda)
       allocate(Job_Info%lambda(n_pat))

       if (allocated(Job_Info%ratio)) deallocate(Job_Info%ratio)
       allocate(Job_Info%ratio(n_pat))

       if (allocated(Job_Info%dtt1)) deallocate(Job_Info%dtt1)
       allocate(Job_Info%dtt1(n_pat))

       if (allocated(Job_Info%dtt2)) deallocate(Job_Info%dtt2)
       allocate(Job_Info%dtt2(n_pat))

       !---- Initialize all variables
       Job_Info%Patt_typ    =" "
       Job_Info%Phas_nam    =" "
       Job_Info%range_stl%mina=0.0
       Job_Info%range_stl%maxb=0.0
       Job_Info%range_q%mina=0.0
       Job_Info%range_q%maxb=0.0
       Job_Info%range_d%mina=0.0
       Job_Info%range_d%maxb=0.0
       Job_Info%range_2theta%mina=0.0
       Job_Info%range_2theta%maxb=0.0
       Job_Info%range_Energy%mina=0.0
       Job_Info%range_Energy%maxb=0.0
       Job_Info%range_tof%mina=0.0
       Job_Info%range_tof%maxb=0.0
       Job_Info%Lambda%mina=0.0
       Job_Info%Lambda%maxb=0.0
       Job_Info%ratio = 0.0
       Job_Info%dtt1 = 0.0
       Job_Info%dtt2 = 0.0
       if (ncmd > 0) then
          if (allocated(Job_Info%cmd)) deallocate(Job_Info%cmd)
          allocate(Job_Info%cmd(ncmd))
          Job_Info%cmd=" "
       end if

       !---- Fill the different fields of Job_Info
       !---- Start with patterns
       fmtfields = "9fffff"

       !---- First asks if there is a PATT_ card, if not a standard is taken
       if (ipt(1) /= 0) then
          do n_pat=1, Job_info%Num_Patterns
             i=ipt(n_pat)
             line=u_case(adjustl(file_dat(i)))
             line=line(8:)
             call findfmt(0,line,fmtfields,fmtformat)
             read(unit=line,fmt=fmtformat) Job_Info%Patt_typ(n_pat), a1,a2,a3,a4,a5
             if (ierr_fmt /= 0) return
             line=u_case(Job_Info%Patt_typ(n_pat))

             select case(line(1:9))
                case("XRAY_2THE","NEUT_2THE","XRAY_SXTA","NEUT_SXTA")
                   if ( a1 <= 0.000001) a1=1.5405
                   if ( a2 <= 0.000001) then
                      a2=a1
                      a3=0.0
                   end if
                   if (a5 <= a4) a5=120.0
                   Job_Info%Lambda(n_pat)%mina=a1
                   Job_Info%Lambda(n_pat)%maxb=a2
                   Job_Info%ratio(n_pat)=a3
                   Job_Info%range_2theta(n_pat)%mina=a4
                   Job_Info%range_2theta(n_pat)%maxb=a5
                   a4=sind(0.5*a4)/a1
                   a5=sind(0.5*a5)/a2
                   Job_Info%range_stl(n_pat)%mina=a4
                   Job_Info%range_stl(n_pat)%maxb=a5
                   Job_Info%range_q(n_pat)%mina=a4*4.0*pi
                   Job_Info%range_q(n_pat)%maxb=a5*4.0*pi
                   Job_Info%range_d(n_pat)%mina=0.5/a5
                   Job_Info%range_d(n_pat)%maxb=0.5/a4

                case("NEUT_TOF ")
                   if (a1 <= 0.000001) a1=1000.0
                   if (a4 <= a3) a4=2.0*abs(a3)
                   Job_Info%dtt1(n_pat)=a1
                   Job_Info%dtt2(n_pat)=a2
                   Job_Info%range_tof(n_pat)%mina=a3
                   Job_Info%range_tof(n_pat)%maxb=a4
                   Job_Info%range_d(n_pat)%mina=0.5*(-1.0+sqrt(1.0+4.0*a2*a3/a1/a1))
                   Job_Info%range_d(n_pat)%maxb=0.5*(-1.0+sqrt(1.0+4.0*a2*a4/a1/a1))
                   Job_Info%range_stl(n_pat)%mina=0.5/Job_Info%range_d(n_pat)%maxb
                   Job_Info%range_stl(n_pat)%maxb=0.5/Job_Info%range_d(n_pat)%mina
                   Job_Info%range_q(n_pat)%mina=Job_Info%range_stl(n_pat)%mina*4.0*pi
                   Job_Info%range_q(n_pat)%maxb=Job_Info%range_stl(n_pat)%maxb*4.0*pi

                case("XRAY_ENER")
                   if (a1 <= 0.000001) a1=12.4 !(=hc(keV.Angstr.)
                   Job_Info%dtt1(n_pat)=a1
                   Job_Info%dtt2(n_pat)=0.0
                   Job_Info%range_energy(n_pat)%mina=a3
                   Job_Info%range_energy(n_pat)%maxb=a4
                   if (a3 <= 0.00001) a3=0.01
                   if (a4 <= 0.00001) a4=2.00
                   Job_Info%range_d(n_pat)%mina=a1/a4
                   Job_Info%range_d(n_pat)%maxb=a1/a3
                   Job_Info%range_stl(n_pat)%mina=0.5/Job_Info%range_d(n_pat)%maxb
                   Job_Info%range_stl(n_pat)%maxb=0.5/Job_Info%range_d(n_pat)%mina
                   Job_Info%range_q(n_pat)%mina=Job_Info%range_stl(n_pat)%mina*4.0*pi
                   Job_Info%range_q(n_pat)%maxb=Job_Info%range_stl(n_pat)%maxb*4.0*pi

             end select
          end do

       else
          n_pat=1
          a1=1.5405
          a2=a1
          a3=0.0
          a4=0.0
          a5=120.0
          Job_Info%Patt_typ(n_pat)="XRAY_2THE"
          Job_Info%Lambda(n_pat)%mina=a1
          Job_Info%Lambda(n_pat)%maxb=a2
          Job_Info%ratio(n_pat)=a3
          Job_Info%range_2theta(n_pat)%mina=a4
          Job_Info%range_2theta(n_pat)%maxb=a5
          a4=sind(0.5*a4)/a1
          a5=sind(0.5*a5)/a2
          Job_Info%range_stl(n_pat)%mina=a4
          Job_Info%range_stl(n_pat)%maxb=a5
          Job_Info%range_q(n_pat)%mina=a4*4.0*pi
          Job_Info%range_q(n_pat)%maxb=a5*4.0*pi
          Job_Info%range_d(n_pat)%mina=0.5/a5
          Job_Info%range_d(n_pat)%maxb=0.5/a4
       end if

       !---- Phase names
       if (ip(1) /= 0) then
          do i=1,nphas
             j=ip(i)
             line=adjustl(file_dat(j))
             Job_Info%Phas_nam(i)=line(8:)
          end do
       else
          Job_Info%Phas_nam(1)= Job_info%title
       end if

       !---- Command Lines, stored but not analysed here
       do i=1,ncmd
          j=ic(i)
          line=adjustl(file_dat(j))
          Job_Info%cmd(i)=line(8:)
       end do

       return
    End Subroutine Get_Job_Info




    !!----
    !!---- Subroutine Readn_Set_Magnetic_Structure_MCIF(file_mcif,mCell,MGp,Am)
    !!----    character(len=*),               intent (in)  :: file_mcif
    !!----    type(Crystal_Cell_type),        intent (out) :: mCell
    !!----    type(Magnetic_Space_Group_Type),intent (out) :: MGp
    !!----    type(Atom_List_Type),           intent (out) :: Am
    !!----
    !!----    Subroutine for reading and construct a magnetic structure.
    !!----    The atom list and the unit cell reading an mCIF file.
    !!----
    !!----  Created: January-2014 (JRC)
    !!----  Updated: August-2014 (JRC), January 2020
    !!
    Subroutine Readn_Set_Magnetic_Structure_MCIF(file_mcif,mCell,MGp,Am)
       character(len=*),               intent (in)  :: file_mcif
       type(Crystal_Cell_type),        intent (out) :: mCell
       type(Magnetic_Space_Group_Type),intent (out) :: MGp
       type(Atom_List_Type),           intent (out) :: Am

       !---- Local Variables ----!
       integer :: i,num_sym, num_constr, num_kvs,num_matom, num_mom, num_magscat, ier, j, m, n, k, L,   &
                  ncar,mult,nitems,iv, num_irreps, nitems_irreps, num_rsym, num_centering,det,kfin
       integer,          dimension(10)     :: lugar
       integer,          dimension(7)      :: irrep_pos
       integer,          dimension(5)      :: pos
       integer,          dimension(3,3)    :: Rot
       real(kind=cp),    dimension(3)      :: cel,ang,cel_std,ang_std,tr,v
       real(kind=cp),    dimension(6)      :: values,std
       real(kind=cp),    dimension(3,3)    :: matr
       real(kind=cp),    dimension(3,384)  :: orb
       character(len=132)                  :: lowline,keyword,line, mxmymz_op,linat
       character(len=132),dimension(384)   :: sym_strings, cent_strings
       character(len=132),dimension(384)   :: atm_strings
       character(len=132),dimension(384)   :: mom_strings
       character(len=132),dimension(30)    :: constr_strings, mag_scatt_string
       character(len=132),dimension(30)    :: irreps_strings
       character(len=132),dimension(30)    :: kv_strings
       character(len=20), dimension(15)    :: lab_items
       character(len=50)                   :: shubk
       character(len=2)                    :: chars
       character(len=10)                   :: label
       character(len=4)                    :: symbcar
       logical                             :: ktag,no_symop_mxmymz,no_cent_mxmymz,mom_symmform

       !type(Magnetic_Group_Type)  :: SG
       type(file_list_type)       :: mcif

       call init_err_Form()
       call File_To_FileList(file_mcif,mcif)
       !Remove all possible tabs and non-ASCII characters in the CIF
       do i=1,mcif%nlines
         do j=1,len_trim(mcif%line(i))
           if(mcif%line(i)(j:j) == char(9)) mcif%line(i)(j:j)=" "
         end do
       end do
       num_constr=0; num_kvs=0; num_matom=0; num_mom=0; num_sym=0; num_magscat=0; num_rsym=0; num_centering=0
       cel=0.0; ang=0.0; num_irreps=0; nitems_irreps=0
       i=0
       call Init_Magnetic_Space_Group_Type(MGp)
       ktag=.false.
       no_symop_mxmymz=.false.
       no_cent_mxmymz=.false.
       mom_symmform=.false.

       do
          i=i+1
          if(i > mcif%nlines) exit
          if (index(mcif%line(i)(1:1),"!")/=0 .or. index(mcif%line(i)(1:1),"#")/=0 .or. len_trim(mcif%line(i)) == 0) cycle
          line=adjustl(mcif%line(i))
          lowline=l_case(line)
          j=index(lowline," ")
          keyword=lowline(1:j-1)
          !write(*,"(a)") " Keyword: "//trim(keyword)

          Select Case (trim(keyword))

             Case("_magnetic_space_group_standard_setting","_magnetic_space_group.standard_setting")
                chars=adjustl(line(j+1:))
                if(chars(2:2) == "y" .or. chars(2:2) == "Y") MGp%standard_setting=.true.
                !write(unit=*,fmt="(a)") "  Treating item: _magnetic_space_group_standard_setting -> "//trim(chars)

             Case("_parent_space_group.name_h-m", "_parent_space_group_name_h-m","_parent_space_group.name_h-m_alt")
                shubk=adjustl(line(j+1:))
                m=len_trim(shubk)
                MGp%Parent_spg=shubk(2:m-1)
                !write(unit=*,fmt="(a)") "  Treating item: _parent_space_group_name_h-m -> "// MGp%Parent_spg

             Case("_parent_space_group.it_number","_parent_space_group_it_number")
                read(unit=lowline(j:),fmt=*,iostat=ier) m
                if(ier /= 0) then
                  Err_Form=.true.
                  Err_Form_Mess=" Error reading the number of the parent space group"
                  return
                end if
                MGp%Parent_num=m
                !write(unit=*,fmt="(a,i4)") "  Treating item: _parent_space_group_it_number -> ", MGp%Parent_num

             Case("_magnetic_space_group_bns_number","_space_group.magn_number_bns","_space_group_magn.number_bns")
                shubk=adjustl(line(j+1:))
                k=len_trim(shubk)
                if(shubk(1:1) == '"' .or. shubk(1:1) == "'") shubk=adjustl(shubk(2:k-1))
                MGp%BNS_number=shubk
                !write(unit=*,fmt="(a)") "  Treating item: _space_group.magn_number_bns -> "//trim(MGp%BNS_number)

             Case("_magnetic_space_group_bns_name","_space_group_magn.name_bns","_space_group.magn_name_bns")
                shubk=adjustl(line(j+1:))
                k=len_trim(shubk)
                if(shubk(1:1) == '"' .or. shubk(1:1) == "'") shubk=adjustl(shubk(2:k-1))
                MGp%BNS_symbol=pack_string(shubk)
                !write(unit=*,fmt="(a)") "  Treating item: _space_group.magn_name_bns -> "//trim(MGp%BNS_symbol)

             Case("_magnetic_space_group_og_number","_space_group_magn.number_og","_space_group.magn_number_og")
                shubk=adjustl(line(j+1:))
                k=len_trim(shubk)
                if(shubk(1:1) == '"' .or. shubk(1:1) == "'") shubk=adjustl(shubk(2:k-1))
                MGp%OG_number=shubk
                !write(unit=*,fmt="(a)") "  Treating item: _space_group.magn_number_og -> "//trim(MGp%OG_number)

             Case("_magnetic_space_group_point_group","_space_group_magn.point_group","_space_group.magn_point_group")
                shubk=adjustl(line(j+1:))
                k=len_trim(shubk)
                if(shubk(1:1) == '"' .or. shubk(1:1) == "'") shubk=adjustl(shubk(2:k-1))
                MGp%PG_symbol=pack_string(shubk)
                !write(unit=*,fmt="(a)") "  Treating item: _space_group_magn.point_group -> "//trim(MGp%PG_symbol)

             Case("_magnetic_space_group_og_name","_space_group_magn.name_og","_space_group.magn_name_og")
                shubk=adjustl(line(j+1:))
                k=len_trim(shubk)
                if(shubk(1:1) == '"' .or. shubk(1:1) == "'") shubk=adjustl(shubk(2:k-1))
                MGp%OG_symbol=pack_string(shubk)
                !write(unit=*,fmt="(a)") "  Treating item: _space_group.magn_name_og -> "//trim(MGp%OG_symbol)

             Case("_magnetic_space_group.transform_from_parent_pp_abc","_magnetic_space_group_transform_from_parent_pp_abc", &
                   "_parent_space_group.child_transform_pp_abc")
                shubk=adjustl(line(j+1:))
                k=len_trim(shubk)
                if(shubk(1:1) == '"' .or. shubk(1:1) == "'") shubk=adjustl(shubk(2:k-1))
                MGp%trn_from_parent=pack_string(shubk)
                !write(unit=*,fmt="(a)") "  Treating item: _magnetic_space_group_transform_from_parent_pp_abc -> "//trim(MGp%trn_from_parent)

             Case("_parent_space_group.transform_pp_abc")
                shubk=adjustl(line(j+1:))
                k=len_trim(shubk)
                if(shubk(1:1) == '"' .or. shubk(1:1) == "'") shubk=adjustl(shubk(2:k-1))
                MGp%trn_to_parent=pack_string(shubk)
                !write(unit=*,fmt="(a)") "  Treating item: _magnetic_space_group_transform_from_parent_pp_abc -> "//trim(MGp%trn_from_parent)

             Case("_magnetic_space_group.transform_to_standard_pp_abc","_magnetic_space_group_transform_to_standard_pp_abc")
                shubk=adjustl(line(j+1:))
                k=len_trim(shubk)
                if(shubk(1:1) == '"' .or. shubk(1:1) == "'") shubk=adjustl(shubk(2:k-1))
                MGp%trn_to_standard=pack_string(shubk)

             Case("_space_group_magn.transform_bns_pp_abc")
                shubk=adjustl(line(j+1:))
                k=len_trim(shubk)
                if(shubk(1:1) == '"' .or. shubk(1:1) == "'") shubk=adjustl(shubk(2:k-1))
                MGp%trn_to_standard=pack_string(shubk)

             Case("_magnetic_cell_length_a","_cell_length_a")
                call getnum_std(lowline(j:),values,std,iv)
                if(err_string) then
                  Err_Form=.true.
                  Err_Form_Mess=" Error reading the magnetic unit cell parameter 'a' -> "//trim(err_string_mess)
                  return
                end if
                cel(1)=values(1)
                cel_std(1)=std(1)
                MGp%m_cell=.true.
                !write(unit=*,fmt="(a)") "  Treating item: _cell_length_a"

             Case("_magnetic_cell_length_b","_cell_length_b")
                call getnum_std(lowline(j:),values,std,iv)
                if(err_string) then
                  Err_Form=.true.
                  Err_Form_Mess=" Error reading the magnetic unit cell parameter 'b' -> "//trim(err_string_mess)
                  return
                end if
                cel(2)=values(1)
                cel_std(2)=std(1)
                !write(unit=*,fmt="(a)") "  Treating item: _cell_length_b"

             Case("_magnetic_cell_length_c","_cell_length_c")
                call getnum_std(lowline(j:),values,std,iv)
                if(err_string) then
                  Err_Form=.true.
                  Err_Form_Mess=" Error reading the magnetic unit cell parameter 'c' -> "//trim(err_string_mess)
                  return
                end if
                cel(3)=values(1)
                cel_std(3)=std(1)
                !write(unit=*,fmt="(a)") "  Treating item: _cell_length_c"

             Case("_magnetic_cell_angle_alpha","_cell_angle_alpha")
                call getnum_std(lowline(j:),values,std,iv)
                if(err_string) then
                  Err_Form=.true.
                  Err_Form_Mess=" Error reading the magnetic unit cell parameter 'alpha' -> "//trim(err_string_mess)
                  return
                end if
                ang(1)=values(1)
                ang_std(1)=std(1)
                !write(unit=*,fmt="(a)") "  Treating item: _cell_angle_alpha"

             Case("_magnetic_cell_angle_beta","_cell_angle_beta")
                call getnum_std(lowline(j:),values,std,iv)
                if(err_string) then
                  Err_Form=.true.
                  Err_Form_Mess=" Error reading the magnetic unit cell parameter 'beta' -> "//trim(err_string_mess)
                  return
                end if
                ang(2)=values(1)
                ang_std(2)=std(1)
                !write(unit=*,fmt="(a)") "  Treating item: _cell_angle_beta"

             Case("_magnetic_cell_angle_gamma","_cell_angle_gamma")
                call getnum_std(lowline(j:),values,std,iv)
                if(err_string) then
                  Err_Form=.true.
                  Err_Form_Mess=" Error reading the magnetic unit cell parameter 'gamma' -> "//trim(err_string_mess)
                  return
                end if
                ang(3)=values(1)
                ang_std(3)=std(1)
                !write(unit=*,fmt="(a)") "  Treating item: _cell_angle_gamma"

             Case("loop_")
                 i=i+1
                 line=adjustl(mcif%line(i))
                 lowline=l_case(line)
                 j=index(lowline," ")
                 keyword=lowline(1:j-1)
                 !write(*,"(a)") "         Loop_Keyword: "//trim(keyword)
                 Select Case(trim(keyword))

                   Case("_space_group_magn_transforms.id")
                      !write(*,"(a)") "         Loop_Keyword: "//trim(keyword)

                      do k=1,2
                        i=i+1
                        if(index(mcif%line(i),"_space_group_magn_transforms") == 0) then
                          Err_Form=.true.
                          Err_Form_Mess=" Error reading _space_group_magn_transforms in loop"
                          return
                        end if
                      end do
                      i=i+1
                      call getword(mcif%line(i),lab_items,iv)
                      !write(unit=*,fmt="(3a)")  (lab_items(k),k=1,3)
                      if(lab_items(3)(1:3) == "BNS") then
                        MGp%trn_to_standard=lab_items(2)
                      end if
                      i=i+1
                      call getword(mcif%line(i),lab_items,iv)
                      !write(unit=*,fmt="(3a)")  (lab_items(k),k=1,3)
                      if(lab_items(3)(1:2) == "OG") then
                        !nothing to do
                      end if

                   Case("_irrep_id")
                      irrep_pos=0
                      irrep_pos(1)=1
                      j=1
                      do k=1,6
                         i=i+1
                         if(index(mcif%line(i),"_irrep_dimension") /= 0) then
                            j=j+1
                            irrep_pos(2)=j
                            cycle
                         end if
                         if(index(mcif%line(i),"_small_irrep_dimension") /= 0 .or.  &
                            index(mcif%line(i),"_irrep_small_dimension") /= 0) then
                            j=j+1
                            irrep_pos(3)=j
                            cycle
                         end if
                         if(index(mcif%line(i),"_irrep_direction_type") /= 0) then
                            j=j+1
                            irrep_pos(4)=j
                            cycle
                         end if
                         if(index(mcif%line(i),"_irrep_action") /= 0) then
                            j=j+1
                            irrep_pos(5)=j
                            cycle
                         end if
                         if(index(mcif%line(i),"_irrep_modes_number") /= 0) then
                            j=j+1
                            irrep_pos(6)=j
                            cycle
                         end if
                         if(index(mcif%line(i),"_irrep_presence") /= 0) then
                            j=j+1
                            irrep_pos(7)=j
                            cycle
                         end if
                         exit
                      end do

                      i=i-1
                      nitems_irreps=count(irrep_pos > 0)

                      k=0
                      do
                        i=i+1
                        if(i > mcif%nlines) exit
                        if(len_trim(mcif%line(i)) == 0) exit
                        k=k+1
                        irreps_strings(k)=mcif%line(i)
                      end do
                      num_irreps=k
                      !Treat later the list of irreps

                   Case("_magnetic_propagation_vector_seq_id")
                      do k=1,3
                        i=i+1
                        if(index(mcif%line(i),"_magnetic_propagation_vector") == 0) then
                          Err_Form=.true.
                          Err_Form_Mess=" Error reading the propagation vector loop"
                          return
                        end if
                        if(index(mcif%line(i),"_magnetic_propagation_vector_kxkykz") /= 0) then
                          ktag=.true.  !new format for k-vector klabel '0,1/2,0'
                          exit
                        end if
                      end do
                      k=0
                      do
                        i=i+1
                        if(len_trim(mcif%line(i)) == 0) exit
                        k=k+1
                        kv_strings(k)=mcif%line(i)
                      end do
                      num_kvs=k
                      MGp%n_kv=k
                      if(allocated(Mgp%kv)) deallocate(Mgp%kv)
                      allocate(Mgp%kv(3,k))
                      if(allocated(Mgp%kv_label)) deallocate(Mgp%kv_label)
                      allocate(Mgp%kv_label(k))
                      !Treat later the propagation vectors

                   Case("_atom_type_symbol")
                      !write(unit=*,fmt="(a)") "  Treating item: _atom_type_symbol"
                      do k=1,3
                        i=i+1
                        if(index(mcif%line(i),"_atom_type_symbol") == 0) then
                          Err_Form=.true.
                          Err_Form_Mess=" Error reading the _atom_type_symbol in loop"
                          return
                        end if
                      end do
                      k=0
                      do
                        i=i+1
                        if(len_trim(mcif%line(i)) == 0) exit
                        k=k+1
                        mag_scatt_string(k)=mcif%line(i)
                      end do
                      num_magscat=k
                      !Treat later the scattering factor

                   Case("_magnetic_atom_site_moment_symmetry_constraints_label")
                      !write(unit=*,fmt="(a)") "  Treating item: _magnetic_atom_site_moment_symmetry_constraints_label"
                      i=i+1
                      if(index(mcif%line(i),"_atom_site_magnetic_moment_symmetry_constraints_mxmymz") == 0) then
                        Err_Form=.true.
                        Err_Form_Mess=" Error reading the magnetic_atom_site_moment_symmetry_constraints loop"
                        return
                      end if
                      k=0
                      do
                        i=i+1
                        if(len_trim(mcif%line(i)) == 0) exit
                        k=k+1
                        constr_strings(k)=mcif%line(i)
                      end do
                      num_constr=k
                      MGp%m_constr=.true.
                      !Treat later the constraints

                   Case("_magnetic_space_group_symop_id")
                      !write(unit=*,fmt="(a)") "  Treating item: _space_group_symop_magn_operation.id"
                      do k=1,3
                        i=i+1
                        j=index(mcif%line(i),"_magnetic_space_group_symop_operation")
                        if(j == 0 ) then
                          Err_Form=.true.
                          Err_Form_Mess=" Error reading the _magnetic_space_group_symop_operation loop"
                          return
                        end if
                      end do
                      k=0
                      do
                        i=i+1
                        if(len_trim(mcif%line(i)) == 0) exit
                        k=k+1
                        sym_strings(k)=mcif%line(i)
                      end do
                      !now allocate the list of symmetry operators
                      num_sym=k
                      MGp%Multip=k

                   Case("_space_group_symop_magn_operation.id","_space_group_symop_magn.id") !The second item is added to be compatible with BCS error
                      !write(unit=*,fmt="(a)") "  Treating item: _space_group_symop_magn_operation.id"

                      i=i+1
                      j=index(mcif%line(i),"_space_group_symop_magn_operation.xyz")
                      if(j == 0 ) then
                        Err_Form=.true.
                        Err_Form_Mess=" Error reading the _space_group_symop_magn_operation loop"
                        return
                      end if

                      k=0
                      do
                        i=i+1
                        if(len_trim(mcif%line(i)) == 0) exit
                        k=k+1
                        sym_strings(k)=mcif%line(i)
                      end do
                      !now allocate the list of symmetry operators
                      num_sym=k
                      MGp%Multip=k

                   Case("_space_group_symop.magn_id")
                      !write(unit=*,fmt="(a)") "  Treating item: _space_group_symop_magn_id"
                      do k=1,2
                        i=i+1
                        if(index(mcif%line(i),"_space_group_symop.magn_operation") == 0 .and. &
                           index(mcif%line(i),"_space_group_symop_magn_operation") == 0) then
                          Err_Form=.true.
                          Err_Form_Mess=" Error reading the _space_group_symop_magn_operation loop"
                          return
                        end if
                      end do
                      if(index(mcif%line(i),"_space_group_symop.magn_operation_mxmymz") == 0 .and. &
                         index(mcif%line(i),"_space_group_symop_magn_operation_mxmymz") == 0) then
                         i=i-1
                         no_symop_mxmymz=.true.
                      end if
                      k=0
                      do
                        i=i+1
                        if(len_trim(mcif%line(i)) == 0) exit
                        k=k+1
                        sym_strings(k)=mcif%line(i)
                      end do

                      num_rsym=k

                   Case("_space_group_symop_magn_id")   !here the symmetry operators are separated from the translations
                      !write(unit=*,fmt="(a)") "  Treating item: _space_group_symop_magn_id"
                      i=i+1
                      if(index(mcif%line(i),"_space_group_symop.magn_operation") == 0 .and. &
                         index(mcif%line(i),"_space_group_symop_magn_operation") == 0) then
                        Err_Form=.true.
                        Err_Form_Mess=" Error reading the _space_group_symop.magn_operation loop"
                        return
                      end if
                      if(index(mcif%line(i),"_space_group_symop.magn_operation_mxmymz") == 0 .and. &
                         index(mcif%line(i),"_space_group_symop_magn_operation_mxmymz") == 0) then
                         no_symop_mxmymz=.true.
                      end if
                      k=0
                      do
                        i=i+1
                        if(len_trim(mcif%line(i)) == 0) exit
                        k=k+1
                        sym_strings(k)=mcif%line(i)
                      end do

                      num_rsym=k

                   Case("_space_group_symop.magn_centering_id")   !here we read the translations and anti-translations
                      !write(unit=*,fmt="(a)") "  Treating item: _space_group_symop_magn_centering_id"
                      do k=1,2
                        i=i+1
                        if(index(mcif%line(i),"_space_group_symop.magn_centering") == 0 .and. &
                           index(mcif%line(i),"_space_group_symop_magn_centering") == 0 ) then
                          Err_Form=.true.
                          Err_Form_Mess=" Error reading the _space_group_symop_magn_centering loop"
                          return
                        end if
                      end do
                      if(index(mcif%line(i),"_space_group_symop.magn_centering_mxmymz") == 0 .and. &
                         index(mcif%line(i),"_space_group_symop_magn_centering_mxmymz") == 0 ) then
                         i=i-1
                         no_cent_mxmymz=.true.
                      end if
                      k=0
                      do
                        i=i+1
                        if(len_trim(mcif%line(i)) == 0) exit
                        k=k+1
                        cent_strings(k)=mcif%line(i)
                      end do
                      num_centering=k

                   Case("_space_group_symop_magn_centering.id")   !here we read the translations and anti-translations
                      !write(unit=*,fmt="(a)") "  Treating item: _space_group_symop_magn_centering_id"
                      i=i+1
                      if(index(mcif%line(i),"_space_group_symop_magn_centering.xyz") == 0) then
                        Err_Form=.true.
                        Err_Form_Mess=" Error reading the _space_group_symop_magn_centering.xyz loop"
                        return
                      end if
                      k=0
                      do
                        i=i+1
                        if(len_trim(mcif%line(i)) == 0) exit
                        k=k+1
                        cent_strings(k)=mcif%line(i)
                      end do
                      num_centering=k

                   Case("_magnetic_atom_site_label","_atom_site_label")
                      !write(unit=*,fmt="(a)") "  Treating item: _atom_site_label"
                      !Count the number of keywords following the _loop
                      do k=1,10
                        linat=adjustl(mcif%line(i+k))
                        if(linat(1:1) /=  "_") then
                          kfin=k+1
                          iv=i+k
                          exit
                        end if
                      end do
                      lugar=0
                      lugar(1)=1
                      j=1
                      do k=1,kfin
                         i=i+1
                         if(index(mcif%line(i),"_atom_site_type_symbol") /= 0) then
                            j=j+1
                            lugar(2)=j
                            cycle
                         end if
                         if(index(mcif%line(i),"_atom_site_fract_x") /= 0) then
                            j=j+1
                            lugar(3)=j
                            cycle
                         end if
                         if(index(mcif%line(i),"_atom_site_fract_y") /= 0) then
                            j=j+1
                            lugar(4)=j
                            cycle
                         end if
                         if(index(mcif%line(i),"_atom_site_fract_z") /= 0) then
                            j=j+1
                            lugar(5)=j
                            cycle
                         end if
                         if (index(mcif%line(i),"_atom_site_U_iso_or_equiv") /= 0) then
                            j=j+1
                            lugar(6)=j
                            cycle
                         end if
                         if (index(mcif%line(i),"_atom_site_B_iso_or_equiv") /= 0) then
                            j=j+1
                            lugar(10)=j
                            cycle
                         end if
                         if (index(mcif%line(i),"_atom_site_occupancy") /= 0) then
                            j=j+1
                            lugar(7)=j
                            cycle
                         end if
                         if (index(mcif%line(i),"_atom_site_symmetry_multiplicity") /= 0) then
                            j=j+1
                            lugar(8)=j
                            cycle
                         end if
                         if (index(mcif%line(i),"_atom_site_Wyckoff_label") /= 0) then
                            j=j+1
                            lugar(9)=j
                            cycle
                         end if
                         exit
                      end do

                      if (any(lugar(3:5) == 0)) then
                          Err_Form=.true.
                          Err_Form_Mess=" Error reading the asymmetric unit of magnetic atoms"
                          return
                      end if

                      i=iv-1
                      nitems=count(lugar > 0)

                      k=0
                      do
                        i=i+1
                        if(i > mcif%nlines) exit
                        if(len_trim(mcif%line(i)) == 0) exit
                        k=k+1
                        atm_strings(k)=adjustl(mcif%line(i))
                      end do
                      num_matom=k
                      !Treat late the list atoms

                   Case("_magnetic_atom_site_moment_label","_atom_site_moment_label","_atom_site_moment.label")
                      !write(unit=*,fmt="(a)") "  Treating item: _atom_site_moment_label"
                      do k=1,3
                        i=i+1
                        if(index(mcif%line(i),"_atom_site_moment_crystalaxis") == 0 .and. &
                           index(mcif%line(i),"_atom_site_moment.crystalaxis") == 0) then
                          Err_Form=.true.
                          Err_Form_Mess=" Error reading the magnetic_atom_site_moment loop"
                          return
                        end if
                      end do
                      i=i+1
                      if(index(mcif%line(i),"_atom_site_moment.symmform") /= 0) then
                        !write(*,*) " _atom_site_moment.symmform FOUND"
                        mom_symmform=.true.
                      else
                        i=i-1
                      end if
                      k=0
                      do
                        i=i+1
                        if(i > mcif%nlines) exit
                        if(len_trim(mcif%line(i)) == 0) exit
                        k=k+1
                        mom_strings(k)=mcif%line(i)
                      end do
                      num_mom=k
                      !Treat later the magnetic moment of the atoms
                 End Select
          End Select
       end do

       if(MGp%m_cell) then
         call Set_Crystal_Cell(cel,ang,mCell)
         mCell%cell_std=cel_std
         mCell%ang_std=ang_std
       end if

       !Treat symmetry operators
       !write(unit=*,fmt="(a,2i4)") " num_sym, num_rsym :",num_sym,num_rsym
       if(num_sym == 0 .and. num_rsym == 0) then
          Err_Form=.true.
          Err_Form_Mess=" No symmetry operators have been provided in the MCIF file "//trim(file_mcif)
          return
       else
          if(no_cent_mxmymz) then  !Full number of symmetry operators is not separated from the centering

            if(allocated(Mgp%SymopSymb)) deallocate(Mgp%SymopSymb)
            allocate(Mgp%SymopSymb(num_sym))
            if(allocated(Mgp%Symop)) deallocate(Mgp%Symop)
            allocate(Mgp%Symop(num_sym))
            if(allocated(Mgp%MSymopSymb)) deallocate(Mgp%MSymopSymb)
            allocate(Mgp%MSymopSymb(num_sym))
            if(allocated(Mgp%MSymop)) deallocate(Mgp%MSymop)
            allocate(Mgp%MSymop(num_sym))
            !write(unit=*,fmt="(a)") "  Decoding symmetry operators 1"

            ! Decode the symmetry operators
            do i=1,num_sym
              line=adjustl(sym_strings(i))
              j=index(line," ")
              line=adjustl(line(j+1:))
              j=index(line," ")
              MGp%SymopSymb(i)=line(1:j-1)
              line=adjustl(line(j+1:))
              j=index(line," ")
              MGp%MSymopSymb(i)=line(1:j-1)
              read(unit=line(j:),fmt=*,iostat=ier) n
              if(ier /= 0) then
                 Err_Form=.true.
                 Err_Form_Mess=" Error reading the time inversion in line: "//trim(sym_strings(i))
                 return
              else
                 MGp%MSymOp(i)%phas=real(n)
              end if
              call Read_Xsym(MGp%SymopSymb(i),1,MGp%Symop(i)%Rot,MGp%Symop(i)%tr)
              line=MGp%MSymopSymb(i)
              do k=1,len_trim(line)
                if(line(k:k) == "m") line(k:k)=" "
              end do
              line=Pack_String(line)
              call Read_Xsym(line,1,MGp%MSymop(i)%Rot)
            end do

          else

            if( num_rsym == 0) num_rsym=num_sym
            ! First allocate the full number of symmetry operators after decoding if centering lattice
            ! have been provided and if the group is centred or not
            if(num_centering == 0) then
               MGp%Multip=num_rsym
            else
               MGp%Multip=num_rsym*num_centering
            end if

            num_sym=MGp%Multip
            if(allocated(Mgp%SymopSymb)) deallocate(Mgp%SymopSymb)
            allocate(Mgp%SymopSymb(num_sym))
            if(allocated(Mgp%Symop)) deallocate(Mgp%Symop)
            allocate(Mgp%Symop(num_sym))
            if(allocated(Mgp%MSymopSymb)) deallocate(Mgp%MSymopSymb)
            allocate(Mgp%MSymopSymb(num_sym))
            if(allocated(Mgp%MSymop)) deallocate(Mgp%MSymop)
            allocate(Mgp%MSymop(num_sym))
            ! Decode the symmetry operators
            !write(unit=*,fmt="(a)") "  Decoding symmetry operators 2"
            do i=1,num_rsym
              line=adjustl(sym_strings(i))
              j=index(line," ")
              line=adjustl(line(j+1:))
              j=index(line," ")
              MGp%SymopSymb(i)=line(1:j-1)
              k=index(MGp%SymopSymb(i),",",back=.true.)
              read(unit=MGp%SymopSymb(i)(k+1:),fmt=*,iostat=ier) n
              if(ier /= 0) then
                 Err_Form=.true.
                 Err_Form_Mess=" Error reading the time inversion in line: "//trim(sym_strings(i))
                 return
              else
                 MGp%MSymOp(i)%phas=real(n)
              end if
              MGp%SymopSymb(i)=MGp%SymopSymb(i)(1:k-1)
              call Read_Xsym(MGp%SymopSymb(i),1,MGp%Symop(i)%Rot,MGp%Symop(i)%tr)

              !Now construc the magnetic rotation symbols
              line=adjustl(line(j+1:))
              if(len_trim(line) /= 0) then
                j=index(line," ")
                MGp%MSymopSymb(i)=line(1:j-1)
                line=MGp%MSymopSymb(i)
                do k=1,len_trim(line)
                  if(line(k:k) == "m") line(k:k)=" "
                end do
                line=Pack_String(line)
                call Read_Xsym(line,1,MGp%MSymop(i)%Rot)
              else
                det=determ_a(MGp%Symop(i)%Rot)
                MGp%MSymop(i)%Rot=MGp%Symop(i)%Rot*det*nint(MGp%MSymOp(i)%phas)
                call Get_Symsymb(MGp%MSymOp(i)%Rot,(/0.0,0.0,0.0/),line)
                !Expand the operator "line" to convert it to mx,my,mz like
                mxmymz_op=" "
                do j=1,len_trim(line)
                  Select Case(line(j:j))
                    case("x")
                       mxmymz_op=trim(mxmymz_op)//"mx"
                    case("y")
                       mxmymz_op=trim(mxmymz_op)//"my"
                    case("z")
                       mxmymz_op=trim(mxmymz_op)//"mz"
                    case default
                       mxmymz_op=trim(mxmymz_op)//line(j:j)
                  End Select
                end do
                MGp%MSymopSymb(i)=trim(mxmymz_op)
              end if


            end do
            !Decode lattice translations and anti-translations

            !write(unit=*,fmt="(a)") "  Decoding lattice translations and anti-translations"
            m=num_rsym
            do L=2,num_centering
              line=adjustl(cent_strings(L))
              j=index(line," ")
              line=adjustl(line(j+1:))
              j=index(line," ")
              line=line(1:j-1)
              k=index(line,",",back=.true.)
              read(unit=line(k+1:),fmt=*,iostat=ier) n
              if(ier /= 0) then
                 Err_Form=.true.
                 Err_Form_Mess=" Error reading the time inversion in line: "//trim(cent_strings(i))
                 return
              end if
              line=line(1:k-1)
              call Read_Xsym(line,1,Rot,tr)

              do j=1,num_rsym
                m=m+1
                v=MGp%SymOp(j)%tr(:) + tr
                MGp%SymOp(m)%Rot  = MGp%SymOp(j)%Rot
                MGp%SymOp(m)%tr   = modulo_lat(v)
                MGp%MSymOp(m)%Rot = n*MGp%MSymOp(j)%Rot
                MGp%MSymOp(m)%phas= n*MGp%MSymOp(j)%phas
                call Get_Symsymb(MGp%SymOp(m)%Rot,MGp%SymOp(m)%tr,MGp%SymopSymb(m))
                call Get_Symsymb(MGp%MSymOp(m)%Rot,(/0.0,0.0,0.0/),line)
                !Expand the operator "line" to convert it to mx,my,mz like
                mxmymz_op=" "
                do i=1,len_trim(line)
                  Select Case(line(i:i))
                    case("x")
                       mxmymz_op=trim(mxmymz_op)//"mx"
                    case("y")
                       mxmymz_op=trim(mxmymz_op)//"my"
                    case("z")
                       mxmymz_op=trim(mxmymz_op)//"mz"
                    case default
                       mxmymz_op=trim(mxmymz_op)//line(i:i)
                  End Select
                end do
                MGp%MSymopSymb(m)=trim(mxmymz_op)
              end do
            end do
          end if
       end if
       ! Symmetry operators treatment done
       Call cleanup_symmetry_operators(MGp)
       if(Err_Form) then
          return
          !write(unit=*,fmt="(a)") " => "//trim(Err_Form)
       end if

       !Treating irreps

       if(num_irreps == 0) then

          MGp%n_irreps=0

       else
          !write(*,"(a,i3)") " Treating irreps: ",num_irreps
          MGp%n_irreps=num_irreps
          if(allocated(MGp%irrep_dim))          deallocate(MGp%irrep_dim)
          if(allocated(MGp%small_irrep_dim))    deallocate(MGp%small_irrep_dim)
          if(allocated(MGp%irrep_id))           deallocate(MGp%irrep_id)
          if(allocated(MGp%irrep_direction))    deallocate(MGp%irrep_direction)
          if(allocated(MGp%irrep_action))       deallocate(MGp%irrep_action)
          if(allocated(MGp%irrep_modes_number)) deallocate(MGp%irrep_modes_number)
          allocate(MGp%irrep_dim(num_irreps),MGp%small_irrep_dim(num_irreps),MGp%irrep_id(num_irreps), &
                   MGp%irrep_direction(num_irreps),MGp%irrep_action(num_irreps),MGp%irrep_modes_number(num_irreps))

          MGp%irrep_dim=0; MGp%small_irrep_dim=0; MGp%irrep_id=" "; MGp%irrep_direction=" "; MGp%irrep_action=" "
          MGp%irrep_modes_number=0

          do i=1,MGp%n_irreps

            call getword(irreps_strings(i),lab_items,iv)

            !if(iv /= nitems_irreps) write(*,"(2(a,i2))") " => Warning irreps_nitems=",nitems_irreps," /= items read=",iv

            MGp%irrep_id(i)=lab_items(irrep_pos(1))
            if(MGp%irrep_id(i) == "?") then
               MGp%n_irreps=0
               exit
            end if

            if (irrep_pos(2) /= 0) then
               read(unit=lab_items(irrep_pos(2)),fmt=*,iostat=ier) MGp%irrep_dim(i)
               if(ier /= 0) MGp%irrep_dim(i)=0
            end if

            if (irrep_pos(3) /= 0) then
               read(unit=lab_items(irrep_pos(3)),fmt=*,iostat=ier) MGp%small_irrep_dim(i)
               if(ier /= 0) MGp%small_irrep_dim(i)=0
            end if

            if (irrep_pos(4) /= 0) then
               MGp%irrep_direction(i)=lab_items(irrep_pos(4))
            end if

            if (irrep_pos(5) /= 0) then
               MGp%irrep_action(i)=lab_items(irrep_pos(5))
            end if

            if (irrep_pos(6) /= 0) then
               read(unit=lab_items(irrep_pos(6)),fmt=*,iostat=ier) MGp%irrep_modes_number(i)
               if(ier /= 0) MGp%irrep_modes_number(i)=0
            end if

          end do
       end if
       ! End treatment of irreps

       ! Treating propagation vectors
       if(num_kvs == 0) then
         MGp%n_kv=0
       else
         !write(*,"(a,i3)") " Treating propagation vectors: ",num_kvs
         do i=1,MGp%n_kv
            line=adjustl(kv_strings(i))
            j=index(line," ")
            MGp%kv_label(i)=line(1:j-1)
            line=adjustl(line(j+1:))
            n=len_trim(line)
            if(ktag) then
              line=adjustl(line(2:n-1))
              n=n-2
              Call Get_Separator_Pos(line,",",pos,ncar)
            else
              Call Get_Separator_Pos(line," ",pos,ncar)
            end if
            keyword=line(1:pos(1)-1)//"a,"//line(pos(1)+1:pos(2)-1)//"b,"//trim(line(pos(2)+1:))//"c"
            keyword=Pack_String(keyword)
            call Get_Mat_From_Symb(keyword,Matr, (/"a","b","c"/) )
            do k=1,3
               MGp%kv(k,i)=Matr(k,k)
            end do
         end do
       end if
       ! Propagation vectors treatment done!

       !Treating magnetic atoms
       if(num_matom == 0) then
          Am%natoms = 0
          return
       else
          !write(*,"(a,i4)") " Treating magnetic atoms:  ",num_matom
          Call Allocate_Atom_list(num_matom,Am)

          do i=1,Am%natoms

            call getword(atm_strings(i),lab_items,iv)
            !if(iv /= nitems) write(*,"(2(a,i2))") " => Warning nitems=",nitems," /= items read=",iv
            Am%atom(i)%lab=lab_items(lugar(1))
            if (lugar(2) /= 0) then
               Am%atom(i)%SfacSymb=lab_items(lugar(2))(1:4)
               if(index("1234567890+-",lab_items(lugar(2))(2:2)) /= 0 ) then
                  Am%atom(i)%chemSymb=U_case(lab_items(lugar(2))(1:1))
               else
                  Am%atom(i)%chemSymb=U_case(lab_items(lugar(2))(1:1))//L_case(lab_items(lugar(2))(2:2))
               end if
            else
               if(index("1234567890+-",lab_items(lugar(1))(2:2)) /= 0 ) then
                  Am%atom(i)%chemSymb=U_case(lab_items(lugar(1))(1:1))
               else
                  Am%atom(i)%chemSymb=U_case(lab_items(lugar(1))(1:1))//L_case(lab_items(lugar(1))(2:2))
               end if
               Am%atom(i)%SfacSymb=Am%atom(i)%chemSymb
            end if
            call getnum_std(lab_items(lugar(3)),values,std,iv)    ! _atom_site_fract_x
            Am%atom(i)%x(1)=values(1)
            Am%atom(i)%x_std(1)=std(1)
            call getnum_std(lab_items(lugar(4)),values,std,iv)    ! _atom_site_fract_y
            Am%atom(i)%x(2)=values(1)
            Am%atom(i)%x_std(2)=std(1)
            call getnum_std(lab_items(lugar(5)),values,std,iv)    ! _atom_site_fract_z
            Am%atom(i)%x(3)=values(1)
            Am%atom(i)%x_std(3)=std(1)

            if (lugar(6) /= 0) then  ! _atom_site_U_iso_or_equiv
               call getnum_std(lab_items(lugar(6)),values,std,iv)
               Am%atom(i)%ueq=values(1)
               Am%atom(i)%Biso=values(1)*78.95683521     !If anisotropic they
               Am%atom(i)%Biso_std=std(1)*78.95683521    !will be put to zero
            else if (lugar(10) /= 0) then    ! _atom_site_B_iso_or_equiv
               call getnum_std(lab_items(lugar(10)),values,std,iv)
               Am%atom(i)%ueq=values(1)/78.95683521
               Am%atom(i)%Biso=values(1)     !If anisotropic they
               Am%atom(i)%Biso_std=std(1)    !will be put to zero
            else
               Am%atom(i)%ueq=0.0
               Am%atom(i)%Biso=0.0
               Am%atom(i)%Biso_std=0.0
            end if
            Am%atom(i)%utype="u_ij"

            if (lugar(7) /= 0) then ! _atom_site_occupancy
               call getnum_std(lab_items(lugar(7)),values,std,iv)
            else
               values=1.0
               std=0.0
            end if
            Am%atom(i)%occ=values(1)
            Am%atom(i)%occ_std=std(1)

            if(lugar(8) /= 0) then
              read(unit=lab_items(lugar(8)),fmt=*) Mult
              Am%atom(i)%mult=Mult
            else
              Call Get_mOrbit(Am%atom(i)%x,MGp,Mult,orb)
              Am%atom(i)%mult=Mult
            end if
            !Conversion from occupancy to occupation factor
            Am%atom(i)%occ=Am%atom(i)%occ*real(Mult)/real(MGp%Multip)

            if(lugar(9) /= 0) then
               Am%atom(i)%wyck=adjustl(trim(lab_items(lugar(9))))
            end if

          end do
       end if

       !Treating moments of magnetic atoms
       if(num_mom /= 0) then
          !write(*,"(a,i4)") " Treating magnetic moments:  ",num_mom
          m=4
          if(mom_symmform) m=5
          do i=1,num_mom
            call getword(mom_strings(i),lab_items,iv)
            !write(*,"(2i6,tr4,5(a,tr3))") k,iv,lab_items(1:iv)
            if(iv /= m) then
               Err_Form=.true.
               write(unit=Err_Form_Mess,fmt="(a,i4)")" Error reading magnetic moment #",i
               Err_Form_Mess=trim(Err_Form_Mess)//" -> 4-5 items expected in this line: 'Label mx my mz', read: "// &
                                                      trim(mom_strings(i))
               return
            end if
            label=Lab_items(1)
            do j=1,Am%natoms
               if(label == Am%Atom(j)%lab) then
                 do k=1,3
                     call getnum_std(lab_items(1+k),values,std,iv)
                     Am%Atom(j)%M_xyz(k)=values(1)
                     Am%Atom(j)%sM_xyz(k)=std(1)
                 end do
                 Am%Atom(j)%moment=99.0  !used for indicating that this atom is susceptible to bring a magnetic moment
               end if
            end do
          end do
       end if

       if(num_constr /= 0) then

         !write(*,"(a,i4)") " Treating constraints:  ",num_constr
         do i=1,num_constr
           line=adjustl(constr_strings(i))
           j=index(line," ")
           label=line(1:j-1)
           keyword=adjustl(line(j+1:))
           Call Get_Separator_Pos(keyword,",",pos,ncar)
           if(ncar == 0) then !There are no ","
             j=index(keyword," ")
             shubk=keyword(1:j-1)//","
             keyword=adjustl(keyword(j+1:))
             j=index(keyword," ")
             shubk=trim(shubk)//keyword(1:j-1)//","
             keyword=trim(shubk)//trim(adjustl(keyword(j+1:)))
           end if
           do j=1,len_trim(keyword)
             if(keyword(j:j) == "m") keyword(j:j) = " "
           end do
           keyword=Pack_String(keyword)
           !write(*,"(a)") "  constr_string: "//trim(line)
           !write(*,"(a)") "        keyword: "//trim(keyword)
           call Get_Mat_From_Symb(keyword,Matr, (/"x","y","z"/) )
           !write(*,"(9f10.3)") Matr
           do j=1,Am%natoms
             if(label == Am%Atom(j)%lab) then
                Am%Atom(j)%M_xyz=matmul(Matr,Am%Atom(j)%M_xyz)
                Am%Atom(j)%AtmInfo=constr_strings(i)
                Am%Atom(j)%moment=99.0  !used for indicating that this atom is susceptible to bring a magnetic moment
                exit
             end if
           end do
           !The treatment of the codes will be done in the future
         end do
       end if

       if(num_magscat > 0) then !Reading the valence for determining the magnetic form factor
         do i=1,num_magscat
           call getword(mag_scatt_string(i),lab_items,iv)
           do j=1,Am%natoms
             if(Am%atom(j)%chemSymb == lab_items(1)) then
               Am%atom(j)%SfacSymb=lab_items(2)
               if(lab_items(2) /= ".") then !magnetic atoms
                  Am%Atom(j)%moment=99.0  !used for indicating that this atom is susceptible to bring a magnetic moment
               end if
             end if
           end do
         end do
       end if

       !Get pointers to the magnetic form factors
       !Stored for each atom in the component ind(1)
       call Set_Magnetic_Form()

       !---- Find Species in Magnetic_Form ----!
       do i=1,Am%natoms
          symbcar=get_magnetic_form_factor(Am%atom(i)%SfacSymb)
          do j=1,num_mag_form
             if (symbcar /= Magnetic_Form(j)%Symb) cycle
             Am%atom(i)%ind(2)=j
             Am%atom(i)%SfacSymb=symbcar
             exit
          end do
       end do

       return
    End Subroutine Readn_Set_Magnetic_Structure_MCIF

    !!--++
    !!--++ Subroutine Readn_Set_XTal_CFL(file_dat,nlines,Cell,SpG,A,CFrame,NPhase,Job_Info)
    !!--++    character(len=*),dimension(:),intent(in)   :: file_dat
    !!--++    integer,                      intent(in)   :: nlines
    !!--++    Type (Crystal_Cell_Type),     intent(out)  :: Cell
    !!--++    Type (Space_Group_Type),      intent(out)  :: SpG
    !!--++    Type (atom_list_type),        intent(out)  :: A
    !!--++    character(len=*),    optional,intent(in)   :: CFrame
    !!--++    Integer,             optional,intent( in)  :: Nphase
    !!--++    Type(Job_Info_type), optional,intent(out)  :: Job_Info
    !!--++
    !!--++ (Private)
    !!--++ Read and Set Crystal Information in a CFL File
    !!--++
    !!--++ Update: April - 2005
    !!
    Subroutine Readn_Set_XTal_CFL(file_dat,nlines,Cell,SpG,A,CFrame,NPhase,Job_Info)
       !---- Arguments ----!
       character(len=*),dimension(:),intent(in)   :: file_dat
       integer,                      intent(in)   :: nlines
       Type (Crystal_Cell_Type),     intent(out)  :: Cell
       Type (Space_Group_Type),      intent(out)  :: SpG
       Type (atom_list_type),        intent(out)  :: A
       character(len=*),    optional,intent(in)   :: CFrame
       Integer,             optional,intent( in)  :: Nphase
       Type(Job_Info_type), optional,intent(out)  :: Job_Info

       !---- Local variables ----!
       character(len=132)               :: line
       character(len= 20)               :: Spp
       character(len= 40),dimension(192):: gen
       integer                          :: i, nauas, ndata, iph, n_ini,n_end,ngen,k,nsym
       integer, parameter               :: maxph=21  !Maximum number of phases "maxph-1"
       integer, dimension(maxph)        :: ip

       real(kind=cp),dimension(3):: vet

       !---- Standard CrysFML file *.CFL ----!
       nauas=0
       ndata=0
       ip=nlines
       ip(1)=1

       !---- Calculating number of Phases ----!
       do i=1,nlines
          line=adjustl(file_dat(i))
          if (l_case(line(1:6)) == "phase_")  then
             ndata=ndata+1
             ip(ndata)=i
          end if
       end do

       !---- Reading Phase Information ----!
       iph=1
       if (present(nphase)) iph=nphase
       if (present(Job_Info)) then
          n_ini=ip(iph)           !Updated values to handle non-conventional order
          n_end=ip(iph+1)
          call Get_Job_Info(file_dat,n_ini,n_end,Job_info)
       end if

       !---- Reading Cell Parameters ----!
       n_ini=ip(iph)           !Updated values to handle non-conventional order
       n_end=ip(iph+1)
       if(present(CFrame)) then
         call read_File_Cell(file_dat,n_ini,n_end,Cell,CFrame) !Read and construct Cell
       else
         call read_File_Cell(file_dat,n_ini,n_end,Cell) !Read and construct Cell
       end if
       if (err_form) return

       !---- Reading Space Group Information ----!
       n_ini=ip(iph)           !Updated values to handle non-conventional order
       n_end=ip(iph+1)
       call read_File_Spg (file_dat,n_ini,n_end,Spp)
       if (err_form) then !Try to read symmetry operators or generators
         ngen=0
         nsym=0
         do i=n_ini, n_end
           line=l_case(adjustl(file_dat(i)))
           if(line(1:4) == "symm") nsym=nsym+1
           if(line(1:3) == "gen")  ngen=ngen+1
         end do
         if(ngen > 0) then
           k=0
           do i=n_ini, n_end
             line=l_case(adjustl(file_dat(i)))
             if(line(1:3) == "gen")  then
              k=k+1
              gen(k)=adjustl(line(5:))
             end if
           end do
           call Set_SpaceGroup(" ",SpG,gen,ngen,"gen")   !Construct the space group from generators
         else if (nsym > 0) then
           k=0
           do i=n_ini, n_end
             line=l_case(adjustl(file_dat(i)))
             if(line(1:4) == "symm")  then
              k=k+1
              gen(k)=adjustl(line(6:))
             end if
           end do
           call Set_SpaceGroup(" ",SpG,gen,nsym,"fix")  !Construct the space group from fixed symmetry elements
         else
           return
         end if
       else
          call Set_SpaceGroup(Spp,SpG) !Construct the space group
       end if
       !---- Read Atoms Information ----!
       n_ini=ip(iph)           !Updated values to handle non-conventional order
       n_end=ip(iph+1)

       !---- Calculating number of Atoms in the Phase ----!
       do i=n_ini,n_end
          line=adjustl(file_dat(i))
          if (l_case(line(1:4)) == "atom")  nauas=nauas+1
       end do

       if (nauas > 0) then
          call Allocate_atom_list(nauas,A)  !allocation space for Atom list
          call read_File_Atom(file_dat,n_ini,n_end,A)
          if (err_form) return

          do i=1,A%natoms
             vet=A%atom(i)%x
             A%atom(i)%Mult=Get_Multip_Pos(vet,SpG)
             if(A%atom(i)%occ < epsv) A%atom(i)%occ=real(A%atom(i)%Mult)/real(SpG%Multip)
             if (A%atom(i)%thtype == "aniso") then
                select case (A%atom(i)%Utype)
                   case ("u_ij")
                      A%atom(i)%u(1:6) =  Convert_U_Betas(A%atom(i)%u(1:6),Cell)
                   case ("b_ij")
                      A%atom(i)%u(1:6) =  Convert_B_Betas(A%atom(i)%u(1:6),Cell)
                end select
                A%atom(i)%Utype="beta"
             end if
          end do
       end if

       return
    End Subroutine Readn_Set_XTal_CFL

    !!--++
    !!--++ Subroutine Readn_Set_XTal_CFL_Molec(file_dat, nlines, Molcrys)
    !!--++    character(len=*),dimension(:),  intent(in)     :: file_dat
    !!--++    integer,                        intent(in)     :: nlines
    !!--++    Type (Molecular_Crystal_Type),  intent(in out) :: Molcrys
    !!--++
    !!--++ (Private)
    !!--++ Read Molecule Information in a CFL
    !!--++
    !!--++ Update: April - 2005
    !!
    Subroutine Readn_Set_XTal_CFL_Molec(file_dat, nlines, Molcrys, Nphase)
       !---- Arguments ----!
       character(len=*),dimension(:),  intent(in)     :: file_dat
       integer,                        intent(in)     :: nlines
       type (Molecular_Crystal_Type),  intent(in out) :: Molcrys
       Integer, optional,              intent(in)     :: Nphase

       !---- Local variables ----!
       character(len=132)            :: line
       integer                       :: i,n,nmol,npos,n_ini,n_end,ierr,nauas, iph, ndata
       integer, parameter               :: maxph=21  !Maximum number of phases "maxph-1"
       integer, dimension(maxph)        :: ip
       real(kind=cp)                 :: theta,phi,chi
       real(kind=cp), dimension(3)   :: x1f,x2f,x3f
       real(kind=cp), dimension(3,3) :: EuM

       !---- Standard CrysFML file *.CFL ----!
       nauas=0
       ndata=0
       ip=nlines
       ip(1)=1

       !---- Calculating number of Phases ----!
       do i=1,nlines
          line=adjustl(file_dat(i))
          if (l_case(line(1:6)) == "phase_")  then
             ndata=ndata+1
             ip(ndata)=i
          end if
       end do

       !---- Reading Phase Information ----!

       if (present(nphase)) then
           iph=nphase
       else
           iph=1
       end if

       n_ini=ip(iph)
       n_end=ip(iph+1)

       !---- Detecting the Molecules defined in the file ----!
       nmol=0
       do i=n_ini,n_end
          line=u_case(adjustl(file_dat(i)))
          if (line(1:1) == " ") cycle
          if (line(1:1) == "!") cycle
          npos=index(line,"MOLE")
          if (npos /= 0) nmol=nmol+1
       end do
       if (nmol==0) return

       !---- Allocating Memory for all molecules ----!
       if (allocated(molcrys%mol)) deallocate(molcrys%mol)
       molcrys%n_mol=nmol
       allocate(molcrys%mol(nmol))

       !---- Reading Molecules ----!

       do n=1,nmol
          !---- Read ----!
          do i=n_ini,n_end
             line=u_case(adjustl(file_dat(i)))
             if (line(1:1) == " ") cycle
             if (line(1:1) == "!") cycle
             npos=index(line,"MOLE")
             if (npos == 0) cycle
             call read_molecule(file_dat,n_ini,n_end,molcrys%mol(n))
             err_form=err_molec
             ERR_Form_Mess=err_molec_mess
             if (err_form) then
                molcrys%n_mol=n-1
                return
             end if
             exit
          end do

          !---- Search for three points (fractional coordinates) ----!
          !---- defining a Cartesian frame                       ----!
          do
             if (n_ini > n_end) exit
             line=adjustl(file_dat(n_ini))
             if (u_case(line(1:9)) == "XYZ_FRAME") then
                read(unit=line(10:),fmt=*,iostat=ierr) x1f,x2f,x3f
                if (ierr == 0) then
                   call get_euler_from_fract(x1f,x2f,x3f,molcrys%Cell%Cr_Orth_cel,phi,theta,chi,EuM, Code="D")
                   molcrys%mol(n)%orient(1)= phi
                   molcrys%mol(n)%orient(2)= theta
                   molcrys%mol(n)%orient(3)= chi
                   molcrys%mol(n)%xcentre= x3f
                   call Set_euler_matrix(molcrys%mol(n)%rot_type, phi,theta,chi,EuM)
                   molcrys%mol(n)%Euler=EuM
                   molcrys%mol(n)%is_EulerMat=.true.
                   molcrys%mol(n)%in_Xtal=.true.
                end if
                n_ini=n_ini+1
                exit
             else
                if (u_case(line(1:4)) =="MOLE") exit
                n_ini=n_ini+1
             end if
          end do

       end do

       return
    End Subroutine Readn_Set_XTal_CFL_Molec

    !!--++
    !!--++ Subroutine Readn_Set_XTal_CFL_Shub(file_dat,nlines,Cell,SpG,A,CFrame,NPhase,Job_Info)
    !!--++    character(len=*),dimension(:),intent(in)   :: file_dat
    !!--++    integer,                      intent(in)   :: nlines
    !!--++    Type (Crystal_Cell_Type),     intent(out)  :: Cell
    !!--++    Type (Magnetic_Space_Group_Type), intent(out)  :: SpG
    !!--++    Type (atom_list_type),        intent(out)  :: A
    !!--++    character(len=*),    optional,intent(in)   :: CFrame
    !!--++    Integer,             optional,intent( in)  :: Nphase
    !!--++    Type(Job_Info_type), optional,intent(out)  :: Job_Info
    !!--++
    !!--++ (Private)
    !!--++ Read and Set Crystal Information in a CFL File
    !!--++
    !!--++ Update: April - 2005
    !!
    Subroutine Readn_Set_XTal_CFL_Shub(file_dat,nlines,Cell,SpG,A,CFrame,NPhase,Job_Info)
       !---- Arguments ----!
       character(len=*),dimension(:),    intent(in)   :: file_dat
       integer,                          intent(in)   :: nlines
       Type (Crystal_Cell_Type),         intent(out)  :: Cell
       Type (Magnetic_Space_Group_Type), intent(out)  :: SpG
       Type (atom_list_type),            intent(out)  :: A
       character(len=*),        optional,intent(in)   :: CFrame
       Integer,                 optional,intent( in)  :: Nphase
       Type(Job_Info_type),     optional,intent(out)  :: Job_Info

       !---- Local variables ----!
       character(len=132)               :: line
       character(len= 50)               :: Spp,setting
       !character(len= 40),dimension(192):: gen
       integer                          :: i, nauas, ndata, iph, n_ini,n_end,k !,ngen
       integer, parameter               :: maxph=21  !Maximum number of phases "maxph-1"
       integer, dimension(maxph)        :: ip

       real(kind=cp),dimension(3):: vet

       !---- Standard CrysFML file *.CFL ----!
       nauas=0
       ndata=0
       ip=nlines
       ip(1)=1

       !---- Calculating number of Phases ----!
       do i=1,nlines
          line=adjustl(file_dat(i))
          if (l_case(line(1:6)) == "phase_")  then
             ndata=ndata+1
             ip(ndata)=i
          end if
       end do

       !---- Reading Phase Information ----!
       iph=1
       if (present(nphase)) iph=nphase
       if (present(Job_Info)) then
          n_ini=ip(iph)           !Updated values to handle non-conventional order
          n_end=ip(iph+1)
          call Get_Job_Info(file_dat,n_ini,n_end,Job_info)
       end if

       !---- Reading Cell Parameters ----!
       n_ini=ip(iph)           !Updated values to handle non-conventional order
       n_end=ip(iph+1)
       if(present(CFrame)) then
         call read_File_Cell(file_dat,n_ini,n_end,Cell,CFrame) !Read and construct Cell
       else
         call read_File_Cell(file_dat,n_ini,n_end,Cell) !Read and construct Cell
       end if
       if (err_form) return

       !---- Reading Space Group Information ----!
       n_ini=ip(iph)           !Updated values to handle non-conventional order
       n_end=ip(iph+1)
       call read_File_Spg (file_dat,n_ini,n_end,Spp)
       i=index(Spp,"{")
       k=len_trim(Spp)
       setting=" "
       if(i /= 0) then
         setting=Spp(i+1:k-1)
         Spp=Spp(1:i-1)
       end if
       call Set_Magnetic_Space_Group(Spp,setting,Spg) !Construct the magnetic space group

       !---- Read Atoms Information ----!
       n_ini=ip(iph)           !Updated values to handle non-conventional order
       n_end=ip(iph+1)

       !---- Calculating number of Atoms in the Phase ----!
       do i=n_ini,n_end
          line=adjustl(file_dat(i))
          if (l_case(line(1:4)) == "atom")  nauas=nauas+1
       end do

       if (nauas > 0) then
          call Allocate_atom_list(nauas,A)  !allocation space for Atom list
          call read_File_Atom(file_dat,n_ini,n_end,A)
          if (err_form) return

          do i=1,A%natoms
             vet=A%atom(i)%x
             A%atom(i)%Mult=Get_Multip_Pos(vet,SpG)
             if(A%atom(i)%occ < epsv) A%atom(i)%occ=real(A%atom(i)%Mult)/real(SpG%Multip)
             if (A%atom(i)%thtype == "aniso") then
                select case (A%atom(i)%Utype)
                   case ("u_ij")
                      A%atom(i)%u(1:6) =  Convert_U_Betas(A%atom(i)%u(1:6),Cell)
                   case ("b_ij")
                      A%atom(i)%u(1:6) =  Convert_B_Betas(A%atom(i)%u(1:6),Cell)
                end select
                A%atom(i)%Utype="beta"
             end if
          end do
       end if

       return
    End Subroutine Readn_Set_XTal_CFL_Shub
    !!--++
    !!--++ Subroutine Readn_Set_XTal_CIF(file_dat, nlines, Cell, Spg, A, CFrame, NPhase)
    !!--++    character(len=*),dimension(:),intent(in)   :: file_dat
    !!--++    integer,                      intent(in)   :: nlines
    !!--++    Type (Crystal_Cell_Type),     intent(out)  :: Cell
    !!--++    Type (Space_Group_Type),      intent(out)  :: SpG
    !!--++    Type (atom_list_type),        intent(out)  :: A
    !!--++    Character(len=*),    optional,intent( in)  :: CFrame
    !!--++    Integer,             optional,intent( in)  :: Nphase
    !!--++
    !!--++ (Private)
    !!--++ Read and Set Crystal Information in a CIF File
    !!--++
    !!--++ Update: April - 2005
    !!
    Subroutine Readn_Set_XTal_CIF(file_dat, nlines, Cell, Spg, A, CFrame, NPhase)
       !---- Arguments ----!
       character(len=*),dimension(:),intent(in)   :: file_dat
       integer,                      intent(in)   :: nlines
       Type (Crystal_Cell_Type),     intent(out)  :: Cell
       Type (Space_Group_Type),      intent(out)  :: SpG
       Type (atom_list_type),        intent(out)  :: A
       Character(len=*),    optional,intent( in)  :: CFrame
       Integer,             optional,intent( in)  :: Nphase

       !---- Local Variables ----!
       character(len=132)                :: line
       character(len= 20)                :: Spp
       character(len=60), dimension(192) :: symm_car

       integer                   :: i, nauas, ndata, iph, n_ini,n_end,noper
       integer, parameter        :: maxph=250  !Maximum number of phases "maxph-1"
       integer, dimension(maxph) :: ip

       real(kind=cp),dimension(6):: vet,vet2

       ip=nlines
       ip(1)=1

       !---- First determine if there is more than one structure ----!
       do i=1,nlines
          line=adjustl(file_dat(i))
          if (l_case(line(1:5)) == "data_" .and. l_case(line(1:11)) /= "data_global" )  then
             n_ini=i
             ip(1)=i
             exit
          end if
       end do

       ndata=0
       do i=n_ini,nlines
          line=adjustl(file_dat(i))
          if (l_case(line(1:5)) == "data_")  then
             ndata=ndata+1
             if (ndata > maxph-1) then
                err_form=.true.
                ERR_Form_Mess=" => Too many phases in this file "
                return
             end if
             ip(ndata)=i   !Pointer to the number of the line starting a single phase
          end if
       end do

       iph=1
       if (present(nphase)) iph=nphase

       !---- Read Cell Parameters ----!
       n_ini=ip(iph)           !Updated values to handle non-conventional order
       n_end=ip(iph+1)
       call Read_Cif_Cell(file_dat,n_ini,n_end,vet,vet2)
       if (err_form) return
       if(present(CFrame)) then
         call Set_Crystal_Cell(vet(1:3),vet(4:6),Cell,CFrame,vet2(1:3),vet2(4:6))
       else
         call Set_Crystal_Cell(vet(1:3),vet(4:6),Cell,"A",vet2(1:3),vet2(4:6))
       end if
       !---- Read Atoms Information ----!
       n_ini=ip(iph)           !Updated values to handle non-conventional order
       n_end=ip(iph+1)
       call Read_Cif_Atom(file_dat,n_ini,n_end,nauas,A)
       if (err_form) return

       !---- SpaceGroup Information ----!
       n_ini=ip(iph)           !Updated values to handle non-conventional order
       n_end=ip(iph+1)
       call Read_Cif_Hm(file_dat,n_ini,n_end,Spp)

       n_ini=ip(iph)           !Updated values to handle non-conventional order
       n_end=ip(iph+1)
       if (len_trim(Spp) == 0) call Read_Cif_Hall(file_dat,n_ini,n_end,Spp)

       if (len_trim(Spp) == 0) then
          n_ini=ip(iph)           !Updated values to handle non-conventional order
          n_end=ip(iph+1)
          call Read_Cif_Symm(file_dat,n_ini,n_end,noper,symm_car)

          if (noper ==0) then
             err_form=.true.
             ERR_Form_Mess=" => No Space Group/No Symmetry information in this file "
             return
          else
             call Set_SpaceGroup("  ",SpG,symm_car,noper,"GEN")
          end if
       else
          call Set_SpaceGroup(Spp,SpG) !Construct the space group
       end if

       !---- Modify occupation factors and set multiplicity of atoms
       !---- in order to be in agreement with the definitions of Sfac in CrysFML
       !---- Convert Us to Betas and Uiso to Biso
       do i=1,A%natoms
          vet(1:3)=A%atom(i)%x
          A%atom(i)%Mult=Get_Multip_Pos(vet(1:3),SpG)
          A%atom(i)%Occ=A%atom(i)%Occ*real(A%atom(i)%Mult)/max(1.0,real(SpG%Multip))
          if(A%atom(i)%occ < epsv) A%atom(i)%occ=real(A%atom(i)%Mult)/max(1.0,real(SpG%Multip))

          select case (A%atom(i)%thtype)
             case ("isotr")
                A%atom(i)%biso= A%atom(i)%ueq*78.95683521

             case ("aniso")
                select case (A%atom(i)%Utype)
                   case ("u_ij")
                      A%atom(i)%u(1:6) =  Convert_U_Betas(A%atom(i)%u(1:6),Cell)
                   case ("b_ij")
                      A%atom(i)%u(1:6) = Convert_B_Betas(A%atom(i)%u(1:6),Cell)
                end select
                A%atom(i)%Utype="beta"

             case default
                A%atom(i)%biso = A%atom(i)%ueq*78.95683521
                A%atom(i)%thtype = "isotr"
          end select
       end do

       return
    End Subroutine Readn_Set_XTal_CIF

    !!--++
    !!--++ Subroutine Readn_Set_XTal_PCR(file_dat, nlines, Cell, Spg, A, CFrame, NPhase)
    !!--++    character(len=*),dimension(:),intent(in)   :: file_dat
    !!--++    integer,                      intent(in)   :: nlines
    !!--++    Type (Crystal_Cell_Type),     intent(out)  :: Cell
    !!--++    Type (Space_Group_Type),      intent(out)  :: SpG
    !!--++    Type (atom_list_type),        intent(out)  :: A
    !!--++    character(len=*),    optional,intent(in)   :: CFrame
    !!--++    Integer,             optional,intent( in)  :: Nphase
    !!--++
    !!--++ (Private)
    !!--++ Read and Set Crystal Information in a PCR File
    !!--++
    !!--++ Update: 17/05/2010
    !!
    Subroutine Readn_Set_XTal_PCR(file_dat, nlines, Cell, Spg, A, CFrame, NPhase)
       !---- Arguments ----!
       character(len=*),dimension(:),intent(in)   :: file_dat
       integer,                      intent(in)   :: nlines
       Type (Crystal_Cell_Type),     intent(out)  :: Cell
       Type (Space_Group_Type),      intent(out)  :: SpG
       Type (atom_list_type),        intent(out)  :: A
       character(len=*),    optional,intent(in)   :: CFrame
       Integer,             optional,intent(in)   :: Nphase

       !---- Local Variables ----!
       logical                           :: multi,ask_phase,is_codewords
       character(len=132)                :: line
       character(len= 20)                :: Spp, label
       integer                           :: i,j, k,iv, nauas, ndata, iph, n_ini,n_end, nlong1
       integer, parameter                :: maxph=21  !Maximum number of phases "maxph-1"
       integer, dimension(maxph)         :: ip
       integer, dimension(30)            :: ivet

       real(kind=cp),dimension(30)       :: vet

       ip=nlines
       ip(1)=1

       !> Simple / Multi format
       multi=.false.
       do i=1,nlines
          line=adjustl(file_dat(i))
          if (line(1:1) =='!' .or. line(1:1)==' ') cycle
          if (index(line,'NPATT ') <=0) cycle
          multi=.true.
       end do

       !> Number of Phases
       if (.not. multi) then
          do i=2,nlines
             line=adjustl(file_dat(i))
             if (line(1:1) =='!' .or. line(1:1)==' ') cycle
             call getnum(line,vet,ivet,iv)
             if (iv > 3) then
                iph=ivet(3)
                exit
             end if
          end do

       else
          do i=1,nlines
             line=adjustl(file_dat(i))
             if (line(1:4) /='!Nph') cycle

             line=adjustl(file_dat(i+1))
             call getnum(line,vet,ivet,iv)
             if (iv > 1) then
                iph=ivet(1)
                exit
             end if
          end do
       end if
       if (iph == 0) then
          err_form=.true.
          ERR_Form_Mess=" No Phase information was found in this PCR file. Please, check it! "
          return
       end if

       !> Locate where begin each Phase
       k=0
       ask_phase=.true.

       do i=1,nlines
          line=adjustl(file_dat(i))
          if (ask_phase) then
             if (index(line,'Data for PHASE') <= 0) cycle
          else
             if (line(1:1) /='!') then
                k=k+1
                ip(k)=i
                if (k == iph) exit

                ask_phase=.true.
             end if
             cycle
          end if
          ask_phase=.false.
       end do
       if (iph /= k) then
          err_form=.true.
          ERR_Form_Mess=" Locating Phases failed in this PCR. Please, check it!"
          return
       end if

       !> Select the Phase
       iph=1
       if (present(nphase)) iph=nphase
       n_ini=ip(iph)
       n_end=ip(iph+1)

       !---- Read Cell Parameters ----!
       do i=n_ini,n_end
          if (index(file_dat(i),'alpha') /=0 .and. index(file_dat(i),'gamma') /=0) then
             do j=i+1,n_end
                line=adjustl(file_dat(j))
                if (line(1:1) == '!' .or. line(1:1) == ' ') cycle
                iv=index(line,'#')
                if (iv > 1) line=line(1:iv-1)

                call getnum(line, vet, ivet,iv)
                if (iv /= 6) then
                   err_form=.true.
                   ERR_Form_Mess=" => Problems reading Cell Parameters on PCR file "
                   return
                end if
                if(present(CFrame)) then
                  call Set_Crystal_Cell(vet(1:3),vet(4:6),Cell,CFrame)
                else
                  call Set_Crystal_Cell(vet(1:3),vet(4:6),Cell)
                end if
                exit
             end do
             exit
          end if
       end do

       !---- SpaceGroup Information ----!
       Spp=' '
       do i=n_ini,n_end
          line=adjustl(file_dat(i))
          if (line(1:1) == '!' .or. line(1:1)==' ') cycle
          if (index(file_dat(i),'<--Space') /=0) then
             j=index(file_dat(i),'<--Space')
             Spp=adjustl(file_dat(i)(1:j-1))
             if (len_trim(Spp) <= 0) then
                err_form=.true.
                ERR_Form_Mess=" => Problems reading Space group on PCR file "
                return
             end if
             call Set_SpaceGroup(Spp,SpG) !Construct the space group
             exit
          end if
       end do

       !---- Read Atoms Information ----!
       do i=n_ini,n_end
          line=adjustl(file_dat(i))
          if (line(1:4) /= '!Nat') cycle
          do j=i+1,n_end
             line=adjustl(file_dat(j))
             if (line(1:1) == '!' .or. line(1:1)==' ') cycle
             call getnum(line(1:5),vet,ivet,iv)
             ndata=ivet(1)
             exit
          end do
          exit
       end do

       if (ndata > 0) then
          call allocate_atom_list(ndata,A)

          is_codewords=.false.
          nauas=0

          do i=n_ini,n_end
             line=adjustl(file_dat(i))
             if (index(line,'!Atom') == 0 .or. index(line,'Typ') == 0) cycle

             do j=i+1,n_end
                line=adjustl(file_dat(j))
                if (line(1:1) == '!' .or. line(1:1)==' ') cycle
                if (is_codewords) then
                   is_codewords=.false.
                   cycle
                end if

                iv=index(line,'#')
                if (iv > 1) line=line(1:iv-1)

                nauas=nauas+1
                ! Atom Label
                call cutst(line,nlong1,label)
                A%atom(nauas)%lab=trim(label)

                ! Atom Type
                call cutst(line,nlong1,label)
                A%Atom(nauas)%chemsymb=U_case(label(1:1))//L_case(label(2:2))

                ! Atom Coordinates,Biso and Occ
                call getnum(line,vet,ivet,iv)
                if (iv < 5) then    !Line reading for the second time anisotropic temperature factors
                   nauas = nauas -1 !see below
                   is_codewords=.true.
                   cycle
                end if

                A%atom(nauas)%x=vet(1:3)
                A%atom(nauas)%Mult=Get_Multip_Pos(vet(1:3),SpG)
                A%atom(nauas)%biso=vet(4)
                A%atom(nauas)%occ=vet(5)
                A%atom(nauas)%thtype='isotr'
                A%atom(nauas)%Utype="beta"
                if (ivet(8) == 2) then    ! Anisotropic reading
                   A%atom(nauas)%thtype='aniso'
                   call getnum(file_dat(j+2),vet,ivet,iv)
                   A%atom(nauas)%u(1:6)=vet(1:6)
                end if
                is_codewords=.true.
                if (nauas == ndata) exit
             end do
             exit
          end do
       end if

       return
    End Subroutine Readn_Set_XTal_PCR



    !!--++
    !!--++ Subroutine Readn_Set_Xtal_Structure_Molcr(filenam,Molcrys,Mode,Iphase, Job_Info, file_list,CFrame)
    !!--++    character(len=*),              intent( in)     :: filenam  ! In -> Name of the file
    !!--++    Type (Molecular_Crystal_Type), intent(out)     :: Molcrys  ! Molecular crytal
    !!--++    Character(len=*),    optional, intent( in)     :: Mode     ! In -> if Mode="CIF" filenam
    !!--++                                                                       is of CIF type format
    !!--++    Integer,             optional, intent( in)     :: Iphase   ! Number of the phase.
    !!--++    Type(Job_Info_type), optional, intent(out)     :: Job_Info ! Diffaction conditions
    !!--++    Type(file_list_type),optional, intent(in out)  :: file_list! Complete file to be used by
    !!--++                                                              the calling program or other procedures
    !!--++    Character(len=*),    optional, intent(in)      :: CFrame
    !!--++    Overloaded
    !!--++    Subroutine to read and input file and construct the crystal structure
    !!--++    in terms of the ofjects Cell, SpG and A. The optional argument Iphase is an integer
    !!--++    telling to the program to read the phase number Iphase in the case of the presence
    !!--++    of more than one phase. If absent only the first phase is read.
    !!--++
    !!--++ Update: April - 2005
    !!
    Subroutine Readn_Set_Xtal_Structure_Molcr(filenam,Molcrys,Mode,Iphase,Job_Info,file_list,CFrame)
       !---- Arguments ----!
       character(len=*),              intent( in)     :: filenam
       Type (Molecular_Crystal_Type), intent(out)     :: Molcrys
       Character(len=*),     optional,intent( in)     :: Mode
       Integer,              optional,intent( in)     :: Iphase
       Type(Job_Info_type),  optional,intent(out)     :: Job_Info
       Type(file_list_type), optional,intent(in out)  :: file_list
       Character(len=*),     optional,intent(in)      :: CFrame
       !---- Local variables -----!
       Type (Atom_list_type)                         :: A
       character(len=132), allocatable, dimension(:) :: file_dat
       character(len=3)                              :: modec
       integer                                       :: i,nlines


       call init_err_form()

       nlines=0
       if (present(file_list)) nlines=file_list%nlines

       !---- Number of Lines in the input file ----!
       if(nlines == 0) then
           call Number_Lines(trim(filenam), nlines)
           if (nlines==0) then
              err_form=.true.
              ERR_Form_Mess="The file "//trim(filenam)//" contains nothing"
              return
           else
              if (allocated(file_dat)) deallocate( file_dat)
              allocate( file_dat(nlines))
              call reading_Lines(trim(filenam),nlines,file_dat)
           end if
           if (present(file_list)) then
              file_list%nlines=nlines
              if (allocated(file_list%line)) deallocate(file_list%line)
              allocate(file_list%line(nlines))
              file_list%line=file_dat
           end if
       else
           if (allocated(file_dat)) deallocate( file_dat)
           allocate( file_dat(nlines))
           file_dat=file_list%line
       end if


       !---- Define the type of file: CIF, CFL, RES,... ----!
       modec=" "
       if (present(mode)) modec=l_case(mode(1:3))

       select case(modec)
           case("cif")
              if (present(iphase)) then
                 if(present(CFrame)) then
                   call readn_set_xtal_cif(file_dat,nlines,molcrys%Cell,molcrys%Spg, A,CFrame,NPhase=IPhase)
                 else
                   call readn_set_xtal_cif(file_dat,nlines,molcrys%Cell,molcrys%Spg, A,NPhase=IPhase)
                 end if
              else
                 if(present(CFrame)) then
                   call readn_set_xtal_cif(file_dat,nlines,molcrys%Cell,molcrys%Spg,A,CFrame)
                 else
                   call readn_set_xtal_cif(file_dat,nlines,molcrys%Cell,molcrys%Spg,A)
                 end if
              end if

           case("pcr")
              if (present(iphase)) then
                 if(present(CFrame)) then
                   call readn_set_xtal_pcr(file_dat,nlines,molcrys%Cell,molcrys%Spg, A,CFrame,NPhase=IPhase)
                 else
                   call readn_set_xtal_pcr(file_dat,nlines,molcrys%Cell,molcrys%Spg, A,NPhase=IPhase)
                 end if
              else
                 if(present(CFrame)) then
                   call readn_set_xtal_pcr(file_dat,nlines,molcrys%Cell,molcrys%Spg,A,CFrame)
                 else
                   call readn_set_xtal_pcr(file_dat,nlines,molcrys%Cell,molcrys%Spg,A)
                 end if
              end if

           case("shx")
              if(present(CFrame)) then
                call readn_set_xtal_shx(file_dat,nlines,molcrys%Cell,molcrys%Spg,A,CFrame)
              else
                call readn_set_xtal_shx(file_dat,nlines,molcrys%Cell,molcrys%Spg,A)
              end if
           case default
              !---- CFL Format ----!
              if (present(Job_Info)) then
                 if (present(iphase)) then
                    if(present(CFrame)) then
                      call readn_set_xtal_cfl(file_dat,nlines,molcrys%Cell,molcrys%Spg,A,CFrame,NPhase=IPhase,Job_Info=Job_Info)
                    else
                      call readn_set_xtal_cfl(file_dat,nlines,molcrys%Cell,molcrys%Spg,A,NPhase=IPhase,Job_Info=Job_Info)
                    end if
                 else
                    if(present(CFrame)) then
                      call readn_set_xtal_cfl(file_dat,nlines,molcrys%Cell,molcrys%Spg,A,CFrame,Job_Info=Job_Info)
                    else
                      call readn_set_xtal_cfl(file_dat,nlines,molcrys%Cell,molcrys%Spg,A,Job_Info=Job_Info)
                    end if
                 end if
              else
                 if (present(iphase)) then
                    if(present(CFrame)) then
                      call readn_set_xtal_cfl(file_dat,nlines,molcrys%Cell,molcrys%Spg,A,CFrame,NPhase=IPhase)
                    else
                      call readn_set_xtal_cfl(file_dat,nlines,molcrys%Cell,molcrys%Spg,A,NPhase=IPhase)
                    end if
                 else
                    if(present(CFrame)) then
                      call readn_set_xtal_cfl(file_dat,nlines,molcrys%Cell,molcrys%Spg,A,CFrame)
                    else
                      call readn_set_xtal_cfl(file_dat,nlines,molcrys%Cell,molcrys%Spg,A)
                    end if
                 end if
              end if
              !---- Reading molecules ----!
              if (present(iphase)) then
                call readn_set_xtal_cfl_molec(file_dat,nlines,molcrys,NPhase=IPhase)
              else
                call readn_set_xtal_cfl_molec(file_dat,nlines,molcrys)
              end if

       end select
       if (err_form) return

       !---- Passing from Atom_List_Type -> Molcrys ----!
       molcrys%n_free=A%natoms
       if (A%natoms > 0) then
          if (allocated(molcrys%Atm)) deallocate(molcrys%Atm)
          allocate(molcrys%Atm(A%natoms))
          molcrys%Atm=A%Atom
       end if

       call deallocate_atom_list(A)

       !---- Testing if Xtal was defined ----!
       if (all(molcrys%cell%cell > 0.0)) then
          do i=1,molcrys%n_mol
             if (.not. molcrys%mol(i)%in_xtal) then
                 molcrys%mol(i)%in_xtal=.true.
             end if
          end do
       end if

       return
    End Subroutine Readn_Set_Xtal_Structure_Molcr



    Subroutine Readn_Set_Xtal_Structure_Magn(filenam,Cell,SpG,A,Mode,Iphase,Job_Info,file_list,CFrame)
       !---- Arguments ----!
       character(len=*),                 intent( in)     :: filenam
       Type (Crystal_Cell_Type),         intent(out)     :: Cell
       Type (Magnetic_Space_Group_Type), intent(out)     :: SpG
       Type (atom_list_type),            intent(out)     :: A
       Character(len=*),    optional,    intent( in)     :: Mode
       Integer,             optional,    intent( in)     :: Iphase
       Type(Job_Info_type), optional,    intent(out)     :: Job_Info
       Type(file_list_type),optional,    intent(in out)  :: file_list
       Character(len=*),    optional,    intent( in)     :: CFrame
       !
       character(len=132), allocatable, dimension(:) :: file_dat
       character(len=3)                              :: modec
       integer                                       :: nlines

       call init_err_form()

       nlines=0
       if (present(file_list)) nlines=file_list%nlines

       !---- Number of Lines in the input file ----!
       if(nlines == 0) then
           call Number_Lines(trim(filenam), nlines)
           if (nlines==0) then
              err_form=.true.
              ERR_Form_Mess="The file "//trim(filenam)//" contains nothing"
              return
           else
              if (allocated(file_dat)) deallocate( file_dat)
              allocate( file_dat(nlines))
              call reading_Lines(trim(filenam),nlines,file_dat)
           end if
           if (present(file_list)) then
              file_list%nlines=nlines
              if (allocated(file_list%line)) deallocate(file_list%line)
              allocate(file_list%line(nlines))
              file_list%line=file_dat
           end if
       else
           if (allocated(file_dat)) deallocate( file_dat)
           allocate( file_dat(nlines))
           file_dat=file_list%line
       end if

       !---- Define the type of file: CIF, CFL, RES,... ----!
       modec=" "
       if (present(mode)) modec=l_case(mode(1:3))

       select case(modec)
           case("cif")

              call Readn_Set_Magnetic_Structure_MCIF(filenam,Cell,Spg,A)

           case default
              !---- CFL Format ----!
              if (present(Job_Info)) then
                 if (present(iphase)) then
                    if(present(CFrame)) then
                      call Readn_Set_XTal_CFL_Shub(file_dat,nlines,Cell,SpG,A,CFrame,NPhase=IPhase,Job_Info=Job_Info)
                    else
                      call Readn_Set_XTal_CFL_Shub(file_dat,nlines,Cell,Spg,A,NPhase=IPhase,Job_Info=Job_Info)
                    end if
                 else
                    if(present(CFrame)) then
                      call Readn_Set_XTal_CFL_Shub(file_dat,nlines,Cell,Spg,A,CFrame,Job_Info=Job_Info)
                    else
                      call Readn_Set_XTal_CFL_Shub(file_dat,nlines,Cell,Spg,A,Job_Info=Job_Info)
                    end if
                 end if
              else
                 if (present(iphase)) then
                    if(present(CFrame)) then
                      call Readn_Set_XTal_CFL_Shub(file_dat,nlines,Cell,Spg,A,CFrame,NPhase=IPhase)
                    else
                      call Readn_Set_XTal_CFL_Shub(file_dat,nlines,Cell,Spg,A,NPhase=IPhase)
                    end if
                 else
                    if(present(CFrame)) then
                      call Readn_Set_XTal_CFL_Shub(file_dat,nlines,Cell,Spg,A,CFrame)
                    else
                      call Readn_Set_XTal_CFL_Shub(file_dat,nlines,Cell,Spg,A)
                    end if
                 end if
              end if

       end select
    End Subroutine Readn_Set_Xtal_Structure_Magn



    !!----
    !!---- Subroutine Get_Phases_File(filecode, Nphas, PhasesName,ILines)
    !!----    character(len=*),                intent(in)   :: filecode
    !!----    Integer,                         intent(out)  :: Nphas
    !!----    Character(len=80), dimension(:), intent(out)  :: PhasesName
    !!----    Integer,dimension(2,:),          intent(out)  :: ILines
    !!----
    !!---- Determine how many phases there are in a CIF or PCR file and
    !!---- give the lines to locate
    !!----
    !!---- Update: 01/05/2013
    !!
    Subroutine Get_Phases_File(filecode, NPhas, PhasesName,ILines)
       !---- Arguments ----!
       character(len=*),             intent(in)   :: filecode
       integer,                      intent(out)  :: Nphas
       character(len=*),dimension(:),intent(out)  :: PhasesName
       integer,dimension(:,:),       intent(out)  :: ILines

       !---- Local Variables ----!
       character(len=3) :: ext
       integer          :: npos

       !> Error
       call init_err_form()

       !> Init
       Nphas=0
       PhasesName=' '
       Ilines=0

       !> PCR or CIF file
       npos=index(filecode,'.',back=.true.)
       if (npos <=0) then
          err_form=.true.
          err_form_mess='No extension was found in the name of the file!'
          return
       end if

       ext=filecode(npos+1:)
       ext=u_case(ext)
       select case (ext)
          case ('CIF')
             call get_nphases_ciffile(filecode, NPhas, PhasesName,ILines)
          case ('PCR')
             call get_nphases_pcrfile(filecode, NPhas, PhasesName,ILines)
          case default
             err_form=.true.
             err_form_mess='Extension for this file not valid!'
       end select

       return
    End Subroutine Get_Phases_File

    !!--++
    !!--++ Subroutine Get_NPhases_CIFFile(Filecode,NPhas,PhasesName,ILines)
    !!--++    character(len=*),                 intent(in)  :: Filecode    ! Filename
    !!--++    integer,                          intent(out) :: NPhas       ! Number of Phases in the file
    !!--++    character(len=*), dimension(:),   intent(out) :: PhasesName     ! Name of Phases in the file
    !!--++    integer,          dimension(:,:), intent(out) :: ILines        ! Index for lines for each Phase
    !!--++
    !!--++ Determine the number of phases are included into the file
    !!--++
    !!--++ Date: 01/05/2013
    !!
    Subroutine Get_NPhases_CIFFile(Filecode,NPhas,PhasesName,ILines)
       !---- Arguments ----!
       character(len=*),                 intent(in)  :: Filecode    ! Filename
       integer,                          intent(out) :: NPhas       ! Number of Phases in the file
       character(len=*), dimension(:),   intent(out) :: PhasesName     ! Name of Phases in the file
       integer,          dimension(:,:), intent(out) :: ILines        ! Index for lines for each Phase

       !---- Local Variables ----!
       character(len=150), dimension(:), allocatable :: filen
       character(len=150)                            :: line
       integer                                       :: i,j,nl

       !> Error
       call init_err_form()

       !> Initialize
       NPhas=0
       PhasesName=' '
       ILines=0

       !> Reading file
       nl=0
       call number_lines(trim(filecode),nl)
       if (nl <=0) then
          err_form=.true.
          err_form_mess='No lines were read for '//trim(filecode)//' !!'
          return
       end if
       allocate(filen(nl))
       call reading_lines(trim(filecode),nl,filen)

       !> Number of Phases
       do i=1,nl
          line=adjustl(filen(i))

          !> empty line
          if (len_trim(line) <= 0) cycle

          !> comment line
          if (line(1:1) =='#') cycle

          !> No data_global
          j=index(line,'data_global')
          if (j > 0) cycle

          !> Just only lines beginning with data...
          j=index(line,'data_')
          if (j /= 1) cycle

          nphas=nphas+1
          ILines(1,Nphas)=i
          PhasesName(nphas)=trim(line(j+5:))
          if (nphas > 1) ILines(2,nphas-1)=i-1
       end do
       if (nphas > 0 .and. ILines(2,nphas)==0) ILines(2,nphas)=nl

       if (allocated(filen)) deallocate(filen)

       return
    End Subroutine  Get_NPhases_CIFFile

    !!--++
    !!--++ Subroutine Get_NPhases_PCRFile(filecode, Nphas,PhasesName,ILines)
    !!--++    character(len=*),                intent(in)   :: filecode
    !!--++    Integer,                         intent(out)  :: Nphas
    !!--++    Character(len=80), dimension(:), intent(out)  :: PhasesName
    !!--++    Integer,dimension(2,:),          intent(out)  :: ILines
    !!--++
    !!--++ Determine how many phases and where there in a PCR file
    !!--++
    !!--++ Update: 01/05/2013
    !!
    Subroutine Get_NPhases_PCRFile(filecode, NPhas, PhasesName,ILines)
       !---- Arguments ----!
       character(len=*),             intent(in)   :: filecode
       integer,                      intent(out)  :: Nphas
       character(len=*),dimension(:),intent(out)  :: PhasesName
       integer,dimension(:,:),       intent(out)  :: ILines

       !---- Local Variables ----!
       logical                                      :: multi, ask_phase
       character(len=80), dimension(:), allocatable :: file_dat
       character(len=80)                            :: line
       integer                                      :: i,k,iv,nlines
       integer, dimension(30)                       :: ivet
       real(kind=cp), dimension(30)                 :: vet

       !> Err
       call init_err_form()

       !> Init
       NPhas=0
       PhasesName=' '
       ILines=0

       !> Reading file
       nlines=0
       call number_lines(trim(filecode),nlines)
       if (nlines <=0) then
          err_form=.true.
          err_form_mess='No lines were read for '//trim(filecode)//' !!'
          return
       end if
       allocate(file_dat(nlines))
       call reading_lines(trim(filecode),nlines,file_dat)

       ILines(1,:)=1
       ILines(2,:)=nlines

       !> Simple / Multi format
       multi=.false.
       do i=1,nlines
          line=adjustl(file_dat(i))
          if (line(1:1) =='!' .or. line(1:1)==' ') cycle
          if (index(line,'NPATT ') <=0) cycle
          multi=.true.
       end do

       !> Number of Phases
       if (.not. multi) then
          do i=2,nlines
             line=adjustl(file_dat(i))
             if (line(1:1) =='!' .or. line(1:1)==' ') cycle
             call getnum(line,vet,ivet,iv)
             if (iv > 3) then
                NPhas=ivet(3)
                exit
             end if
          end do

       else
          do i=1,nlines
             line=adjustl(file_dat(i))
             if (line(1:4) /='!Nph') cycle

             line=adjustl(file_dat(i+1))
             call getnum(line,vet,ivet,iv)
             if (iv > 1) then
                NPhas=ivet(1)
                exit
             end if
          end do
       end if

       if (NPhas == 0) then
          err_form=.true.
          err_form_mess=" No Phase information was found in this PCR file. Please, check it! "
          return
       end if

       !> Locate where begin each Phase
       k=0
       ask_phase=.true.

       do i=1,nlines
          line=adjustl(file_dat(i))
          if (ask_phase) then
             if (index(line,'Data for PHASE') <= 0) cycle
          else
             if (line(1:1) /='!') then
                k=k+1
                ILines(1,k)=i
                PhasesName(k)=trim(adjustl(line))
                if (k == NPhas) exit

                ask_phase=.true.
             end if
             cycle
          end if
          ask_phase=.false.
       end do

       if (NPhas /= k) then
          err_form=.true.
          err_form_mess=" Locating Phases failed in this PCR. Please, check it!"
          return
       end if

       do i=1,Nphas
          if (nphas > 1) then
             ilines(2,i)=ilines(1,i+1)-1
          end if
       end do

       return
    End Subroutine Get_NPhases_PCRFile


    !!----
    !!---- Subroutine Write_CFL(lun,Molx,comment)
    !!----    integer,                       intent(in) :: lun
    !!----    type (Molecular_Crystal_Type), intent(in) :: Molx
    !!----    character(len=*),optional,     intent(in) :: comment
    !!----
    !!----    (OVERLOADED)
    !!----
    !!----    Write a CFL-file with molecular_crystal_type
    !!----
    !!---- Update: July - 2014
    !!
    Subroutine Write_CFL_Molcrys(lun,Molx,comment)
       !---- Arguments ----!
       integer,                       intent(in) :: lun
       type (Molecular_Crystal_Type), intent(in) :: Molx
       character(len=*),optional,     intent(in) :: comment

       !----- Local variables -----!
       integer                         :: j !,loc
       real(kind=cp), dimension(6)     :: a,sa
       character(len=30), dimension(6) :: text

       if(present(comment)) write(unit=lun,fmt="(a)") "TITLE "//trim(comment)
       write(unit=lun,fmt="(a)") "!  Automatically generated CFL file (Write_CFL)"

       a(1:3)=molx%cell%Cell
       a(4:6)=molx%cell%ang
       sa(1:3)=molx%cell%Cell_std
       sa(4:6)=molx%cell%ang_std
       do j=1,6
          call SetNum_Std(a(j), sa(j), text(j))
       end do
       write(unit=lun,fmt="(a)") "!         a               b               c            alpha           beta            gamma"
       write(unit=lun,fmt="(a,6a16)") "Cell ",text
       write(unit=lun,fmt="(a,i3)")"!     Space Group # ",molx%spg%NumSpg
       write(unit=lun,fmt="(a,a)") "Spgr  ",molx%spg%SPG_Symb
       call Write_Atoms_CFL(Molx,Lun)

       return
    End Subroutine Write_CFL_Molcrys



    !!----
    !!---- Subroutine Write_Atoms_CFL(Ats,Lun,Cell)
    !!----    Type (atom_list_type),dimension(:),  intent(in) :: Ats     !  In -> Atom List
    !!----    integer, optional,                   intent(in) :: lun     !  In -> Unit to write
    !!----    Type(Crystal_Cell_Type), optional,   intent(in) :: Cell    !  In -> Transform to thermal parameters
    !!----
    !!----    Write the atoms in the asymmetric unit for a CFL file
    !!----
    !!---- Update: February - 2003
    !!
    Subroutine Write_Atoms_CFL_MOLX(Molx,Lun)
        !---- Arguments ----!
        type (Molecular_Crystal_Type), intent(in) :: Molx
        integer, optional,             intent(in) :: Lun

        !---- Local Variables ----!
        character(len=30),dimension(6) :: text
        character(len=36)              :: forma,fom
        integer                        :: i, j, iunit, leng, maxl,ish
        real(kind=cp), dimension(6)    :: u,bet,sb

        iunit=6
        if (present(lun)) iunit=lun

        if(molx%n_free > 0) then
            !Determine the maximum length of the atom labels
            maxl=0
            do i=1,molx%n_free
                leng=len_trim(molx%atm(i)%lab)
                if(leng > maxl) maxl=leng
            end do
            maxl=max(maxl,4)+1
            ish=maxl-4
            fom   ="(a,tr  ,a)"
            Select Case(ish)
                Case(:9)
                    write(unit=fom(6:6),fmt="(i1)") ish
                Case(10:)
                    write(unit=fom(6:7),fmt="(i2)") ish
            End Select
            forma="(a,a  ,tr2,a,tr3,5a14,2f8.2,tr3,a)"
            Select Case(maxl)
                Case(:9)
                    write(unit=forma(5:5),fmt="(i1)") maxl
                Case(10:)
                    write(unit=forma(5:6),fmt="(i2)") maxl
            End Select
            write (unit=iunit,fmt=fom) "!     ", &
                  "Atom  Type     x/a           y/b           z/c           Biso          Occ           Spin    Charge    Info"
            do i=1,molx%n_free

                do j=1,3
                   call SetNum_Std(molx%atm(i)%x(j), molx%atm(i)%x_std(j), text(j))
                end do
                call SetNum_Std(molx%atm(i)%Biso, molx%atm(i)%Biso_std, text(4))
                call SetNum_Std(molx%atm(i)%Occ, molx%atm(i)%Occ_std, text(5))

                write (unit=iunit,fmt=forma) &
                      "Atom   ",trim(molx%atm(i)%lab),molx%atm(i)%chemsymb, (text(j),j=1,5), &
                       molx%atm(i)%moment,molx%atm(i)%charge,"# "//molx%atm(i)%AtmInfo

                if (molx%atm(i)%thtype == "aniso") then

                    if (molx%atm(i)%utype == "beta") then
                        bet=molx%atm(i)%u(1:6)
                        sb=molx%atm(i)%u_std(1:6)
                        do j=1,6
                            call SetNum_Std(bet(j), sb(j), text(j))
                        end do
                        write (unit=iunit,fmt="(a,tr1,6a14)") "Beta  ", text
                        u=convert_betas_u(bet,molx%cell)
                        sb=convert_betas_u(molx%atm(i)%u_std,molx%cell)
                        do j=1,6
                            call SetNum_Std(u(j), sb(j), text(j))
                        end do
                        write(unit=iunit,fmt="(a,6a14)") "!U_ij  ", text
                    else if(molx%atm(i)%thtype == "u_ij") then
                        u=molx%atm(i)%u(1:6)
                        sb=molx%atm(i)%u_std(1:6)
                        do j=1,6
                            call SetNum_Std(u(j), sb(j), text(j))
                        end do
                        write(unit=iunit,fmt="(a,6a14)") "U_ij  ", text
                        bet=convert_u_betas(u,molx%cell)
                        sb=convert_u_betas(molx%atm(i)%u_std,molx%cell)
                        do j=1,6
                            call SetNum_Std(bet(j), sb(j), text(j))
                        end do
                        write(unit=iunit,fmt="(a,6a14)") "!Beta  ", text
                    end if
                end if
            end do ! i=1,molx%n_free
        end if ! molx%n_free > 0

        if (molx%n_mol > 0) then
            do i=1,molx%n_mol
                write(unit=iunit,fmt="(/,a,tr2,i3,tr2,a,tr2,a)") &
                     "MOLEX",molx%mol(i)%natoms,trim(molx%mol(i)%Name_mol),molx%mol(i)%coor_type
                write(unit=iunit,fmt="(a)") &
                     "!    Xc         Yc          Zc        Phi        Theta      Chi     TypeAngles TypeThermal"
                write(unit=iunit,fmt="(6f11.5,tr6,a,tr10,a)") &
                     molx%mol(i)%xcentre,molx%mol(i)%orient,molx%mol(i)%rot_type,molx%mol(i)%therm_type
                write(unit=iunit,fmt="(t1,6i10,tr2,a)") &
                     molx%mol(i)%lxcentre,molx%mol(i)%lorient," ! Refinemencodes"

                select case (molx%mol(i)%coor_type)
                    case ("C","c")
                        write(unit=iunit,fmt="(a)") &
                        "!Atom   Type        XC          YC          ZC    N1  N2  N3      Biso        Occ "
                    case ("F","f")
                        write(unit=iunit,fmt="(a)") &
                        "!Atom   Type        X           Y           Z     N1  N2  N3      Biso        Occ "
                    case ("S","s")
                        write(unit=iunit,fmt="(a)") &
                        "!Atom   Type    distance      Theta       Phi     N1  N2  N3      Biso        Occ "
                    case ("Z","z")
                        write(unit=iunit,fmt="(a)") &
                        "!Atom   Type    distance  Bond-Angle Torsion-Ang  N1  N2  N3      Biso        Occ "
                    case default
                        write(unit=iunit,fmt="(a)") &
                        "!Atom   Type      Coor1       Coor2       Coor3   N1  N2  N3      Biso        Occ "
                end select ! molx%mol(i)%coor_type

                do j=1,molx%mol(i)%natoms
                    write(unit=iunit,fmt="(a,tr2,a,3f12.5,3i4,2f12.5)")  &
                          molx%mol(i)%AtName(j), molx%mol(i)%AtSymb(j),molx%mol(i)%I_Coor(:,j),  &
                          molx%mol(i)%Conn(:,j), molx%mol(i)%Biso(j),  molx%mol(i)%Occ(j)
                end do ! j = molx%mol(i)%natoms
            end do ! i = 1,molx%n_mol
        end if ! molx%n_mol > 0
        return
    End Subroutine Write_Atoms_CFL_MOLX
    !!----
    !!---- Subroutine Write_Atoms_CFL(Ats,Lun,Cell)
    !!----    Type (atom_list_type),dimension(:),  intent(in) :: Ats     !  In -> Atom List
    !!----    integer, optional,                   intent(in) :: lun     !  In -> Unit to write
    !!----    Type(Crystal_Cell_Type), optional,   intent(in) :: Cell    !  In -> Transform to thermal parameters
    !!----
    !!----    Write the atoms in the asymmetric unit for a CFL file
    !!----
    !!---- Update: February - 2003
    !!
    Subroutine Write_Atoms_CFL_MOLX_orig(Molx,Lun)
       !---- Arguments ----!
       type (Molecular_Crystal_Type), intent(in) :: Molx
       integer, optional,             intent(in) :: Lun

       !---- Local Variables ----!
       character(len=30),dimension(6) :: text
       character(len=36)              :: forma,fom
       integer                        :: i, j, iunit, leng, maxl,ish
       real(kind=cp), dimension(6)    :: u,bet,sb

       iunit=6
       if (present(lun)) iunit=lun

       if(molx%n_free == 0) then
         write (unit=iunit,fmt="(a)") "!  No atoms ..."
         return
       end if
       !Determine the maximum length of the atom labels
       maxl=0
       do i=1,molx%n_free
         leng=len_trim(molx%atm(i)%lab)
         if(leng > maxl) maxl=leng
       end do
       maxl=max(maxl,4)+1
       ish=maxl-4
       fom   ="(a,tr  ,a)"
       Select Case(ish)
          Case(:9)
            write(unit=fom(6:6),fmt="(i1)") ish
          Case(10:)
            write(unit=fom(6:7),fmt="(i2)") ish
       End Select
       forma="(a,a  ,tr2,a,tr3,5a14,2f8.2,tr3,a)"
       Select Case(maxl)
         Case(:9)
             write(unit=forma(5:5),fmt="(i1)") maxl
         Case(10:)
             write(unit=forma(5:6),fmt="(i2)") maxl
       End Select
       write (unit=iunit,fmt=fom) "!     ", &
             "Atom  Type     x/a           y/b           z/c           Biso          Occ           Spin    Charge    Info"
       do i=1,molx%n_free

          do j=1,3
             call SetNum_Std(molx%atm(i)%x(j), molx%atm(i)%x_std(j), text(j))
          end do
          call SetNum_Std(molx%atm(i)%Biso, molx%atm(i)%Biso_std, text(4))
          call SetNum_Std(molx%atm(i)%Occ, molx%atm(i)%Occ_std, text(5))

          write (unit=iunit,fmt=forma) &
                "Atom   ",trim(molx%atm(i)%lab),molx%atm(i)%chemsymb, (text(j),j=1,5), &
                 molx%atm(i)%moment,molx%atm(i)%charge,"# "//molx%atm(i)%AtmInfo

          if (molx%atm(i)%thtype == "aniso") then

             if (molx%atm(i)%utype == "beta") then
                bet=molx%atm(i)%u(1:6)
                sb=molx%atm(i)%u_std(1:6)
                do j=1,6
                   call SetNum_Std(bet(j), sb(j), text(j))
                end do
                write (unit=iunit,fmt="(a,tr1,6a14)") "Beta  ", text
                u=convert_betas_u(bet,molx%cell)
                sb=convert_betas_u(molx%atm(i)%u_std,molx%cell)
                do j=1,6
                    call SetNum_Std(u(j), sb(j), text(j))
                end do
                write(unit=iunit,fmt="(a,6a14)") "!U_ij  ", text
             else if(molx%atm(i)%thtype == "u_ij") then
                u=molx%atm(i)%u(1:6)
                sb=molx%atm(i)%u_std(1:6)
                do j=1,6
                   call SetNum_Std(u(j), sb(j), text(j))
                end do
                write(unit=iunit,fmt="(a,6a14)") "U_ij  ", text
                bet=convert_u_betas(u,molx%cell)
                sb=convert_u_betas(molx%atm(i)%u_std,molx%cell)
                do j=1,6
                    call SetNum_Std(bet(j), sb(j), text(j))
                end do
                write(unit=iunit,fmt="(a,6a14)") "!Beta  ", text
             end if
          end if
       end do

       return
    End Subroutine Write_Atoms_CFL_MOLX_orig

    Subroutine Write_MCIF(Ipr,mCell,MSGp,Am,Cell)
       Integer,                         intent(in)           :: Ipr
       type(Magnetic_Space_Group_Type), intent(in)           :: MSGp
       type(Crystal_Cell_Type),         intent(in)           :: mCell
       type(Atom_List_Type),            intent(in)           :: Am
       type(Crystal_Cell_Type),optional,intent(in)           :: Cell
       !
       Character(len=132)             :: line
       character(len=80),dimension(6) :: text
       character(len=2)               :: invc
       real(kind=cp)                  :: occ,occ_std,uiso,uiso_std
       integer :: i,j

       write(unit=Ipr,fmt="(a)") "#  --------------------------------------"
       write(unit=Ipr,fmt="(a)") "#  Magnetic CIF file generated by CrysFML"
       write(unit=Ipr,fmt="(a)") "#  --------------------------------------"
       write(unit=Ipr,fmt="(a)") "# https://forge.epn-campus.eu/projects/crysfml/repository"
       call Write_Date_Time(dtim=line)
       write(unit=Ipr,fmt="(a)") trim(line)
       write(unit=Ipr,fmt="(a)") " "

       write(unit=Ipr,fmt="(a)") "data_"
       write(unit=Ipr,fmt="(a)") "_citation_journal_abbrev ?"
       write(unit=Ipr,fmt="(a)") "_citation_journal_volume ?"
       write(unit=Ipr,fmt="(a)") "_citation_page_first     ?"
       write(unit=Ipr,fmt="(a)") "_citation_page_last      ?"
       write(unit=Ipr,fmt="(a)") "_citation_article_id     ?"
       write(unit=Ipr,fmt="(a)") "_citation_year           ?"
       write(unit=Ipr,fmt="(a)") "_loop "
       write(unit=Ipr,fmt="(a)") "_citation_author_name"
       write(unit=Ipr,fmt="(a)") "?"
       write(unit=Ipr,fmt="(a)")
       write(unit=Ipr,fmt="(a)") "_atomic_positions_source_database_code_ICSD  ?"
       write(unit=Ipr,fmt="(a)") "_atomic_positions_source_other    .  "
       write(unit=Ipr,fmt="(a)")
       write(unit=Ipr,fmt="(a)") "_Neel_temperature  ?"
       write(unit=Ipr,fmt="(a)") "_magn_diffrn_temperature  ?"
       write(unit=Ipr,fmt="(a)") "_exptl_crystal_magnetic_properties_details"
       write(unit=Ipr,fmt="(a)") ";"
       write(unit=Ipr,fmt="(a)") ";"
       write(unit=Ipr,fmt="(a)") "_active_magnetic_irreps_details"
       write(unit=Ipr,fmt="(a)") ";"
       write(unit=Ipr,fmt="(a)") ";"
       write(unit=Ipr,fmt="(a)") " "
       if(MSGp%standard_setting) then
          write(unit=Ipr,fmt="(a)") "_magnetic_space_group_standard_setting  'yes'"
       else
          write(unit=Ipr,fmt="(a)") "_magnetic_space_group_standard_setting  'no'"
       end if
       write(unit=Ipr,fmt="(a)")    '_parent_space_group.name_H-M  "'//trim(MSGp%Parent_spg)//'"'
       write(unit=Ipr,fmt="(a,i3)") "_parent_space_group.IT_number  ",MSGp%Parent_num
       write(unit=Ipr,fmt="(a)")    "_magnetic_space_group.transform_from_parent_Pp_abc  '"//trim(MSGp%trn_from_parent)//"'"
       write(unit=Ipr,fmt="(a)")    "_magnetic_space_group.transform_to_standard_Pp_abc  '"//trim(MSGp%trn_to_standard)//"'"
       write(unit=Ipr,fmt="(a)")
       if(len_trim(MSGp%BNS_number) /= 0) &
       write(unit=Ipr,fmt="(a)") "_space_group.magn_number_BNS  "//trim(MSGp%BNS_number)
       if(len_trim(MSGp%BNS_symbol) /= 0) &
       write(unit=Ipr,fmt="(a)") '_space_group.magn_name_BNS  "'//trim(MSGp%BNS_symbol)//'"'
       if(len_trim(MSGp%OG_number) /= 0) &
       write(unit=Ipr,fmt="(a)") '_space_group.magn_number_OG '//trim(MSGp%OG_number)
       if(len_trim(MSGp%OG_symbol) /= 0) &
       write(unit=Ipr,fmt="(a)") '_space_group.magn_name_OG  "'//trim(MSGp%OG_symbol)//'"'
       write(unit=Ipr,fmt="(a)")

       if(MSGp%n_irreps /= 0) then
          write(unit=Ipr,fmt="(a)") "loop_"
          write(unit=Ipr,fmt="(a)") "_irrep_id"
          write(unit=Ipr,fmt="(a)") "_irrep_dimension"
          if( any(MSGp%small_irrep_dim > 0) ) write(unit=Ipr,fmt="(a)") "_small_irrep_dimension"
          write(unit=Ipr,fmt="(a)") "_irrep_direction_type"
          write(unit=Ipr,fmt="(a)") "_irrep_action"
          if( any(MSGp%irrep_modes_number > 0) ) write(unit=Ipr,fmt="(a)") "_irrep_modes_number"
          do i=1,MSGp%n_irreps
            if(MSGp%small_irrep_dim(i) > 0) then
               write(unit=line,fmt=("(2i4)"))  MSGp%irrep_dim(i), MSGp%small_irrep_dim(i)
            else
               write(unit=line,fmt=("(i4)"))  MSGp%irrep_dim(i)
            end if
            line= trim(MSGp%irrep_id(i))//"  "//trim(line)//"   "// &
                                      trim(MSGp%irrep_direction(i))//"  "//trim(MSGp%irrep_action(i))
            if( MSGp%irrep_modes_number(i) > 0) then
               j=len_trim(line)
              write(unit=line(j+1:),fmt="(i4)") MSGp%irrep_modes_number(i)
            end if
            write(unit=Ipr,fmt="(a)") trim(line)
          end do
          write(unit=Ipr,fmt="(a)")
       else
          write(unit=Ipr,fmt="(a)") "loop_"
          write(unit=Ipr,fmt="(a)") "_irrep_id"
          write(unit=Ipr,fmt="(a)") "_irrep_dimension"
          write(unit=Ipr,fmt="(a)") "_small_irrep_dimension"
          write(unit=Ipr,fmt="(a)") "_irrep_direction_type"
          write(unit=Ipr,fmt="(a)") "_irrep_action"
          write(unit=Ipr,fmt="(a)") "_irrep_modes_number"
          write(unit=Ipr,fmt="(a)") " ?  ?  ?  ?  ?  ?"
          write(unit=Ipr,fmt="(a)")
       end if

       if(MSGp%m_cell) then
          do i=1,3
            call setnum_std(mCell%Cell(i),mCell%cell_std(i),text(i))
            call setnum_std(mCell%ang(i),mCell%ang_std(i),text(i+3))
          end do
          write(unit=Ipr,fmt="(a)") "_cell_length_a    "//trim(text(1))
          write(unit=Ipr,fmt="(a)") "_cell_length_b    "//trim(text(2))
          write(unit=Ipr,fmt="(a)") "_cell_length_c    "//trim(text(3))
          write(unit=Ipr,fmt="(a)") "_cell_angle_alpha "//trim(text(4))
          write(unit=Ipr,fmt="(a)") "_cell_angle_beta  "//trim(text(5))
          write(unit=Ipr,fmt="(a)") "_cell_angle_gamma "//trim(text(6))
          write(unit=Ipr,fmt="(a)")
       else
          if(present(Cell)) then
             do i=1,3
               call setnum_std(Cell%Cell(i),Cell%cell_std(i),text(i))
               call setnum_std(Cell%ang(i),Cell%ang_std(i),text(i+3))
             end do
             write(unit=Ipr,fmt="(a)") "_cell_length_a    "//trim(text(1))
             write(unit=Ipr,fmt="(a)") "_cell_length_b    "//trim(text(2))
             write(unit=Ipr,fmt="(a)") "_cell_length_c    "//trim(text(3))
             write(unit=Ipr,fmt="(a)") "_cell_angle_alpha "//trim(text(4))
             write(unit=Ipr,fmt="(a)") "_cell_angle_beta  "//trim(text(5))
             write(unit=Ipr,fmt="(a)") "_cell_angle_gamma "//trim(text(6))
             write(unit=Ipr,fmt="(a)")
          end if
       end if
       if(MSGp%n_kv > 0) then
          write(unit=Ipr,fmt="(a)") "loop_"
          write(unit=Ipr,fmt="(a)") "_magnetic_propagation_vector_seq_id"
          write(unit=Ipr,fmt="(a)") "_magnetic_propagation_vector_kxkykz"
          do i=1,MSGp%n_kv
            call Frac_Trans_2Dig(MSGp%kv(:,i),line)
            line=adjustl(line(2:len_trim(line)-1))
            write(unit=Ipr,fmt="(a)") trim(MSGp%kv_label(i))//"  '"//trim(line)//"'"
          end do
       end if
       if(MSGp%m_constr) then
          write(unit=Ipr,fmt="(a)")
          write(unit=Ipr,fmt="(a)") "loop_"
          write(unit=Ipr,fmt="(a)") "_magnetic_atom_site_moment_symmetry_constraints_label"
          write(unit=Ipr,fmt="(a)") "_atom_site_magnetic_moment_symmetry_constraints_mxmymz"
          do i=1,Am%natoms
            line=Am%Atom(i)%AtmInfo
            if(len_trim(line) < 8) cycle
            write(unit=Ipr,fmt="(a)")trim(line)
          end do
       end if
       write(unit=Ipr,fmt="(a)")
       write(unit=Ipr,fmt="(a)")  "loop_"
       write(unit=Ipr,fmt="(a)")  "_space_group_magn_symop_operation.id"
       write(unit=Ipr,fmt="(a)")  "_space_group_magn_symop_operation.xyz"
       write(unit=Ipr,fmt="(a)")  "_space_group_magn_symop_operation.mxmymz"
       do i=1,MSGp%Multip            !New mCIF format
          write(unit=invc,fmt="(i2)") nint(MSgp%MSymop(i)%Phas)
          if(invc(1:1) == " ") invc(1:1)="+"
          write(unit=Ipr,fmt="(i3,a)") i," "//trim(MSgp%SymopSymb(i))//","//invc//" "//trim(MSgp%MSymopSymb(i))
       end do
       write(unit=Ipr,fmt="(a)")
       write(unit=Ipr,fmt="(a)") "loop_"
       write(unit=Ipr,fmt="(a)") "_atom_site_label"
       write(unit=Ipr,fmt="(a)") "_atom_site_type_symbol"
       write(unit=Ipr,fmt="(a)") "_atom_site_fract_x"
       write(unit=Ipr,fmt="(a)") "_atom_site_fract_y"
       write(unit=Ipr,fmt="(a)") "_atom_site_fract_z"
       write(unit=Ipr,fmt="(a)") "_atom_site_U_iso_or_equiv"
       write(unit=Ipr,fmt="(a)") "_atom_site_occupancy"
       write(unit=Ipr,fmt="(a)") "_atom_site_symmetry_multiplicity"
       write(unit=Ipr,fmt="(a)") "_atom_site_Wyckoff_label"
       line=" "
       do i=1,Am%natoms
          do j=1,3
            call setnum_std(Am%atom(i)%x(j),Am%atom(i)%x_std(j),text(j))
          end do
          occ=real(MSgp%Multip)/real(Am%atom(i)%Mult)*Am%atom(i)%occ
          occ_std=real(MSgp%Multip)/real(Am%atom(i)%Mult)*Am%atom(i)%occ_std
          call setnum_std(occ,occ_std,text(5))
          uiso=Am%atom(i)%biso/78.95683521
          uiso_std=Am%atom(i)%biso_std/78.95683521
          call setnum_std(uiso,uiso_std,text(4))
          write(unit=Ipr,fmt="(a6,a6,3a13,2a11,i4,a)") Am%Atom(i)%lab, Am%atom(i)%SfacSymb,(text(j),j=1,5),&
                                                       Am%atom(i)%Mult," "//Am%atom(i)%wyck
       end do
       write(unit=Ipr,fmt="(a)")
       write(unit=Ipr,fmt="(a)") "loop_"
       write(unit=Ipr,fmt="(a)") "_atom_site_moment_label"
       write(unit=Ipr,fmt="(a)") "_atom_site_moment_crystalaxis_x"
       write(unit=Ipr,fmt="(a)") "_atom_site_moment_crystalaxis_y"
       write(unit=Ipr,fmt="(a)") "_atom_site_moment_crystalaxis_z"
       do i=1,Am%natoms
          !if(sum(abs(Am%Atom(i)%Skr(:,1))) < 0.0001) cycle
          if(Am%Atom(i)%moment < 0.01) cycle
          do j=1,3
            call setnum_std(Am%atom(i)%M_xyz(j),Am%atom(i)%sM_xyz(j),text(j))
          end do
          write(unit=Ipr,fmt="(a8,3a12)") Am%Atom(i)%lab,(text(j),j=1,3)
       end do
       write(unit=Ipr,fmt="(a)")
       return
    End Subroutine Write_MCIF

    Subroutine Check_Symmetry_Constraints(SpG,Atm)
       class(SpG_Type),   intent(in)     :: SpG
       type(AtList_Type), intent(in out) :: Atm

       !--- Local variables ---!
       integer :: i,codini
       real(kind=cp), dimension(3)   :: codes
       real(kind=cp), dimension(6,8) :: codeT
       codini=1
       do i=1,Atm%natoms
          if(Atm%Atom(i)%Magnetic) then
            codes=1.0
            call Get_moment_ctr(Atm%Atom(i)%x,Atm%Atom(i)%moment,SpG,codini,codes)
          end if
       end do
       Select Type (SpG)
         type is (SuperSpaceGroup_Type)
            do i=1,Atm%natoms
              Select Type(at => Atm%Atom(i))
                class is (MAtm_Std_Type)
                  if(at%n_mc > 0) then
                    codeT=1.0
                    call Get_TFourier_ctr(at%x,at%Mcs(:,1:at%n_mc),codeT(:,1:at%n_mc),SpG,codini,"M")
                  end if
                  if(at%n_dc > 0) then
                    codeT=1.0
                    call Get_TFourier_ctr(at%x,at%Dcs(:,1:at%n_dc),codeT(:,1:at%n_dc),SpG,codini,"D")
                  end if
              End Select
            end do
       End Select
       Atm%symm_checked=.true.
     End Subroutine Check_Symmetry_Constraints

 End Module CFML_IO_Formats

