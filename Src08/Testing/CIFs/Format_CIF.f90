!!----
!!----
!!----
  Module test_CIF
    !---- Use modules ----!
    Use CFML_GlobalDeps,        only: CP, PI, EPS, TAB, Err_CFML, Clear_Error
    Use CFML_Strings,           only: l_case, u_case, get_num, cut_string, get_words, &
                                      get_numstd, Read_Key_Value, Read_Key_ValueSTD,  &
                                      String_NumStd, Number_Lines, Reading_Lines, FindFMT, &
                                      Init_FindFMT, String_Array_Type, File_type
    Use CFML_Atoms,             only: Atm_Type, Atm_Std_Type, Matm_std_type, Atm_Ref_Type, &
                                      AtList_Type, Allocate_Atom_List,Init_Atom_Type
    Use CFML_Metrics,           only: Cell_Type, Cell_G_Type, Set_Crystal_Cell, U_equiv, &
                                      get_U_from_Betas, get_Betas_from_U, get_Betas_from_B
    Use CFML_gSpaceGroups,      only: SpG_Type, SuperSpaceGroup_Type, Kvect_Info_Type,   &
                                      Change_Setting_SpaceG, Set_SpaceGroup, Get_Multip_Pos,&
                                      Get_Orbit, Get_Moment_Ctr, Get_TFourier_Ctr
    Use CFML_Maths,             only: Get_Eps_Math

    Use CFML_Rational

   public

   Contains


   !!----
   !!---- Write_Cif_Template
   !!----    Write a Cif File
   !!----
   !!---- 28/06/2019
   !!
   Subroutine Write_Cif_Template(filename, Cell, SpG, At_list, Type_data, Code)
      !---- Arguments ----!
      character(len=*),           intent(in) :: filename     ! Filename
      class(Cell_G_Type),         intent(in) :: Cell         ! Cell parameters
      class(SpG_Type),            intent(in) :: SpG          ! Space group information
      Type (AtList_Type),         intent(in) :: At_List      ! Atoms
      integer,                    intent(in) :: Type_data    ! 0,2:Single crystal diffraction; 1:Powder
      character(len=*),           intent(in) :: Code         ! Code or name of the structure

      !---- Local Variables ----!
      logical                                 :: info, aniso
      character(len=1), parameter             :: QMARK='?'
      character(len=:), allocatable           :: line
      character(len=:), allocatable           :: comm,adptyp
      character(len=30),dimension(6)          :: text
      real(kind=cp)                           :: u,su, ocf
      real(kind=cp), dimension(6)             :: Ua,sua,aux
      real(kind=cp), dimension(At_List%Natoms):: occup,soccup
      integer                                 :: iunit,i, j
      type(Atm_Std_Type),dimension(At_list%natoms) :: Atm
      !> Init
      info=.false.
      iunit=0

      !> Is this file opened?
      inquire(file=trim(filename),opened=info)
      if (info) then
         inquire(file=trim(filename),number=iunit)
         close(unit=iunit)
      end if

      !> Writing
      open(newunit=iunit, file=trim(filename),status="unknown",action="write")
      rewind(unit=iunit)

      !> Head Information
      select case (type_data)
         case (0:1)
            write(unit=iunit,fmt="(a)") "##############################################################################"
            write(unit=iunit,fmt="(a)") "###    CIF submission form for molecular structure report (Acta Cryst. C)  ###"
            write(unit=iunit,fmt="(a)") "##############################################################################"
            write(unit=iunit,fmt="(a)") " "
            write(unit=iunit,fmt="(a)") "#============================================================================="
            write(unit=iunit,fmt="(a)") "data_global"
            write(unit=iunit,fmt="(a)") "#============================================================================="
            write(unit=iunit,fmt="(a)") " "

         case (2:)
            write(unit=iunit,fmt="(a)") "##################################################################"
            write(unit=iunit,fmt="(a)") "###    CIF file from CrysFML, contains only structural data    ###"
            write(unit=iunit,fmt="(a)") "##################################################################"
      end select

      !> Processing Summary
      if (type_data < 2) then
         write(unit=iunit,fmt="(a)") "# PROCESSING SUMMARY (IUCr Office Use Only)"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "_journal_data_validation_number      ?"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "_journal_date_recd_electronic        ?"
         write(unit=iunit,fmt="(a)") "_journal_date_to_coeditor            ?"
         write(unit=iunit,fmt="(a)") "_journal_date_from_coeditor          ?"
         write(unit=iunit,fmt="(a)") "_journal_date_accepted               ?"
         write(unit=iunit,fmt="(a)") "_journal_date_printers_first         ?"
         write(unit=iunit,fmt="(a)") "_journal_date_printers_final         ?"
         write(unit=iunit,fmt="(a)") "_journal_date_proofs_out             ?"
         write(unit=iunit,fmt="(a)") "_journal_date_proofs_in              ?"
         write(unit=iunit,fmt="(a)") "_journal_coeditor_name               ?"
         write(unit=iunit,fmt="(a)") "_journal_coeditor_code               ?"
         write(unit=iunit,fmt="(a)") "_journal_coeditor_notes"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_journal_techeditor_code             ?"
         write(unit=iunit,fmt="(a)") "_journal_techeditor_notes"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_journal_coden_ASTM                  ?"
         write(unit=iunit,fmt="(a)") "_journal_name_full                   ?"
         write(unit=iunit,fmt="(a)") "_journal_year                        ?"
         write(unit=iunit,fmt="(a)") "_journal_volume                      ?"
         write(unit=iunit,fmt="(a)") "_journal_issue                       ?"
         write(unit=iunit,fmt="(a)") "_journal_page_first                  ?"
         write(unit=iunit,fmt="(a)") "_journal_page_last                   ?"
         write(unit=iunit,fmt="(a)") "_journal_paper_category              ?"
         write(unit=iunit,fmt="(a)") "_journal_suppl_publ_number           ?"
         write(unit=iunit,fmt="(a)") "_journal_suppl_publ_pages            ?"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "#============================================================================="
         write(unit=iunit,fmt="(a)") " "

         !> Submission details
         write(unit=iunit,fmt="(a)") "# 1. SUBMISSION DETAILS"
         write(unit=iunit,fmt="(a)") " "

         write(unit=iunit,fmt="(a)") "_publ_contact_author_name            ?   # Name of author for correspondence"
         write(unit=iunit,fmt="(a)") "_publ_contact_author_address             # Address of author for correspondence"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_publ_contact_author_email           ?"
         write(unit=iunit,fmt="(a)") "_publ_contact_author_fax             ?"
         write(unit=iunit,fmt="(a)") "_publ_contact_author_phone           ?"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "_publ_contact_letter"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "_publ_requested_journal              ?"
         write(unit=iunit,fmt="(a)") "_publ_requested_coeditor_name        ?"
         write(unit=iunit,fmt="(a)") "_publ_requested_category             ?   # Acta C: one of CI/CM/CO/FI/FM/FO"

         write(unit=iunit,fmt="(a)") "#=============================================================================="
         write(unit=iunit,fmt="(a)") " "

         !> Title  and Author List
         write(unit=iunit,fmt="(a)") "# 3. TITLE AND AUTHOR LIST"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "_publ_section_title"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_publ_section_title_footnote"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") ";"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "# The loop structure below should contain the names and addresses of all "
         write(unit=iunit,fmt="(a)") "# authors, in the required order of publication. Repeat as necessary."

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "loop_"
         write(unit=iunit,fmt="(a)") "    _publ_author_name"
         write(unit=iunit,fmt="(a)") "    _publ_author_footnote"
         write(unit=iunit,fmt="(a)") "    _publ_author_address"
         write(unit=iunit,fmt="(a)") "?                                   #<--'Last name, first name' "
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "#============================================================================="
         write(unit=iunit,fmt="(a)") " "

         !> Text
         write(unit=iunit,fmt="(a)") "# 4. TEXT"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "_publ_section_synopsis"
         write(unit=iunit,fmt="(a)") ";  ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_publ_section_abstract"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";          "
         write(unit=iunit,fmt="(a)") "_publ_section_comment"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_publ_section_exptl_prep      # Details of the preparation of the sample(s)"
         write(unit=iunit,fmt="(a)") "                              # should be given here. "
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_publ_section_exptl_refinement"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_publ_section_references"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_publ_section_figure_captions"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_publ_section_acknowledgements"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "#============================================================================="
         write(unit=iunit,fmt="(a)") " "

         !> Identifier
         write(unit=iunit,fmt="(a)") "#============================================================================="
         write(unit=iunit,fmt="(a)") "# If more than one structure is reported, the remaining sections should be "
         write(unit=iunit,fmt="(a)") "# completed per structure. For each data set, replace the '?' in the"
         write(unit=iunit,fmt="(a)") "# data_? line below by a unique identifier."
      end if !type_data < 2
      write(unit=iunit,fmt="(a)") " "

      if (len_trim(code) == 0) then
         write(unit=iunit,fmt="(a)") "data_?"
      else
         write(unit=iunit,fmt="(a)") "data_"//code(1:len_trim(code))
      end if
      write(unit=iunit,fmt="(a)") " "

      if (type_data < 2) then
         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "#============================================================================="
         write(unit=iunit,fmt="(a)") " "

         !> Chemical Data
         write(unit=iunit,fmt="(a)") "# 5. CHEMICAL DATA"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "_chemical_name_systematic"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_chemical_name_common             ?"
         write(unit=iunit,fmt="(a)") "_chemical_formula_moiety          ?"
         write(unit=iunit,fmt="(a)") "_chemical_formula_structural      ?"
         write(unit=iunit,fmt="(a)") "_chemical_formula_analytical      ?"
         write(unit=iunit,fmt="(a)") "_chemical_formula_iupac           ?"
         write(unit=iunit,fmt="(a)") "_chemical_formula_sum             ?"
         write(unit=iunit,fmt="(a)") "_chemical_formula_weight          ?"
         write(unit=iunit,fmt="(a)") "_chemical_melting_point           ?"
         write(unit=iunit,fmt="(a)") "_chemical_compound_source         ?       # for minerals and "
         write(unit=iunit,fmt="(a)") "                                          # natural products"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "loop_"
         write(unit=iunit,fmt="(a)") "    _atom_type_symbol               "
         write(unit=iunit,fmt="(a)") "    _atom_type_description          "
         write(unit=iunit,fmt="(a)") "    _atom_type_scat_dispersion_real "
         write(unit=iunit,fmt="(a)") "    _atom_type_scat_dispersion_imag "
         write(unit=iunit,fmt="(a)") "    _atom_type_scat_source          "
         write(unit=iunit,fmt="(a)") "    _atom_type_scat_length_neutron       # include if applicable"
         write(unit=iunit,fmt="(a)") "    ?    ?    ?    ?    ?      ?    "
      end if !type_data < 2
      write(unit=iunit,fmt="(a)") " "
      write(unit=iunit,fmt="(a)") "#============================================================================="
      write(unit=iunit,fmt="(a)") " "

      !> Crystal Data
      select case (type_data)
         case (0,2) ! Single Crystal or structural data only
            write(unit=iunit,fmt="(a)") "# 6. CRYSTAL DATA"
         case (1) ! Powder Data + Crystal Data
            write(unit=iunit,fmt="(a)") "# 6. POWDER SPECIMEN AND CRYSTAL DATA"
      end select
      write(unit=iunit,fmt="(a)") " "

      write(unit=iunit,fmt="(a)") "_symmetry_cell_setting               ?"
      line=SpG%SPG_Symb
      write(unit=iunit,fmt="(a)") "_symmetry_space_group_name_H-M       '"//trim(line)//"'"
      line=SpG%Hall
      write(unit=iunit,fmt="(a)") "_symmetry_space_group_name_Hall      '"//trim(line)//"'"

      write(unit=iunit,fmt="(a)") " "
      write(unit=iunit,fmt="(a)") "loop_"
      write(unit=iunit,fmt="(a)") "    _symmetry_equiv_pos_as_xyz"
      do i=1,SpG%multip
         line="'"//trim(SpG%Symb_Op(i))//"'"
         write(iunit,'(a)') trim(line)
      end do
      write(unit=iunit,fmt="(a)") " "

      do i=1,3
         text(i)  =String_NumStd(cell%cell(i),cell%scell(i))
         text(i+3)=String_NumStd(cell%ang(i),cell%sang(i))
      end do
      write(unit=iunit,fmt='(a)')       "_cell_length_a                       "//trim(adjustl(text(1)))
      write(unit=iunit,fmt='(a)')       "_cell_length_b                       "//trim(adjustl(text(2)))
      write(unit=iunit,fmt='(a)')       "_cell_length_c                       "//trim(adjustl(text(3)))
      write(unit=iunit,fmt='(a)')       "_cell_angle_alpha                    "//trim(adjustl(text(4)))
      write(unit=iunit,fmt='(a)')       "_cell_angle_beta                     "//trim(adjustl(text(5)))
      write(unit=iunit,fmt='(a)')       "_cell_angle_gamma                    "//trim(adjustl(text(6)))
      write(unit=iunit,fmt="(a,f14.4)") "_cell_volume                   ",Cell%Vol
      if (type_data < 2) then
         write(unit=iunit,fmt="(a)") "_cell_formula_units_Z                ?"
         write(unit=iunit,fmt="(a)") "_cell_measurement_temperature        ?"
         write(unit=iunit,fmt="(a)") "_cell_special_details"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
      end if

      select case (type_data)
         case (0) ! Single Crystal
            write(unit=iunit,fmt="(a)") "_cell_measurement_reflns_used        ?"
            write(unit=iunit,fmt="(a)") "_cell_measurement_theta_min          ?"
            write(unit=iunit,fmt="(a)") "_cell_measurement_theta_max          ?"

            write(unit=iunit,fmt="(a)") " "
            write(unit=iunit,fmt="(a)") "_exptl_crystal_description           ?"
            write(unit=iunit,fmt="(a)") "_exptl_crystal_colour                ?"
            write(unit=iunit,fmt="(a)") "_exptl_crystal_size_max              ?"
            write(unit=iunit,fmt="(a)") "_exptl_crystal_size_mid              ?"
            write(unit=iunit,fmt="(a)") "_exptl_crystal_size_min              ?"
            write(unit=iunit,fmt="(a)") "_exptl_crystal_size_rad              ?"
            write(unit=iunit,fmt="(a)") "_exptl_crystal_density_diffrn        ?"
            write(unit=iunit,fmt="(a)") "_exptl_crystal_density_meas          ?"
            write(unit=iunit,fmt="(a)") "_exptl_crystal_density_method        ?"
            write(unit=iunit,fmt="(a)") "_exptl_crystal_F_000                 ?"

         case (1) ! Powder Data
            write(unit=iunit,fmt="(a)") "# The next three fields give the specimen dimensions in mm.  The equatorial"
            write(unit=iunit,fmt="(a)") "# plane contains the incident and diffracted beam."

            write(unit=iunit,fmt="(a)") " "
            write(unit=iunit,fmt="(a)") "_pd_spec_size_axial               ?       # perpendicular to "
            write(unit=iunit,fmt="(a)") "                                          # equatorial plane"

            write(unit=iunit,fmt="(a)") "_pd_spec_size_equat               ?       # parallel to "
            write(unit=iunit,fmt="(a)") "                                          # scattering vector"
            write(unit=iunit,fmt="(a)") "                                          # in transmission"
            write(unit=iunit,fmt="(a)") "_pd_spec_size_thick               ?       # parallel to "
            write(unit=iunit,fmt="(a)") "                                          # scattering vector"
            write(unit=iunit,fmt="(a)") "                                          # in reflection"

            write(unit=iunit,fmt="(a)") " "
            write(unit=iunit,fmt="(a)") "# The next five fields are character fields that describe the specimen."

            write(unit=iunit,fmt="(a)") " "
            write(unit=iunit,fmt="(a)") "_pd_spec_mounting                         # This field should be"
            write(unit=iunit,fmt="(a)") "                                          # used to give details of the "
            write(unit=iunit,fmt="(a)") "                                          # container."
            write(unit=iunit,fmt="(a)") "; ?"
            write(unit=iunit,fmt="(a)") ";"
            write(unit=iunit,fmt="(a)") "_pd_spec_mount_mode               ?       # options are 'reflection'"
            write(unit=iunit,fmt="(a)") "                                          # or 'transmission'"
            write(unit=iunit,fmt="(a)") "_pd_spec_shape                    ?       # options are 'cylinder' "
            write(unit=iunit,fmt="(a)") "                                          # 'flat_sheet' or 'irregular'"
            write(unit=iunit,fmt="(a)") "_pd_char_particle_morphology      ?"
            write(unit=iunit,fmt="(a)") "_pd_char_colour                   ?       # use ICDD colour descriptions"

            write(unit=iunit,fmt="(a)") " "
            write(unit=iunit,fmt="(a)") "# The following three fields describe the preparation of the specimen."
            write(unit=iunit,fmt="(a)") "# The cooling rate is in K/min.  The pressure at which the sample was "
            write(unit=iunit,fmt="(a)") "# prepared is in kPa.  The temperature of preparation is in K.        "

            write(unit=iunit,fmt="(a)") " "
            write(unit=iunit,fmt="(a)") "_pd_prep_cool_rate                ?"
            write(unit=iunit,fmt="(a)") "_pd_prep_pressure                 ?"
            write(unit=iunit,fmt="(a)") "_pd_prep_temperature              ?"
      end select
      if (type_data < 2) then
         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "# The next four fields are normally only needed for transmission experiments."
         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "_exptl_absorpt_coefficient_mu        ?"
         write(unit=iunit,fmt="(a)") "_exptl_absorpt_correction_type       ?"
         write(unit=iunit,fmt="(a)") "_exptl_absorpt_process_details       ?"
         write(unit=iunit,fmt="(a)") "_exptl_absorpt_correction_T_min      ?"
         write(unit=iunit,fmt="(a)") "_exptl_absorpt_correction_T_max      ?"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "#============================================================================="
         write(unit=iunit,fmt="(a)") " "

         !> Experimental Data
         write(unit=iunit,fmt="(a)") "# 7. EXPERIMENTAL DATA"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "_exptl_special_details"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"

        if (type_data == 1) then
           write(unit=iunit,fmt="(a)") " "
           write(unit=iunit,fmt="(a)") "# The following item is used to identify the equipment used to record "
           write(unit=iunit,fmt="(a)") "# the powder pattern when the diffractogram was measured at a laboratory "
           write(unit=iunit,fmt="(a)") "# other than the authors' home institution, e.g. when neutron or synchrotron"
           write(unit=iunit,fmt="(a)") "# radiation is used."

           write(unit=iunit,fmt="(a)") " "
           write(unit=iunit,fmt="(a)") "_pd_instr_location"
           write(unit=iunit,fmt="(a)") "; ?"
           write(unit=iunit,fmt="(a)") ";"
           write(unit=iunit,fmt="(a)") "_pd_calibration_special_details           # description of the method used"
           write(unit=iunit,fmt="(a)") "                                          # to calibrate the instrument"
           write(unit=iunit,fmt="(a)") "; ?"
           write(unit=iunit,fmt="(a)") ";"
        end if

        write(unit=iunit,fmt="(a)") " "
        write(unit=iunit,fmt="(a)") "_diffrn_ambient_temperature          ?"
        write(unit=iunit,fmt="(a)") "_diffrn_radiation_type               ?"
        write(unit=iunit,fmt="(a)") "_diffrn_radiation_wavelength         ?"
        write(unit=iunit,fmt="(a)") "_diffrn_radiation_source             ?"
        write(unit=iunit,fmt="(a)") "_diffrn_source                       ?"
        write(unit=iunit,fmt="(a)") "_diffrn_source_target                ?"
        write(unit=iunit,fmt="(a)") "_diffrn_source_type                  ?"

        write(unit=iunit,fmt="(a)") " "
        write(unit=iunit,fmt="(a)") "_diffrn_radiation_monochromator      ?"
        write(unit=iunit,fmt="(a)") "_diffrn_measurement_device_type      ?"
        write(unit=iunit,fmt="(a)") "_diffrn_measurement_method           ?"
        write(unit=iunit,fmt="(a)") "_diffrn_detector_area_resol_mean     ?   # Not in version 2.0.1"
        write(unit=iunit,fmt="(a)") "_diffrn_detector                     ?"
        write(unit=iunit,fmt="(a)") "_diffrn_detector_type                ?   # make or model of detector"
        if (type_data == 1) then
           write(unit=iunit,fmt="(a)") "_pd_meas_scan_method                 ?   # options are 'step', 'cont',"
           write(unit=iunit,fmt="(a)") "                                         # 'tof', 'fixed' or"
           write(unit=iunit,fmt="(a)") "                                         # 'disp' (= dispersive)"
           write(unit=iunit,fmt="(a)") "_pd_meas_special_details"
           write(unit=iunit,fmt="(a)") ";  ?"
           write(unit=iunit,fmt="(a)") ";"
        end if

        select case (type_data)
           case (0)
              write(unit=iunit,fmt="(a)") " "
              write(unit=iunit,fmt="(a)") "_diffrn_reflns_number                ?"
              write(unit=iunit,fmt="(a)") "_diffrn_reflns_av_R_equivalents      ?"
              write(unit=iunit,fmt="(a)") "_diffrn_reflns_av_sigmaI/netI        ?"
              write(unit=iunit,fmt="(a)") "_diffrn_reflns_theta_min             ?"
              write(unit=iunit,fmt="(a)") "_diffrn_reflns_theta_max             ?"
              write(unit=iunit,fmt="(a)") "_diffrn_reflns_theta_full            ?"
              write(unit=iunit,fmt="(a)") "_diffrn_measured_fraction_theta_max  ?"
              write(unit=iunit,fmt="(a)") "_diffrn_measured_fraction_theta_full ?"
              write(unit=iunit,fmt="(a)") "_diffrn_reflns_limit_h_min           ?"
              write(unit=iunit,fmt="(a)") "_diffrn_reflns_limit_h_max           ?"
              write(unit=iunit,fmt="(a)") "_diffrn_reflns_limit_k_min           ?"
              write(unit=iunit,fmt="(a)") "_diffrn_reflns_limit_k_max           ?"
              write(unit=iunit,fmt="(a)") "_diffrn_reflns_limit_l_min           ?"
              write(unit=iunit,fmt="(a)") "_diffrn_reflns_limit_l_max           ?"
              write(unit=iunit,fmt="(a)") "_diffrn_reflns_reduction_process     ?"

              write(unit=iunit,fmt="(a)") " "
              write(unit=iunit,fmt="(a)") "_diffrn_standards_number             ?"
              write(unit=iunit,fmt="(a)") "_diffrn_standards_interval_count     ?"
              write(unit=iunit,fmt="(a)") "_diffrn_standards_interval_time      ?"
              write(unit=iunit,fmt="(a)") "_diffrn_standards_decay_%            ?"
              write(unit=iunit,fmt="(a)") "loop_"
              write(unit=iunit,fmt="(a)") "    _diffrn_standard_refln_index_h"
              write(unit=iunit,fmt="(a)") "    _diffrn_standard_refln_index_k"
              write(unit=iunit,fmt="(a)") "    _diffrn_standard_refln_index_l"
              write(unit=iunit,fmt="(a)") "?   ?   ?"

           case (1)
              write(unit=iunit,fmt="(a)") " "
              write(unit=iunit,fmt="(a)") "#  The following four items give details of the measured (not processed)"
              write(unit=iunit,fmt="(a)") "#  powder pattern.  Angles are in degrees."

              write(unit=iunit,fmt="(a)") " "
              write(unit=iunit,fmt="(a)") "_pd_meas_number_of_points         ?"
              write(unit=iunit,fmt="(a)") "_pd_meas_2theta_range_min         ?"
              write(unit=iunit,fmt="(a)") "_pd_meas_2theta_range_max         ?"
              write(unit=iunit,fmt="(a)") "_pd_meas_2theta_range_inc         ?"

              write(unit=iunit,fmt="(a)") " "
              write(unit=iunit,fmt="(a)") "# The following three items are used for time-of-flight measurements only."

              write(unit=iunit,fmt="(a)") " "
              write(unit=iunit,fmt="(a)") "_pd_instr_dist_src/spec           ?"
              write(unit=iunit,fmt="(a)") "_pd_instr_dist_spec/detc          ?"
              write(unit=iunit,fmt="(a)") "_pd_meas_2theta_fixed             ?"

        end select

        write(unit=iunit,fmt="(a)") " "
        write(unit=iunit,fmt="(a)") "#============================================================================="
        write(unit=iunit,fmt="(a)") " "

        !> Refinement Data
        write(unit=iunit,fmt="(a)") "# 8. REFINEMENT DATA"

        write(unit=iunit,fmt="(a)") " "

        write(unit=iunit,fmt="(a)") "_refine_special_details"
        write(unit=iunit,fmt="(a)") "; ?"
        write(unit=iunit,fmt="(a)") ";"

        if (type_data == 1) then
           write(unit=iunit,fmt="(a)") " "
           write(unit=iunit,fmt="(a)") "# Use the next field to give any special details about the fitting of the"
           write(unit=iunit,fmt="(a)") "# powder pattern."

           write(unit=iunit,fmt="(a)") " "
           write(unit=iunit,fmt="(a)") "_pd_proc_ls_special_details"
           write(unit=iunit,fmt="(a)") "; ?"
           write(unit=iunit,fmt="(a)") ";"

           write(unit=iunit,fmt="(a)") " "
           write(unit=iunit,fmt="(a)") "# The next three items are given as text."
           write(unit=iunit,fmt="(a)") " "

           write(unit=iunit,fmt="(a)") "_pd_proc_ls_profile_function      ?"
           write(unit=iunit,fmt="(a)") "_pd_proc_ls_background_function   ?"
           write(unit=iunit,fmt="(a)") "_pd_proc_ls_pref_orient_corr"
           write(unit=iunit,fmt="(a)") "; ?"
           write(unit=iunit,fmt="(a)") ";"
        end if

        select case (type_data)
           case (0)
              write(unit=iunit,fmt="(a)") " "
              write(unit=iunit,fmt="(a)") "_reflns_number_total                 ?"
              write(unit=iunit,fmt="(a)") "_reflns_number_gt                    ?"
              write(unit=iunit,fmt="(a)") "_reflns_threshold_expression         ?"

           case (1)
              write(unit=iunit,fmt="(a)") " "
              write(unit=iunit,fmt="(a)") "_pd_proc_ls_prof_R_factor         ?"
              write(unit=iunit,fmt="(a)") "_pd_proc_ls_prof_wR_factor        ?"
              write(unit=iunit,fmt="(a)") "_pd_proc_ls_prof_wR_expected      ?"

              write(unit=iunit,fmt="(a)") " "
              write(unit=iunit,fmt="(a)") "# The following four items apply to angular dispersive measurements."
              write(unit=iunit,fmt="(a)") "# 2theta minimum, maximum and increment (in degrees) are for the "
              write(unit=iunit,fmt="(a)") "# intensities used in the refinement."

              write(unit=iunit,fmt="(a)") " "
              write(unit=iunit,fmt="(a)") "_pd_proc_2theta_range_min         ?"
              write(unit=iunit,fmt="(a)") "_pd_proc_2theta_range_max         ?"
              write(unit=iunit,fmt="(a)") "_pd_proc_2theta_range_inc         ?"
              write(unit=iunit,fmt="(a)") "_pd_proc_wavelength               ?"

              write(unit=iunit,fmt="(a)") " "
              write(unit=iunit,fmt="(a)") "_pd_block_diffractogram_id        ?  # The id used for the block containing"
              write(unit=iunit,fmt="(a)") "                                     # the powder pattern profile (section 11)."

              write(unit=iunit,fmt="(a)") " "
              write(unit=iunit,fmt="(a)") "# Give appropriate details in the next two text fields."
              write(unit=iunit,fmt="(a)") " "
              write(unit=iunit,fmt="(a)") "_pd_proc_info_excluded_regions    ?"
              write(unit=iunit,fmt="(a)") "_pd_proc_info_data_reduction      ?"
        end select

        write(unit=iunit,fmt="(a)") " "
        write(unit=iunit,fmt="(a)") "_refine_ls_structure_factor_coef     ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_matrix_type               ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_R_I_factor                ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_R_Fsqd_factor             ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_R_factor_all              ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_R_factor_gt               ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_wR_factor_all             ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_wR_factor_ref             ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_goodness_of_fit_all       ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_goodness_of_fit_ref       ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_restrained_S_all          ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_restrained_S_obs          ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_number_reflns             ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_number_parameters         ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_number_restraints         ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_number_constraints        ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_hydrogen_treatment        ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_weighting_scheme          ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_weighting_details         ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_shift/su_max              ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_shift/su_mean             ?"
        write(unit=iunit,fmt="(a)") "_refine_diff_density_max             ?"
        write(unit=iunit,fmt="(a)") "_refine_diff_density_min             ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_extinction_method         ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_extinction_coef           ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_abs_structure_details     ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_abs_structure_Flack       ?"
        write(unit=iunit,fmt="(a)") "_refine_ls_abs_structure_Rogers      ?"

        write(unit=iunit,fmt="(a)") " "
        write(unit=iunit,fmt="(a)") "# The following items are used to identify the programs used."
        write(unit=iunit,fmt="(a)") " "

        write(unit=iunit,fmt="(a)") "_computing_data_collection           ?"
        write(unit=iunit,fmt="(a)") "_computing_cell_refinement           ?"
        write(unit=iunit,fmt="(a)") "_computing_data_reduction            ?"
        write(unit=iunit,fmt="(a)") "_computing_structure_solution        ?"
        write(unit=iunit,fmt="(a)") "_computing_structure_refinement      ?"
        write(unit=iunit,fmt="(a)") "_computing_molecular_graphics        ?"
        write(unit=iunit,fmt="(a)") "_computing_publication_material      ?"
      end if  !(type_data < 2) then
      write(unit=iunit,fmt="(a)") " "
      write(unit=iunit,fmt="(a)") "#============================================================================="
      write(unit=iunit,fmt="(a)") " "

      !> Atomic Coordinates and Displacement Parameters
      write(unit=iunit,fmt="(a)") "# 9. ATOMIC COORDINATES AND DISPLACEMENT PARAMETERS"

      write(unit=iunit,fmt="(a)") " "

      write(unit=iunit,fmt="(a)") "loop_"
      write(unit=iunit,fmt='(a)') "    _atom_site_label"
      write(unit=iunit,fmt='(a)') "    _atom_site_type_symbol"
      write(unit=iunit,fmt='(a)') "    _atom_site_fract_x"
      write(unit=iunit,fmt='(a)') "    _atom_site_fract_y"
      write(unit=iunit,fmt='(a)') "    _atom_site_fract_z"
      write(unit=iunit,fmt='(a)') "    _atom_site_U_iso_or_equiv"
      write(unit=iunit,fmt='(a)') "    _atom_site_occupancy"
      write(unit=iunit,fmt='(a)') "    _atom_site_adp_type"
      write(unit=iunit,fmt='(a)') "    _atom_site_type_symbol"

      !> Calculation of the factor corresponding to the occupation factor provided in A
      Select Type(atms => At_list%Atom)
        Type is (Atm_Std_Type)
          atm=atms
        Class default
          do i=1,At_List%natoms
            call Init_Atom_Type(Atm(i),0)
            Atm(i)%lab     =Atms(i)%lab
            Atm(i)%ChemSymb=Atms(i)%ChemSymb
            Atm(i)%SfacSymb=Atms(i)%SfacSymb
            Atm(i)%Z       =Atms(i)%Z
            Atm(i)%mult    =Atms(i)%mult
            Atm(i)%mult    =Atms(i)%charge
            Atm(i)%x       =Atms(i)%x
            Atm(i)%x       =Atms(i)%U_iso
            Atm(i)%occ     =Atms(i)%occ
            Atm(i)%UType   =Atms(i)%UType
            Atm(i)%ThType  =Atms(i)%ThType
            Atm(i)%U       =Atms(i)%U
            Atm(i)%magnetic=Atms(i)%magnetic
            Atm(i)%mom     =Atms(i)%mom
            Atm(i)%Ind_ff  =Atms(i)%Ind_ff
            Atm(i)%AtmInfo =Atms(i)%AtmInfo
            Atm(i)%wyck    =Atms(i)%wyck
          end do

      End Select

      do i=1,At_List%natoms
         occup(i)=Atm(i)%occ/(real(Atm(i)%mult)/real(SpG%multip))
        soccup(i)=Atm(i)%occ_std/(real(Atm(i)%mult)/real(SpG%multip))
      end do
      ocf=sum(abs(Atm(1)%x-Atm(2)%x))
      if ( ocf < 0.001) then
         ocf=occup(1)+occup(2)
      else
         ocf=occup(1)
      end if
      occup=occup/ocf; soccup=soccup/ocf
      aniso=.false.
      do i=1,At_List%natoms
         line=" "
         line(2:)= Atm(i)%Lab//"  "//Atm(i)%SfacSymb

         do j=1,3
            comm=String_NumStd(Atm(i)%x(j),Atm(i)%x_std(j))
            line=trim(line)//" "//trim(comm)
         end do

         comm=" "
         select case (Atm(i)%Thtype)
            case ('iso')
               adptyp='Uiso'
               select case (trim(Atm(i)%UType))
                  case ("U")
                     u=Atm(i)%U_iso
                     su=Atm(i)%U_iso_std

                  case ("B")
                     u=Atm(i)%U_iso/(8.0*pi*pi)
                     su=Atm(i)%U_iso_std/(8.0*pi*pi)

                  case ("beta")
                     u=Atm(i)%U_iso
                     su=Atm(i)%U_iso_std
               end select
               comm=String_NumStd(u,su)

            case ('ani')
               aniso=.true.
               adptyp='Uani'
               select case (trim(Atm(i)%UType))
                  case ("U")
                     ua=Atm(i)%u
                     sua=Atm(i)%u_std

                  case ("B")
                     ua=Atm(i)%u/(8.0*pi*pi)
                     sua=Atm(i)%u_std/(8.0*pi*pi)

                  case ("beta")
                     aux=Atm(i)%u
                     ua=get_U_from_Betas(aux,cell)
                     aux=Atm(i)%u_std
                     sua=get_U_from_Betas(aux,cell)
               end select
               u=(ua(1)+ua(2)+ua(3))/3.0
               su=(ua(1)+ua(2)+ua(3))/3.0
               comm=String_NumStd(u,su)

            case default
               adptyp='.'
         end select
         line=trim(line)//" "//trim(comm)

         comm=String_NumStd(occup(i),soccup(i))
         line=trim(line)//" "//trim(comm)
         write(unit=iunit,fmt="(a)") trim(line)//" "//trim(adptyp)//" "//Atm(i)%SfacSymb
      end do

      if (aniso) then
         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "loop_"
         write(unit=iunit,fmt="(a)") "    _atom_site_aniso_label "
         write(unit=iunit,fmt="(a)") "    _atom_site_aniso_U_11  "
         write(unit=iunit,fmt="(a)") "    _atom_site_aniso_U_22  "
         write(unit=iunit,fmt="(a)") "    _atom_site_aniso_U_33  "
         write(unit=iunit,fmt="(a)") "    _atom_site_aniso_U_12  "
         write(unit=iunit,fmt="(a)") "    _atom_site_aniso_U_13  "
         write(unit=iunit,fmt="(a)") "    _atom_site_aniso_U_23  "
         write(unit=iunit,fmt="(a)") "    _atom_site_aniso_type_symbol"

         do i=1,At_List%natoms
            if (Atm(i)%thtype /= "ani") cycle

            line=" "
            line(2:)= Atm(i)%Lab

            select case (trim(Atm(i)%UType))
               case ("U")
                  ua=Atm(i)%u
                  sua=Atm(i)%u_std

               case ("B")
                  ua=Atm(i)%u/(8.0*pi*pi)
                  sua=Atm(i)%u_std/(8.0*pi*pi)

               case ("beta")
                  aux=Atm(i)%u
                  ua=get_U_from_Betas(aux,cell)
                  aux=Atm(i)%u_std
                  sua=get_U_from_Betas(aux,cell)
            end select

            do j=1,6
              comm=" "
              comm=String_NumStd(ua(j),sua(j))
              line=trim(line)//" "//trim(comm)
            end do
            write(iunit,"(a)") trim(line)//"  "//Atm(i)%SfacSymb
         end do
      end if

      if(type_data < 2) then
         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "# Note: if the displacement parameters were refined anisotropically"
         write(unit=iunit,fmt="(a)") "# the U matrices should be given as for single-crystal studies."

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "#============================================================================="
         write(unit=iunit,fmt="(a)") " "

         !---- Molecular Geometry ----!
         write(unit=iunit,fmt="(a)") "# 10. MOLECULAR GEOMETRY"

         write(unit=iunit,fmt="(a)") " "


         write(unit=iunit,fmt="(a)") "_geom_special_details                ?"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "loop_"
         write(unit=iunit,fmt="(a)") "    _geom_bond_atom_site_label_1  "
         write(unit=iunit,fmt="(a)") "    _geom_bond_atom_site_label_2  "
         write(unit=iunit,fmt="(a)") "    _geom_bond_site_symmetry_1    "
         write(unit=iunit,fmt="(a)") "    _geom_bond_site_symmetry_2    "
         write(unit=iunit,fmt="(a)") "    _geom_bond_distance           "
         write(unit=iunit,fmt="(a)") "    _geom_bond_publ_flag          "
         write(unit=iunit,fmt="(a)") "    ?   ?   ?   ?   ?   ?"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "loop_"
         write(unit=iunit,fmt="(a)") "    _geom_contact_atom_site_label_1 "
         write(unit=iunit,fmt="(a)") "    _geom_contact_atom_site_label_2 "
         write(unit=iunit,fmt="(a)") "    _geom_contact_distance          "
         write(unit=iunit,fmt="(a)") "    _geom_contact_site_symmetry_1   "
         write(unit=iunit,fmt="(a)") "    _geom_contact_site_symmetry_2   "
         write(unit=iunit,fmt="(a)") "    _geom_contact_publ_flag         "
         write(unit=iunit,fmt="(a)") "    ?   ?   ?   ?   ?   ?"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "loop_"
         write(unit=iunit,fmt="(a)") "_geom_angle_atom_site_label_1 "
         write(unit=iunit,fmt="(a)") "_geom_angle_atom_site_label_2 "
         write(unit=iunit,fmt="(a)") "_geom_angle_atom_site_label_3 "
         write(unit=iunit,fmt="(a)") "_geom_angle_site_symmetry_1   "
         write(unit=iunit,fmt="(a)") "_geom_angle_site_symmetry_2   "
         write(unit=iunit,fmt="(a)") "_geom_angle_site_symmetry_3   "
         write(unit=iunit,fmt="(a)") "_geom_angle                   "
         write(unit=iunit,fmt="(a)") "_geom_angle_publ_flag         "
         write(unit=iunit,fmt="(a)") "?   ?   ?   ?   ?   ?   ?   ?"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "loop_"
         write(unit=iunit,fmt="(a)") "_geom_torsion_atom_site_label_1 "
         write(unit=iunit,fmt="(a)") "_geom_torsion_atom_site_label_2 "
         write(unit=iunit,fmt="(a)") "_geom_torsion_atom_site_label_3 "
         write(unit=iunit,fmt="(a)") "_geom_torsion_atom_site_label_4 "
         write(unit=iunit,fmt="(a)") "_geom_torsion_site_symmetry_1   "
         write(unit=iunit,fmt="(a)") "_geom_torsion_site_symmetry_2   "
         write(unit=iunit,fmt="(a)") "_geom_torsion_site_symmetry_3   "
         write(unit=iunit,fmt="(a)") "_geom_torsion_site_symmetry_4   "
         write(unit=iunit,fmt="(a)") "_geom_torsion                   "
         write(unit=iunit,fmt="(a)") "_geom_torsion_publ_flag         "
         write(unit=iunit,fmt="(a)") "?   ?   ?   ?   ?   ?   ?   ?   ?   ?"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "loop_"
         write(unit=iunit,fmt="(a)") "_geom_hbond_atom_site_label_D "
         write(unit=iunit,fmt="(a)") "_geom_hbond_atom_site_label_H "
         write(unit=iunit,fmt="(a)") "_geom_hbond_atom_site_label_A "
         write(unit=iunit,fmt="(a)") "_geom_hbond_site_symmetry_D   "
         write(unit=iunit,fmt="(a)") "_geom_hbond_site_symmetry_H   "
         write(unit=iunit,fmt="(a)") "_geom_hbond_site_symmetry_A   "
         write(unit=iunit,fmt="(a)") "_geom_hbond_distance_DH       "
         write(unit=iunit,fmt="(a)") "_geom_hbond_distance_HA       "
         write(unit=iunit,fmt="(a)") "_geom_hbond_distance_DA       "
         write(unit=iunit,fmt="(a)") "_geom_hbond_angle_DHA         "
         write(unit=iunit,fmt="(a)") "_geom_hbond_publ_flag         "
         write(unit=iunit,fmt="(a)") "?   ?   ?   ?   ?   ?   ?   ?   ?   ?   ?"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "#============================================================================="
         write(unit=iunit,fmt="(a)") " "


         !---- Final Informations ----!
         write(unit=iunit,fmt="(a)") "#============================================================================="
         write(unit=iunit,fmt="(a)") "# Additional structures (last six sections and associated data_? identifiers) "
         write(unit=iunit,fmt="(a)") "# may be added at this point.                                                 "
         write(unit=iunit,fmt="(a)") "#============================================================================="

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "# The following lines are used to test the character set of files sent by     "
         write(unit=iunit,fmt="(a)") "# network email or other means. They are not part of the CIF data set.        "
         write(unit=iunit,fmt="(a)") "# abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789              "
         write(unit=iunit,fmt="(a)") "# !@#$%^&*()_+{}:"//""""//"~<>?|\-=[];'`,./ "
      end if

      close(unit=iunit)

      return
   End Subroutine Write_Cif_Template



 End Module test_CIF
