program XRMS_ref
!
!  This program refines a magnetic structure model using the integrated intensities 
!  obtained in the experiment described in the draft (Mardegan, 2017).
!  Sample: EuPtIn4
!  Polarization: sigma -> pi'
!

   !  Modules
   use CFML_Math_General
   use CFML_crystal_metrics,              only: write_crystal_cell
   use CFML_crystallographic_symmetry,    only: write_spacegroup
   use CFML_atom_typedef,                 only: write_atom_list
   use CFML_IO_Formats,                   only: readn_set_xtal_structure
   use CFML_Magnetic_Symmetry,            only: readn_set_magnetic_structure, write_magnetic_structure
   use CFML_Geometry_Calc
   use CFML_Optimization_LSQ

   use global_data
   use XRMS_procedures

   implicit none
   
   !  Variables
   integer                                      :: i, n_ini, n_end
   real                                         :: fcalc2, chi2
   complex                                      :: m_str_f
   type(magh_type_xrms)                         :: m_refl
   character(len=256)                           :: inform
   
   real,dimension(3,3)           :: matr_man, matr_cfml
   
   real,dimension(3)    :: orig = (/0.0, 0.0, 0.0/)
   real,dimension(3)    :: m_mom

   !  Procedures

   !  Getting the name of the file with the structure information
   call get_filenames(filename, dataname, outname)

   !  Open the output file.
   open(unit=luo,file=outname,status="replace",action="write")
   write(unit=*,fmt="(a)") "The output file is " // outname

   ! Reading the information of the structure
   call readn_set_xtal_structure(filename,unit_cell,space_gr,atom_list,file_list=str_file)
   n_ini = 1
   n_end = str_file%nlines
   call readn_set_magnetic_structure(str_file,n_ini,n_end,magn_model,magn_atom_list)

   call read_hkl_file(dataname,sample_refl_list)
   call read_hkl_output('hkl_matrices_output.txt',sample_refl_list)
   
   matrix_mom_star = 0.0
   do i=1,3
      matrix_mom_star(i,i) = unit_cell%cell(i)/tpi
   end do

   !  Writting the information of the structure in the output file.
   call write_crystal_cell(unit_cell,luo)
   call write_spacegroup(space_gr,luo,full=.true.)
   call write_atom_list(atom_list,lun=luo)
   call write_magnetic_structure(luo,magn_model,magn_atom_list)
   
   call setup_LSQ_ref()
   
   write(unit=*,fmt='(a)') 'Starting the refinement...'
   write(unit=*,fmt='(a)') ''
   write(unit=luo,fmt='(a)') 'Starting the refinement...'
   write(unit=luo,fmt='(a)') ''
   write(unit=*,fmt='(a)') 'Initial parameters -----'
   write(unit=*,fmt='(a,i5)') 'Number of parameters: ', vs%np
   do i=1,vs%np
      write(unit=*,fmt='(a,i1,a,f14.6,a,f14.6,a)') 'Parameter (', i, ') = ', vs%pv(i),' (+/-',vs%spv(i),')'
   end do
   
   call Levenberg_Marquardt_Fit(intensities_noder,dat%nobs,cond,vs,chi2,inform) ! Numerical derivatives
   !call Levenberg_Marquardt_Fit(intensities_der,dat%nobs,cond,vs,chi2,.true.,inform) ! Analytical derivatives
   write(unit=*,fmt='(a)') inform
   
   write(unit=*,fmt='(a)') 'Final parameters -----'
   write(unit=*,fmt='(a,i5)') 'Number of parameters: ', vs%np
   do i=1,vs%np
      write(unit=*,fmt='(a,i1,a,f14.6,a,f14.6,a)') 'Parameter (', i, ') = ', vs%pv(i),' (+/-',vs%spv(i),')'
   end do
   write(unit=*,fmt='(a)') 'Final resctrictions -----'
   !write(unit=*,fmt=*) dat%yc(dat%nobs-1)
   write(unit=*,fmt=*) dat%yc(dat%nobs)
   
   write(unit=*,fmt='(/3a8,2a12)') 'h', 'k', 'l', 'Iobs', 'Icalc'
   do i=1,sample_refl_list%nref
      m_refl = sample_refl_list%mh(i)
      m_str_f = magn_str_factor(m_refl)
      fcalc2 = scale_fact*m_str_f*conjg(m_str_f)
      write(unit=*,fmt='(3f8.3,2f12.5)') m_refl%h, m_refl%fobs2, fcalc2
   end do

   call write_intensities()
   
   close(unit=luo)
   
end program XRMS_ref
