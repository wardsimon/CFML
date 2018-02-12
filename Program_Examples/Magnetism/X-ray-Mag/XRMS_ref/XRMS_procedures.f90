module XRMS_procedures

   !  Modules
   use CFML_globaldeps,                   only: tpi
   use CFML_math_general
   use CFML_IO_formats,                   only: file_list_type, file_to_filelist
   use CFML_reflections_utilities,        only: hkl_s, hkl_equal
   use CFML_atom_typedef,                 only: matom_type
   use CFML_crystallographic_symmetry,    only: applyso
   use CFML_geometry_calc
   use CFML_String_Utilities,             only: getword

   use global_data

   implicit none
   
   integer           :: ierr
  
   contains
   
   subroutine get_filenames(f_str,f_data,f_out)
      character(len=*),    intent(out)    :: f_str
      character(len=*),    intent(out)    :: f_data
      character(len=*),    intent(out)    :: f_out

      if (command_argument_count() > 0) then
         call get_command_argument(1,f_str)
      else
         write(unit=*,fmt="(/,a)",advance="no") "File with structure information (extension .cfl): "
         read(unit=*,fmt="(a)",iostat=ierr) f_str
      end if
      
      if (command_argument_count() > 1) then
         call get_command_argument(2,f_data)
      else
         write(unit=*,fmt="(/,a)",advance="no") "File with integrated intensities (extension .hkl): "
         read(unit=*,fmt="(a)",iostat=ierr) f_data
      end if
      
      allocate(character(len=len_trim(f_str)-4) :: codename)
      codename = f_str(1:len_trim(f_str)-4)
      f_out = codename // ".out"
      return
   end subroutine get_filenames
   
   !  Subroutine read_hkl_file(f_hkl,refl_list)
   !     !  Arguments
   !     character(len=*),             intent(in)  :: f_hkl
   !     type(magh_list_type_xrms),    intent(out) :: refl_list
   !
   !  This subroutine reads the experimental data from the hkl file.
   !  The data are integrated intensities and the file format is:
   !  h k l Fo2 [sigma(Fo2)]  (sigma is optional)
   !
   subroutine read_hkl_file(f_hkl,refl_list)
      !  Arguments
      character(len=*),             intent(in)  :: f_hkl
      type(magh_list_type_xrms),    intent(out) :: refl_list
      
      !  Local variables
      type(file_list_type)    :: hkl_list
      integer                 :: i, n_head, n_refl, j, n_col
      character(len=132)      :: new_line
      character(len=20), dimension(50) :: dire
      integer, dimension(1000):: point_to_ref
      
      call file_to_filelist(f_hkl,hkl_list)
      
      n_head = 0
      n_refl = 0
      do i=1,hkl_list%nlines         
         new_line = adjustl(hkl_list%line(i))
         if (new_line(1:1) == '!' .or. len_trim(new_line) < 1 ) cycle
         n_refl = n_refl + 1
         point_to_ref(n_refl)=i
      end do
      
      if (n_refl < 1) then
         write(unit=*,fmt='(a)') 'There is some problem with the file '//trim(f_hkl)//'.'
         write(unit=*,fmt='(a)') 'No reflections have been detected.'
         stop
      end if
      call Getword(hkl_list%line(n_refl),dire,n_col)
      
      if (n_col < 4) then
         write(unit=*,fmt='(a)') 'Please, check the data file '//trim(f_hkl)//'.'
         write(unit=*,fmt='(a)') 'There is some problem with the file format.'
         stop
      end if
      
      refl_list%nref = n_refl
      if (allocated(refl_list%mh)) deallocate(refl_list%mh)
      allocate(refl_list%mh(n_refl))
      
      do i=1,n_refl
         new_line = hkl_list%line(point_to_ref(i))
         write(*,"(a)") trim(new_line)
         if (n_col == 4) then
            read(new_line,fmt=*) refl_list%mh(i)%h(:), refl_list%mh(i)%fobs2
            refl_list%mh(i)%sfobs2 = sqrt(refl_list%mh(i)%fobs2)
         else if (n_col == 5) then
            read(new_line,fmt=*) refl_list%mh(i)%h(:), refl_list%mh(i)%fobs2, &
                                 refl_list%mh(i)%sfobs2
         end if
      end do
      return
   end subroutine read_hkl_file
   
   subroutine read_hkl_output(filename,refl_list)
      character(len=*),             intent(in)        :: filename
      type(magh_list_type_xrms),    intent(inout)     :: refl_list
      
      !  Local variables
      integer           :: i, j
      
      open(unit=8,file=filename,status='old',action='read')
      
      do i=1,refl_list%nref
         do j=1,3
            read(unit=8,fmt=*)
         end do
         read(unit=8,fmt=*) refl_list%mh(i)%kf(:)
         read(unit=8,fmt=*)
         do j=1,3
            read(unit=8,fmt=*) refl_list%mh(i)%tmatrix(j,:)
         end do
         read(unit=8,fmt=*)
      end do
      return
   end subroutine read_hkl_output
   
   subroutine setup_LSQ_ref()
      
      !  Local variables
      integer        :: nobs, i, m, l
      
      !  Data values
      !nobs = sample_refl_list%nref ! No restrictions 
      nobs = sample_refl_list%nref + 1 ! With restriction: |m|^2 = 1
      !nobs = sample_refl_list%nref + 2
      dat%nobs = nobs
      allocate(dat%x(nobs),dat%y(nobs),dat%sw(nobs),dat%yc(nobs))
      dat%iw = 0
      do i=1,sample_refl_list%nref
         dat%x(i) = i
         dat%y(i) = sample_refl_list%mh(i)%fobs2
         dat%sw(i) = sample_refl_list%mh(i)%sfobs2
      end do
      
      !  Restriction: m(perp) = 0
      !dat%x(nobs-1) = nobs - 1
      !dat%y(nobs-1) = 0.0
      !dat%sw(nobs-1) = 1e-4
      
      !  Restriction: |m|^2 = 1
      dat%x(nobs) = nobs
      dat%y(nobs) = 1.0
      dat%sw(nobs) = 1e-4
      
      !  Initial values for state vector
      vs%np = 4
      vs%pv(1) = scale_fact
      vs%pv(2) = magn_atom_list%atom(1)%cbas(1,1)
      vs%pv(3) = magn_atom_list%atom(1)%cbas(2,1)
      vs%pv(4) = magn_atom_list%atom(1)%cbas(3,1)
      vs%code(1:vs%np) = 1
   
      !  Conditions for running LSQ refinement
      cond%npvar = vs%np
      cond%icyc = 5000
      cond%tol = 1e-6
      !cond%constr = .true.
      !cond%percent = 25
      return
   end subroutine setup_LSQ_ref

   !  Function theta(refl) result(th_value)
   !     !  Arguments
   !     type(magh_type_desy),    intent(in)  :: refl
   !     real,                    intent(out) :: th_value
   !
   !  This function calculates theta angle (Bragg angle/2) given the reflection
   !  and the unit cell.
   !
   function theta(m_refl) result(th_value)
      !  Arguments
      type(magh_type_xrms),      intent(in)  :: m_refl
      real                                   :: th_value
      
      th_value = asind(lambda*hkl_s(m_refl%h,unit_cell))
      return
   end function theta
   
   !  Function alpha(m_refl) result(alpha_value)
   !     !  Arguments
   !     type(magh_type_desy),      intent(in)  :: refl
   !     real                       intent(out) :: alpha_value
   !
   !  This function calculates alpha angle between the specular axis
   !  and the wave vector transfer.
   !
   function alpha(m_refl) result(alpha_val)
      !  Arguments
      type(magh_type_xrms),      intent(in)  :: m_refl
      real                                   :: alpha_val
      
      !  Local variables
      integer                     :: i
      
      alpha_val = acosd(dot_product(specular_axis,m_refl%h)/&
                     (2*norm(specular_axis,unit_cell%gd)*hkl_s(m_refl%h,unit_cell)))
      return
   end function alpha
   
   function factor_abs_lor(m_refl) result(factor_val)
      !  Arguments
      type(magh_type_xrms),      intent(in)  :: m_refl
      real                                   :: factor_val
      
      !  Local variables
      real                 :: theta, norm_h, norm_n, alpha
      real,dimension(3)    :: n_surf = (/ 0, 1, 0 /)
      
      theta = asind(lambda*hkl_s(m_refl%h,unit_cell))
      norm_h = norm(m_refl%h,unit_cell%gr)
      norm_n = norm(n_surf,unit_cell%gd)
      alpha = acosd(dot_product(m_refl%h,n_surf)/(norm_h*norm_n))
      
      factor_val = (sind(theta + alpha)*sind(theta - alpha))/&
                     (sind(theta)*cosd(alpha)*sind(2*theta))
      return
   end function factor_abs_lor
   
   !  Function magn_moment(m_atom,n_atom) result(mag_mom)
   !     !  Arguments
   !     type(matom_type),       intent(in)  :: m_atom
   !     integer,                intent(in)  :: n_so
   !     real,dimension(3)       intent(out) :: mag_mom
   !
   !  This function calculates the magnetic moment of the atom
   !  corresponding to magn_model%symop(n_so).
   !
   function magn_moment(m_atom,n_so) result(m_mom)
      !  Arguments
      type(matom_type),       intent(in)  :: m_atom
      integer,                intent(in)  :: n_so
      real,dimension(3)                   :: m_mom
      
      !  Local variables
      integer                 :: c,l,m
      complex                 :: ci
      complex,dimension(3)    :: Sk
      
      m = 1 ! There is only 1 Irrep
      Sk = cmplx(0,0)
      do l=1,abs(magn_model%Nbas(m))
         c = magn_model%icomp(l,m)
         ci = cmplx(1-c,c)
         Sk(:) = Sk(:) + ci*m_atom%cbas(l,m)*magn_model%basf(:,l,n_so,m)
      end do
      m_mom = real(Sk)
      return
   end function magn_moment
   
   function magn_str_factor(m_refl) result(m_str_f)
      !  Arguments
      type(magh_type_xrms),      intent(in)     :: m_refl
      complex                                   :: m_str_f
      
      !  Local variables
      integer                 :: i, j
      real,dimension(3)       :: x_atom, m_mom, m_mom_star, m_mom_lab, z_mom_lab 
      complex                 :: scatt_ampl
      real                    :: phase
      
      m_str_f = cmplx(0.0)
      do i=1,magn_atom_list%natoms
         do j=1,magn_model%numops
            m_mom = magn_moment(magn_atom_list%atom(i),j)
            m_mom_star = matmul(matrix_mom_star,m_mom)
            m_mom_lab = matmul(m_refl%tmatrix,m_mom_star)
            z_mom_lab = m_mom_lab/euclidean_norm(3,m_mom_lab)
            scatt_ampl = cmplx(0,1)*dot_product(m_refl%kf,z_mom_lab)
            x_atom = applyso(magn_model%symop(j),magn_atom_list%atom(i)%x)
            phase = tpi*dot_product(m_refl%h,x_atom)
            m_str_f = m_str_f + scatt_ampl*exp(cmplx(0,1)*phase)
         end do
      end do
      return
   end function

   function jacobian() result(jacob)
      real(kind=cp),dimension(:,:),allocatable     :: jacob
      
      !  Local variables
      integer                 :: i, j, k, l, r, c
      complex                 :: m_str_f, DfmiDxj, ci, DfxrmsiDxj
      complex,dimension(3)    :: DmomDxj
      real,dimension(3)       :: m_mom, z_mom, DznDxj, DznDxj_lab, x_atom
      real                    :: norm_m, phase
      type(magh_type_xrms)    :: m_refl
      real,dimension(3,3)     :: matr_mom_lab
      !  The variables that contain the different derivatives:
      !  DfmiDxj: derivative of F_m(i) (magnetic str factor) w.r.t. x(j)
      !  DfxrmsiDxj: derivative of f^XRMS(i) (XRMS amplitude) w.r.t. x(j)
      !  DmomDxj: derivative of m_mom (magnetic moment) w.r.t. x(j)
      !  DznDxj: derivative of z_mom (unitary magnetic moment) w.r.t. x(j)
      !  DznDxj_lab: DznDxj in laboratory reference frame
      
      allocate(jacob(dat%nobs,vs%np))
      
      !  Jacobian for yc of the reflections
      do i=1,sample_refl_list%nref
         m_refl = sample_refl_list%mh(i)
         m_str_f = magn_str_factor(sample_refl_list%mh(i))
         jacob(i,1) = m_str_f*conjg(m_str_f)
         
         do j=2,vs%np
            DfmiDxj = 0.0
            r = 1 ! There is only 1 Irrep
            do k=1,magn_atom_list%natoms
               do l=1,magn_model%numops
                  c = magn_model%icomp(j,r)
                  ci = cmplx(1-c,c)
                  DmomDxj(:) = real(ci*magn_model%basf(:,j,l,r))
                  m_mom = magn_moment(magn_atom_list%atom(k),l)
                  norm_m = euclidean_norm(size(m_mom),m_mom)
                  z_mom = m_mom/norm_m
                  DznDxj(:) = (DmomDxj(:) - z_mom(:)*dot_product(z_mom,DmomDxj))/norm_m
                  matr_mom_lab = matmul(m_refl%tmatrix,matrix_mom_star)
                  DznDxj_lab = matmul(matr_mom_lab,DznDxj)
                  DfxrmsiDxj = cmplx(0,1)*dot_product(m_refl%kf,DznDxj_lab)
                  x_atom = applyso(magn_model%symop(l),magn_atom_list%atom(k)%x)
                  phase = dot_product(m_refl%h,x_atom)
                  DfmiDxj = DfmiDxj + DfxrmsiDxj*exp(tpi*cmplx(0,1)*phase)
               end do
            end do
            jacob(i,j) = 2*scale_fact*real(conjg(m_str_f)*DfmiDxj)
         end do
      end do
      
      !  Jacobian for the restriction
      jacob(dat%nobs,1) = 0.0
      do j=2,vs%np
         m_mom = magn_moment(magn_atom_list%atom(1),1)
         DmomDxj(:) = real(ci*magn_model%basf(:,j,1,r))
         jacob(dat%nobs,j) = 2*dot_product(m_mom,DmomDxj)
      end do
      return
   end function jacobian
   
   subroutine intensities_der(m,n,x,fvec,fjac,iflag)
      integer,                      intent(in)       :: m, n
      real(kind=cp),dimension(:),   intent(in)       :: x
      real(kind=cp),dimension(:),   intent(in out)   :: fvec
      real(kind=cp),dimension(:,:), intent(out)      :: fjac
      integer,                      intent(in out)   :: iflag
      
      !  Local variables
      integer                          :: i,j,l
      real(kind=cp)                    :: chi, chi_old = 1.0e20
      real,dimension(3)                :: m_mom
      complex                          :: m_str_f
      real,dimension(m,n)              :: DyciDxj
      !  The variable that contains the Jacobian:
      !  DyciDxj: derivative of yc(i) w.r.t. x(j)

      !  Update the state vector with the values in array x
      do i=1,vs%np
         vs%pv(i) = x(i)
      end do
      
      !  Update scale factor
      scale_fact = vs%pv(1)
      
      !  Update the magnetic model parameters (mixing coefficients)
      do i=2,vs%np
         magn_atom_list%atom(1)%cbas(i-1,1) = vs%pv(i)
      end do

      if (iflag == 1) then
         chi = 0.0
         do i=1,sample_refl_list%nref
            m_str_f = magn_str_factor(sample_refl_list%mh(i))
            dat%yc(i) = scale_fact*m_str_f*conjg(m_str_f)
            fvec(i) = (dat%y(i)-dat%yc(i))/dat%sw(i)
            chi = chi + fvec(i)*fvec(i)
         end do
         
         !  Value of the restriction: |m_mom|^2 = 1
         !m_mom = magn_moment(magn_atom_list%atom(1),1)
         !dat%yc(m) = dot_product(m_mom,m_mom)
         !fvec(m) = (dat%y(m)-dat%yc(m))/dat%sw(m)
         !chi = chi + fvec(m)*fvec(m)
         
         chi = chi/real(m-n)
         if (chi <= chi_old) then
            cond%nfev = cond%nfev + 1
            write(unit=*,fmt='(a,i4,a,g16.8)') '=> Iteration number: ', cond%nfev, '  Chi2 = ', chi
            write(unit=luo,fmt='(a,i4,a,g16.8)') '=> Iteration number: ', cond%nfev, '  Chi2 = ', chi
            chi_old = chi
            write(unit=*,fmt='(a)') 'Current parameters:'
            write(unit=*,fmt='(4(f12.6,2x))') vs%pv(1:vs%np)
            write(unit=*,fmt='(a)') ''
         end if
      elseif (iflag == 2) then
         DyciDxj = jacobian()
         do i=1,m
            fjac(i,:) = -DyciDxj(i,:)/dat%sw(i)
         end do
         cond%njev = cond%njev + 1
      end if
      return
   end subroutine intensities_der

   subroutine intensities_noder(m,n,x,fvec,iflag)
      integer,                      intent(in)       :: m, n
      real(kind=cp),dimension(:),   intent(in)       :: x
      real(kind=cp),dimension(:),   intent(in out)   :: fvec
      integer,                      intent(in out)   :: iflag
      
      !  Local variables
      integer                          :: i,l,j
      real(kind=cp)                    :: chi, chi_old = 1.0e20, abs_lor
      real,dimension(3)                :: m_mom, m_mom_lab, bstar_lab
      complex                          :: m_str_f
      real,dimension(3,3)              :: matr_mom_lab, matr_lab, matr_det

      !  Update the state vector with the values in array x
      do i=1,vs%np
         vs%pv(i) = x(i)
      end do
      
      !  Update scale factor
      scale_fact = vs%pv(1)
      
      !  Update the magnetic model (mixing coefficients)
      do i=2,vs%np
         magn_atom_list%atom(1)%cbas(i-1,1) = vs%pv(i)
      end do

      chi = 0.0
      do i=1,sample_refl_list%nref
         m_str_f = magn_str_factor(sample_refl_list%mh(i))
         abs_lor = factor_abs_lor(sample_refl_list%mh(i))
         dat%yc(i) = scale_fact*m_str_f*conjg(m_str_f)
         fvec(i) = (dat%y(i)-dat%yc(i))/dat%sw(i)
         chi = chi + fvec(i)*fvec(i)
      end do
      
      !  Value of the first restriction: m(perp) = 0
      !m_mom = magn_moment(magn_atom_list%atom(1),1)
      !matr_lab = sample_refl_list%mh(1)%tmatrix
      !matr_mom_lab = matmul(matr_lab,matrix_mom_star)
      !m_mom_lab = matmul(matr_mom_lab,m_mom)
      !bstar_lab = matmul(matr_lab,(/ 0, 1, 0 /))
      !matr_det(1,:) = m_mom_lab(:)
      !matr_det(2,:) = sample_refl_list%mh(1)%kf(:)
      !matr_det(3,:) = bstar_lab(:)
      !call determinant(matr_det,3,dat%yc(m-1))
      !fvec(m-1) = (dat%y(m-1)-dat%yc(m-1))/dat%sw(m-1)
      !chi = chi + fvec(m-1)*fvec(m-1)
      
      !  Value of the second restriction: |m_mom|^2 = 1
      m_mom = magn_moment(magn_atom_list%atom(1),1)
      matr_lab = sample_refl_list%mh(1)%tmatrix
      matr_mom_lab = matmul(matr_lab,matrix_mom_star)
      m_mom_lab = matmul(matr_mom_lab,m_mom)
      dat%yc(m) = dot_product(m_mom_lab,m_mom_lab)
      fvec(m) = (dat%y(m)-dat%yc(m))/dat%sw(m)
      chi = chi + fvec(m)*fvec(m)
      
      chi = chi/real(m-n)
      if (chi <= chi_old) then
         cond%nfev = cond%nfev + 1
         write(unit=*,fmt='(a,i4,a,g16.8)') '=> Iteration number: ', cond%nfev, '  Chi2 = ', chi
         write(unit=luo,fmt='(a,i4,a,g16.8)') '=> Iteration number: ', cond%nfev, '  Chi2 = ', chi
         chi_old = chi
         write(unit=*,fmt='(a)') 'Current parameters:'
         write(unit=*,fmt='(4f12.6,2x)') vs%pv(1:vs%np)
         write(unit=*,fmt='(a)') ''
      end if
      return
   end subroutine intensities_noder

   subroutine write_intensities()
      
      !  Local variables
      integer        :: i
      complex        :: m_str_f
      
      do i=1,sample_refl_list%nref
         m_str_f = magn_str_factor(sample_refl_list%mh(i))
         dat%yc(i) = vs%pv(1)*m_str_f*conjg(m_str_f)
      end do
      
      write(unit=luo,fmt='(a)') ''
      write(unit=luo,fmt='(a)') '===================================='
      write(unit=luo,fmt='(a)') '======= Magnetic Intensities ======='
      write(unit=luo,fmt='(a)') '===================================='
      write(unit=luo,fmt='(a)') ''
      
      write(unit=luo,fmt='(2a6,3a12)') 'data', 'h', 'y_obs', 'sigma(obs)', 'y_calc'
      do i=1,sample_refl_list%nref
         write(unit=luo,fmt='(6x,f6.2,3f12.4)') sample_refl_list%mh(i)%h(2), dat%y(i), dat%sw(i), dat%yc(i)
      end do
      write(unit=luo,fmt='(a12)') ''
      write(unit=luo,fmt='(a12)') 'restrictions'
      do i=dat%nobs-2,dat%nobs
         write(unit=luo,fmt='(6x,f6.2,3f12.6)') dat%x(i), dat%y(i), dat%sw(i), dat%yc(i)
      end do
      write(unit=luo,fmt='(a12)') ''
      write(unit=luo,fmt='(a12)') 'parameters'
      write(unit=luo,fmt='(12x,f12.6)') vs%pv(1:vs%np)
      close(unit=luo)
      return
   end subroutine write_intensities

end module XRMS_procedures
