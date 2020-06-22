SubModule (CFML_HDF5) HDF5_D19
    Contains

    !!----
    !!---- Read_NXS_Numor_D19
    !!----
    !!---- Read a nexus file of d19 and build the corresponding numor.
    !!----
    !!---- 04/05/20
    !!
    Module Subroutine Read_Nexus_D19(filename,numor,counts)
        !---- Arguments ----!
        character(len=*),                       intent(in)  :: filename
        type(SXTAL_NUMOR_type),                 intent(out) :: numor
        integer, dimension(:,:,:), allocatable, intent(out) :: counts

        !---- Local variables ----!
        integer :: i,j,k,i_om,i_ti,i_mo,i_mc
        integer :: hdferr
        integer :: inum,nfrm,nrows,ncols
        real :: wave,start,step,width,&
                phi,chi,gamma,omega
        real,    dimension(9)                  :: ub
        integer, dimension(:),     allocatable :: axis
        real,    dimension(:,:),   allocatable :: scan
        integer, dimension(:,:,:), allocatable :: cnts

        !---- Local variables with hdf5 types ----!
        integer(SIZE_T), PARAMETER :: sdim = 20 ! maximum string length
        integer(HID_T) :: file_id,inum_data_id,nfrm_data_id,wave_data_id,&
                          scan_data_id,axis_data_id,cnts_data_id,&
                          start_data_id,step_data_id,width_data_id,&
                          ub_data_id,nrows_data_id,ncols_data_id,&
                          phi_data_id,chi_data_id,gamma_data_id,&
                          omega_data_id,name_data_id,memtype,filetype
        integer(HID_T) :: axis_space_id,scan_space_id,cnts_space_id,name_space_id
        integer(HSIZE_T), dimension(1) :: scalar,axis_dim,axis_maxdim,&
                                          ub_dim,ub_maxdim,name_dim,name_maxdim
        integer(HSIZE_T), dimension(2) :: scan_dims,scan_maxdims
        integer(HSIZE_T), dimension(3) :: cnts_dims,cnts_maxdims

        character(len=sdim), dimension(:), allocatable :: name_
        integer(HSIZE_T), dimension(2) :: name_dims
        integer(SIZE_T), dimension(:), allocatable :: str_len

        !---- Initialize variables
        ub_dim(1) = 9
        ub_maxdim(1) = 9
        if (allocated(cnts)) deallocate(cnts)

        !---- Initialize fortran interface
        call h5open_f(hdferr)

        !---- Open NEXUS file
        write(*,'(2A)') " => Opening Nexus file ",trim(filename)
        write(*,'(A)')  "    Reading data..."
        call h5fopen_f(trim(filename),H5F_ACC_RdoNLY_F,file_id,hdferr)

        !---- Open datasets
        call h5dopen_f(file_id,'entry0/run_number',inum_data_id,hdferr)
        call h5dopen_f(file_id,'entry0/data_scan/actual_step',nfrm_data_id,hdferr)
        call h5dopen_f(file_id,'entry0/wavelength',wave_data_id,hdferr)
        call h5dopen_f(file_id,'entry0/instrument/ScanInfo/first_start_value',start_data_id,hdferr)
        call h5dopen_f(file_id,'entry0/instrument/ScanInfo/first_step_value',step_data_id,hdferr)
        call h5dopen_f(file_id,'entry0/instrument/ScanInfo/first_width_value',width_data_id,hdferr)
        call h5dopen_f(file_id,'entry0/instrument/SingleCrystalSettings/orientation_matrix',ub_data_id,hdferr)
        call h5dopen_f(file_id,'entry0/instrument/Det1/nrows',nrows_data_id,hdferr)
        call h5dopen_f(file_id,'entry0/instrument/Det1/ncols',ncols_data_id,hdferr)
        call h5dopen_f(file_id,'entry0/instrument/gamma/value',gamma_data_id,hdferr)
        call h5dopen_f(file_id,'entry0/instrument/chi/value',chi_data_id,hdferr)
        call h5dopen_f(file_id,'entry0/instrument/phi/value',phi_data_id,hdferr)
        call h5dopen_f(file_id,'entry0/instrument/omega/value',omega_data_id,hdferr)
        call h5dopen_f(file_id,'entry0/data_scan/scanned_variables/variables_names/axis',axis_data_id,hdferr)
        call h5dopen_f(file_id,'entry0/data_scan/scanned_variables/variables_names/name',name_data_id,hdferr)
        call h5dopen_f(file_id,'entry0/data_scan/scanned_variables/data',scan_data_id,hdferr)
        call h5dopen_f(file_id,'entry0/data_scan/detector_data/data',cnts_data_id,hdferr)

        !---- Get type of string datasets
        call h5dget_type_f(name_data_id, filetype, hdferr)

        !---- Get dimensions of datasets
        call h5dget_space_f(axis_data_id,axis_space_id,hdferr)
        call h5dget_space_f(name_data_id,name_space_id,hdferr)
        call h5dget_space_f(scan_data_id,scan_space_id,hdferr)
        call h5dget_space_f(cnts_data_id,cnts_space_id,hdferr)
        call h5sget_simple_extent_dims_f(axis_space_id,axis_dim,axis_maxdim,hdferr)
        call h5sget_simple_extent_dims_f(name_space_id,name_dim,name_maxdim,hdferr)
        call h5sget_simple_extent_dims_f(scan_space_id,scan_dims,scan_maxdims,hdferr)
        call h5sget_simple_extent_dims_f(cnts_space_id,cnts_dims,cnts_maxdims,hdferr)

        !---- Assign memory to arrays
        allocate(axis(axis_dim(1)))
        allocate(name_(name_dim(1)))
        allocate(scan(scan_dims(1),scan_dims(2)))
        allocate(cnts(cnts_dims(1),cnts_dims(2),cnts_dims(3)))
        allocate(str_len(name_dim(1)))
        str_len(:) = sdim
        name_dims(:) = (/ sdim, name_dim(1) /)

        !---- Read datasets
        call h5dread_f(inum_data_id,H5T_NATIVE_INTEGER,inum,scalar,hdferr)
        call h5dread_f(nfrm_data_id,H5T_NATIVE_INTEGER,nfrm,scalar,hdferr)
        call h5dread_f(wave_data_id,H5T_NATIVE_REAL,wave,scalar,hdferr)
        call h5dread_f(start_data_id,H5T_NATIVE_REAL,start,scalar,hdferr)
        call h5dread_f(step_data_id,H5T_NATIVE_REAL,step,scalar,hdferr)
        call h5dread_f(width_data_id,H5T_NATIVE_REAL,width,scalar,hdferr)
        call h5dread_f(ub_data_id,H5T_NATIVE_REAL,ub,ub_dim,hdferr)
        call h5dread_f(gamma_data_id,H5T_NATIVE_REAL,gamma,scalar,hdferr)
        call h5dread_f(chi_data_id,H5T_NATIVE_REAL,chi,scalar,hdferr)
        call h5dread_f(phi_data_id,H5T_NATIVE_REAL,phi,scalar,hdferr)
        call h5dread_f(omega_data_id,H5T_NATIVE_REAL,omega,scalar,hdferr)
        call h5dread_f(nrows_data_id,H5T_NATIVE_INTEGER,nrows,scalar,hdferr)
        call h5dread_f(ncols_data_id,H5T_NATIVE_INTEGER,ncols,scalar,hdferr)
        call h5dread_f(axis_data_id,H5T_NATIVE_INTEGER,axis,axis_dim,hdferr)
        call h5dread_f(scan_data_id,H5T_NATIVE_REAL,scan,scan_dims,hdferr)
        call h5dread_f(cnts_data_id,H5T_NATIVE_INTEGER,cnts,cnts_dims,hdferr)

        !---- Read keywords of the data scan
        call h5dread_vl_f(name_data_id,filetype,name_,name_dims,str_len,hdferr,name_space_id)

        !---- Find keywords
        i_om = -1
        i_ti = -1
        i_mo = -1
        i_mc = -1
        do i = 1 , scan_dims(2)
            select case(l_case(trim(name_(i))))
                case("omega")
                    i_om = i
                case("acquisitionspy") ! time
                    i_ti = i
                case("monitor1")
                    i_mo = i
                case("multicalib") ! total counts
                    i_mc = i
            end select
        end do
        if (i_om == -1) then
            err_CFML%Ierr = 1
            err_CFML%Msg = "Error reading numor. Omega not found in scanned variables"
            return
        end if

        !---- Close datasets
        call h5dclose_f(inum_data_id,hdferr)
        call h5dclose_f(nfrm_data_id,hdferr)
        call h5dclose_f(wave_data_id,hdferr)
        call h5dclose_f(start_data_id,hdferr)
        call h5dclose_f(step_data_id,hdferr)
        call h5dclose_f(width_data_id,hdferr)
        call h5dclose_f(ub_data_id,hdferr)
        call h5dclose_f(nrows_data_id,hdferr)
        call h5dclose_f(ncols_data_id,hdferr)
        call h5dclose_f(axis_data_id,hdferr)
        call h5dclose_f(name_data_id,hdferr)
        call h5dclose_f(scan_data_id,hdferr)
        call h5dclose_f(cnts_data_id,hdferr)
        call h5dclose_f(gamma_data_id,hdferr)
        call h5dclose_f(chi_data_id,hdferr)
        call h5dclose_f(phi_data_id,hdferr)
        call h5dclose_f(omega_data_id,hdferr)

        !---- Close NEXUS file.
        write(*,'(2A)') "    Closing Nexus file ",trim(filename)
        call h5fclose_f(file_id,hdferr)

        !---- Close FORTRAN interface.
        call h5close_f(hdferr)

        !---- Build the Numor object
        !     Only omega scans can be processed for the moment.
        !     Raise an error otherwise
        write(*,'(A)') " => Building Numor..."
        call Initialize_Numor(numor)

        numor%numor     = inum
        numor%instrm    = 'D19'
        numor%scantype  = 'omega'
        numor%manip     = 2
        numor%nbang     = 1
        numor%wave      = wave
        numor%nframes   = nfrm
        numor%nbdata    = nrows * ncols
        numor%angles(1) = phi
        numor%angles(2) = chi
        numor%angles(3) = omega
        numor%angles(4) = gamma
        numor%scans(:)  = [ start,step,width ]
        numor%ub        = transpose(RESHAPE(ub,[3,3]))
        allocate(numor%tmc_ang(numor%nbang+3,numor%nframes))
        if (i_ti == -1) then
            write(*,'(2a)')  " => WARNING: time (acquisitionspy) not found",&
            " in scanned variables..."
            numor%tmc_ang(1,:) = 0
        else
            numor%tmc_ang(1,:) = scan(:,i_ti)
        end if
        if (i_mo == -1) then
            write(*,'(2a)')  " => WARNING: monitor (monitor1) not found",&
            " in scanned variables..."
            numor%tmc_ang(2,:) = 0
        else
            numor%tmc_ang(2,:) = scan(:,i_mo)
        end if
        if (i_mc == -1) then
            write(*,'(2a)')  " => WARNING: total counts (multicalib)",&
            " not found in scanned variables..."
            numor%tmc_ang(3,:) = 0
        else
            numor%tmc_ang(3,:) = scan(:,i_mc)
        end if
        numor%tmc_ang(4,:) = scan(:,i_om) ! angle
        allocate(counts(nrows,ncols,nfrm))
        do k = 1 , nfrm
            do j = 1 , ncols
                do i = 1 , nrows
                    counts(i,j,k) = cnts(nrows-i+1,j,k) ! Flip
                end do
            end do
        end do
        deallocate(cnts)

    End Subroutine Read_Nexus_D19

End SubModule HDF5_D19