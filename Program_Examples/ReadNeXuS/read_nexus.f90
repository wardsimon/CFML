program test_nexus_opening
    !---- Use Modules ----!
    Use CFML_ILL_Instrm_Data

    !---- Local Variables ----!
    character(len=256)     :: args
    character(len=256)      :: fname
    character(len=4)       :: instrm
    integer                :: i,j, n, nsep, iv, iattr
    type(powder_numor_type):: numpwd
    type(SXTAL_numor_type) :: numsxt
    logical success
    integer, dimension(:), allocatable    :: frames
    integer, dimension(:,:,:), allocatable    :: counts

    if (iargc()/=2) then
      write(*,*) "this program requires two arguments (file_path, instrument_name)"
      stop 9
    end if
    call getarg(1, args)
    fname=trim(args)
    call getarg(2, args)
    instrm=trim(args)

	write(*,*) "======= Reading header"
    call read_numor(trim(fname),trim(instrm),numsxt,.true.)
    if (abs(numsxt%wave - 1.45000005) > 0.00000001) then
      write(*,*) "Error in reading wavelength"
      stop 9
    endif

    write(*,*) "======= Reading header and data"
    call read_numor(trim(fname),trim(instrm),numsxt,.false.,counts=counts)
    numrows=size(counts, 1)
    numcols=size(counts, 2)
    numframes=size(counts, 3)
    do i=1,numframes
      if (abs(numsxt%tmc_ang(3,i) - sum(counts(:,:,i))) > 0) then
        write(*,*) "Error in reading counts"
        stop 9
    endif
    end do

    stop 0


end program  test_nexus_opening
