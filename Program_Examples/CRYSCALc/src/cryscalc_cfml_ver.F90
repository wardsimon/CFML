

Program create_cryscalc_version

 ! creation d'un fichier fortran cryscalc_version.F90 contenant
 ! la definition des variables suivantes :
 !  . version de CRYSCALC
 !  . compilateur utiise
 !  . version de CRYSFML utilisee

 implicit none

  character (len=32), dimension(3)      :: file_F90
  character (len=32)                    :: F90_file
  CHARACTER (LEN=10)                    :: date, time
  character (len=256)                   :: CRYSCALC_date, CFML_env, CFML_html, CFML_version
  character (len=256)                   :: archi
  character (len=256)                   :: read_line
  logical                               :: file_exist, archi64
  integer                               :: i, i1, i2


  file_F90(1) = "cryscalc_ver_lf95.f90"
  file_F90(2) = "cryscalc_ver_gfort.f90"
  file_F90(3) = "cryscalc_ver_ifort.f90"
  F90_file    = "cc_cfml_ver.f90"

  call date_and_time(date, time)
  WRITE(CRYSCALC_date, '(11a)') date(7:8),'-', date(5:6),'-', date(1:4), '  at ', &
                                time(1:2),':', time(3:4),':', time(5:6)

  call getenv('CRYSFML', CFML_env)

  write(CFML_html, '(2a)') trim(CFML_env), '\html\files\Welcome.html'
  call test_file_exist(trim(CFML_html), file_exist)

  if(.not. file_exist) then
   CFML_version = "?"

  else
   open(unit=1, file=trim(CFML_html))
   do
    read(unit=1, fmt='(a)') read_line
    if(index(read_line, "Version:") /=0) then
     i1 = index(read_line, ":")
     i2 = index(read_line, "</span>")
     write(CFML_version, '(a)') read_line(i1+1: i2-1)
     exit
    end if
   end do
  end if
  close (unit=1)
  CFML_version = adjustl(CFML_version)

  archi64 = .false.
  call getenv('TARGET_ARCH', archi)
  if(len_trim(archi) == 7) then
   if(archi(1:7) == 'intel64') archi64 = .true.
  endif

  open (unit=1, file = trim(F90_file))
    write(unit=1, fmt='(a)') ''
    write(unit=1, fmt='(a)') 'subroutine Cryscalc_CFML_Versions'
    write(unit=1, fmt='(a)') ' use cryscalc_module, only : cryscalc, CFML_version, Architecture'
    write(unit=1, fmt='(a)') ' implicit none'
    write(unit=1, fmt='(a)') '  '
    write(unit=1, fmt='(4a)') '   cryscalc%date       = ', '"',trim(cryscalc_date),'"'
    write(unit=1, fmt='(4a)') '   CFML_version        = ', '"',trim(CFML_version),' (JRC, JGP)"'
    if(archi64) then
     write(unit=1, fmt='(4a)') '   Architecture        = "64 bits"'
    else
     write(unit=1, fmt='(4a)') '   Architecture        = "32 bits"'
    end if
    write(unit=1, fmt='(a)') '  '
    write(unit=1, fmt='(a)') ' return '
    write(unit=1, fmt='(a)') 'end subroutine Cryscalc_CFML_Versions'
    write(unit=1, fmt='(a)') ''
  close (unit=1)

 stop
end program create_cryscalc_version



!-------------------------------------------------------------------
subroutine test_file_exist(file_name, file_exist)
  CHARACTER(LEN=*), INTENT(IN)          :: file_name
  LOGICAL,          INTENT(INOUT)       :: file_exist
   character(len=256)                    :: text_info

  inquire (FILE=TRIM(file_name), EXIST=file_exist)

  if(.not. file_exist) then
   write(*,*)
   write(*, '(3a)') ' >>> ', trim(file_name), ' does not exist !'
   write(*,*)
  endif

  return
 end subroutine test_file_exist
!---------------------------------------------------------------
