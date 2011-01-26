!     Last change:  TR    5 Dec 2006    5:26 pm

 ! lecture fichier final.y cree par EVALCCD

!------------------------------------------------------------------------------
subroutine get_cell_parameters_from_FINAL_Y_file(file_unit, cell_param)
 use IO_module
 implicit none
  INTEGER, INTENT(IN)                  :: file_unit
  REAL,    INTENT(OUT),   dimension(6) :: cell_param
  character (len=256)                  :: line, new_line
  integer                              :: ier, i1, i2
  CHARACTER (LEN=32)                   :: required_string
  integer                              :: long_string


  do
   read(unit=file_unit, '(a)', iostat=ier) line
   if(ier<0) then
    call write_info('')
    call write_info(' ... No cell parameters in the final.y file !!')
    call write_info('')
    stop
   ELSEIF(ier >0) then
    call write_info('')
    call write_info(' ... Problem reading final.y file !!')
    call write_info('')
    stop
   endif

   required_string = '<DIMENSIONS TYPE="standard"'
   long_string = LEN_TRIM(required_string)
   IF(TRIM(line(1:long_string)) == required_string(1:long_string)) then
    line = line(long_string+1:)
    line = ADJUSTL(line)

    call get_real_value_from_FINAL_Y('A="', cell_param(1), line)
    call get_real_value_from_FINAL_Y('B="', cell_param(2), line)
    call get_real_value_from_FINAL_Y('C="', cell_param(3), line)
    call get_real_value_from_FINAL_Y('ALPHA="', cell_param(4), line)
    call get_real_value_from_FINAL_Y('BETA="', cell_param(5), line)
    call get_real_value_from_FINAL_Y('GAMMA="', cell_param(6), line)

    exit
   endif

  end do



 RETURN
end subroutine get_cell_parameters_from_FINAL_Y_file

!----------------------------------------------------------------------
subroutine get_modulation_vect(file_unit, qvec_type, QVEC)
 use IO_module
 
 
 implicit none
  INTEGER, INTENT(IN)                  :: file_unit
  LOGICAL, INTENT(OUT)                 :: QVEC_type
  real, dimension(3), INTENT(inout)    :: QVEC
  
  character (len=256)                  :: line, new_line
  integer                              :: ier, i1, i2
  CHARACTER (LEN=32)                   :: required_string
  integer                              :: long_string


  do
   read(unit=file_unit, '(a)', iostat=ier) line
   if(ier<0) then
    QVEC_type = .false.
   ELSEIF(ier >0) then
    call write_info('')
    call write_info(' ... Problem reading final.y file !!')
    call write_info('')
    stop
   endif

   required_string = '<QVEC TYPE="MODULATED"'
   long_string = LEN_TRIM(required_string)
   IF(TRIM(line(1:long_string)) == required_string(1:long_string)) then
    line = line(long_string+1:)
    line = ADJUSTL(line)

    call get_real_value_from_FINAL_Y('DH="', QVEC(1), line)
    call get_real_value_from_FINAL_Y('DK="', QVEC(2), line)
    call get_real_value_from_FINAL_Y('DL="', QVEC(3), line)
    QVEC_type = .true.
    exit
   endif

  end do



 RETURN


end subroutine get_modulation_vect


!----------------------------------------------------------------------
subroutine get_real_value_from_FINAL_Y(required_string, real_value, line)
 CHARACTER (LEN=*), INTENT(IN)    :: required_string
 REAL,              INTENT(OUT)   :: real_value
 CHARACTER (LEN=*), INTENT(INOUT) :: line
 !local variables
 INTEGER                  :: long_string
 INTEGER                  :: i1, i2
 CHARACTER (LEN=256)      :: new_line


  long_string = LEN_TRIM(required_string)
  i1 = INDEX(line, required_string(1:long_string))
  new_line = line(long_string+1:)
  i2 = INDEX(new_line, '"')
  READ(line(i1+long_string:i2-1+long_string), *) real_value
  line = new_line(i2+1:)
  line = ADJUSTL(line)

  return
end subroutine get_real_value_from_FINAL_Y
!----------------------------------------------------------------------
subroutine get_integer_value_from_FINAL_Y(required_string, integer_value, line)
 CHARACTER (LEN=*), INTENT(IN)    :: required_string
 integer,           INTENT(OUT)   :: integer_value
 CHARACTER (LEN=*), INTENT(INOUT) :: line
 !local variables
 INTEGER                  :: long_string
 INTEGER                  :: i1, i2
 CHARACTER (LEN=256)      :: new_line


  long_string = LEN_TRIM(required_string)
  i1 = INDEX(line, required_string(1:long_string))
  new_line = line(long_string+1:)
  i2 = INDEX(new_line, '"')
  READ(line(i1+long_string:i2-1+long_string), *) integer_value
  line = new_line(i2+1:)
  line = ADJUSTL(line)

  return
end subroutine get_integer_value_from_FINAL_Y

!----------------------------------------------------------------------

subroutine read_reflexion_from_FINAL_Y(file_unit, H_h, H_k, H_l, H_F2, H_sig, H_ok)
 INTEGER, INTENT(IN)  :: file_unit
 INTEGER, INTENT(OUT) :: H_h, H_k, H_l
 REAL,    INTENT(OUT) :: H_F2, H_sig
 LOGICAL, INTENT(OUT) :: H_ok
 !local variables
  character (len=256)                  :: line, new_line
  integer                              :: ier, i1, i2
  CHARACTER (LEN=32)                   :: required_string
  integer                              :: long_string

 do
   read(unit=file_unit, '(a)', iostat=ier) line
   IF(ier <0) then
    H_ok = .false.
    return
   endif

   required_string = '<INDEX'
   long_string = LEN_TRIM(required_string)
   IF(TRIM(line(1:long_string)) == required_string(1:long_string)) then
    line = line(long_string+1:)
    line = ADJUSTL(line)

    call get_integer_value_from_FINAL_Y('H="', H_h, line)
    call get_integer_value_from_FINAL_Y('K="', H_k, line)
    call get_integer_value_from_FINAL_Y('L="', H_l, line)

   endif

   required_string = '<INTENSITY'
   long_string = LEN_TRIM(required_string)
   IF(TRIM(line(1:long_string)) == required_string(1:long_string)) then
    line = line(long_string+1:)
    line = ADJUSTL(line)

    call get_real_value_from_FINAL_Y('I="',     H_F2,  line)
    call get_real_value_from_FINAL_Y('SIGMA="', H_sig, line)

    H_ok = .true.
    return

   endif
END do

end subroutine read_reflexion_from_FINAL_Y
!----------------------------------------------------------------------

subroutine read_mod_reflexion_from_FINAL_Y(file_unit, H_h, H_k, H_l, H_m, H_F2, H_sig, H_ok)
 INTEGER, INTENT(IN)  :: file_unit
 INTEGER, INTENT(OUT) :: H_h, H_k, H_l, H_m
 REAL,    INTENT(OUT) :: H_F2, H_sig
 LOGICAL, INTENT(OUT) :: H_ok
 !local variables
  character (len=256)                  :: line, new_line
  integer                              :: ier, i1, i2
  CHARACTER (LEN=32)                   :: required_string
  integer                              :: long_string

 do
   read(unit=file_unit, '(a)', iostat=ier) line
   IF(ier <0) then
    H_ok = .false.
    return
   endif

   required_string = '<INDEX'
   long_string = LEN_TRIM(required_string)
   IF(TRIM(line(1:long_string)) == required_string(1:long_string)) then
    line = line(long_string+1:)
    line = ADJUSTL(line)

    call get_integer_value_from_FINAL_Y('H="',  H_h, line)
    call get_integer_value_from_FINAL_Y('K="',  H_k, line)
    call get_integer_value_from_FINAL_Y('L="',  H_l, line)
    call get_integer_value_from_FINAL_Y('M1="', H_m, line)

   endif

   required_string = '<INTENSITY'
   long_string = LEN_TRIM(required_string)
   IF(TRIM(line(1:long_string)) == required_string(1:long_string)) then
    line = line(long_string+1:)
    line = ADJUSTL(line)

    call get_real_value_from_FINAL_Y('I="',     H_F2,  line)
    call get_real_value_from_FINAL_Y('SIGMA="', H_sig, line)

    H_ok = .true.
    return

   endif
END do


end subroutine read_mod_reflexion_from_FINAL_Y
!----------------------------------------------------------------------


subroutine get_wavelength_from_FINAL_Y_file(file_unit, wave)
 USE wavelength_module
 USE IO_module,       only : write_info
 implicit none
  INTEGER, INTENT(IN)                  :: file_unit
  REAL,    INTENT(OUT)                 :: wave
  character (len=256)                  :: line
  CHARACTER (LEN=256)                  :: champ_text
  REAL                                 :: champ_value
  integer                              :: i, ier


  do
   read(unit=file_unit, '(a)', iostat=ier) line
   if(ier<0) then
    call write_info('')
    call write_info(' ... No wavelength  in the final.y file !!')
    call write_info('')
    stop
   ELSEIF(ier >0) then
    call write_info('')
    call write_info(' ... Problem reading wavelength value in final.y file !!')
    call write_info('')
    stop
   endif

   line = ADJUSTL(line)


   if(trim(line(1:13))=='<SOURCE TYPE=') then
    IF(INDEX(line, 'TARGET="AG"') /=0) then
     wave = X_target(1)%wave(1)
    elseIF(INDEX(line, 'TARGET="MO"') /=0) then
     wave = X_target(2)%wave(1)
    ELSEIF(INDEX(line, 'TARGET="CU"') /=0) then
     wave = X_target(3)%wave(1)
    ELSEIF(INDEX(line, 'TARGET="NI"') /=0) then
     wave = X_target(4)%wave(1)
    ELSEIF(INDEX(line, 'TARGET="CO"') /=0) then
     wave = X_target(5)%wave(1)
    ELSEIF(INDEX(line, 'TARGET="FE"') /=0) then
     wave = X_target(6)%wave(1)
    ELSEIF(INDEX(line, 'TARGET="CR"') /=0) then
     wave = X_target(7)%wave(1)

    else
     call write_info('')
     call write_info(' ... No available wavelength  in the final.y file !!')
     call write_info('')
     stop
    endif
    exit
   endif
  END do

 return
end subroutine get_wavelength_from_FINAL_Y_file

!--------------------------------------------------------------------------


