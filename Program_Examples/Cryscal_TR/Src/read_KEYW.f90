!     Last change:  TR    5 Sep 2007    6:21 pm
subroutine read_input_file_KEYWORDS()
 USE cryscal_module, ONLY : input_unit, keyword_BEAM, keyword_SFAC_UNIT, keyword_CONT, keyword_CHEM, keyword_ZUNIT
 USE macros_module,  ONLY : u_case
 USE IO_module,      ONLY : write_info
 implicit none
  CHARACTER (LEN=256)                      :: read_line
  INTEGER                                  :: i_error


  REWIND(UNIT=input_unit)

 do        ! lecture du fichier d'entree
  READ(UNIT=input_unit, '(a)', IOSTAT=i_error) read_line
  IF(i_error < 0) EXIT   ! fin du fichier
  read_line = ADJUSTL(read_line)
  read_line = u_case(TRIM(read_line))

  IF (LEN_TRIM(read_line) == 0)                           CYCLE  ! ligne vide
  IF (read_line(1:1) == '! ' .or. read_line(1:1) == '#' ) CYCLE  ! commentaire
  !IF (read_line(1:1) == '_')                              CYCLE  ! cas d'un fichier.CIF

  
  call identification_keywords(read_line)
 ! call run_keyword_interactive(current_keyword)

 END do


 IF(.NOT. keyword_BEAM .AND. ( keyword_SFAC_UNIT .or. keyword_CONT .or. keyword_CHEM)) then
  call incident_beam()
 endif

 IF(keyword_CHEM) then
  IF(.NOT. keyword_ZUNIT) then
   call write_info('')
   call write_info('  >> ZUNIT keyword is missing ! ')
   call write_info('')
   stop
  endif
 endif


  RETURN
end subroutine read_input_file_KEYWORDS


!----------------------------------------------------------------------
subroutine identification_keywords(read_line)
 USE cryscal_module
 USE hkl_module,                ONLY : n_sig, threshold, ratio_criteria, requested_H, requested_H_string, search_H_string,    &
                                       HKL_list, HKL_list_NEG, HKL_list_POS,   HKL_file,                           &
                                       HKL_data_known, HKL_list_ABSENT, search_equiv, search_friedel, HKL_rule_nb

 use CFML_crystallographic_symmetry, only : set_spacegroup
 USE MATRIX_list_module
 USE macros_module
 USE IO_module,                 ONLY : write_info
 USE matrix_module,             ONLY : get_mat_from_setting
 USE wavelength_module
 USE text_module ,              ONLY : CIF_lines_nb, CIF_title_line

 implicit none
  CHARACTER (LEN=*), INTENT(IN)              :: read_line
  CHARACTER (LEN=32)                         :: input_keyword
  INTEGER                                    :: long_kw, long, long1, long2
  CHARACTER (LEN=256)                        :: arg_line, new_line, arg1
  CHARACTER (LEN=12), DIMENSION(nb_help_max) :: temp_string
  INTEGER                                    :: i, j, i1, i2, i_error, nb_arg
  INTEGER                                    :: current_i, i_ok
  REAL,               DIMENSION(10)          :: var
  CHARACTER (LEN=64), DIMENSION(20)          :: arg_string
  CHARACTER (LEN=64)                         :: arg
  LOGICAL                                    :: file_exist

  arg_string = ''

  READ(read_line, *) input_keyword
  input_keyword = ADJUSTL(input_keyword)
  long_kw = LEN_TRIM(input_keyword)
  if (input_keyword(long_kw:long_kw) == '=') then
   input_keyword = input_keyword(1:long_kw)
   long_kw = LEN_TRIM(input_keyword)
  end if

  READ(read_line(long_kw+1:),'(a)') arg_line
  arg_line = ADJUSTL(arg_line)

  long=LEN_TRIM(arg_line)
  !IF(arg_line(1:1) == "'" .AND. arg_line(long:long)=="'") arg_line = arg_line(2:long-1) ! chaine entre ''
  !call remove_car(arg_line, "'")    ! caractere "'" dans la chaine
  arg_line = remove_car(arg_line, "'")
  IF(arg_line(1:1) == '=') then
   arg_line = arg_line(2:)
   arg_line = ADJUSTL(arg_line)
  endif

  call nombre_de_colonnes(arg_line, nb_arg)
  IF(nb_arg/=0)  then
   READ(arg_line, *, IOSTAT=i_error) arg_string(1:nb_arg)

   IF(i_error /=0) then
    call write_info('')
    WRITE(message_text, '(2a)') '  input line: ', TRIM(read_line)
    call write_info(TRIM(message_text))
    WRITE(message_text, '(3a)') '  ... Error reading ', TRIM(input_keyword), ' arguments ...'
    call write_info(TRIM(message_text))
    return
   endif
   do i=1, nb_arg
    arg_string(i) = u_case(arg_string(i))
   end do
  endif


  select case (TRIM(input_keyword))

   case ('ACTA', 'CIF', 'CREATE_CIF')
    keyword_create_CIF = .true.
    OPEN(UNIT=CIF_unit, FILE='cryscal.cif')
    do i=1, CIF_lines_nb
     WRITE(UNIT=CIF_unit, '(a)') trim(CIF_title_line(i))
    end do
    WRITE(UNIT=CIF_unit, '(a)') ''
    WRITE(UNIT=CIF_unit, '(a)') 'data_cryscal'
    WRITE(UNIT=CIF_unit, '(a)') ''
    
   case ('CREATE_ACE')
    keyword_create_ACE = .true.  

   case ('CREATE_CEL')
    keyword_create_CEL = .true.  
    
   case ('CREATE_CFL')
    keyword_create_CFL = .true.

   case ('CREATE_FST')
    keyword_create_FST = .true.
        
   case ('CREATE_INS')
    keyword_create_INS = .true.  

   CASE ('LST_ATOMS', 'ATOM_LIST', 'LIST_ATOM_LIST', 'LIST_ATOMS', 'WRITE_ATOMS', 'WRITE_ATMS')
    keyword_atom_list = .true.

   case ('EDIT')
    IF(nb_arg == 0) then
     call check_arg('0', 'EDIT')
     return
    endif
    file_to_edit = arg_string(1)
    call test_file_exist(file_to_edit, file_exist)
    IF(.not. file_exist) then
     file_to_edit = ''
     return
    endif
    keyword_EDIT = .true.

   CASE ('FILE')
    HKL_file%NAME    = ''
    HKL_file%plot     = .false.
    HKL_file%read_NEG = .true.
    HKL_data_known    = .false.
    HKL_file%cif      = .false.
    HKL_file%SHELX    = .false.
    HKL_file%final_y  = .false.
    HKL_file%RAW      = .false.
    HKL_file%M91      = .false.
    HKL_file%M95      = .false.
    nb_shell = 0
    nb_sort  = 0

    IF(nb_arg == 0) then
     call check_arg('0', 'FILE')
     return
    endif
    READ(arg_string(1), *) HKL_file%NAME
    IF(nb_arg >1) then
     READ(arg_string(2), *) arg1
     arg1 = ADJUSTL(arg1)
    endif

    i = INDEX(HKL_file%NAME, '.')
    if (i==0) then
     HKL_file%NAME = TRIM(HKL_file%NAME)//'.hkl'
     HKL_file%SHELX = .true.
    else
     long = LEN_TRIM(HKL_file%NAME)
     IF(HKL_file%NAME(long-2:long) == 'HKL')     HKL_file%SHELX   = .true.
     IF(HKL_file%NAME(long-2:long) == 'CIF')     HKL_file%CIF     = .true.
     IF(HKL_file%NAME(long-2:long) == 'RAW')     HKL_file%RAW     = .true.
     IF(HKL_file%NAME(long-2:long) == 'M91')     HKL_file%M91     = .true.
     IF(HKL_file%NAME(long-2:long) == 'M95')     HKL_file%M95     = .true.
     IF(HKL_file%NAME(1:long)      == 'FINAL.Y') HKL_file%final_y = .true.
    endif
    call test_file_exist(trim(HKL_file%NAME), file_exist)
    IF(.NOT. file_exist) then
     HKL_file%name = ''
     return
    endif
    keyword_FILE = .true.


    do i=2, nb_arg
     IF(arg_string(i)(1:4) == 'PLOT')     then
      HKL_file%plot     = .true.
     elseIF(arg_string(i)(1:3) == 'NEG')  then
      HKL_file%read_NEG = .false.
     endif
    end do

   case ('P4P', 'READ_P4P')
    keyword_P4P = .true.
    if(nb_arg/=0) then
     i1 = INDEX(arg_string(1), '.')
     if (i1==0) then
      P4P_file_name = TRIM(arg_string(1))//'.P4P'
     else
      READ(arg_string(1), *, iostat = i_error) P4P_file_name
     end if
    else
     P4P_file_name = ""
    endif
    keyword_read_INS = .true.

   case ('SEARCH_EXTI', 'FIND_EXTI')
    !n_sig = 0.
    ratio_criteria = 0.03
    if(nb_arg/=0) then
     READ(arg_string(1),*, IOSTAT=i_error) var(1)
     IF(i_error ==0) n_sig = var(1)
    end if 
    IF(nb_arg >1) then
     READ(arg_string(2),*, IOSTAT=i_error) var(2)
     IF(i_error ==0) ratio_criteria = var(2)
    endif
    keyword_search_exti = .true.
    
    case ("SEARCH_SPGR", "SEARCH_SPACE_GROUP", "SEARCH_GROUP",   &
          "CHECK_SPGR",  "CHECK_SPACE_GROUP",  "CHECK_GROUP")
     
     if (nb_arg /=0) then
     
      current_i = 0
      do i=1, nb_arg
       IF(arg_string(i)(1:4) == 'TRIC') then
        crystal_system = "TRIC"
        current_i = i
        exit
       elseIF(arg_string(i)(1:4) == 'MONO'  .or. (len_trim(arg_string(i))==1 .and. arg_string(i)(1:1) == 'M')) then
        crystal_system = "MONO"
        current_i = i
        exit
       elseIF(arg_string(i)(1:5) == 'ORTHO' .or. (len_trim(arg_string(i))==1 .and. arg_string(i)(1:1) == 'O')) then
        crystal_system = "ORTHO"
        current_i = i
        exit
       elseIF(arg_string(i)(1:5) == 'TETRA' .or. (len_trim(arg_string(i))==1 .and. arg_string(i)(1:1) == 'T')) then
        crystal_system = "TETRA"
        current_i = i
        exit
       elseIF(arg_string(i)(1:4) == 'TRIG'  .or. arg_string(i)(1:5) == 'RHOMB' &
                                            .or. (len_trim(arg_string(i))==1 .and. arg_string(i)(1:1) == 'R')) then
        crystal_system = "TRIG"
        current_i = i
        exit
       elseIF(arg_string(i)(1:4) == 'HEXA'  .or. (len_trim(arg_string(i))==1 .and. arg_string(i)(1:1) == 'H')) then
        crystal_system = "HEXA"
        current_i = i
        exit
       elseIF(arg_string(i)(1:3) == 'CUB'   .or. (len_trim(arg_string(i))==1 .and. arg_string(i)(1:1) == 'C')) then
        crystal_system = "CUB"
        current_i = i 
        exit
       endif
      end do  
      
      if(current_i /=0) then
       if(nb_arg > 1) then
        do i=1, nb_arg
         if(i==current_i) cycle
         read(arg_string(i), *, iostat=i_error) var(1)
         if(i_error ==0) then
          n_sig = var(1)
          exit
         endif
        end do
       endif
       
       if(nb_arg > 2) then
        i_ok = 0
        do i=1, nb_arg
         if(i==current_i) cycle
         i_ok = i_ok + 1
         read(arg_string(i), *, iostat=i_error) var(i)
         if(i_error /=0) exit
         if(i_ok == 1) then
          n_sig = var(i)
         elseif(i_ok == 2) then
          threshold = var(i)
          exit
         endif
        end do        
       endif
      
      else ! le systeme n'est pas precise
       if(nb_arg == 1) then
        read(arg_string(1), *, iostat=i_error) var(1)
        if (i_error==0) n_sig = var(1)
       else
        read(arg_string(1), *, iostat=i_error) var(1)
        if (i_error==0) n_sig = var(1)
        read(arg_string(2), *, iostat=i_error) var(2)
        if (i_error==0) threshold = var(2)
       endif
      endif
      
     end if  ! fin de la condition if(n_arg /=0
     keyword_search_spgr = .true.
   

   case ('FIND_HKL', 'SEARCH_HKL')
    IF(nb_arg < 3) then
     call check_arg('3', 'FIND_HKL')
     return
    endif
    !READ(arg_string(1), * ) requested_H(1)
    !READ(arg_string(2), * ) requested_H(2)
    !READ(arg_string(3), * ) requested_H(3)
 
    requested_H_string(1:3) = arg_string(1:3)
    if(requested_H_string(1)(1:1) == 'H'  .OR.  &
       requested_H_string(2)(1:1) == 'K'  .OR.  &
       requested_H_string(3)(1:1) == 'L')  then
     search_H_string = .true.            
    else
     READ(arg_string(1), * ) requested_H(1)
     READ(arg_string(2), * ) requested_H(2)
     READ(arg_string(3), * ) requested_H(3)
     search_H_string = .false.       
    endif
    
    search_equiv   = .false.
    search_friedel = .false.
    keyword_FIND_HKL = . true.
    IF(nb_arg == 3) return
    do i=1,nb_arg-3
     IF(arg_string(i+3)(1:5) == 'EQUIV'   ) search_equiv   = .true.
     IF(arg_string(i+3)(1:7) == 'FRIEDEL' ) search_friedel = .true.
    end do

   case ('EQUIV', 'EQUIV_HKL', 'SEARCH_EQUIV', 'SEARCH_EQUIV_HKL', 'FIND_EQUIV', 'FIND_EQUIV_HKL')
    IF(nb_arg < 3) then
     call check_arg('3', 'EQUIV')
     return
    endif
    READ(arg_string(1), * ) requested_H(1)
    READ(arg_string(2), * ) requested_H(2)
    READ(arg_string(3), * ) requested_H(3)
    search_friedel = .false.
    keyword_FIND_HKL_EQUIV = . true.
    IF(nb_arg == 3) return
    do i=1,nb_arg-3
     IF(arg_string(i+3)(1:7) == 'FRIEDEL' ) search_friedel = .true.
    end do


   case ('ABSENT_HKL', 'HKL_ABSENT')
    HKL_list_ABSENT%OUT   = .false.
    HKL_list_ABSENT%ALL   = .false.
    HKL_list_ABSENT%WRITE = .false.
    !n_sig = 0.
    !IF(nb_arg ==0) then
    ! call check_arg('0', 'HKL_ABSENT')
    ! return
    !endif
    IF(nb_arg /=0) then
     do i=1, nb_arg
      IF(arg_string(i)(1:3) == 'ALL')  then
       HKL_list_ABSENT%ALL   = .true.
      elseIF(arg_string(i)(1:3) == 'OUT') then
       HKL_list_ABSENT%OUT   = .true.
      elseIF(arg_string(i)(1:5) == 'WRITE')   then
       HKL_list_ABSENT%WRITE = .true.
      else
       READ(arg_string(i),*) n_sig
       HKL_list_ABSENT%ALL = .false.
      endif
     end do
    END if
    CLOSE(UNIT=31)
    IF(HKL_list_ABSENT%WRITE) then
     i1 = INDEX(HKL_file%NAME,'.')
     WRITE(HKL_file%OUTPUT, '(a,a)') HKL_file%NAME(1:i1-1),'_ABSENT.hkl'
     OPEN(UNIT=31, FILE=TRIM(HKL_file%OUTPUT))
    endif
    keyword_FIND_HKL_ABSENT    = .true.


   CASE ('HKL_NEG', 'HKL_NEGATIVE', 'NEG_HKL', 'NEGATIVE_HKL')
    HKL_list_NEG%OUT   = .false.
    HKL_list_NEG%ALL   = .false.
    HKL_list_NEG%WRITE = .false.
    !n_sig = 0.
    !IF(nb_arg ==0) then
    ! call check_arg('0', 'HKL_NEG')
    ! return
    !endif
    IF(nb_arg /=0) then
     do i=1, nb_arg
      IF(arg_string(i)(1:3) == 'ALL')  then
       HKL_list_NEG%ALL   = .true.
      elseIF(arg_string(i)(1:3) == 'OUT') then
       HKL_list_NEG%OUT   = .true.
      elseIF(arg_string(i)(1:5) == 'WRITE')   then
       HKL_list_NEG%WRITE = .true.
      !else
      ! READ(arg_string(i),*) n_sig
      endif
     end do
    endif
    CLOSE(UNIT=31)
    IF(HKL_list_NEG%WRITE) then
     i1 = INDEX(HKL_file%NAME,'.')
     WRITE(HKL_file%OUTPUT, '(a,a)') HKL_file%NAME(1:i1-1),'_NEG.hkl'
     OPEN(UNIT=31, FILE=TRIM(HKL_file%OUTPUT))
    endif
    keyword_FIND_HKL_NEG = .true.


   CASE ('HKL_POS', 'HKL_POSITIVE', 'POS_HKL', 'POSITIVE_HKL')
    HKL_list_POS%OUT   = .false.
    HKL_list_POS%ALL   = .false.
    HKL_list_POS%WRITE = .false.
    !n_sig = 0.
    !IF(nb_arg ==0) then
    ! call check_arg('0', 'HKL_POS')
    ! return
    !endif
    if(nb_arg /=0) then
     do i=1, nb_arg
      IF(arg_string(i)(1:3) == 'ALL')  then
       HKL_list_POS%ALL   = .true.
      elseIF(arg_string(i)(1:3) == 'OUT') then
       HKL_list_POS%OUT   = .true.
      elseIF(arg_string(i)(1:5) == 'WRITE')   then
       HKL_list_POS%WRITE = .true.
      else
       READ(arg_string(i),*) n_sig
      endif
     end do
    END if
    CLOSE(UNIT=31)
    IF(HKL_list_POS%WRITE) then
     i1 = INDEX(HKL_file%NAME,'.')
     WRITE(HKL_file%OUTPUT, '(a,a)') HKL_file%NAME(1:i1-1),'_POS.hkl'
     OPEN(UNIT=31, FILE=TRIM(HKL_file%OUTPUT))
    endif
    keyword_FIND_HKL_POS = .true.


   case ('STRUCT', 'STRUCTURE_FACTOR', 'FSTR')
    keyword_STRUCT_factor = .true.
    if (nb_arg ==0) then
     call check_arg('0', 'STRUCTURE_FACTOR')
     return
    end if


   case ('RINT', 'R_INT')
    keyword_RINT = .true.
    if (nb_arg /=0) then
     call get_SPG(trim(arg_string(1)))
    end if

   case ("MERGE")
    keyword_MERGE_hkl = .true.
    if (nb_arg /=0) then
     IF(arg_string(1)(1:4) == 'TRIC') then
      space_group_symbol = 'P 1'
      call set_spacegroup(space_group_symbol, SPG)
     elseIF(arg_string(1)(1:4) == 'MONO')  then
      space_group_symbol = 'P 2'
      call set_spacegroup(space_group_symbol, SPG)
     elseIF(arg_string(1)(1:5) == 'ORTHO') then
      space_group_symbol = 'P 2 2 2'
      call set_spacegroup(space_group_symbol, SPG)
     elseIF(arg_string(1)(1:5) == 'TETRA') then
      space_group_symbol = 'P 4'
      call set_spacegroup(space_group_symbol, SPG)
     elseIF(arg_string(1)(1:4) == 'TRIG')  then
      space_group_symbol = 'P 3'
      call set_spacegroup(space_group_symbol, SPG)
     elseIF(arg_string(1)(1:4) == 'HEXA')  then
      space_group_symbol = 'P 6'
      call set_spacegroup(space_group_symbol, SPG)
     elseIF(arg_string(1)(1:3) == 'CUB')   then
      space_group_symbol = 'P 2 3'
      call set_spacegroup(space_group_symbol, SPG)
     endif
    end if



   !case ('LIST_HKL', 'LST_HKL', 'LIST_HKL_EXTI', 'LST_HKL_EXTI')
   !CASE ("WRITE_HKL", "WRITE_HKL_LIST", "WRITE_HKL_LST")
    CASE ("FIND_HKL_LIST", "FIND_HKL_LST", "EXTRACT_HKL_LIST", "EXTRACT_HKL_LST") 
    if (nb_arg == 0) then
     call check_arg('0', 'FIND_HKL_LIST')
     return
    endif
    
    if(arg_string(1)(1:1)     == 'A') then
     HKL_list%EXTI_number = 13
    elseif(arg_string(1)(1:1) == 'B') then
     HKL_list%EXTI_number = 14
    elseif(arg_string(1)(1:1) == 'C') then
     HKL_list%EXTI_number = 15
    elseif(arg_string(1)(1:1) == 'F') then
     HKL_list%EXTI_number = 16
    elseif(arg_string(1)(1:1) == 'I') then
     HKL_list%EXTI_number = 17
    else  
     READ(arg_string(1),*) var(1)
     IF(abs(var(1)) < 1. .or. abs(var(1)) > HKL_rule_nb) then
      call write_info('... Wrong selection rule number ...')
      return
     endif
     HKL_list%EXTI_number = INT(var(1))
    endif
    
    HKL_list%OUT      = .false.
    HKL_list%ALL      = .false.
    HKL_list%WRITE    = .false.
    HKL_list%SUPPRESS = .false.
    IF(nb_arg > 1) then
     do i=2, nb_arg
      IF(arg_string(i)(1:3) == 'OUT') then
       HKL_list%OUT = .TRUE.
      elseIF(arg_string(i)(1:3) == 'ALL') then
       HKL_list%ALL = .TRUE.
      ELSEIF(arg_string(i)(1:5) == 'WRITE') then
       HKL_list%WRITE = .true.
      ELSEIF(arg_string(i)(1:8) == 'SUPPRESS' .or. arg_string(i)(1:6) == 'REMOVE') then
       HKL_list%SUPPRESS = .true.
      else
       READ(arg_string(i),*, iostat=i_error) n_sig
      endif
     end do
    endif
    keyword_FIND_HKL_list = .true.
    IF(HKL_list%WRITE) then
     CLOSE(UNIT=HKL_list_out1_unit)
     i1 = INDEX(HKL_file%NAME,'.')
     IF(i1/=0) then
      if(HKL_list%EXTI_number > 0) then
       IF(HKL_list%EXTI_number < 10) then
        WRITE(HKL_file%OUTPUT, '(a,a,i1,a)') HKL_file%NAME(1:i1-1),'_exti_',HKL_list%EXTI_number,'.hkl'
       else
        WRITE(HKL_file%OUTPUT, '(a,a,i2,a)') HKL_file%NAME(1:i1-1),'_exti_',HKL_list%EXTI_number,'.hkl'
       endif
      else
       IF(abs(HKL_list%EXTI_number) < 10) then
        WRITE(HKL_file%OUTPUT, '(a,a,i1,a)') HKL_file%NAME(1:i1-1),'_exti_m',ABS(HKL_list%EXTI_number),'.hkl'
       else
        WRITE(HKL_file%OUTPUT, '(a,a,i2,a)') HKL_file%NAME(1:i1-1),'_exti_m',ABS(HKL_list%EXTI_number),'.hkl'
       endif       
      endif
      
      OPEN(UNIT=HKL_list_out1_unit, FILE=TRIM(HKL_file%OUTPUT))
     endif
    endif
    
    IF(HKL_list%SUPPRESS) then
     CLOSE(UNIT=HKL_list_out2_unit)
     i1 = INDEX(HKL_file%NAME,'.')
     IF(i1/=0) then
      if(HKL_list%EXTI_number > 0) then
       IF(HKL_list%EXTI_number < 10) then
        WRITE(HKL_file%OUTPUT2, '(a,a,i1,a)') HKL_file%NAME(1:i1-1),'_suppress_exti_',HKL_list%EXTI_number,'.hkl'
       else
        WRITE(HKL_file%OUTPUT2, '(a,a,i2,a)') HKL_file%NAME(1:i1-1),'_suppress_exti_',HKL_list%EXTI_number,'.hkl'
       endif
      else      
       IF(abs(HKL_list%EXTI_number) < 10) then
        WRITE(HKL_file%OUTPUT2, '(a,a,i1,a)') HKL_file%NAME(1:i1-1),'_suppress_exti_',HKL_list%EXTI_number,'.hkl'
       else
        WRITE(HKL_file%OUTPUT2, '(a,a,i2,a)') HKL_file%NAME(1:i1-1),'_suppress_exti_',HKL_list%EXTI_number,'.hkl'
       endif
      endif
      OPEN(UNIT=HKL_list_out2_unit, FILE=TRIM(HKL_file%OUTPUT2))
     endif
    endif
    

   case ('LIST_EXTI', 'LIST_EXTI_RULE', 'LST_EXTI', 'LST_EXTI_RULE')
    keyword_LIST_EXTI_RULE = .true.

   case ('MATMUL')
    keyword_MATMUL = .true.

   CASE ('MAT', 'MATR', 'MATRIX')
    ! matrice de transformation
    keyword_MATR = .false.
    matrix_num = 0
    IF(nb_arg == 1) then
     IF(arg_string(1)(1:1) == 'I' .OR. arg_string(1)(1:2) == '+I') then
      matrix_num = 1
     ELSEIF(arg_string(1)(1:2) == '-I') then
      matrix_num = 2
     ELSEIF(arg_string(1)(1:1) == '#') then
      READ(arg_string(1)(2:),*) matrix_num
     else
      call check_arg('9', 'MATRIX')
      return
     endif

    ELSEIF(nb_arg ==3)  then   ! entree de la matrice sous la forme d'un repere abc
     do i=1,3
      call get_mat_from_setting(arg_STRING(i), i, 'A',  Mat)
      call get_mat_from_setting(arg_STRING(i), i, 'B',  Mat)
      call get_mat_from_setting(arg_STRING(i), i, 'C',  Mat)
     end do

    ELSEIF(nb_arg == 2) then
     IF    (arg_string(1)(1:1) == 'P' .AND. arg_string(2)(1:1) == 'R') then
      matrix_num = 8
     elseif(arg_string(1)(1:1) == 'R' .AND. arg_string(2)(1:1) == 'P') then
      matrix_num = 9     
     elseif(arg_string(1)(1:5) == 'R_REV' .AND. arg_string(2)(1:5) == 'R_OBV') then
      matrix_num = 10
     elseif(arg_string(1)(1:5) == 'R_OBV' .AND. arg_string(2)(1:5) == 'R_REV') then
      matrix_num = 11
     elseif(arg_string(1)(1:1) == 'P' .AND. arg_string(2)(1:1) == 'A') then
      matrix_num = 12
     elseif(arg_string(1)(1:1) == 'A' .AND. arg_string(2)(1:1) == 'P') then
      matrix_num = 13
     elseif(arg_string(1)(1:1) == 'P' .AND. arg_string(2)(1:1) == 'B') then
      matrix_num = 14
     elseif(arg_string(1)(1:1) == 'B' .AND. arg_string(2)(1:1) == 'P') then
      matrix_num = 15
     elseif(arg_string(1)(1:1) == 'P' .AND. arg_string(2)(1:1) == 'C') then
      matrix_num = 16
     elseif(arg_string(1)(1:1) == 'C' .AND. arg_string(2)(1:1) == 'P') then
      matrix_num = 17
     elseif(arg_string(1)(1:1) == 'P' .AND. arg_string(2)(1:1) == 'F') then
      matrix_num = 18
     elseif(arg_string(1)(1:1) == 'F' .AND. arg_string(2)(1:1) == 'P') then
      matrix_num = 19
     elseif(arg_string(1)(1:1) == 'P' .AND. arg_string(2)(1:1) == 'I') then
      matrix_num = 20
     elseif(arg_string(1)(1:1) == 'I' .AND. arg_string(2)(1:1) == 'P') then
      matrix_num = 21
     else
      matrix_num = 0
      WRITE(message_text, '(a)') '  > Unknown matrix !! '
      call write_info(TRIM(message_text))
      return
     endif

    elseif(nb_arg /=9) then
     call check_arg('9', 'MATRIX')
     return
    endif

    IF(matrix_num /=0 .and. (matrix_num <1 .OR. matrix_num > max_mat_nb)) then
     call write_info('')
     WRITE(message_text, '(a,i2,a)') ' ... Numor of matrix has to be in the 1-', max_mat_nb, ' range...'
     call write_info(trim(message_text))
     call write_info('')
     return
    endif
    IF(matrix_num /=0 .and. matrix_num <=max_mat_nb) then
     call def_transformation_matrix()
      Mat(1:3,1:3) = transf_mat(1:3,1:3,matrix_num)
    ELSEIF(nb_arg /=3) then
     i1 = index(arg_line, '/')
     if(i1 /=0) then ! elements fractionnaires
      call Get_matrix_coord(arg_line)
     else
      READ(arg_string(1:3), *) Mat(1,1:3)
      READ(arg_string(4:6), *) Mat(2,1:3)
      READ(arg_string(7:9), *) Mat(3,1:3)
     endif 
    endif
    keyword_MATR = .true.
    call check_matrice(Mat)    
    WRITE(message_text, '(a,F6.3)') '  > Matrix determinant: ', Mat_det
    call write_info(TRIM(message_text))
    call write_info('')

   CASE ('LST_MAT', 'LST_MATR', 'LST_MATRIX',  'LIST_MAT',  'LIST_MATR',  'LIST_MATRIX', 'LIST_TRANSFORMATION_MATRIX')
    keyword_LST_MAT  = .true.



   CASE ('TRANSLATION', 'TRANS', 'TRANSLATE', 'MOVE')
    ! translation
    if (nb_arg/=3) then
     call check_arg('3', 'TRANSLATION')
     return
    end if
    READ(arg_string(1:3), *) translat(1:3)
    keyword_TRANSL = .true.


   case ('SET')
    IF(nb_arg/=2) then
     call check_arg('2', 'SET')
     return
    endif
    IF(arg_string(1)(1:7) == 'BROWSER') then
     my_browser%name = arg_string(2)
    elseIF(arg_string(1)(1:6) == 'EDITOR')  then
     my_editor%name  = arg_string(2)
    else
     call check_arg('unknown', 'SET')
     return
    endif
    IF(LEN_TRIM(my_browser%name) == 0 .AND. LEN_TRIM(my_editor%name) == 0) return
    keyword_set    = .true.

   case ('REPORT', 'CREATE_REPORT')
    long_report = .false.
    if (nb_arg ==0) then
     archive_CIF = 'archive.cif'
    ELSEIF(nb_arg==1) then
     IF(arg_string(1)(1:4)=='LONG' .OR. arg_string(1)(1:3)=='EXT') then
      long_report = .TRUE.
      archive_CIF = 'archive.cif'
     endif
    ELSEIF(nb_arg==2) then
     IF(arg_string(1)(1:4)=='LONG' .OR. arg_string(1)(1:3)=='EXT') then
      long_report = .true.
      archive_CIF = arg_string(2)
     ELSEIF(arg_string(2)(1:4)=='LONG' .OR. arg_string(2)(1:3)=='EXT') then
      long_report = .true.
      archive_CIF = arg_string(1)
     endif
    end if
    keyword_create_REPORT = .true.

   CASE ('SETTING')
    keyword_setting = .true.

   CASE ('SHELL')
    shell_arg_min_max = .false.
    IF(nb_arg ==0) then
     call check_arg('0', 'SHELL')
     return
    endif
    IF(mode_interactif) then
     nb_shell = 1
    else
     nb_shell = nb_shell + 1
    endif

    IF(arg_string(1)(1:1) == 'D')    then
     shell_type(nb_shell) = 'd'
    elseif(arg_string(1)(1:3) == 'STL')   then
     shell_type(nb_shell) = 'stl'
    elseif(arg_string(1)(1:5) == 'THETA') then
     shell_type(nb_shell) = 'theta'
    elseif(arg_string(1)(1:3) == 'INT')   then
     shell_type(nb_shell) = 'int'
    elseif(arg_string(1)(1:4) == 'ISIG')  then
     shell_type(nb_shell) = 'isig'
    else
     call check_arg('unknown', 'SHELL')
     return
    endif

    IF(nb_arg >=3) then
     READ(arg_string(2), *) shell_min(nb_shell)
     READ(arg_string(3), *) shell_max(nb_shell)
     shell_arg_min_max = .true.
    endif

    IF(nb_arg > 1 .AND. arg_string(2)(1:4) == 'PLOT') shell_plot(nb_shell) = .true.


   CASE ('SORT')
    IF(nb_arg ==0) then
     call check_arg('0', 'SORT')
     return
    endif
    IF(mode_interactif) then
     nb_sort = 1
    else
     nb_sort = nb_sort + 1
    endif

    !IF(nb_arg > 1 .AND. arg_string(2)(1:4) == 'PLOT') sort_plot(nb_sort) = .true.
    sort_plot = .false.
    sort_out  = .false.
    IF(nb_arg > 1) then
     do i=2, nb_arg
      IF(arg_string(i)(1:4) == 'PLOT') then
       sort_plot(nb_sort) = .TRUE.
      elseIF(arg_string(i)(1:3) == 'OUT') then
       sort_out(nb_sort) = .true.
       IF(arg_string(i)(4:4) == '_')  READ(arg_string(i)(5:),*) sort_out_n(nb_sort)
      endif
     END do
    endif

    IF(arg_string(1)(1:1) == 'D')     then
     sort_type(nb_sort) = 'd'
    ELSEIF(arg_string(1)(1:3) == 'STL')   then
     sort_type(nb_sort) = 'stl'
    ELSEIF(arg_string(1)(1:5) == 'THETA') then
     sort_type(nb_sort) = 'theta'
    ELSEIF(arg_string(1)(1:3) == 'INT')   then
     sort_type(nb_sort) = 'int'
    ELSEIF(arg_string(1)(1:4) == 'ISIG')  then
     sort_type(nb_sort) = 'isig'
    ELSE
     call check_arg('unknown', 'SORT')
     return
    ENDIF


   CASE ('GEN_HKL', 'GENERATE_HKL', 'GENERATE_HKL_LIST')
    i1 = INDEX(read_line, '=')
    if (i1==0) then
     call check_arg('', 'GEN_HKL')
     return
    endif

    if (arg_string(1)(1:3) == '2TH') then
     HKL_2THETA = .true.
    ELSEIF(arg_string(1)(1:2) == 'TH') then
     HKL_THETA = .true.
    ELSEIF(arg_string(1)(1:3) == 'STL') then
     HKl_STL = .true.
    ELSEIF(arg_string(1)(1:1) == 'D') then
     HKL_D = .true.
    ELSEIF(arg_string(1)(1:1) == 'Q') then
     HKL_Q = .true. 
    end if

    READ(read_line(i1+1:), *) X_min

    i2 = INDEX(read_line, '=', back=.TRUE.)
    IF(i2==0) return
    new_line = read_line(i2+1:)
    call nombre_de_colonnes(new_line, nb_arg)
    if(nb_arg == 1) then
     read(new_line, *) X_max
     arg_line = ''
    else
     read(new_line, *) X_max, arg_line
    endif


    !READ(read_line(i2+1:), *) X_max, arg_line
    if(len_trim(arg_line) /=0 ) then
     write_HKL = .false.
     arg_line = adjustl(arg_line)
     arg_line = u_case(arg_line)
     if(arg_line(1:3) =='OUT') write_HKL = .true.
    else
     WRITE_HKL = .false.
    endif
    keyword_GENHKL = .true.

    call write_info('')
    call write_info('  > GEN_HKL ')
    call write_info('')

   CASE ('GEN_SAT')
    keyword_GENSAT = .true.

   CASE ('INSIDE')
    keyword_INSIDE = .true.

   CASE ('DIST', 'DISTANCE', 'ATOMIC_DISTANCE')
    IF(nb_arg < 2) then
     call check_arg('2', 'DIST')
     return
    endif
    IF(.NOT. keyword_CELL) then
     call write_info('')
     call write_info('  CELL keyword is mandatory for interatomic distance calculations')
     call write_info('')
     return
    endif
    nb_dist_calc = nb_dist_calc + 1
    READ(arg_string(1:2), *) atom1_dist(nb_dist_calc), atom2_dist(nb_dist_calc)
    atom1_dist(nb_dist_calc) = ADJUSTL(atom1_dist(nb_dist_calc))
    atom2_dist(nb_dist_calc) = ADJUSTL(atom2_dist(nb_dist_calc))

   CASE ('CONN', 'CONNECT', 'CONNECTIVITY')
    !if(nb_arg < 1) then
    ! call check_arg('1', 'CONN')
    ! return
    !end if
    CONN_all = .false.
    CIF_dist%n_text = 0
    if(nb_arg == 0) CONN_all = .true.
    if(.not. keyword_CELL) then
     call write_info('')
     call write_info(' CELL keyword is mandatory for connectivity calculations')
     call write_info('')
     return
    endif
    if(.not. keyword_SPGR) then
     call write_info('')
     call write_info(' SPGR keyword is mandatory for connectivity calculations')
     call write_info('')
     return
    end if
    if(nb_atom == 0) then
     call write_info('')
     call write_info(' ATOM keyword is mandatory for connectivity calculations')
     call write_info('')
     return
    end if
    
    if (nb_arg /=0) then
     if(arg_string(1)(1:3) == 'ALL') then
      CONN_all = .true.
     else
      read(arg_string(1), *) atom_CONN%label
      atom_CONN%label = adjustl(atom_CONN%label)
     endif
    endif 
    
    if(nb_arg > 1) then
     read(arg_string(nb_arg), *, iostat=i_error) CONN_dmax
     if(i_error /=0) CONN_dmax = CONN_dmax_ini
     if(CONN_dmax < 0.1 .or. CONN_dmax > 15.) CONN_dmax = CONN_dmax_ini
    endif
    keyword_CONN = .true.

   CASE ('DIST_')
    keyword_DIST_ = .false.
    IF(nb_arg < 3) then
     call check_arg('3', 'DIST_')
     return
    endif
    IF(.NOT. keyword_CELL) then
     call write_info('')
     call write_info('  CELL keyword is mandatory for interatomic distance calculations')
     call write_info('')
     return
    endif
    nb_dist_calc = nb_dist_calc + 1
    READ(arg_string(1)  , *) DIST_coef
    READ(arg_string(2:3), *) atom1_dist(nb_dist_calc), atom2_dist(nb_dist_calc)
    atom1_dist(nb_dist_calc) = ADJUSTL(atom1_dist(nb_dist_calc))
    atom2_dist(nb_dist_calc) = ADJUSTL(atom2_dist(nb_dist_calc))
    keyword_DIST_ = .true.


   CASE ('ANG', 'ANGLE')
    IF(nb_arg < 2) then
     call check_arg('', 'ANG')
     return
    endif
    IF(.NOT. keyword_CELL) then
     call write_info('')
     call write_info('  CELL keyword is mandatory for interatomic angle calculation')
     call write_info('')
     return
    endif

    nb_ang_calc = nb_ang_calc + 1

    if (nb_arg == 3) then
     READ(arg_string(1:3),*) atom1_ang(nb_ang_calc),  atom2_ang(nb_ang_calc), atom3_ang(nb_ang_calc)
     atom1_ang(nb_ang_calc) = ADJUSTL(atom1_ang(nb_ang_calc))
     atom2_ang(nb_ang_calc) = ADJUSTL(atom2_ang(nb_ang_calc))
     atom3_ang(nb_ang_calc) = ADJUSTL(atom3_ang(nb_ang_calc))
     atom4_ang(nb_ang_calc) = 'ZZZ'

    ELSEIF (nb_arg == 4) then
     READ(arg_string(1:4),*) atom1_ang(nb_ang_calc),  atom2_ang(nb_ang_calc), atom3_ang(nb_ang_calc), atom4_ang(nb_ang_calc)
     atom1_ang(nb_ang_calc) = ADJUSTL(atom1_ang(nb_ang_calc))
     atom2_ang(nb_ang_calc) = ADJUSTL(atom2_ang(nb_ang_calc))
     atom3_ang(nb_ang_calc) = ADJUSTL(atom3_ang(nb_ang_calc))
     atom4_ang(nb_ang_calc) = ADJUSTL(atom4_ang(nb_ang_calc))

    else
     call write_info('')
     call write_info('    ... uncorrect list of atoms ...')
     call write_info('')
     nb_ang_calc = nb_ang_calc - 1
    end if


   case ('DIR_ANG', 'DIRANG', 'DIRECT_ANGLE')
    IF(nb_arg <6 ) then
     call check_arg('6', 'DIR_ANG')
    endif
    READ(arg_string(1:6), *, IOSTAT=i_error) var(1:6)
    IF(i_error /=0) then
     call check_arg('uncorrect', 'DIR_ANG')
    endif
    IF(mode_interactif) then
     nb_da = 1
    else
     nb_da = nb_da + 1
    endif
    U1_da(1:3, nb_da) = var(1:3)
    U2_da(1:3, nb_da) = var(4:6)
    keyword_DIR_ANG = .true.



   case ('REC_ANG', 'RECANG', 'RECIPROCAL_ANGLE')
    IF(nb_arg <6 ) then
     call check_arg('6', 'REC_ANG')
    endif
    READ(arg_string(1:6), *, IOSTAT=i_error) var(1:6)
    IF(i_error /=0) then
     call check_arg('uncorrect', 'REC_ANG')
    endif
    IF(mode_interactif) then
     nb_ra = 1
    else
     nb_ra = nb_ra + 1
    endif
    U1_ra(1:3, nb_ra) = var(1:3)
    U2_ra(1:3, nb_ra) = var(4:6)
    keyword_REC_ANG = .true.

   CASE ('BARY', 'CENTROID')
    IF(nb_arg < 3) then
     call check_arg('at least 3', 'BARY')
     return
    endif
    nb_bary_calc = nb_bary_calc + 1
    nb_atom_bary(nb_bary_calc) = nb_arg
    read(arg_string(1:nb_arg), *) atom_bary(nb_bary_calc,1:nb_arg)

   case ('STL', 'STL_HKL', 'SINTHETA/WAVE', 'SINTHETA/LAMBDA')
    IF(nb_arg == 0) then
     call check_arg('0', 'STL')
     return
    endif
    nb_STL_value = nb_arg
    READ(arg_string(1:nb_arg), *) STL_value(1:nb_arg)
    keyword_STL = .true.

   case ('D_HKL', 'DHKL')
    IF(nb_arg == 0) then
     call check_arg('0', 'D_HKL')
     return
    endif
    nb_dhkl_value = nb_arg
    READ(read_line(long_kw+1:), *) (dhkl_value(i), i=1, nb_dhkl_value)
    keyword_DHKL = .true.

   case ('D_STAR', 'D_STAR_HKL', 'DSTAR', 'DSTARHKL', 'DSTAR_HKL')
    IF(nb_arg == 0) then
     call check_arg('0', 'D_STAR')
     return
    endif
    nb_dstar_value = nb_arg
    READ(arg_string(1:nb_arg), *) dstar_value(1:nb_arg)
    keyword_DSTAR = .true.

   case ('Q_HKL', 'QHKL')
    IF(nb_arg == 0) then
     call check_arg('0', 'Q_HKL')
     return
    endif
    nb_Qhkl_value = nb_arg
    READ(arg_string(1:nb_arg), *) Qhkl_value(1:nb_arg)
    keyword_QHKL = .true.

   case ('THETA', 'TH', 'TH_HKL', 'THETAHKL', 'THETA_HKL')
    IF(nb_arg== 0) then
     call check_arg('0', 'THETA_HKL')
     return
    endif
    nb_theta_value = nb_arg
    READ(arg_string(1:nb_arg), *) theta_value(1:nb_arg)
    keyword_theta = .true.

   case ('2THETA', '2TH', '2TH_HKL', '2THETA_HKL', 'TWO_THETA', 'TWO_THETA_HKL')
    IF(nb_arg == 0) then
     call check_arg('0', '2THETA_HKL')
     return
    endif
    nb_2theta_value = nb_arg
    READ(arg_string(1:nb_arg), *) TWO_theta_value(1:nb_arg)
    keyword_2theta = .true.


   case ('RESET', 'RAZ', 'INITIALIZE', 'INIT')
    keyword_RESET = .true.

   CASE ('SYST', 'CMD', 'COMMAND', 'DOS', 'DOS_COMMAND')
    IF(nb_arg == 0) then
     call check_arg('0', 'SYST')
     return
    endif
    READ(arg_line,'(a)') SYST_command
    SYST_command = ADJUSTL(SYST_command)
    IF(LEN_TRIM(SYST_command) /=0) keyword_SYST =.true.

   case ('DIR')
    IF(nb_arg == 0) then
     SYST_command = 'dir'
    else
     WRITE(SYST_command, '(2a)') 'dir ', TRIM(arg_line)
    endif
    keyword_SYST = .true.



   CASE ('THERM', 'THERMAL', 'ADP')
    THERM_Uiso              = .false.
    THERM_Biso              = .false.
    THERM_Uij               = .false.
    THERM_Bij               = .false.
    THERM_Beta              = .false.
    THERM_aniso             = .false.
    keyword_THERM           = .false.

    IF(nb_arg <2) then
     call check_arg('at least 2', 'THERM')
     return
    endif

    IF(arg_string(1)(1:4)     == 'BISO'    .or. arg_string(1)(1:5) == 'B_ISO') then
     THERM_Biso = .true.
    ELSEIF(arg_string(1)(1:4) == 'UISO'    .or. arg_string(1)(1:5) == 'U_ISO') then
     THERM_Uiso = .true.
    ELSEIF(arg_string(1)(1:4) == 'U_IJ'    .OR. arg_string(1)(1:3) == 'UIJ') then
     THERM_Uij   = .true.
     THERM_aniso = .true.
    ELSEIF(arg_string(1)(1:4) == 'B_IJ'    .OR. arg_string(1)(1:3) == 'BIJ') then
     THERM_Bij   = .true.
     THERM_aniso = .true.
    ELSEIF(arg_string(1)(1:7) == 'BETA_IJ' .OR. arg_string(1)(1:6) == 'BETAIJ' .or. arg_string(1)(1:4) == 'BETA') then
     THERM_Beta  = .true.
     THERM_aniso = .true.
    endif
    nb_therm_values = nb_arg - 1
    IF(THERM_aniso) then
     IF(.NOT. keyword_CELL) then
      call write_info('')
      call write_info('  Cell parameters has to be known for anisotropic ADP calculations')
      call write_info('')
      return
     endif
     if (nb_therm_values <6) then
      call check_arg('6', 'THERM')
      return
     endif
    !else
    ! nb_therm_values = 6
    endif
    READ(arg_string(2:nb_arg),*) therm_values(1:nb_therm_values)

    keyword_THERM = .true.

   CASE ('TRANSMISSION')
    IF(nb_arg == 0) then
     call check_arg('0', 'TRANSMISSION')
     return
    endif
    keyword_transmission = .true.
    nb_dim_transmission = nb_arg
    READ(arg_string(1:nb_arg),*)  dim_transmission(1:nb_arg)


   case ('SG_ALL',  'SP_ALL')
    write_SPG_info = .true.
    write_SPG_exti = .true.
    write_SPG_subgroups = .true.
    write_SPG_info_all  = .true.

   CASE ('SG_INFO', 'SP_INFO', 'SPACE_GROUP_INFO', 'LIST_SPACE_GROUP_INFO')
    WRITE_SPG_info     = .true.
    write_SPG_info_all = .true.
    write_spg_exti     = .false.
    IF(nb_arg /=0 .and. arg_string(1)(1:3)=='RED') write_SPG_info_all = .false.

   CASE ('SG_EXTI', 'SP_EXTI', 'SG_EXTINCTIONS', 'SPACE_GROUP_EXTI', 'SPACE_GROUP_EXTINCTIONS')
    WRITE_SPG_exti = .true.
    WRITE_SPG_info = .false.

   CASE ('SG_SUB', 'SG_SUBGROUP')
    write_SPG_subgroups = .true.
    write_SPG_info      = .false.
    write_SPG_info_all  = .false.


   !CASE ('SYM_OP', 'SYMM_OP', 'SYM_OPERATOR', 'SYMM_OPERATOR')
   CASE ('WRITE_SYM_OP', 'WRITE_SYMM_OP', 'WRITE_SYM_OP', 'WRITE_SYMM_OP', 'WRITE_SYMMETRY_OPERATORS')
    WRITE_symm_op  = .true.


   CASE ('APPLY_OP', 'APPLY_SYMMETRY_OPERATOR')
    WRITE_APPLY_symm = .true.


   CASE ('SITE_INFO', 'LIST_SITE_INFO')
    WRITE_site_info = .true.

    IF(nb_arg==0 .or. arg_string(1)(1:3) == 'ALL') then
     site_info_all_atoms = .true.
    else
     site_info_all_atoms = .false.
     nb_atom_site_info = nb_arg
     READ(arg_string(1:nb_arg),*) site_info_label(1:nb_arg)
     site_info_label(1:nb_arg) = ADJUSTL(site_info_label(1:nb_arg))
    endif


   CASE ('SHIFT_2TH', 'SHIFT_2THETA', '2TH_SHIFT', '2THETA_SHIFT')
    IF(nb_arg == 0) then
     call check_arg('0', 'SHIFT_2TH')
     return
    endif

    keyword_SH_2th = .true.
    READ(arg_string(1),*, iostat=i_error) var(1)
    IF(i_error /=0) then
     call error_message('SHIFT_2TH')
     return
    endif

    shift_2theta = var(1)
    call write_info('')
    WRITE(message_text,'( a,F10.5)') '  > Shift 2theta (deg):', shift_2theta
    call write_info(TRIM(message_text))
    keyword_SH_2th = .true.


   CASE ('HELP', 'MAN')
    keyword_HELP = .true.
    nb_help = nb_help_max
    !IF(nb_col /=0) then
    ! nb_help = nb_col
    ! read(read_line(long_kw+1:), *) HELP_string(1:nb_help)
    !endif
    IF(nb_arg >1) then
     nb_help = nb_arg
     READ(arg_string(1:nb_arg),*) HELP_arg(1:nb_arg)
    ELSEIF(nb_arg == 1) then
     i1 = INDEX(arg_string(1), '*')
     IF(i1==0) then
      nb_help = 1
      READ(arg_string(1), *) HELP_arg(1)
      return
     else
      IF(i1==1) then
       READ(arg_string(1)(i1+1:),*) arg1
      else
       READ(arg_string(1)(1:i1-1),*) arg1
      endif
      nb_help = 0
      temp_string = help_string
      do i=1, nb_help_max
       i1=INDEX(temp_string(i),TRIM(arg1))
       if (i1/=0) then
        nb_help= nb_help+1
        help_arg(nb_help) = temp_string(i)
       end if
      end do
     endif
    END if

   case ('MAN_HTML', 'HTML_MAN', 'HTML')
    keyword_create_CRYSCAL_HTML = .true.
    browse_cryscal_HTML         = .false.
    IF(arg_string(1)(1:6) == 'BROWSE') browse_cryscal_HTML = .true.

   case ('NEWS')
    keyword_NEWS = .true.
    if(nb_arg /=0) then
     read(arg_string(1), *) news_year
    else
     news_year = 'all'
    end if

   case ('HEADER', 'HEAD')
    keyword_HEADER = .true.

   CASE ('KEY', 'KEYS', 'LST_KEYS', 'LIST_KEYS', 'LST_KEYWORDS', 'LIST_KEYWORDS')
    keyword_KEY = .true.
    IF(nb_arg == 0) then
     write_keys(1: nb_help_max) = .true.
    ELSE
     READ(read_line(long_kw+1:), *) arg_line
     i1=INDEX(arg_string(1), '*')
     if (i1==0) then
      write_keys(1:nb_help_max) = .false.
      do i=1, nb_help_max
       do j=1, nb_arg
        if (arg_string(j) == help_string(i)) write_keys(i) = .true.
       END do
      end do
      return
     elseif(i1 > 1) then
      READ(arg_string(1)(1:i1-1),*) arg1
      temp_string = help_string
      do i=1, nb_help_max
       i1 = INDEX(temp_string(i), TRIM(arg1))
       if (i1/=0) then
        write_keys(i) = .true.
       else
        write_keys(i) = .false.
       end if
      end do
     end if
    endif


   CASE ('LIST_SG', 'LST_SG', 'LIST_SPACE_GROUPS')
    list_sg(1:7)         = .false.
    list_sg_Bravais(1:7) = .true.
    list_sg_centric(1:2) = .false.
    list_sg_laue(1:14)   = .true.
    list_sg_multip       = .false.

    IF(nb_arg == 0) then
     list_sg(1:7)  = .true.
     keyword_LSPGR = .true.
     return
    endif


    do i=1, nb_arg
     arg = l_case(arg_string(i))

     select case (arg)
       case ('all')
         list_sg(1:7) = .true.

       case ('tri', 'tric', 'tricl', 'tricli', 'triclini', 'triclinic')
         list_sg(1) = .true.

       case ('mono', 'monoc', 'monocl', 'monocli', 'monoclini', 'monoclinic')
         list_sg(2) = .true.

       case ('ortho', 'orthor', 'orthorh', 'orthorho', 'orthorhom', 'orthorhomb', 'orthorhombic', 'orthorhombic')
         list_sg(3) = .true.

       case ('tetra',  'tetrag',  'tetrago',  'tetragon',  'tetragona',  'tetragonal', &
             'quadra', 'quadrat', 'quadrati', 'quadratiq', 'quadratiqu', 'quadratique')
         list_sg(4) =  .true.

       case ('trig', 'trigo', 'trigon', 'trigona', 'trigonal')
         list_sg(5) = .true.

       case ('hex', 'hexa', 'hexag', 'hexago', 'hexagon', 'hexagona', 'hexagonal')
         list_sg(6) = .true.

       case ('cub', 'cubi', 'cubic')
         list_sg(7) = .true.

       case ('centric', 'centro')
         list_sg_centric(1:2) = .false.
         list_sg_centric(1) = .true.
         IF(nb_arg == 1) list_sg(1:7) = .true.

       case ('acentric', 'non-centro')
        list_sg_centric(1:2) = .false.
        list_sg_centric(2) = .true.
        IF(nb_arg == 1) list_sg(1:7) = .true.

       case ('mult', 'multip', 'multipl', 'multiplicity')
        list_sg_multip = .true.

       case ('p')
        list_sg_bravais(1:7) = .false.
        list_sg_bravais(1)   = .true.
       case ('a')
        list_sg_bravais(1:7) = .false.
        list_sg_bravais(2)   = .true.
       case ('b')
        list_sg_bravais(1:7) = .false.
        list_sg_bravais(3)   = .true.
       case ('c')
        list_sg_bravais(1:7) = .false.
        list_sg_bravais(4)   = .true.
       case ('i')
        list_sg_bravais(1:7) = .false.
        list_sg_bravais(5)   = .true.
       case ('f')
        list_sg_bravais(1:7) = .false.
        list_sg_bravais(6)   = .true.
       case ('r')
        list_sg_bravais(1:7) = .false.
        list_sg_bravais(7)   = .true.

       case ('laue_1')    ! -1
        list_sg_laue(1:14) = .false.
        list_sg_laue(1)    = .true.
        list_sg(1)         = .true.
       case ('laue_2')    ! 2/m
        list_sg_laue(1:14) = .false.
        list_sg_laue(2)    = .true.
        list_sg(2)         = .true.
       case ('laue_3')    ! mmm
        list_sg_laue(1:14) = .false.
        list_sg_laue(3)    = .true.
        list_sg(3)         = .true.
       case ('laue_4')    ! 4/m
        list_sg_laue(1:14) = .false.
        list_sg_laue(4)    = .true.
        list_sg(4)         = .true.
       case ('laue_5')   ! 4/mmm
        list_sg_laue(1:14) = .false.
        list_sg_laue(5)    = .true.
        list_sg(4)         = .true.
       case ('laue_6')   ! -3 (rhomb. axes)
        list_sg_laue(1:14) = .false.
        list_sg_laue(6)    = .true.
        list_sg(5)         = .true.
       case ('laue_7')   ! -3 (hex. axes)
        list_sg_laue(1:14) = .false.
        list_sg_laue(7)    = .true.
        list_sg(5)         = .true.
       case ('laue_8')  ! -3m (rhomb. axes)
        list_sg_laue(1:14) = .false.
        list_sg_laue(8)    = .true.
        list_sg(5)         = .true.
       case ('laue_9')  ! -31m (hex. axes)
        list_sg_laue(1:14) = .false.
        list_sg_laue(9)    = .true.
        list_sg(5)         = .true.
       case ('laue_10') ! -3m1 (hex. axes)
        list_sg_laue(1:14) = .false.
        list_sg_laue(10)   = .true.
        list_sg(5)         = .true.
       case ('laue_11') ! 6/m
        list_sg_laue(1:14) = .false.
        list_sg_laue(11)    = .true.
        list_sg(6)         = .true.
       case ('laue_12') ! 6/mmm
        list_sg_laue(1:14) = .false.
        list_sg_laue(12)    = .true.
        list_sg(6)         = .true.
       case ('laue_13') ! m3
        list_sg_laue(1:14) = .false.
        list_sg_laue(13)    = .true.
        list_sg(7)         = .true.
       case ('laue_14') ! m3m
        list_sg_laue(1:14) = .false.
        list_sg_laue(14)    = .true.
        list_sg(7)         = .true.

     end select

    end do
    keyword_LSPGR = .true.

   case ('LST_LAUE', 'LIST_LAUE', 'LST_LAUE_CLASS', 'LIST_LAUE_CLASS')
    keyword_laue = .true.

   case ('PAUSE')
    keyword_PAUSE = .true.

   case ('PERMUT', 'PERMUTATION', 'PERMUT_ABC', 'PERMUTATION_ABC')
    write_permutation_abc= .true.

   CASE ('TRICLINIC', 'TRICL')
    write_triclinic_transf = .true.
    
   CASE ('MONOCLINIC', 'MONOC', 'MONOCL')
    WRITE_monoclinic_transf = .true.
    
   CASE ('HEXA_TWIN', 'HEXA_TWINNING',  'HEXAGONAL_TWIN', 'HEXAGONAL_TWINNING', &
         'TWIN_HEXA', 'TWIN_HEXAGONAL', 'TWINNING_HEXA',  'TWINNING_HEXAGONAL')
    WRITE_twin_hexa = .true. 
    
   CASE ('TWIN_PSEUDO_HEXA')
    WRITE_twin_pseudo_hexa = .true.


   CASE ('OBV_REV', 'OBVERSE_REVERSE', 'TWIN_OBVERSE_REVERSE', 'TWIN_OBV_REV', &
            'TWINNING_OBVERSE_REVERSE', 'TWINNING_OBV_REV')
    keyword_OBV_REV = .true.
    if (nb_arg/=0) then
     READ(arg_string(1), *) OBV_REV_twin_matrix
     IF(OBV_REV_twin_matrix /=2)  OBV_REV_twin_matrix = 1
    end if

   CASE ('RHOMB_HEX', 'RHOMB_HEXA', 'RHOMB_TO_HEX', 'RHOMB_TO_HEXA')
    WRITE_rhomb_hex_transf = .true.

   CASE ('HEX_RHOMB', 'HEXA_RHOMB', 'HEX_TO_RHOMB', 'HEXA_TO_RHOMB')
    WRITE_hex_rhomb_transf = .true.

   case ('SHAN', 'SHANNON')
    shannon_atom_label = ''
    !IF(nb_arg == 0) then
    ! call check_arg('0', 'SHANNON')
    ! return
    !endif
    IF(nb_arg /=0 ) READ(arg_string(1),*) shannon_atom_label
    keyword_SHANNON = .true.

   case ('MAG', 'MAGNETISM', 'MAGNETIC')
    mag_atom_label = ''
    IF(nb_arg /=0) READ(arg_string(1), *) mag_atom_label
    keyword_MAG = .true.

   case ('MENDEL')
    mendel_atom_nb   = 0
    mendel_plot      = .false.
    IF(nb_arg == 0 .or. arg_string(1) == 'PLOT')  then
     call check_arg('0', 'MENDEL')
     return
    endif

    if (arg_string(nb_arg) == 'PLOT') then
     mendel_plot = .true.
     mendel_atom_nb = nb_arg - 1
    else
     mendel_atom_nb = nb_arg
    endif
    READ(arg_string(1:nb_arg),*) mendel_atom_label(1:nb_arg)
    keyword_MENDEL = .true.

   case ('DATA_NEUTRONS', 'DATA_NEUTRON', 'NEUTRON_DATA', 'NEUTRONS_DATA')
    data_neutrons_PLOT = .false.
    if(nb_arg/=0 .and. arg_string(1) == 'PLOT') data_neutrons_PLOT = .true.
    keyword_DATA_NEUTRONS = .true.

   case ('DATA_XRAYS', 'DATA_XRAY', 'XRAYS_DATA', 'XRAY_DATA')
    data_xrays_PLOT = .false.
    if(nb_arg/=0 .and. arg_string(1) == 'PLOT') data_Xrays_PLOT = .true.
    keyword_DATA_XRAYS = .true.

   case ('DATA_DENSITY', 'DENSITY_DATA', 'DATA_ATOMIC_DENSITY', 'ATOMIC_DENSITY')
    data_atomic_density_PLOT = .false.
    if(nb_arg/=0 .and. arg_string(1) == 'PLOT') data_atomic_density_PLOT = .true.
    keyword_DATA_ATOMIC_DENSITY = .true.

   case ('DATA_RADIUS', 'RADIUS_DATA', 'DATA_ATOMIC_RADIUS', 'ATOMIC_RADIUS')
    data_atomic_radius_PLOT = .false.
    if(nb_arg/=0 .and. arg_string(1) == 'PLOT') data_atomic_radius_PLOT = .true.
    keyword_DATA_ATOMIC_RADIUS  = .true.

   case ('DATA_WEIGHT', 'WEIGHT_DATA', 'DATA_ATOMIC_WEIGHT', 'ATOMIC_WEIGHT')
    data_atomic_weight_PLOT = .false.
    if(nb_arg/=0 .and. arg_string(1) == 'PLOT') data_atomic_weight_PLOT = .true.
    keyword_DATA_ATOMIC_WEIGHT = .true.


   case ('WEB', 'INTERNET')
    IF(nb_arg == 0)  then
     call check_arg('0', 'WEB')
     return
    endif
    URL_address = ''
    if (arg_string(1)(1:5) == 'CDIFX') then
     URL_address = 'www.cdifx.univ-rennes1.fr'
    else
     long1 = LEN_TRIM(arg_string(1))
     do i = 1, WEB%num_site
      long2 = LEN_TRIM(WEB%NAME(i))
      if (arg_string(1)(1:long1) == WEB%NAME(i)(1:long2)) then
       URL_address = WEB%address(i)
       exit
      end if
     end do
     IF(LEN_TRIM(URL_address) == 0) URL_address = arg_string(1)
    end if
    keyword_WEB = .true.

   case ('WRITE_CELL', 'OUTPUT_CELL')
    IF(.NOT. keyword_CELL) then
     call write_info('')
     call write_info('  Cell parameters has to be known for WRITE_CELL keyword')
     call write_info('')
     return
    endif
    keyword_WRITE_CELL = .true.

 
   case ('WRITE_CHEM', 'WRITE_CHEMICAL_FORMULA', 'OUTPUT_CHEM')
    IF(len_trim(molecule%formula) == 0) then
     call write_info('')
     call write_info('  Chemical formula not known !')
     call write_info('')
     return
    endif
    keyword_WRITE_CHEM = .true.

   case ('WRITE_QVEC', 'OUTPUT_QVEC')
    if(.not. keyword_QVEC) then
     call write_info('')
     call write_info('  QVEC has to be known for WRITE_QVEC keyword')
     call write_info('')
     return
    endif
    keyword_WRITE_QVEC = .true.
   
   case ('WRITE_SG', 'WRITE_SPACE_GROUP')
    IF(.NOT. keyword_SPGR) then
     call write_info('')
     call write_info('  SPGR keyword is mandatory for WRITE_SG keyword')
     call write_info('')
     return
    endif
    list_sg(1:7)  = .true.
    keyword_WRITE_SG = .true.


   case ('WRITE_WAVE', 'OUTPUT_WAVE')
    IF(.NOT. keyword_wave) then
     call write_info('')
     call write_info('  Wave has to be known for WRITE_WAVE keyword')
     call write_info('')
     return
    endif
    keyword_WRITE_WAVE = .true.
    
   case ('WRITE_DEVICE', 'OUTPUT_DEVICE')
    keyword_WRITE_DEVICE =.true. 

   case ('WRITE_BEAM', 'WRITE_INCIDENT_BEAM')  
    keyword_WRITE_BEAM = .true.


   case ('NIGGLI', 'NIGGLI_CELL')
    IF(.NOT. keyword_CELL) then
     call write_info('')
     call write_info(' CELL paramters has to be known')
     call write_info('')
     return
    endif
    keyword_NIGGLI = .true.

   case ('REF_KCCD', 'KCCD')
    keyword_WRITE_REF_KCCD = .true.

   case ('REF_APEX', 'REF_APEXII' , 'WRITE_APEX', 'WRITE_APEXII',  'APEX', 'APEXII')
    keyword_WRITE_REF_APEX = .true.

   case ('REF_EVAL', 'REF_EVALCCD', 'WRITE_EVAL', 'WRITE_EVALCCD', 'EVAL', 'EVALCCD')
    keyword_WRITE_REF_EVAL = .true.

   case ('REF_DENZO', 'WRITE_DENZO', 'DENZO')   
    keyword_WRITE_REF_DENZO = .true.  
    
   case ('REF_SADABS', 'REF_SAD', 'SADABS')
    keyword_WRITE_REF_SADABS = .true.   

   case ('READ_NREPORT', 'READ_NREPORT_HTML', 'READ_HTMLREPORT')
    keyword_read_NREPORT = .true.

   case ('READ_CEL', 'READ_CEL_FILE', 'READ_POWDERCELL')    
    if(nb_arg == 0) then
     call check_arg('0', 'READ_CEL')
     return
    endif
    i1 = index(arg_string(1), '.')
    if (i1 == 0) then
     CEL_file_name = trim(arg_string(1))//'.CEL'
    else
     READ(arg_string(1), *) CEL_file_name
    endif
    call test_file_exist(trim(CEL_file_name), file_exist)
    if(.not. file_exist) then
     CEL_file_name = ''
     return
    endif 
    keyword_read_CEL  = .true.
    
   case ('READ_CIF', 'READ_CIF_FILE', 'CIF_FILE')
    IF(nb_arg == 0)  then
     call check_arg('0', 'READ_CIF')
     return
    endif
    i1 = INDEX(arg_string(1), '.')
    if (i1==0) then
     CIF_file_name = TRIM(arg_string(1))//'.CIF'
    else
     READ(arg_string(1), *) CIF_file_name
    end if
    call test_file_exist(TRIM(CIF_file_name), file_exist)
    if(.not. file_exist) then
     CIF_file_name = ''
     return
    endif 
    keyword_read_CIF = .true.

   case ('READ_INS', 'READ_INS_FILE', 'INS_FILE')
    IF(nb_arg == 0)  then
     call check_arg('0', 'READ_INS')
     return
    endif
    i1 = INDEX(arg_string(1), '.')
    if (i1==0) then
     INS_file_name = TRIM(arg_string(1))//'.INS'
    else
     READ(arg_string(1), *, iostat=i_error) INS_file_name
    end if
    keyword_read_INS = .true.

   case ('READ_PCR', 'READ_PCR_FILE', 'PCR_FILE')
    IF(nb_arg == 0)  then
     call check_arg('0', 'READ_PCR')
     return
    endif
    i1 = INDEX(arg_string(1), '.PCR')
    if (i1==0) then
     PCR_file_name = TRIM(arg_string(1))//'.PCR'
    else
     READ(arg_string(1), *, iostat=i_error) PCR_file_name
    end if
    call test_file_exist(PCR_file_name, file_exist)
    IF(.NOT. file_exist) then
     PCR_file_name = ''
     return
    endif
    keyword_read_PCR = .true.



   case ('SIR', 'SIR97')
    keyword_SIR = .true.

   case ('XRAYS_WAVELENGTH', 'X_WAVE')
    keyword_X_WAVE = .true.
    IF(nb_arg /= 0) then
     X_target(1: tabulated_target_nb)%write = .false.
     do i=1, tabulated_target_nb
      do j=1, nb_arg
       if (arg_string(j)(1:2) == u_case(X_target(i)%label)) X_target(i)%write = .true.
      END do
     end do
    else
     X_target(1: tabulated_target_nb)%write = .true.
    endif

   CASE ('END')
    return

  ! CASE default
  !  call write_info('')
  !  call write_info('  > Unknown '//TRIM(input_keyword)//' keyword !')
  !  call write_info('')

   case default

    unknown_keyword = .true.

  end select


 return
END subroutine identification_keywords

!------------------------------------------------------------------------------


subroutine check_arg(input_string, input_keyword)
 USE IO_module, ONLY : write_info
 implicit none
  CHARACTER (LEN=*),    INTENT(IN) :: input_string
  CHARACTER (LEN=*),    INTENT(IN) :: input_keyword

  select case (input_string)
      case ('0')
     call write_info('')
     call write_info('  Warning: Argument(s) is(are) mandatory for '//TRIM(input_keyword)// ' keyword ...')
     call write_info('')

      case ('-1')
     call write_info('')
     call write_info('  Warning: At least one argument is mandatory for '//TRIM(input_keyword)// ' keyword ...')
     call write_info('')

      case ('9')
     call write_info('')
     call write_info('  Warning: Arguments (components of the 3*3 matrix) are mandatory for MATRIX keyword ...')
     call write_info('')

      case ('unknown')
     call write_info('')
     call write_info(' Warning: Unknown argument for '//TRIM(input_keyword)// ' ...')
     call write_info('')

      case default
     call write_info('')
     call write_info('  Warning: '//TRIM(input_string)//' arguments are mandatory for '//TRIM(input_keyword)// ' keyword ...')
     call write_info('')

  end select

end subroutine check_arg

!------------------------------------------------------------------------------------------------

subroutine get_SPG(input_string)
 use cryscal_module,                 only  : space_group_symbol, SPG
 use CFML_crystallographic_symmetry, only  : set_spacegroup

 implicit none
  character (len=*), intent(in)  :: input_string
  
     IF(input_string(1:4) == 'TRIC') then
       space_group_symbol = 'P 1'
       call set_spacegroup(space_group_symbol, SPG)
     ELSEIF(input_string(1:4) == 'MONO')  then
      space_group_symbol = 'P 2'
      call set_spacegroup(space_group_symbol, SPG)
     ELSEIF(input_string(1:5) == 'ORTHO') then
      space_group_symbol = 'P 2 2 2'
      call set_spacegroup(space_group_symbol, SPG)
     ELSEIF(input_string(1:5) == 'TETRA') then
      space_group_symbol = 'P 4'
      call set_spacegroup(space_group_symbol, SPG)
     ELSEIF(input_string(1:4) == 'TRIG')  then
      space_group_symbol = 'P 3'
      call set_spacegroup(space_group_symbol, SPG)
     ELSEIF(input_string(1:4) == 'HEXA')  then
      space_group_symbol = 'P 6'
      call set_spacegroup(space_group_symbol, SPG)
     ELSEIF(input_string(1:3) == 'CUB')   then
      space_group_symbol = 'P 2 3'
      call set_spacegroup(space_group_symbol, SPG)
     endif

 return
end subroutine get_SPG

!---------------------------------------------------------------------------------------------------------

subroutine get_matrix_coord(input_line)
 USE cryscal_module, ONLY : mat
 use macros_module,  only : nombre_de_colonnes
 USE Definition_fractions
 
 implicit none
  CHARACTER (LEN=*), INTENT(INOUT)            :: input_line
  CHARACTER (LEN=32)                          :: string

  INTEGER                                     :: i1, i,j, k, nb_col

  call def_fractions

  
  do i=1, 3
   do j=1, 3
    input_line = ADJUSTL(input_line)
    i1 = INDEX(input_line, ' ')

    string = input_line(1:i1-1)
    string = ADJUSTL(string)
    if(string(1:1) == '+') string = string(2:)
    do k = 1, nb_fraction
     IF(string(1:3) == ratio_string_pos(k)(1:3)) then
      Mat(i,j) = ratio_real_pos(k)
      exit 

     ELSEIF(string(1:4) == ratio_string_neg(k)(1:4)) then
      Mat(i,j) = ratio_real_neg(k)
      exit 

     else
      READ(string, *) Mat(i,j)
     endif
    end do
    input_line = input_line(i1+1:)
   end do   ! loop_i
  end do    ! loop_i

  




 return
end subroutine get_matrix_coord
