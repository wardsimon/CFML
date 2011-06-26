!     Last change:  TR   23 Feb 2007   11:00 am
subroutine create_CFL_file(input_file, extension)
 USE cryscal_module
 USE SHELX_module,  ONLY : fmt_SFAC, fmt_UNIT
 USE IO_module

 implicit none
  CHARACTER (LEN=*), INTENT(IN)        :: input_file, extension
  ! local variables
  INTEGER                              :: i, n_atom
  CHARACTER(LEN=256)                   :: created_CFL_file

  i= INDEX(input_file, '.')
  !WRITE(created_CFL_file, '(a)') input_file(1:i-1)//'_ins.CFL'
  WRITE(created_CFL_file, '(a)') input_file(1:i-1)//'_'//TRIM(extension)//'.CFL'

  OPEN(UNIT=4, FILE=TRIM(created_CFL_file))

   WRITE(4, '(2a)')                        'TITL  ', TRIM(main_title)
   WRITE(4, '(a,3(1x,F8.4),3(1x,F8.3))')   'CELL  ', unit_cell%param(1:3), unit_cell%param(4:6)
   IF(keyword_WAVE)   WRITE(4, '(a,F8.5)') 'WAVE  ', wavelength
   IF(keyword_ZUNIT)  WRITE(4, '(a,F4.1)') 'ZUNIT ', Z_unit

   IF(keyword_SIZE) then
   WRITE(4, '(a,3(1x,F6.3))')              'SIZE  ', crystal%size(1:3)
   endif

   IF(nb_atoms_type < 10) then
    WRITE(fmt_SFAC, '(a,i1,a)') "(a,", nb_atoms_type,"(1x,a6))"
    WRITE(fmt_UNIT, '(a,i1,a)') "(a,", nb_atoms_type,"(1x,F6.1))"
   ELSE
    WRITE(fmt_SFAC, '(a,i2,a)') "(a,", nb_atoms_type,"(1x,a6))"
    WRITE(fmt_UNIT, '(a,i2,a)') "(a,", nb_atoms_type,"(1x,F6.1))"
   endif
   IF(keyword_SFAC_UNIT) then
    WRITE(4, fmt_SFAC) 'SFAC ', (SFAC_type(i)  ,i=1,nb_atoms_type)
    WRITE(4, fmt_UNIT) 'UNIT ', (SFAC_number(i),i=1,nb_atoms_type)
   END if

   IF(keyword_SPGR) then
    WRITE(4, '(2a)')   'SPGR ', TRIM(space_group_symbol)
   endif

   IF(nb_atom /=0) then
   do n_atom=1, nb_atom
    if((atom_occ_perc(n_atom) -1) .lt. .0001) then
     WRITE(4,'(a,2a6,4(1x,F8.5))') 'ATOM ', trim(atom_label(n_atom)),trim(atom_type(n_atom)),  &
	                               (atom_coord(i,n_atom),i=1,3), atom_occ_perc(n_atom)
    else
     WRITE(4,'(a,2a6,3(1x,F8.5))') 'ATOM ', trim(atom_label(n_atom)),trim(atom_type(n_atom)),  &
	                               (atom_coord(i,n_atom),i=1,3)
    endif
   END do
   endif


   ! CRYSCAL OUTPUTS
   !WRITE(UNIT=4,'(a)') '!----------------------------------------------------------'
   !WRITE(UNIT=4,'(a)') '!SP_INFO       ! get space group features'
   !WRITE(UNIT=4,'(a)') '!SP_EXTI       ! get space group extinctions'
   !WRITE(UNIT=4,'(a)') '!SYM_OP        ! write symmetry operators'
   !WRITE(UNIT=4,'(a)') '!SITE_INFO     ! get informations about the site symmetry'
   !WRITE(UNIT=4,'(a)') '!GEN_HKL 2theta_min=0. 2theta_max=80.'
  close (UNIT=4)


  call write_info('')
  write(message_text, '(3a)') '   >> ', trim(created_CFL_file), ' file has been created.'  
  call write_info(trim(message_text))
  call write_info('')


 RETURN
end subroutine create_CFL_file

!-----------------------------------------------------------------------------------
subroutine create_FST_file(input_file, extension)
 USE cryscal_module
 USE SHELX_module,  ONLY : fmt_SFAC, fmt_UNIT
 USE IO_module

 implicit none
  CHARACTER (LEN=*), INTENT(IN)        :: input_file, extension
  ! local variables
  INTEGER                              :: i, n_atom
  CHARACTER(LEN=256)                   :: created_FST_file
  LOGICAL                              :: atom_Al, atom_Au, atom_B, atom_Ba, atom_Br, atom_C, atom_Ca, atom_Cl, atom_Co, &
                                          atom_Cu, atom_F, atom_Fe, atom_H, &
                                          atom_Hf, atom_I, atom_In, atom_Ir, atom_K, atom_Li, atom_Mg, atom_Mo, atom_N,  &
										  atom_Na, atom_Nd, atom_Ni, &
                                          atom_O, atom_P, atom_Pd, atom_Pt, &
                                          atom_Re, atom_Ru, atom_S, atom_Se, atom_Si, atom_Sn, atom_Sr, atom_Y, atom_Zn, &
										  atom_Zr
 
  atom_Al = .false.
  atom_Au = .false.
  atom_B  = .false.
  atom_Ba = .false.
  atom_Br = .false.
  atom_C  = .false.
  atom_Ca = .false.
  atom_Cl = .false.
  atom_Co = .false.
  atom_Cu = .false.
  atom_F  = .false.
  atom_Fe = .false.
  atom_Hf = .false.
  atom_I  = .false.
  atom_In = .false.
  atom_Ir = .false.
  atom_K  = .false.
  atom_Li = .false.
  atom_Mg = .false.
  atom_Mo = .false.
  atom_N  = .false.
  atom_Na = .false.
  atom_Nd = .false.
  atom_Ni = .false.
  atom_O  = .false.
  atom_P  = .false.
  atom_Pd = .false.
  atom_Pt = .false.
  atom_Re = .false.
  atom_Ru = .false.
  atom_S  = .false.
  atom_Se = .false.
  atom_Si = .false.
  atom_Sn = .false.
  atom_Sr = .false.
  atom_Y  = .false.
  atom_Zn = .false.
  atom_Zr = .false.


  if(.not. keyword_SPGR) then
   call write_info('')
   call write_info(' Space group has be known. Please enter SPGR instruction !')
   call write_info('')
   return
  endif 
  
  if(nb_atom == 0) then
   call write_info('')
   call write_info(' No atoms. Please enter ATOM instruction !')
   call write_info('')
   return
  endif
  
  
  i= INDEX(input_file, '.')
  WRITE(created_FST_file, '(a)') input_file(1:i-1)//'_'//TRIM(extension)//'.FST'

  OPEN(UNIT=4, FILE=TRIM(created_FST_file))
   WRITE(4, '(a)')                         '!  File for FullProf Studio (created by CRYSCAL)'
   WRITE(4, '(2a)')                        'TITLE   ', TRIM(main_title)
   WRITE(4, '(2a)')                        'SPACEG  ', TRIM(space_group_symbol)       
   WRITE(4, '(a,3(1x,F8.4),3(1x,F8.3),a)') 'CELL  ', unit_cell%param(1:3), unit_cell%param(4:6), '   DISPLAY MULTIPLE'
   WRITE(4, '(a)')                         'BOX  -0.15  1.15  -0.15 1.15  -0.15 1.15'
 
   do n_atom=1, nb_atom    
    WRITE(4,'(a,2a6,3(1x,F8.5))') 'ATOM ', trim(atom_label(n_atom)),trim(atom_type(n_atom)),  (atom_coord(i,n_atom),i=1,3)    
    if(len_trim(atom_type(n_atom)) == 2 .and. trim(atom_type(n_atom)(1:2)) == 'AL') atom_Al =.true.
    if(len_trim(atom_type(n_atom)) == 2 .and. trim(atom_type(n_atom)(1:2)) == 'AU') atom_Au =.true.
    if(len_trim(atom_type(n_atom)) == 1 .and. trim(atom_type(n_atom)(1:1)) == 'B')  atom_B  =.true.
    if(len_trim(atom_type(n_atom)) == 2 .and. trim(atom_type(n_atom)(1:2)) == 'BA') atom_Ba =.true.
    if(len_trim(atom_type(n_atom)) == 2 .and. trim(atom_type(n_atom)(1:2)) == 'BR') atom_Br =.true.
    if(len_trim(atom_type(n_atom)) == 1 .and. trim(atom_type(n_atom)(1:1)) == 'C')  atom_C  =.true.
    if(len_trim(atom_type(n_atom)) == 2 .and. trim(atom_type(n_atom)(1:2)) == 'CA') atom_Ca =.true.
    if(len_trim(atom_type(n_atom)) == 2 .and. trim(atom_type(n_atom)(1:2)) == 'CL') atom_Cl =.true.
    if(len_trim(atom_type(n_atom)) == 2 .and. trim(atom_type(n_atom)(1:2)) == 'CO') atom_Co =.true.
    
    if(len_trim(atom_type(n_atom)) == 2 .and. trim(atom_type(n_atom)(1:2)) == 'CU') atom_Cu =.true.
    if(len_trim(atom_type(n_atom)) == 1 .and. trim(atom_type(n_atom)(1:1)) == 'F')  atom_F  =.true.
    
    if(len_trim(atom_type(n_atom)) == 2 .and. trim(atom_type(n_atom)(1:1)) == 'FE') atom_Fe =.true.
    if(len_trim(atom_type(n_atom)) == 2 .and. trim(atom_type(n_atom)(1:1)) == 'HF') atom_Hf =.true.
    if(len_trim(atom_type(n_atom)) == 1 .and. trim(atom_type(n_atom)(1:2)) == 'I')  atom_I  =.true.
    if(len_trim(atom_type(n_atom)) == 2 .and. trim(atom_type(n_atom)(1:2)) == 'IN') atom_In =.true.
    if(len_trim(atom_type(n_atom)) == 2 .and. trim(atom_type(n_atom)(1:2)) == 'IR') atom_Ir =.true.
    if(len_trim(atom_type(n_atom)) == 1 .and. trim(atom_type(n_atom)(1:1)) == 'K')  atom_K  =.true.

    if(len_trim(atom_type(n_atom)) == 2 .and. trim(atom_type(n_atom)(1:2)) == 'LI') atom_Li =.true.
    if(len_trim(atom_type(n_atom)) == 2 .and. trim(atom_type(n_atom)(1:2)) == 'MO') atom_Mo =.true.
    if(len_trim(atom_type(n_atom)) == 2 .and. trim(atom_type(n_atom)(1:2)) == 'MG') atom_Mg =.true.
    if(len_trim(atom_type(n_atom)) == 1 .and. trim(atom_type(n_atom)(1:1)) == 'N')  atom_N  =.true.
    if(len_trim(atom_type(n_atom)) == 2 .and. trim(atom_type(n_atom)(1:2)) == 'NA') atom_Na =.true.
    if(len_trim(atom_type(n_atom)) == 2 .and. trim(atom_type(n_atom)(1:2)) == 'ND') atom_Nd =.true.
    if(len_trim(atom_type(n_atom)) == 2 .and. trim(atom_type(n_atom)(1:2)) == 'NI') atom_Ni =.true.

    if(len_trim(atom_type(n_atom)) == 1 .and. trim(atom_type(n_atom)(1:1)) == 'O')  atom_O  =.true.
    if(len_trim(atom_type(n_atom)) == 1 .and. trim(atom_type(n_atom)(1:1)) == 'P')  atom_P  =.true.
    if(len_trim(atom_type(n_atom)) == 2 .and. trim(atom_type(n_atom)(1:2)) == 'PD') atom_Pd =.true.
    if(len_trim(atom_type(n_atom)) == 2 .and. trim(atom_type(n_atom)(1:2)) == 'PT') atom_Pt =.true.
    if(len_trim(atom_type(n_atom)) == 2 .and. trim(atom_type(n_atom)(1:2)) == 'RE') atom_Re =.true.
    if(len_trim(atom_type(n_atom)) == 2 .and. trim(atom_type(n_atom)(1:2)) == 'RU') atom_Ru =.true.
    if(len_trim(atom_type(n_atom)) == 1 .and. trim(atom_type(n_atom)(1:1)) == 'S')  atom_S  =.true.
    if(len_trim(atom_type(n_atom)) == 2 .and. trim(atom_type(n_atom)(1:2)) == 'SE') atom_Se =.true.
    if(len_trim(atom_type(n_atom)) == 2 .and. trim(atom_type(n_atom)(1:2)) == 'SI') atom_Si =.true.
    if(len_trim(atom_type(n_atom)) == 2 .and. trim(atom_type(n_atom)(1:2)) == 'SN') atom_Sn =.true.
    if(len_trim(atom_type(n_atom)) == 2 .and. trim(atom_type(n_atom)(1:2)) == 'SR') atom_Sr =.true.
    if(len_trim(atom_type(n_atom)) == 1 .and. trim(atom_type(n_atom)(1:1)) == 'Y')  atom_Y  =.true.

    if(len_trim(atom_type(n_atom)) == 2 .and. trim(atom_type(n_atom)(1:2)) == 'ZN') atom_Zn =.true.
    if(len_trim(atom_type(n_atom)) == 2 .and. trim(atom_type(n_atom)(1:2)) == 'ZR') atom_Zr =.true.

    if(trim(atom_type(n_atom)(1:1)) == 'H') atom_H  =.true.    
   END do

   if(atom_C)                WRITE(4, '(a)')   'CONN C   C   0.   2.  RADIUS 1.'
   if(atom_C .and. atom_B)   WRITE(4, '(a)')   'CONN C   B   0.   2.  RADIUS 1.'
   if(atom_C .and. atom_Br)  WRITE(4, '(a)')   'CONN C   BR  0.   2.  RADIUS 1.'
   if(atom_C .and. atom_F)   WRITE(4, '(a)')   'CONN C   F   0.   2.  RADIUS 1.'
   if(atom_C .and. atom_N)   WRITE(4, '(a)')   'CONN C   N   0.   2.  RADIUS 1.'
   if(atom_C .and. atom_O)   WRITE(4, '(a)')   'CONN C   O   0.   2.  RADIUS 1.'
   if(atom_C .and. atom_P)   WRITE(4, '(a)')   'CONN C   P   0.   2.  RADIUS 1.'
   if(atom_C .and. atom_S)   WRITE(4, '(a)')   'CONN C   S   0.   2.  RADIUS 1.'
   if(atom_C .and. atom_Se)  WRITE(4, '(a)')   'CONN C   SE  0.   2.  RADIUS 1.'
   if(atom_C .and. atom_H)   WRITE(4, '(a)')   'CONN C   H   0.   1.4 RADIUS 1.'
   
   if(atom_B)                WRITE(4, '(a)')   'CONN B   B   0.   2.  RADIUS 1.'
   if(atom_B .and. atom_F)   WRITE(4, '(a)')   'CONN B   F   0.   2.  RADIUS 1.'
   if(atom_B .and. atom_N)   WRITE(4, '(a)')   'CONN B   N   0.   2.  RADIUS 1.'
   if(atom_B .and. atom_O)   WRITE(4, '(a)')   'CONN B   O   0.   2.  RADIUS 1.'
   
   if(atom_N)                WRITE(4, '(a)')   'CONN N   N   0.   2.  RADIUS 1.'
   if(atom_N .and. atom_O)   WRITE(4, '(a)')   'CONN N   O   0.   2.  RADIUS 1.'   
   if(atom_N .and. atom_H)   WRITE(4, '(a)')   'CONN N   H   0.   1.4 RADIUS 1.'

   if(atom_Al .and. atom_C)  WRITE(4, '(a)')   'CONN AL  C   0.   2.6 RADIUS 1.'
   if(atom_Al .and. atom_Cl) WRITE(4, '(a)')   'CONN AL  CL  0.   2.6 RADIUS 1.'
   if(atom_Al .and. atom_O)  WRITE(4, '(a)')   'CONN AL  O   0.   2.6 RADIUS 1.'
   if(atom_Al .and. atom_P)  WRITE(4, '(a)')   'CONN AL  P   0.   2.6 RADIUS 1.'
   if(atom_Al .and. atom_S)  WRITE(4, '(a)')   'CONN AL  S   0.   2.6 RADIUS 1.'
   if(atom_Al .and. atom_Se) WRITE(4, '(a)')   'CONN AL  SE  0.   2.6 RADIUS 1.'



   
   if(atom_Au .and. atom_C)  WRITE(4, '(a)')   'CONN AU  C   0.   2.6 RADIUS 1.'
   if(atom_Au .and. atom_Cl) WRITE(4, '(a)')   'CONN AU  CL  0.   2.6 RADIUS 1.'
   if(atom_Au .and. atom_I)  WRITE(4, '(a)')   'CONN AU  I   0.   2.6 RADIUS 1.'
   if(atom_Au .and. atom_O)  WRITE(4, '(a)')   'CONN AU  O   0.   2.6 RADIUS 1.'
   if(atom_Au .and. atom_S)  WRITE(4, '(a)')   'CONN AU  S   0.   2.6 RADIUS 1.'
   if(atom_Au .and. atom_Se) WRITE(4, '(a)')   'CONN AU  SE  0.   2.6 RADIUS 1.'

   if(atom_Ba .and. atom_N)  WRITE(4, '(a)')   'CONN BA  N   0.   2.9 RADIUS 1.'
   if(atom_Ba .and. atom_O)  WRITE(4, '(a)')   'CONN BA  O   0.   2.9 RADIUS 1.'

   if(atom_Ca .and. atom_N)  WRITE(4, '(a)')   'CONN CA  N   0.   2.9 RADIUS 1.'
   if(atom_Ca .and. atom_O)  WRITE(4, '(a)')   'CONN CA  O   0.   2.9 RADIUS 1.'

   if(atom_Cl .and. atom_O)  WRITE(4, '(a)')   'CONN CL  O   0.   2.0 RADIUS 1.'

   if(atom_Cu .and. atom_C)  WRITE(4, '(a)')   'CONN CU  C   0.   2.6 RADIUS 1.'
   if(atom_Cu .and. atom_Cl) WRITE(4, '(a)')   'CONN CU  CL  0.   2.6 RADIUS 1.'
   if(atom_Cu .and. atom_F)  WRITE(4, '(a)')   'CONN CU  F   0.   2.6 RADIUS 1.'
   if(atom_Cu .and. atom_N)  WRITE(4, '(a)')   'CONN CU  N   0.   2.6 RADIUS 1.'
   if(atom_Cu .and. atom_O)  WRITE(4, '(a)')   'CONN CU  O   0.   2.6 RADIUS 1.'
   if(atom_Cu .and. atom_P)  WRITE(4, '(a)')   'CONN CU  P   0.   2.6 RADIUS 1.'
   if(atom_Cu .and. atom_Se) WRITE(4, '(a)')   'CONN CU  SE  0.   2.6 RADIUS 1.'

   if(atom_Co .and. atom_C)  WRITE(4, '(a)')   'CONN CO  C   0.   2.6 RADIUS 1.'
   if(atom_Co .and. atom_Cl) WRITE(4, '(a)')   'CONN CO  CL  0.   2.6 RADIUS 1.'
   if(atom_Co .and. atom_F)  WRITE(4, '(a)')   'CONN CO  F   0.   2.6 RADIUS 1.'
   if(atom_Co .and. atom_N)  WRITE(4, '(a)')   'CONN CO  N   0.   2.6 RADIUS 1.'
   if(atom_Co .and. atom_O)  WRITE(4, '(a)')   'CONN CO  O   0.   2.6 RADIUS 1.'
   if(atom_Co .and. atom_P)  WRITE(4, '(a)')   'CONN CO  P   0.   2.6 RADIUS 1.'
   if(atom_Co .and. atom_Se) WRITE(4, '(a)')   'CONN CO  SE  0.   2.6 RADIUS 1.'


   if(atom_Fe .and. atom_C)  WRITE(4, '(a)')   'CONN FE  C   0.   2.6 RADIUS 1.'
   if(atom_Fe .and. atom_Cl) WRITE(4, '(a)')   'CONN FE  CL  0.   2.6 RADIUS 1.'
   if(atom_Fe .and. atom_N)  WRITE(4, '(a)')   'CONN FE  N   0.   2.6 RADIUS 1.'
   if(atom_Fe .and. atom_O)  WRITE(4, '(a)')   'CONN FE  O   0.   2.6 RADIUS 1.'
   if(atom_Fe .and. atom_P)  WRITE(4, '(a)')   'CONN FE  P   0.   2.6 RADIUS 1.'
   if(atom_Fe .and. atom_Se) WRITE(4, '(a)')   'CONN FE  SE  0.   2.6 RADIUS 1.'

   if(atom_Hf .and. atom_Cl) WRITE(4, '(a)')   'CONN HF  CL  0.   2.9 RADIUS 1.'
   if(atom_Hf .and. atom_N)  WRITE(4, '(a)')   'CONN HF  N   0.   2.9 RADIUS 1.'
   if(atom_Hf .and. atom_O)  WRITE(4, '(a)')   'CONN HF  O   0.   2.9 RADIUS 1.'
   if(atom_Hf .and. atom_P)  WRITE(4, '(a)')   'CONN HF  P   0.   2.6 RADIUS 1.'
   if(atom_Hf .and. atom_S)  WRITE(4, '(a)')   'CONN HF  S   0.   2.6 RADIUS 1.'
   if(atom_Hf .and. atom_Se) WRITE(4, '(a)')   'CONN HF  SE  0.   2.6 RADIUS 1.'

   if(atom_I)                WRITE(4, '(a)')   'CONN I   I   0.   2.6 RADIUS 1.'
   if(atom_I  .and. atom_C)  WRITE(4, '(a)')   'CONN I   C   0.   2.6 RADIUS 1.'


   if(atom_In .and. atom_C)  WRITE(4, '(a)')   'CONN IN  C   0.   2.6 RADIUS 1.'
   if(atom_In .and. atom_Cl) WRITE(4, '(a)')   'CONN IN  CL  0.   2.6 RADIUS 1.'
   if(atom_In .and. atom_O)  WRITE(4, '(a)')   'CONN IN  O   0.   2.6 RADIUS 1.'
   if(atom_In .and. atom_P)  WRITE(4, '(a)')   'CONN IN  P   0.   2.6 RADIUS 1.'
   if(atom_In .and. atom_S)  WRITE(4, '(a)')   'CONN IN  S   0.   2.6 RADIUS 1.'
   if(atom_In .and. atom_Se) WRITE(4, '(a)')   'CONN IN  SE  0.   2.6 RADIUS 1.'

   if(atom_Ir .and. atom_C)  WRITE(4, '(a)')   'CONN IR  C   0.   2.6 RADIUS 1.'
   if(atom_Ir .and. atom_Cl) WRITE(4, '(a)')   'CONN IR  CL  0.   2.6 RADIUS 1.'
   if(atom_Ir .and. atom_O)  WRITE(4, '(a)')   'CONN IR  O   0.   2.6 RADIUS 1.'
   if(atom_Ir .and. atom_P)  WRITE(4, '(a)')   'CONN IR  P   0.   2.6 RADIUS 1.'
   if(atom_Ir .and. atom_S)  WRITE(4, '(a)')   'CONN IR  S   0.   2.6 RADIUS 1.'
   if(atom_Ir .and. atom_Se) WRITE(4, '(a)')   'CONN IR  SE  0.   2.6 RADIUS 1.'

   if(atom_K  .and. atom_C)  WRITE(4, '(a)')   'CONN K   C   0.   2.9 RADIUS 1.'
   if(atom_K  .and. atom_N)  WRITE(4, '(a)')   'CONN K   N   0.   2.9 RADIUS 1.'
   if(atom_K  .and. atom_O)  WRITE(4, '(a)')   'CONN K   O   0.   2.9 RADIUS 1.'

   if(atom_Li .and. atom_N)  WRITE(4, '(a)')   'CONN LI  N   0.   2.5 RADIUS 1.'
   if(atom_Li .and. atom_O)  WRITE(4, '(a)')   'CONN LI  O   0.   2.5 RADIUS 1.'

   if(atom_Mg .and. atom_C)  WRITE(4, '(a)')   'CONN MG  C   0.   2.6 RADIUS 1.'
   if(atom_Mg .and. atom_Cl) WRITE(4, '(a)')   'CONN MG  CL  0.   2.6 RADIUS 1.'
   if(atom_Mg .and. atom_O)  WRITE(4, '(a)')   'CONN MG  O   0.   2.6 RADIUS 1.'
   if(atom_Mg .and. atom_P)  WRITE(4, '(a)')   'CONN MG  P   0.   2.6 RADIUS 1.'
   if(atom_Mg .and. atom_S)  WRITE(4, '(a)')   'CONN MG  S   0.   2.6 RADIUS 1.'
   if(atom_Mg .and. atom_Se) WRITE(4, '(a)')   'CONN MG  SE  0.   2.6 RADIUS 1.'

   if(atom_Mo .and. atom_C)  WRITE(4, '(a)')   'CONN MO  C   0.   2.6 RADIUS 1.'
   if(atom_Mo .and. atom_Cl) WRITE(4, '(a)')   'CONN MO  CL  0.   2.6 RADIUS 1.'
   if(atom_Mo .and. atom_F)  WRITE(4, '(a)')   'CONN MO  F   0.   2.6 RADIUS 1.'
   if(atom_Mo .and. atom_N)  WRITE(4, '(a)')   'CONN MO  N   0.   2.6 RADIUS 1.'
   if(atom_Mo .and. atom_O)  WRITE(4, '(a)')   'CONN MO  O   0.   2.6 RADIUS 1.'
   if(atom_Mo .and. atom_P)  WRITE(4, '(a)')   'CONN MO  P   0.   2.6 RADIUS 1.'
   if(atom_Mo .and. atom_Se) WRITE(4, '(a)')   'CONN MO  SE  0.   2.6 RADIUS 1.'

   if(atom_ND .and. atom_Cl) WRITE(4, '(a)')   'CONN ND  CL  0.   2.9 RADIUS 1.'
   if(atom_ND .and. atom_N)  WRITE(4, '(a)')   'CONN ND  N   0.   2.9 RADIUS 1.'
   if(atom_ND .and. atom_O)  WRITE(4, '(a)')   'CONN ND  O   0.   2.9 RADIUS 1.'
   if(atom_ND .and. atom_P)  WRITE(4, '(a)')   'CONN ND  P   0.   2.6 RADIUS 1.'
   if(atom_ND .and. atom_S)  WRITE(4, '(a)')   'CONN ND  S   0.   2.6 RADIUS 1.'
   if(atom_ND .and. atom_Se) WRITE(4, '(a)')   'CONN ND  SE  0.   2.6 RADIUS 1.'

   if(atom_Ni .and. atom_C)  WRITE(4, '(a)')   'CONN NI  C   0.   2.6 RADIUS 1.'
   if(atom_Ni .and. atom_Cl) WRITE(4, '(a)')   'CONN NI  CL  0.   2.6 RADIUS 1.'
   if(atom_Ni .and. atom_F)  WRITE(4, '(a)')   'CONN NI  F   0.   2.6 RADIUS 1.'
   if(atom_Ni .and. atom_N)  WRITE(4, '(a)')   'CONN NI  N   0.   2.6 RADIUS 1.'
   if(atom_Ni .and. atom_O)  WRITE(4, '(a)')   'CONN NI  O   0.   2.6 RADIUS 1.'
   if(atom_Ni .and. atom_P)  WRITE(4, '(a)')   'CONN NI  P   0.   2.6 RADIUS 1.'
   if(atom_Ni .and. atom_S)  WRITE(4, '(a)')   'CONN NI  S   0.   2.6 RADIUS 1.'
   if(atom_Ni .and. atom_Se) WRITE(4, '(a)')   'CONN NI  SE  0.   2.6 RADIUS 1.'


   if(atom_P  .and. atom_F)  WRITE(4, '(a)')   'CONN P   F   0.   2.0 RADIUS 1.'
   if(atom_P  .and. atom_O)  WRITE(4, '(a)')   'CONN P   O   0.   2.0 RADIUS 1.'

   if(atom_Pd .and. atom_C)  WRITE(4, '(a)')   'CONN PD  C   0.   2.6 RADIUS 1.' 
   if(atom_Pd .and. atom_Cl) WRITE(4, '(a)')   'CONN PD  CL  0.   2.6 RADIUS 1.'
   if(atom_Pd .and. atom_N)  WRITE(4, '(a)')   'CONN PD  N   0.   2.6 RADIUS 1.'   
   if(atom_Pd .and. atom_P)  WRITE(4, '(a)')   'CONN PD  P   0.   2.6 RADIUS 1.'
   if(atom_Pd .and. atom_S)  WRITE(4, '(a)')   'CONN PD  S   0.   2.6 RADIUS 1.'
   if(atom_Pd .and. atom_Se) WRITE(4, '(a)')   'CONN PD  SE  0.   2.6 RADIUS 1.'
  
   if(atom_Pt .and. atom_C)  WRITE(4, '(a)')   'CONN PT  C   0.   2.6 RADIUS 1.' 
   if(atom_Pt .and. atom_Cl) WRITE(4, '(a)')   'CONN PT  CL  0.   2.6 RADIUS 1.'
   if(atom_Pt .and. atom_N)  WRITE(4, '(a)')   'CONN PT  N   0.   2.6 RADIUS 1.'
   if(atom_Pt .and. atom_O)  WRITE(4, '(a)')   'CONN PT  O   0.   2.6 RADIUS 1.'
   if(atom_Pt .and. atom_P)  WRITE(4, '(a)')   'CONN PT  P   0.   2.6 RADIUS 1.'
   if(atom_Pt .and. atom_S)  WRITE(4, '(a)')   'CONN PT  S   0.   2.6 RADIUS 1.'   
   if(atom_Pt .and. atom_Se) WRITE(4, '(a)')   'CONN PT  SE  0.   2.6 RADIUS 1.'

   if(atom_Re .and. atom_C)  WRITE(4, '(a)')   'CONN RE  C   0.   2.6 RADIUS 1.' 
   if(atom_Re .and. atom_Cl) WRITE(4, '(a)')   'CONN RE  CL  0.   2.6 RADIUS 1.'
   if(atom_Re .and. atom_N)  WRITE(4, '(a)')   'CONN RE  N   0.   2.6 RADIUS 1.'
   if(atom_Re .and. atom_O)  WRITE(4, '(a)')   'CONN RE  N   0.   2.6 RADIUS 1.'
   if(atom_Re .and. atom_P)  WRITE(4, '(a)')   'CONN RE  P   0.   2.6 RADIUS 1.'
   if(atom_Re .and. atom_S)  WRITE(4, '(a)')   'CONN RE  S   0.   2.6 RADIUS 1.'   
   if(atom_Re .and. atom_Se) WRITE(4, '(a)')   'CONN RE  SE  0.   2.6 RADIUS 1.'

   if(atom_Ru .and. atom_C)  WRITE(4, '(a)')   'CONN RU  C   0.   2.6 RADIUS 1.'
   if(atom_Ru .and. atom_Cl) WRITE(4, '(a)')   'CONN RU  CL  0.   2.6 RADIUS 1.'
   if(atom_Ru .and. atom_N)  WRITE(4, '(a)')   'CONN RU  N   0.   2.6 RADIUS 1.'
   if(atom_Ru .and. atom_O)  WRITE(4, '(a)')   'CONN RU  O   0.   2.6 RADIUS 1.'
   if(atom_Ru .and. atom_P)  WRITE(4, '(a)')   'CONN RU  P   0.   2.6 RADIUS 1.'
   if(atom_Ru .and. atom_S)  WRITE(4, '(a)')   'CONN RU  S   0.   2.6 RADIUS 1.'
   if(atom_Ru .and. atom_Se) WRITE(4, '(a)')   'CONN RU  SE  0.   2.6 RADIUS 1.'

   if(atom_S  .and. atom_F)  WRITE(4, '(a)')   'CONN S   F   0.   2.0 RADIUS 1.'
   if(atom_S  .and. atom_O)  WRITE(4, '(a)')   'CONN S   O   0.   2.0 RADIUS 1.'

   if(atom_Si)               WRITE(4, '(a)')   'CONN SI  SI  0.   2.  RADIUS 1.'
   if(atom_Si .and. atom_C)  WRITE(4, '(a)')   'CONN SI  C   0.   2.  RADIUS 1.'
   if(atom_Si .and. atom_N)  WRITE(4, '(a)')   'CONN SI  N   0.   2.  RADIUS 1.'
   if(atom_Si .and. atom_O)  WRITE(4, '(a)')   'CONN SI  O   0.   2.  RADIUS 1.'
   
   
   
   if(atom_Sn .and. atom_C)  WRITE(4, '(a)')   'CONN SN  C  0.   2.  RADIUS 1.'
   if(atom_Sn .and. atom_Cl) WRITE(4, '(a)')   'CONN SN  CL 0.   2.  RADIUS 1.'
   if(atom_Sn .and. atom_F)  WRITE(4, '(a)')   'CONN SN  F  0.   2.  RADIUS 1.'
   if(atom_Sn .and. atom_N)  WRITE(4, '(a)')   'CONN SN  N  0.   2.  RADIUS 1.'
   if(atom_Sn .and. atom_O)  WRITE(4, '(a)')   'CONN SN  O  0.   2.  RADIUS 1.'
   if(atom_Sn .and. atom_S)  WRITE(4, '(a)')   'CONN SN  S  0.   2.  RADIUS 1.'
   if(atom_Sn .and. atom_Se) WRITE(4, '(a)')   'CONN SN  SE 0.   2.  RADIUS 1.'

   if(atom_Sr .and. atom_N)  WRITE(4, '(a)')   'CONN SR  N   0.   2.9 RADIUS 1.'
   if(atom_Sr .and. atom_O)  WRITE(4, '(a)')   'CONN SR  O   0.   2.9 RADIUS 1.'
   if(atom_Sr .and. atom_Cl) WRITE(4, '(a)')   'CONN SR  CL  0.   2.9 RADIUS 1.'

   if(atom_Y  .and. atom_Cl) WRITE(4, '(a)')   'CONN Y   CL  0.   2.9 RADIUS 1.'
   if(atom_Y  .and. atom_N)  WRITE(4, '(a)')   'CONN Y   N   0.   2.9 RADIUS 1.'
   if(atom_Y  .and. atom_O)  WRITE(4, '(a)')   'CONN Y   O   0.   2.9 RADIUS 1.'
   if(atom_Y  .and. atom_P)  WRITE(4, '(a)')   'CONN Y   P   0.   2.6 RADIUS 1.'
   if(atom_Y  .and. atom_S)  WRITE(4, '(a)')   'CONN Y   S   0.   2.6 RADIUS 1.'
   if(atom_Y  .and. atom_Se) WRITE(4, '(a)')   'CONN Y   SE  0.   2.6 RADIUS 1.'
   

   if(atom_Zn .and. atom_Cl) WRITE(4, '(a)')   'CONN ZN  CL  0.   2.6 RADIUS 1.'
   if(atom_Zn .and. atom_N)  WRITE(4, '(a)')   'CONN ZN  N   0.   2.6 RADIUS 1.'
   if(atom_Zn .and. atom_O)  WRITE(4, '(a)')   'CONN ZN  O   0.   2.6 RADIUS 1.'
   if(atom_Zn .and. atom_P)  WRITE(4, '(a)')   'CONN ZR  P   0.   2.6 RADIUS 1.'
   if(atom_Zn .and. atom_S)  WRITE(4, '(a)')   'CONN ZN  S   0.   2.6 RADIUS 1.'
   if(atom_Zn .and. atom_Se) WRITE(4, '(a)')   'CONN ZN  SE  0.   2.6 RADIUS 1.'

   if(atom_Zr .and. atom_Cl) WRITE(4, '(a)')   'CONN ZR  CL  0.   2.9 RADIUS 1.'
   if(atom_Zr .and. atom_N)  WRITE(4, '(a)')   'CONN ZR  N   0.   2.9 RADIUS 1.'
   if(atom_Zr .and. atom_O)  WRITE(4, '(a)')   'CONN ZR  O   0.   2.9 RADIUS 1.'
   if(atom_Zr .and. atom_P)  WRITE(4, '(a)')   'CONN ZR  P   0.   2.6 RADIUS 1.'
   if(atom_Zr .and. atom_S)  WRITE(4, '(a)')   'CONN ZR  S   0.   2.6 RADIUS 1.'
   if(atom_Zr .and. atom_Se) WRITE(4, '(a)')   'CONN ZN  SE  0.   2.6 RADIUS 1.'
   
 
   WRITE(4, '(a)') '!MOLECULE'   
   WRITE(4, '(a)') '!CONN Ni O 0. 2.3'
   WRITE(4, '(a)') '!POLY Ni1 COLOR 1 1 0 0.5'
   
  close (UNIT=4)


  call write_info('')
  write(message_text, '(3a)') '   >> ', trim(created_FST_file), ' file has been created.'  
  call write_info(trim(message_text))
  call write_info('')


 RETURN
end subroutine create_FST_file

!------------------------------------------------------------------------
subroutine create_INS_file(input_file, extension)
 USE cryscal_module
! USE CFML_IO_FORMATS, ONLY : Write_Shx_Template
 USE SHELX_module
 USE IO_module

 implicit none
  CHARACTER (LEN=*), INTENT(IN)        :: input_file, extension
  ! local variables
  INTEGER                              :: i, n_atom, atom_order
  CHARACTER(LEN=256)                   :: created_INS_file

  i= INDEX(input_file, '.')
  WRITE(created_INS_file, '(a)') input_file(1:i-1)//'_'//TRIM(extension)//'.INS'

  !call Write_SHX_template(TRIM(created_INS_file), 2, TRIM(main_title), wavelength, INT(Z_unit), crystal_cell, SPG, Atm_list)
  OPEN(UNIT=4, FILE=TRIM(created_INS_file))

   WRITE(4, '(2a)')                                   'TITL ', TRIM(main_title)
   if ((wavelength < 0.01)) then
    WRITE(4, '(a,F8.5, 3(1x,F8.4),3(1x,F8.3))')        'CELL ', 0.71073, unit_cell%param(1:3), unit_cell%param(4:6)
   else
    WRITE(4, '(a,F8.5, 3(1x,F8.4),3(1x,F8.3))')        'CELL ', wavelength, unit_cell%param(1:3), unit_cell%param(4:6)
   endif
   IF(keyword_ZUNIT)  WRITE(4, '(a,F8.2,3(1x,F8.4),3(1x,F8.3))')   'ZERR ', Z_unit, unit_cell%param_ESD(1:3), &
                                                                    unit_cell%param_ESD(4:6)

   IF(keyword_SPGR) then
    call get_SHELX_nlatt()
    WRITE(4, '(a,I6)')   'LATT ', n_latt

   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! operateurs de symetrie !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   do i=1, SPG%multip
    write(4, '(2a)') 'SYMM ', SPG%SymopSymb(i)
   end do

   IF(nb_atoms_type == 0) then
    WRITE(4, '(a)') 'SFAC '
    WRITE(4, '(a)') 'UNIT '
   ELSE
    IF(nb_atoms_type < 10) then
     WRITE(fmt_SFAC, '(a,i1,a)') "(a,", nb_atoms_type,"(1x,a6))"
     WRITE(fmt_UNIT, '(a,i1,a)') "(a,", nb_atoms_type,"(1x,F6.1))"
    ELSE
     WRITE(fmt_SFAC, '(a,i2,a)') "(a,", nb_atoms_type,"(1x,a6))"
     WRITE(fmt_UNIT, '(a,i2,a)') "(a,", nb_atoms_type,"(1x,F6.1))"
    ENDIF 
    WRITE(4, fmt_SFAC) 'SFAC ', (SFAC_type(i)  ,i=1,nb_atoms_type)
    WRITE(4, fmt_UNIT) 'UNIT ', (SFAC_number(i),i=1,nb_atoms_type)
   endif
   
   write(4, "(a)")    'WGHT 0.2'
   write(4, "(a)")    'FVAR 1.0'

             !---- Weight ----!

  do i=1,nb_atom
    !write(UNIT=ins_unit,'(a,1x,a4,5F10.5)')  atom_label(i), atom_type (i) , new_atom_coord(1:3,i), 10.+atom_occ(i) ,  atom_Ueq(i)
    call get_atom_order(TRIM(atom_type(i)), atom_order)

    write(4,'(a,1x,I4,5F10.5)')  atom_label(i), atom_order , atom_coord(1:3,i), 10.+atom_occ(i) ,   atom_adp_equiv(i)
  end do

  WRITE(4, '(a)') 'HKLF  4'
  WRITE(4, '(a)') 'END'


  close (UNIT=4)
  
  call write_info('')
  write(message_text, '(3a)') '   >> ', trim(created_INS_file), ' file has been created.'  
  call write_info(trim(message_text))
  call write_info('')

  return
end subroutine create_INS_file


!------------------------------------------------------------------------
subroutine create_SIR_file()
 USE cryscal_module, ONLY : keyword_CELL, keyword_SPGR, keyword_FILE, keyword_CONT, keyword_SFAC_UNIT, keyword_CHEM,   &
                            keyword_ZUNIT, SPG, molecule, unit_cell
 USE HKL_module,     ONLY : HKL_file
 USE IO_module
 implicit none
  INTEGER :: i1

 IF(.NOT. keyword_CELL) then
  call write_info('')
  call write_info('  >>> CELL keyword is mandatory to create an input file for SIR97 <<<')
  call write_info('')
  return
 endif

 IF(.NOT. keyword_SPGR) then
  call write_info('')
  call write_info('  >>> SPGR keyword is mandatory to create an input file for SIR97 <<<')
  call write_info('')
  return
 endif

 IF(.NOT. keyword_FILE) then
  call write_info('')
  call write_info('  >>> FILE keyword is mandatory to create an input file for SIR97 <<<')
  call write_info('')
  return
 endif

 IF(.NOT. keyword_CONT .and. .not. keyword_SFAC_UNIT ) then
  IF(.NOT. keyword_CHEM .or. (keyword_CHEM .and. .not. keyword_ZUNIT)) then
   call write_info('')
   call write_info('  >>> Atom content (keywords: CONT, SFAC/UNIT, CHEM/ZUNIT) is mandatory to create an input file for SIR97 <<<')
   call write_info('')
   return
  endif
 endif

 i1 = INDEX(HKL_file%name, '.')

  OPEN(UNIT=51, FILE='cryscal_SIR97.in')
   WRITE(51, '(a)' )       '%window'
   WRITE(51, '(2a)')       '%structure ',  HKL_file%name(1:i1-1)
   WRITE(51, '(a)' )       '%init'
   WRITE(51, '(4a)')       '%job ',  HKL_file%name(1:i1-1), '   in ', TRIM(SPG%SPG_symb)
   WRITE(51, '(a)' )       '%data'
   WRITE(51, '(a,6F12.5)') '      Cell         ',  unit_cell%param(1:6)
   WRITE(51, '(a)' )       '      Space        ',  TRIM(SPG%SPG_symb)
   WRITE(51, '(a)' )       '      Content      ',  molecule%content
   WRITE(51, '(2a)')       '      Reflections  ',  TRIM(HKl_file%name)
   WRITE(51, '(a)' )       '      Format (3i4,2F8.2)'
   WRITE(51, '(a)' )       '      Fosquare'
   WRITE(51, '(a)' )       '%normal'
   WRITE(51, '(a)' )       '%invariants'
   WRITE(51, '(a)' )       '%phase'
   WRITE(51, '(a)' )       '%fourier'
   WRITE(51, '(a)' )       '%menu'
   WRITE(51, '(2a)')       '      shelx  ',   HKl_file%name(1:i1-1)//'.res'
   WRITE(51, '(a)' )       '%end'
  close (51)

 RETURN
end subroutine create_SIR_file

!------------------------------------------------------------------------
!-------------------------------------------------------------

subroutine get_SHELX_nlatt()
 USE cryscal_module, ONLY : SPG
 use SHELX_module,   only : n_latt

       !---- Latt ----!
       n_latt=1
       select case (SPG%centred)
          case (2) ! Centric

          case (1) ! Acentric
             n_latt=-1
       end select

       select case (SPG%spg_lat)
          case ("P")

          case ("I")
             n_latt=2*n_latt

          case ("R")
             n_latt=3*n_latt

          case ("F")
             n_latt=4*n_latt

          case ("A")
             n_latt=5*n_latt

          case ("B")
             n_latt=6*n_latt

          case ("C")
             n_latt=7*n_latt

       end select

END subroutine get_SHELX_nlatt

