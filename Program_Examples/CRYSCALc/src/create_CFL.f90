!     Last change:  TR   23 Feb 2007   11:00 am
subroutine create_CFL_file(input_file, extension)
 USE cryscalc_module
 USE SHELX_module,  ONLY : fmt_SFAC, fmt_UNIT
 USE IO_module

 implicit none
  CHARACTER (LEN=*), INTENT(IN)        :: input_file, extension
  ! local variables
  INTEGER                              :: i, i1, n_atom
  CHARACTER(LEN=256)                   :: created_CFL_file

  if(debug_proc%level_2)  call write_debug_proc_level(2, "create_CFL_file")

  i1= INDEX(input_file, '.')
  if(i1 /=0) then
   !WRITE(created_CFL_file, '(a)') input_file(1:i1-1)//'_ins.CFL'
   WRITE(created_CFL_file, '(a)') input_file(1:i1-1)//'_'//TRIM(extension)//'.CFL'
  else
   WRITE(created_CFL_file, '(a)') 'cryscalc_.CFL'
  end if

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
    WRITE(4, fmt_SFAC) 'SFAC ', (SFAC%type(i)  ,i=1,nb_atoms_type)
    WRITE(4, fmt_UNIT) 'UNIT ', (SFAC%number(i),i=1,nb_atoms_type)
   END if

   IF(keyword_SPGR) then
    WRITE(4, '(2a)')   'SPGR ', TRIM(space_group_symbol)
   endif

   IF(nb_atom /=0) then
   do n_atom=1, nb_atom
    if((atom_occ_perc(n_atom) -1) .lt. .0001) then
     WRITE(4,'(a,2a6,4(1x,F8.5))') 'ATOM ', trim(atom_label(n_atom)),trim(atom_typ(n_atom)),  &
                                   (atom_coord(i,n_atom),i=1,3), atom_occ_perc(n_atom)
    else
     WRITE(4,'(a,2a6,3(1x,F8.5))') 'ATOM ', trim(atom_label(n_atom)),trim(atom_typ(n_atom)),  &
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
 USE cryscalc_module
 USE SHELX_module,  ONLY : fmt_SFAC, fmt_UNIT
 USE macros_module, ONLY : u_case
 USE IO_module

 implicit none
  CHARACTER (LEN=*), INTENT(IN)        :: input_file, extension
  ! local variables
  INTEGER                              :: i, i1, n_atom
  INTEGER                              :: long
  CHARACTER(LEN=256)                   :: created_FST_file
  CHARACTER(len=4)                     :: atom_sp
  LOGICAL                              :: atom_Al, atom_Au, atom_B, atom_Ba, atom_Br, atom_C, atom_Ca, atom_Ce, atom_Cl, atom_Co, &
                                          atom_Cr, atom_Cs, atom_Cu, atom_F, atom_Fe, atom_H, &
                                          atom_Hf, atom_Hg, atom_I, atom_In, atom_Ir, atom_K, atom_Li,                            &
                                          atom_Mg, atom_Mn, atom_Mo, atom_N,  &
                                          atom_Na, atom_Nd, atom_Ni, &
                                          atom_O, atom_P, atom_Pd, atom_Pt, &
                                          atom_Re, atom_Ru, atom_S, atom_Sb, atom_Se, atom_Si, atom_Sn, atom_Sr, atom_Ti, atom_U, &
                                          atom_V, atom_Y,  atom_Yb, atom_Zn, atom_Zr

  character (len=16)                    :: conn_radius_string, conn_H_radius_string, atom_radius_string
  character (len=16)                    :: H_radius_string, C_radius_string, N_radius_string, O_radius_string
  real                                  :: min_x, min_y, min_z
  real                                  :: max_x, max_y, max_z

  if(debug_proc%level_2)  call write_debug_proc_level(2, "create_fst_file")


  conn_radius_string   = 'RADIUS 1.'
  conn_H_radius_string = 'RADIUS 0.5'
  atom_radius_string   = 'RADIUS 1.'
  H_radius_string      = 'RADIUS 0.3'
  C_radius_string      = 'RADIUS 0.8'
  N_radius_string      = 'RADIUS 0.8'
  O_radius_string      = 'RADIUS 0.8'
  if(create_FST_poly) then
   conn_radius_string  = 'RADIUS 0.'
   atom_radius_string  = 'RADIUS 0.'
  end if


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


  i1= INDEX(input_file, '.')
  if(i1 /=0) then
   WRITE(created_FST_file, '(a)') input_file(1:i1-1)//'_'//TRIM(extension)//'.FST'
  else
   WRITE(created_FST_file, '(a)') 'cryscalc.FST'
  endif

  OPEN(UNIT=4, FILE=TRIM(created_FST_file))
   WRITE(4, '(a)')                         '!  File for FullProf Studio (created by CRYSCALC)'
   WRITE(4, '(2a)')                        'TITLE   ', TRIM(main_title)
   if(create_FST_MOLE) then
    WRITE(4, '(a)')                         '! next line is commented : only the asymetric unit cell is drawn'
    WRITE(4, '(2a)')                        '!SPACEG  ', TRIM(space_group_symbol)
   else
    WRITE(4, '(2a)')                        'SPACEG  ', TRIM(space_group_symbol)
   end if
   WRITE(4, '(a,3(1x,F8.4),3(1x,F8.3),a)') 'CELL  ', unit_cell%param(1:3), unit_cell%param(4:6), '   DISPLAY MULTIPLE'

   if(create_FST_MOLE) then
    min_x  = MINVAL(atom_coord(1,1:nb_atom))
    min_y  = MINVAL(atom_coord(2,1:nb_atom))
    min_z  = MINVAL(atom_coord(3,1:nb_atom))
    max_x  = MAXVAL(atom_coord(1,1:nb_atom))
    max_y  = MAXVAL(atom_coord(2,1:nb_atom))
    max_z  = MAXVAL(atom_coord(3,1:nb_atom))
    write(4, '(a,2x,6(F6.2,2x))') 'BOX', 0.99*min_x, 1.01*max_x, 0.99*min_y, 1.01*max_y, 0.99*min_z, 1.01*max_z
   else
    WRITE(4, '(a)')                         'BOX  -0.15  1.15  -0.15 1.15  -0.15 1.15'
   end if

  atom_Al = .false.
  atom_Au = .false.
  atom_B  = .false.
  atom_Ba = .false.
  atom_Br = .false.
  atom_C  = .false.
  atom_Ca = .false.
  atom_Ce = .false.
  atom_Cl = .false.
  atom_Co = .false.
  atom_Cr = .false.
  atom_Cs = .false.
  atom_Cu = .false.
  atom_F  = .false.
  atom_Fe = .false.
  atom_H  = .false.
  atom_Hf = .false.
  atom_Hg = .false.
  atom_I  = .false.
  atom_In = .false.
  atom_Ir = .false.
  atom_K  = .false.
  atom_Li = .false.
  atom_Mg = .false.
  atom_Mn = .false.
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
  atom_Sb = .false.
  atom_Se = .false.
  atom_Si = .false.
  atom_Sn = .false.
  atom_Sr = .false.
  atom_Ti = .false.
  atom_U  = .false.
  atom_V  = .false.
  atom_Yb = .false.
  atom_Y  = .false.
  atom_Zn = .false.
  atom_Zr = .false.

   do n_atom=1, nb_atom
    call get_specie_from_type(atom_typ(n_atom), atom_sp)

    long = len_trim(atom_typ(n_atom))
    long = len_trim(atom_sp)

    if(long == 2 .and. trim(u_case(atom_sp(1:long))) == 'AL') atom_Al =.true.
    if(long == 2 .and. trim(u_case(atom_sp(1:long))) == 'AU') atom_Au =.true.
    if(long == 1 .and. trim(u_case(atom_sp(1:long))) == 'B')  atom_B  =.true.
    if(long == 2 .and. trim(u_case(atom_sp(1:long))) == 'BA') atom_Ba =.true.
    if(long == 2 .and. trim(u_case(atom_sp(1:long))) == 'BR') atom_Br =.true.
    if(long == 1 .and. trim(u_case(atom_sp(1:long))) == 'C')  atom_C  =.true.
    if(long == 2 .and. trim(u_case(atom_sp(1:long))) == 'CA') atom_Ca =.true.
    if(long == 2 .and. trim(u_case(atom_sp(1:long))) == 'CE') atom_Ce =.true.

    if(long == 2 .and. trim(u_case(atom_sp(1:long))) == 'CL') atom_Cl =.true.
    if(long == 2 .and. trim(u_case(atom_sp(1:long))) == 'CO') atom_Co =.true.
    if(long == 2 .and. trim(u_case(atom_sp(1:long))) == 'CR') atom_Cr =.true.

    if(long == 2 .and. trim(u_case(atom_sp(1:long))) == 'CS') atom_Cs =.true.

    if(long == 2 .and. trim(u_case(atom_sp(1:long))) == 'CU') atom_Cu =.true.
    if(long == 1 .and. trim(u_case(atom_sp(1:long))) == 'F')  atom_F  =.true.

    if(long == 2 .and. trim(u_case(atom_sp(1:long))) == 'FE') atom_Fe =.true.
    if(long == 1 .and. trim(u_case(atom_sp(1:long))) == 'H')  atom_H  =.true.

    if(long == 2 .and. trim(u_case(atom_sp(1:long))) == 'HF') atom_Hf =.true.
    if(long == 2 .and. trim(u_case(atom_sp(1:long))) == 'HG') atom_Hg =.true.
    if(long == 1 .and. trim(u_case(atom_sp(1:long))) == 'I')  atom_I  =.true.
    if(long == 2 .and. trim(u_case(atom_sp(1:long))) == 'IN') atom_In =.true.
    if(long == 2 .and. trim(u_case(atom_sp(1:long))) == 'IR') atom_Ir =.true.
    if(long == 1 .and. trim(u_case(atom_sp(1:long))) == 'K')  atom_K  =.true.

    if(long == 2 .and. trim(u_case(atom_sp(1:long))) == 'LI') atom_Li =.true.
    if(long == 2 .and. trim(u_case(atom_sp(1:long))) == 'MO') atom_Mo =.true.
    if(long == 2 .and. trim(u_case(atom_sp(1:long))) == 'MG') atom_Mg =.true.
    if(long == 2 .and. trim(u_case(atom_sp(1:long))) == 'MN') atom_Mn =.true.
    if(long == 1 .and. trim(u_case(atom_sp(1:long))) == 'N')  atom_N  =.true.
    if(long == 2 .and. trim(u_case(atom_sp(1:long))) == 'NA') atom_Na =.true.
    if(long == 2 .and. trim(u_case(atom_sp(1:long))) == 'ND') atom_Nd =.true.
    if(long == 2 .and. trim(u_case(atom_sp(1:long))) == 'NI') atom_Ni =.true.

    if(long == 1 .and. trim(u_case(atom_sp(1:long))) == 'O')  atom_O  =.true.
    if(long == 1 .and. trim(u_case(atom_sp(1:long))) == 'P')  atom_P  =.true.
    if(long == 2 .and. trim(u_case(atom_sp(1:long))) == 'PD') atom_Pd =.true.
    if(long == 2 .and. trim(u_case(atom_sp(1:long))) == 'PT') atom_Pt =.true.
    if(long == 2 .and. trim(u_case(atom_sp(1:long))) == 'RE') atom_Re =.true.
    if(long == 2 .and. trim(u_case(atom_sp(1:long))) == 'RU') atom_Ru =.true.
    if(long == 1 .and. trim(u_case(atom_sp(1:long))) == 'S')  atom_S  =.true.
    if(long == 2 .and. trim(u_case(atom_sp(1:long))) == 'SB') atom_Sb =.true.
    if(long == 2 .and. trim(u_case(atom_sp(1:long))) == 'SE') atom_Se =.true.
    if(long == 2 .and. trim(u_case(atom_sp(1:long))) == 'SI') atom_Si =.true.
    if(long == 2 .and. trim(u_case(atom_sp(1:long))) == 'SN') atom_Sn =.true.
    if(long == 2 .and. trim(u_case(atom_sp(1:long))) == 'SR') atom_Sr =.true.
    if(long == 2 .and. trim(u_case(atom_sp(1:long))) == 'TI') atom_Ti =.true.
    if(long == 1 .and. trim(u_case(atom_sp(1:long))) == 'U')  atom_U  =.true.
    if(long == 1 .and. trim(u_case(atom_sp(1:long))) == 'V')  atom_V  =.true.
    if(long == 1 .and. trim(u_case(atom_sp(1:long))) == 'Y')  atom_Y  =.true.
    if(long == 1 .and. trim(u_case(atom_sp(1:long))) == 'Y')  atom_Y  =.true.

    if(long == 2 .and. trim(u_case(atom_sp(1:long))) == 'YB') atom_Yb =.true.

    if(long == 2 .and. trim(u_case(atom_sp(1:long))) == 'ZN') atom_Zn =.true.
    if(long == 2 .and. trim(u_case(atom_sp(1:long))) == 'ZR') atom_Zr =.true.

    if(u_case(atom_label(n_atom)(1:1)) == 'Q') cycle   ! exclude Q peaks

    if(long == 1 .and. trim(u_case(atom_sp(1:long))) == 'H') then
    if(FST_no_H) cycle
     WRITE(4,'(a,2a,3(1x,F8.5), 2a)') 'ATOM ', u_case(atom_label(n_atom)(1:6)), u_case(atom_sp(1:4)),  &
                                  (atom_coord(i,n_atom),i=1,3), '  ', trim(H_radius_string)
    elseif(long == 1 .and. trim(u_case(atom_sp(1:long))) == 'C') then
     WRITE(4,'(a,2a,3(1x,F8.5), 2a)') 'ATOM ', u_case(atom_label(n_atom)(1:6)), u_case(atom_sp(1:4)),  &
                                  (atom_coord(i,n_atom),i=1,3), '  ', trim(C_radius_string)
    elseif(long == 1 .and. trim(u_case(atom_sp(1:long))) == 'N') then
     WRITE(4,'(a,2a,3(1x,F8.5), 2a)') 'ATOM ', u_case(atom_label(n_atom)(1:6)), u_case(atom_sp(1:4)),  &
                                  (atom_coord(i,n_atom),i=1,3), '  ', trim(N_radius_string)
    elseif(long == 1 .and. trim(u_case(atom_sp(1:long))) == 'O') then
     WRITE(4,'(a,2a,3(1x,F8.5), 2a)') 'ATOM ', u_case(atom_label(n_atom)(1:6)), u_case(atom_sp(1:4)),  &
                                  (atom_coord(i,n_atom),i=1,3), '  ', trim(O_radius_string)

    else
     WRITE(4,'(a,2a,3(1x,F8.5), 2a)') 'ATOM ', u_case(atom_label(n_atom)(1:6)), u_case(atom_sp(1:4)),  &
                                  (atom_coord(i,n_atom),i=1,3), '  ', trim(atom_radius_string)
    end if
   END do


   if(atom_C)                WRITE(4, '(2a)')   'CONN C   C   0.   1.7 ', trim(conn_radius_string)
   if(atom_C .and. atom_B)   WRITE(4, '(2a)')   'CONN C   B   0.   1.7 ', trim(conn_radius_string)
   if(atom_C .and. atom_Br)  WRITE(4, '(2a)')   'CONN C   BR  0.   2.  ', trim(conn_radius_string)
   if(atom_C .and. atom_F)   WRITE(4, '(2a)')   'CONN C   F   0.   2.  ', trim(conn_radius_string)
   if(atom_C .and. atom_N)   WRITE(4, '(2a)')   'CONN C   N   0.   1.7 ', trim(conn_radius_string)
   if(atom_C .and. atom_O)   WRITE(4, '(2a)')   'CONN C   O   0.   1.7 ', trim(conn_radius_string)
   if(atom_C .and. atom_P)   WRITE(4, '(2a)')   'CONN C   P   0.   2.  ', trim(conn_radius_string)
   if(atom_C .and. atom_S)   WRITE(4, '(2a)')   'CONN C   S   0.   2.  ', trim(conn_radius_string)
   if(atom_C .and. atom_Se)  WRITE(4, '(2a)')   'CONN C   SE  0.   2.  ', trim(conn_radius_string)
   if(atom_C .and. atom_H .and. .not. FST_no_H)   WRITE(4, '(2a)')   'CONN C   H   0.   1.4 ', trim(conn_H_radius_string)


   if(atom_B)                WRITE(4, '(2a)')   'CONN B   B   0.   2.  ', trim(conn_radius_string)
   if(atom_B .and. atom_F)   WRITE(4, '(2a)')   'CONN B   F   0.   2.  ', trim(conn_radius_string)
   if(atom_B .and. atom_N)   WRITE(4, '(2a)')   'CONN B   N   0.   2.  ', trim(conn_radius_string)
   if(atom_B .and. atom_O)   WRITE(4, '(2a)')   'CONN B   O   0.   2.  ', trim(conn_radius_string)

   if(atom_N)                WRITE(4, '(2a)')   'CONN N   N   0.   1.7 ', trim(conn_radius_string)
   if(atom_N .and. atom_O)   WRITE(4, '(2a)')   'CONN N   O   0.   1.7 ', trim(conn_radius_string)
   if(atom_N .and. atom_H .and. .not. FST_no_H)   WRITE(4, '(2a)')   'CONN N   H   0.   1.4 ', trim(conn_H_radius_string)

   if(atom_Al .and. atom_C)  WRITE(4, '(2a)')   'CONN AL  C   0.   2.6 ', trim(conn_radius_string)
   if(atom_Al .and. atom_Cl) WRITE(4, '(2a)')   'CONN AL  CL  0.   2.6 ', trim(conn_radius_string)
   if(atom_Al .and. atom_O)  WRITE(4, '(2a)')   'CONN AL  O   0.   2.6 ', trim(conn_radius_string)
   if(atom_Al .and. atom_P)  WRITE(4, '(2a)')   'CONN AL  P   0.   2.6 ', trim(conn_radius_string)
   if(atom_Al .and. atom_S)  WRITE(4, '(2a)')   'CONN AL  S   0.   2.6 ', trim(conn_radius_string)
   if(atom_Al .and. atom_Se) WRITE(4, '(2a)')   'CONN AL  SE  0.   2.6 ', trim(conn_radius_string)




   if(atom_Au .and. atom_C)  WRITE(4, '(2a)')   'CONN AU  C   0.   2.6 ', trim(conn_radius_string)
   if(atom_Au .and. atom_Cl) WRITE(4, '(2a)')   'CONN AU  CL  0.   2.6 ', trim(conn_radius_string)
   if(atom_Au .and. atom_I)  WRITE(4, '(2a)')   'CONN AU  I   0.   2.6 ', trim(conn_radius_string)
   if(atom_Au .and. atom_O)  WRITE(4, '(2a)')   'CONN AU  O   0.   2.6 ', trim(conn_radius_string)
   if(atom_Au .and. atom_S)  WRITE(4, '(2a)')   'CONN AU  S   0.   2.6 ', trim(conn_radius_string)
   if(atom_Au .and. atom_Se) WRITE(4, '(2a)')   'CONN AU  SE  0.   2.6 ', trim(conn_radius_string)

   if(atom_Ba .and. atom_N)  WRITE(4, '(2a)')   'CONN BA  N   0.   2.9 ', trim(conn_radius_string)
   if(atom_Ba .and. atom_O)  WRITE(4, '(2a)')   'CONN BA  O   0.   2.9 ', trim(conn_radius_string)

   if(atom_Ca .and. atom_N)  WRITE(4, '(2a)')   'CONN CA  N   0.   2.9 ', trim(conn_radius_string)
   if(atom_Ca .and. atom_O)  WRITE(4, '(2a)')   'CONN CA  O   0.   2.9 ', trim(conn_radius_string)

   if(atom_Ce .and. atom_N)  WRITE(4, '(2a)')   'CONN CE  N   0.   2.9 ', trim(conn_radius_string)
   if(atom_Ce .and. atom_O)  WRITE(4, '(2a)')   'CONN CE  O   0.   2.9 ', trim(conn_radius_string)

   if(atom_Cl .and. atom_O)  WRITE(4, '(2a)')   'CONN CL  O   0.   2.0 ', trim(conn_radius_string)

   if(atom_Co .and. atom_C)  WRITE(4, '(2a)')   'CONN CO  C   0.   2.6 ', trim(conn_radius_string)
   if(atom_Co .and. atom_Cl) WRITE(4, '(2a)')   'CONN CO  CL  0.   2.6 ', trim(conn_radius_string)
   if(atom_Co .and. atom_F)  WRITE(4, '(2a)')   'CONN CO  F   0.   2.6 ', trim(conn_radius_string)
   if(atom_Co .and. atom_N)  WRITE(4, '(2a)')   'CONN CO  N   0.   2.6 ', trim(conn_radius_string)
   if(atom_Co .and. atom_O)  WRITE(4, '(2a)')   'CONN CO  O   0.   2.6 ', trim(conn_radius_string)
   if(atom_Co .and. atom_P)  WRITE(4, '(2a)')   'CONN CO  P   0.   2.6 ', trim(conn_radius_string)
   if(atom_Co .and. atom_Se) WRITE(4, '(2a)')   'CONN CO  SE  0.   2.6 ', trim(conn_radius_string)

   if(atom_Cr .and. atom_C)  WRITE(4, '(2a)')   'CONN CR  C   0.   2.6 ', trim(conn_radius_string)
   if(atom_Cr .and. atom_Cl) WRITE(4, '(2a)')   'CONN CR  CL  0.   2.6 ', trim(conn_radius_string)
   if(atom_Cr .and. atom_N)  WRITE(4, '(2a)')   'CONN CR  N   0.   2.6 ', trim(conn_radius_string)
   if(atom_Cr .and. atom_O)  WRITE(4, '(2a)')   'CONN CR  O   0.   2.6 ', trim(conn_radius_string)
   if(atom_Cr .and. atom_P)  WRITE(4, '(2a)')   'CONN CR  P   0.   2.6 ', trim(conn_radius_string)
   if(atom_Cr .and. atom_Se) WRITE(4, '(2a)')   'CONN CR  SE  0.   2.6 ', trim(conn_radius_string)



   if(atom_Cu .and. atom_C)  WRITE(4, '(2a)')   'CONN CU  C   0.   2.6 ', trim(conn_radius_string)
   if(atom_Cu .and. atom_Cl) WRITE(4, '(2a)')   'CONN CU  CL  0.   2.6 ', trim(conn_radius_string)
   if(atom_Cu .and. atom_F)  WRITE(4, '(2a)')   'CONN CU  F   0.   2.6 ', trim(conn_radius_string)
   if(atom_Cu .and. atom_N)  WRITE(4, '(2a)')   'CONN CU  N   0.   2.6 ', trim(conn_radius_string)
   if(atom_Cu .and. atom_O)  WRITE(4, '(2a)')   'CONN CU  O   0.   2.6 ', trim(conn_radius_string)
   if(atom_Cu .and. atom_P)  WRITE(4, '(2a)')   'CONN CU  P   0.   2.6 ', trim(conn_radius_string)
   if(atom_Cu .and. atom_Se) WRITE(4, '(2a)')   'CONN CU  SE  0.   2.6 ', trim(conn_radius_string)

   if(atom_Cs .and. atom_N)  WRITE(4, '(2a)')   'CONN CS  N   0.   2.9 ', trim(conn_radius_string)
   if(atom_Cs .and. atom_O)  WRITE(4, '(2a)')   'CONN CS  O   0.   2.9 ', trim(conn_radius_string)


   if(atom_Fe .and. atom_C)  WRITE(4, '(2a)')   'CONN FE  C   0.   2.6 ', trim(conn_radius_string)
   if(atom_Fe .and. atom_Cl) WRITE(4, '(2a)')   'CONN FE  CL  0.   2.6 ', trim(conn_radius_string)
   if(atom_Fe .and. atom_N)  WRITE(4, '(2a)')   'CONN FE  N   0.   2.6 ', trim(conn_radius_string)
   if(atom_Fe .and. atom_O)  WRITE(4, '(2a)')   'CONN FE  O   0.   2.6 ', trim(conn_radius_string)
   if(atom_Fe .and. atom_P)  WRITE(4, '(2a)')   'CONN FE  P   0.   2.6 ', trim(conn_radius_string)
   if(atom_Fe .and. atom_Se) WRITE(4, '(2a)')   'CONN FE  SE  0.   2.6 ', trim(conn_radius_string)

   if(atom_Hf .and. atom_Cl) WRITE(4, '(2a)')   'CONN HF  CL  0.   2.9 ', trim(conn_radius_string)
   if(atom_Hf .and. atom_N)  WRITE(4, '(2a)')   'CONN HF  N   0.   2.9 ', trim(conn_radius_string)
   if(atom_Hf .and. atom_O)  WRITE(4, '(2a)')   'CONN HF  O   0.   2.9 ', trim(conn_radius_string)
   if(atom_Hf .and. atom_P)  WRITE(4, '(2a)')   'CONN HF  P   0.   2.6 ', trim(conn_radius_string)
   if(atom_Hf .and. atom_S)  WRITE(4, '(2a)')   'CONN HF  S   0.   2.6 ', trim(conn_radius_string)
   if(atom_Hf .and. atom_Se) WRITE(4, '(2a)')   'CONN HF  SE  0.   2.6 ', trim(conn_radius_string)

   if(atom_Hg)               WRITE(4, '(2a)')   'CONN HG  HG  0.   3.0 ', trim(conn_radius_string)
   if(atom_Hg .and. atom_Sb) WRITE(4, '(2a)')   'CONN HG  SB  0.   3.0 ', trim(conn_radius_string)
   if(atom_Hg .and. atom_I)  WRITE(4, '(2a)')   'CONN HG  I   0.   3.5 ', trim(conn_radius_string)


   if(atom_I)                WRITE(4, '(2a)')   'CONN I   I   0.   2.6 ', trim(conn_radius_string)
   if(atom_I  .and. atom_C)  WRITE(4, '(2a)')   'CONN I   C   0.   2.6 ', trim(conn_radius_string)
   if(atom_I  .and. atom_Hg) WRITE(4, '(2a)')   'CONN I   HG  0.   3.4 ', trim(conn_radius_string)
   if(atom_I  .and. atom_Sb) WRITE(4, '(2a)')   'CONN I   SB  0.   3.4 ', trim(conn_radius_string)


   if(atom_In .and. atom_C)  WRITE(4, '(2a)')   'CONN IN  C   0.   2.6 ', trim(conn_radius_string)
   if(atom_In .and. atom_Cl) WRITE(4, '(2a)')   'CONN IN  CL  0.   2.6 ', trim(conn_radius_string)
   if(atom_In .and. atom_O)  WRITE(4, '(2a)')   'CONN IN  O   0.   2.6 ', trim(conn_radius_string)
   if(atom_In .and. atom_P)  WRITE(4, '(2a)')   'CONN IN  P   0.   2.6 ', trim(conn_radius_string)
   if(atom_In .and. atom_S)  WRITE(4, '(2a)')   'CONN IN  S   0.   2.6 ', trim(conn_radius_string)
   if(atom_In .and. atom_Se) WRITE(4, '(2a)')   'CONN IN  SE  0.   2.6 ', trim(conn_radius_string)

   if(atom_Ir .and. atom_C)  WRITE(4, '(2a)')   'CONN IR  C   0.   2.6 ', trim(conn_radius_string)
   if(atom_Ir .and. atom_Cl) WRITE(4, '(2a)')   'CONN IR  CL  0.   2.6 ', trim(conn_radius_string)
   if(atom_Ir .and. atom_O)  WRITE(4, '(2a)')   'CONN IR  O   0.   2.6 ', trim(conn_radius_string)
   if(atom_Ir .and. atom_P)  WRITE(4, '(2a)')   'CONN IR  P   0.   2.6 ', trim(conn_radius_string)
   if(atom_Ir .and. atom_S)  WRITE(4, '(2a)')   'CONN IR  S   0.   2.6 ', trim(conn_radius_string)
   if(atom_Ir .and. atom_Se) WRITE(4, '(2a)')   'CONN IR  SE  0.   2.6 ', trim(conn_radius_string)

   if(atom_K  .and. atom_C)  WRITE(4, '(2a)')   'CONN K   C   0.   2.9 ', trim(conn_radius_string)
   if(atom_K  .and. atom_N)  WRITE(4, '(2a)')   'CONN K   N   0.   2.9 ', trim(conn_radius_string)
   if(atom_K  .and. atom_O)  WRITE(4, '(2a)')   'CONN K   O   0.   2.9 ', trim(conn_radius_string)

   if(atom_Li .and. atom_N)  WRITE(4, '(2a)')   'CONN LI  N   0.   2.5 ', trim(conn_radius_string)
   if(atom_Li .and. atom_O)  WRITE(4, '(2a)')   'CONN LI  O   0.   2.5 ', trim(conn_radius_string)

   if(atom_Mg .and. atom_C)  WRITE(4, '(2a)')   'CONN MG  C   0.   2.6 ', trim(conn_radius_string)
   if(atom_Mg .and. atom_Cl) WRITE(4, '(2a)')   'CONN MG  CL  0.   2.6 ', trim(conn_radius_string)
   if(atom_Mg .and. atom_O)  WRITE(4, '(2a)')   'CONN MG  O   0.   2.6 ', trim(conn_radius_string)
   if(atom_Mg .and. atom_P)  WRITE(4, '(2a)')   'CONN MG  P   0.   2.6 ', trim(conn_radius_string)
   if(atom_Mg .and. atom_S)  WRITE(4, '(2a)')   'CONN MG  S   0.   2.6 ', trim(conn_radius_string)
   if(atom_Mg .and. atom_Se) WRITE(4, '(2a)')   'CONN MG  SE  0.   2.6 ', trim(conn_radius_string)

   if(atom_Mn .and. atom_C)  WRITE(4, '(2a)')   'CONN MN  C   0.   2.6 ', trim(conn_radius_string)
   if(atom_Mn .and. atom_Cl) WRITE(4, '(2a)')   'CONN MN  CL  0.   2.6 ', trim(conn_radius_string)
   if(atom_Mn .and. atom_O)  WRITE(4, '(2a)')   'CONN MN  O   0.   2.6 ', trim(conn_radius_string)
   if(atom_Mn .and. atom_P)  WRITE(4, '(2a)')   'CONN MN  P   0.   2.6 ', trim(conn_radius_string)
   if(atom_Mn .and. atom_S)  WRITE(4, '(2a)')   'CONN MN  S   0.   2.6 ', trim(conn_radius_string)
   if(atom_Mn .and. atom_Se) WRITE(4, '(2a)')   'CONN MN  SE  0.   2.6 ', trim(conn_radius_string)

   if(atom_Mo .and. atom_C)  WRITE(4, '(2a)')   'CONN MO  C   0.   2.6 ', trim(conn_radius_string)
   if(atom_Mo .and. atom_Cl) WRITE(4, '(2a)')   'CONN MO  CL  0.   2.6 ', trim(conn_radius_string)
   if(atom_Mo .and. atom_F)  WRITE(4, '(2a)')   'CONN MO  F   0.   2.6 ', trim(conn_radius_string)
   if(atom_Mo .and. atom_N)  WRITE(4, '(2a)')   'CONN MO  N   0.   2.6 ', trim(conn_radius_string)
   if(atom_Mo .and. atom_O)  WRITE(4, '(2a)')   'CONN MO  O   0.   2.6 ', trim(conn_radius_string)
   if(atom_Mo .and. atom_P)  WRITE(4, '(2a)')   'CONN MO  P   0.   2.6 ', trim(conn_radius_string)
   if(atom_Mo .and. atom_Se) WRITE(4, '(2a)')   'CONN MO  SE  0.   2.6 ', trim(conn_radius_string)

   if(atom_Na .and. atom_C)  WRITE(4, '(2a)')   'CONN NA  C   0.   2.6 ', trim(conn_radius_string)
   if(atom_Na .and. atom_Cl) WRITE(4, '(2a)')   'CONN NA  CL  0.   2.6 ', trim(conn_radius_string)
   if(atom_Na .and. atom_O)  WRITE(4, '(2a)')   'CONN NA  O   0.   2.6 ', trim(conn_radius_string)
   if(atom_Na .and. atom_P)  WRITE(4, '(2a)')   'CONN NA  P   0.   2.6 ', trim(conn_radius_string)
   if(atom_Na .and. atom_S)  WRITE(4, '(2a)')   'CONN NA  S   0.   2.6 ', trim(conn_radius_string)
   if(atom_Na .and. atom_Se) WRITE(4, '(2a)')   'CONN NA  SE  0.   2.6 ', trim(conn_radius_string)


   if(atom_ND .and. atom_Cl) WRITE(4, '(2a)')   'CONN ND  CL  0.   2.9 ', trim(conn_radius_string)
   if(atom_ND .and. atom_N)  WRITE(4, '(2a)')   'CONN ND  N   0.   2.9 ', trim(conn_radius_string)
   if(atom_ND .and. atom_O)  WRITE(4, '(2a)')   'CONN ND  O   0.   2.9 ', trim(conn_radius_string)
   if(atom_ND .and. atom_P)  WRITE(4, '(2a)')   'CONN ND  P   0.   2.6 ', trim(conn_radius_string)
   if(atom_ND .and. atom_S)  WRITE(4, '(2a)')   'CONN ND  S   0.   2.6 ', trim(conn_radius_string)
   if(atom_ND .and. atom_Se) WRITE(4, '(2a)')   'CONN ND  SE  0.   2.6 ', trim(conn_radius_string)

   if(atom_Ni .and. atom_C)  WRITE(4, '(2a)')   'CONN NI  C   0.   2.6 ', trim(conn_radius_string)
   if(atom_Ni .and. atom_Cl) WRITE(4, '(2a)')   'CONN NI  CL  0.   2.6 ', trim(conn_radius_string)
   if(atom_Ni .and. atom_F)  WRITE(4, '(2a)')   'CONN NI  F   0.   2.6 ', trim(conn_radius_string)
   if(atom_Ni .and. atom_N)  WRITE(4, '(2a)')   'CONN NI  N   0.   2.6 ', trim(conn_radius_string)
   if(atom_Ni .and. atom_O)  WRITE(4, '(2a)')   'CONN NI  O   0.   2.6 ', trim(conn_radius_string)
   if(atom_Ni .and. atom_P)  WRITE(4, '(2a)')   'CONN NI  P   0.   2.6 ', trim(conn_radius_string)
   if(atom_Ni .and. atom_S)  WRITE(4, '(2a)')   'CONN NI  S   0.   2.6 ', trim(conn_radius_string)
   if(atom_Ni .and. atom_Se) WRITE(4, '(2a)')   'CONN NI  SE  0.   2.6 ', trim(conn_radius_string)


   if(atom_O  .and. atom_H .and. .not. FST_no_H)  WRITE(4, '(2a)')   'CONN O   H   0.   1.4 ', trim(conn_H_radius_string)

   if(atom_P  .and. atom_F)  WRITE(4, '(2a)')   'CONN P   F   0.   2.0 ', trim(conn_radius_string)
   if(atom_P  .and. atom_O)  WRITE(4, '(2a)')   'CONN P   O   0.   2.0 ', trim(conn_radius_string)

   if(atom_Pd .and. atom_C)  WRITE(4, '(2a)')   'CONN PD  C   0.   2.6 ', trim(conn_radius_string)
   if(atom_Pd .and. atom_Cl) WRITE(4, '(2a)')   'CONN PD  CL  0.   2.6 ', trim(conn_radius_string)
   if(atom_Pd .and. atom_N)  WRITE(4, '(2a)')   'CONN PD  N   0.   2.6 ', trim(conn_radius_string)
   if(atom_Pd .and. atom_P)  WRITE(4, '(2a)')   'CONN PD  P   0.   2.6 ', trim(conn_radius_string)
   if(atom_Pd .and. atom_S)  WRITE(4, '(2a)')   'CONN PD  S   0.   2.6 ', trim(conn_radius_string)
   if(atom_Pd .and. atom_Se) WRITE(4, '(2a)')   'CONN PD  SE  0.   2.6 ', trim(conn_radius_string)

   if(atom_Pt .and. atom_C)  WRITE(4, '(2a)')   'CONN PT  C   0.   2.6 ', trim(conn_radius_string)
   if(atom_Pt .and. atom_Cl) WRITE(4, '(2a)')   'CONN PT  CL  0.   2.6 ', trim(conn_radius_string)
   if(atom_Pt .and. atom_N)  WRITE(4, '(2a)')   'CONN PT  N   0.   2.6 ', trim(conn_radius_string)
   if(atom_Pt .and. atom_O)  WRITE(4, '(2a)')   'CONN PT  O   0.   2.6 ', trim(conn_radius_string)
   if(atom_Pt .and. atom_P)  WRITE(4, '(2a)')   'CONN PT  P   0.   2.6 ', trim(conn_radius_string)
   if(atom_Pt .and. atom_S)  WRITE(4, '(2a)')   'CONN PT  S   0.   2.6 ', trim(conn_radius_string)
   if(atom_Pt .and. atom_Se) WRITE(4, '(2a)')   'CONN PT  SE  0.   2.6 ', trim(conn_radius_string)

   if(atom_Re .and. atom_C)  WRITE(4, '(2a)')   'CONN RE  C   0.   2.6 ', trim(conn_radius_string)
   if(atom_Re .and. atom_Cl) WRITE(4, '(2a)')   'CONN RE  CL  0.   2.6 ', trim(conn_radius_string)
   if(atom_Re .and. atom_N)  WRITE(4, '(2a)')   'CONN RE  N   0.   2.6 ', trim(conn_radius_string)
   if(atom_Re .and. atom_O)  WRITE(4, '(2a)')   'CONN RE  N   0.   2.6 ', trim(conn_radius_string)
   if(atom_Re .and. atom_P)  WRITE(4, '(2a)')   'CONN RE  P   0.   2.6 ', trim(conn_radius_string)
   if(atom_Re .and. atom_S)  WRITE(4, '(2a)')   'CONN RE  S   0.   2.6 ', trim(conn_radius_string)
   if(atom_Re .and. atom_Se) WRITE(4, '(2a)')   'CONN RE  SE  0.   2.6 ', trim(conn_radius_string)

   if(atom_Ru .and. atom_C)  WRITE(4, '(2a)')   'CONN RU  C   0.   2.6 ', trim(conn_radius_string)
   if(atom_Ru .and. atom_Cl) WRITE(4, '(2a)')   'CONN RU  CL  0.   2.6 ', trim(conn_radius_string)
   if(atom_Ru .and. atom_N)  WRITE(4, '(2a)')   'CONN RU  N   0.   2.6 ', trim(conn_radius_string)
   if(atom_Ru .and. atom_O)  WRITE(4, '(2a)')   'CONN RU  O   0.   2.6 ', trim(conn_radius_string)
   if(atom_Ru .and. atom_P)  WRITE(4, '(2a)')   'CONN RU  P   0.   2.6 ', trim(conn_radius_string)
   if(atom_Ru .and. atom_S)  WRITE(4, '(2a)')   'CONN RU  S   0.   2.6 ', trim(conn_radius_string)
   if(atom_Ru .and. atom_Se) WRITE(4, '(2a)')   'CONN RU  SE  0.   2.6 ', trim(conn_radius_string)

   if(atom_S  .and. atom_F)  WRITE(4, '(2a)')   'CONN S   F   0.   2.0 ', trim(conn_radius_string)
   if(atom_S  .and. atom_O)  WRITE(4, '(2a)')   'CONN S   O   0.   2.0 ', trim(conn_radius_string)

   if(atom_Sb)               WRITE(4, '(2a)')   'CONN SB  SB  0.   3.0 ', trim(conn_radius_string)
   if(atom_Se .and. atom_O)  WRITE(4, '(2a)')   'CONN SE  O   0.   2.0 ', trim(conn_radius_string)

   if(atom_Si)               WRITE(4, '(2a)')   'CONN SI  SI  0.   2.4 ', trim(conn_radius_string)
   if(atom_Si .and. atom_C)  WRITE(4, '(2a)')   'CONN SI  C   0.   2.  ', trim(conn_radius_string)
   if(atom_Si .and. atom_H .and. .not. FST_no_H)  WRITE(4, '(2a)')   'CONN SI  H   0.   2.  ', trim(conn_H_radius_string)
   if(atom_Si .and. atom_N)  WRITE(4, '(2a)')   'CONN SI  N   0.   2.4 ', trim(conn_radius_string)
   if(atom_Si .and. atom_O)  WRITE(4, '(2a)')   'CONN SI  O   0.   2.4 ', trim(conn_radius_string)




   if(atom_Sn .and. atom_C)  WRITE(4, '(2a)')   'CONN SN  C  0.   2.  ', trim(conn_radius_string)
   if(atom_Sn .and. atom_Cl) WRITE(4, '(2a)')   'CONN SN  CL 0.   2.  ', trim(conn_radius_string)
   if(atom_Sn .and. atom_F)  WRITE(4, '(2a)')   'CONN SN  F  0.   2.  ', trim(conn_radius_string)
   if(atom_Sn .and. atom_N)  WRITE(4, '(2a)')   'CONN SN  N  0.   2.  ', trim(conn_radius_string)
   if(atom_Sn .and. atom_O)  WRITE(4, '(2a)')   'CONN SN  O  0.   2.  ', trim(conn_radius_string)
   if(atom_Sn .and. atom_S)  WRITE(4, '(2a)')   'CONN SN  S  0.   2.  ', trim(conn_radius_string)
   if(atom_Sn .and. atom_Se) WRITE(4, '(2a)')   'CONN SN  SE 0.   2.  ', trim(conn_radius_string)

   if(atom_Sr .and. atom_N)  WRITE(4, '(2a)')   'CONN SR  N   0.   2.9 ', trim(conn_radius_string)
   if(atom_Sr .and. atom_O)  WRITE(4, '(2a)')   'CONN SR  O   0.   2.9 ', trim(conn_radius_string)
   if(atom_Sr .and. atom_Cl) WRITE(4, '(2a)')   'CONN SR  CL  0.   2.9 ', trim(conn_radius_string)

   if(atom_Ti .and. atom_N)  WRITE(4, '(2a)')   'CONN TI  N   0.   2.6 ', trim(conn_radius_string)
   if(atom_Ti .and. atom_O)  WRITE(4, '(2a)')   'CONN TI  O   0.   2.6 ', trim(conn_radius_string)
   if(atom_Ti .and. atom_Cl) WRITE(4, '(2a)')   'CONN TI  CL  0.   2.6 ', trim(conn_radius_string)

   if(atom_U  .and. atom_O)  WRITE(4, '(2a)')   'CONN U   O   0.   2.9 ', trim(conn_radius_string)
   if(atom_V  .and. atom_O)  WRITE(4, '(2a)')   'CONN V   O   0.   2.9 ', trim(conn_radius_string)

   if(atom_Y  .and. atom_Cl) WRITE(4, '(2a)')   'CONN Y   CL  0.   2.9 ', trim(conn_radius_string)
   if(atom_Y  .and. atom_N)  WRITE(4, '(2a)')   'CONN Y   N   0.   2.9 ', trim(conn_radius_string)
   if(atom_Y  .and. atom_O)  WRITE(4, '(2a)')   'CONN Y   O   0.   2.9 ', trim(conn_radius_string)
   if(atom_Y  .and. atom_P)  WRITE(4, '(2a)')   'CONN Y   P   0.   2.6 ', trim(conn_radius_string)
   if(atom_Y  .and. atom_S)  WRITE(4, '(2a)')   'CONN Y   S   0.   2.6 ', trim(conn_radius_string)
   if(atom_Y  .and. atom_Se) WRITE(4, '(2a)')   'CONN Y   SE  0.   2.6 ', trim(conn_radius_string)

   if(atom_Yb .and. atom_Cl) WRITE(4, '(2a)')   'CONN YB  CL  0.   2.9 ', trim(conn_radius_string)
   if(atom_Yb .and. atom_N)  WRITE(4, '(2a)')   'CONN YB  N   0.   2.9 ', trim(conn_radius_string)
   if(atom_Yb .and. atom_O)  WRITE(4, '(2a)')   'CONN YB  O   0.   2.9 ', trim(conn_radius_string)
   if(atom_Yb .and. atom_P)  WRITE(4, '(2a)')   'CONN YB  P   0.   2.6 ', trim(conn_radius_string)
   if(atom_Yb .and. atom_S)  WRITE(4, '(2a)')   'CONN YB  S   0.   2.6 ', trim(conn_radius_string)
   if(atom_Yb .and. atom_Se) WRITE(4, '(2a)')   'CONN YB  SE  0.   2.6 ', trim(conn_radius_string)

   if(atom_Zn .and. atom_Cl) WRITE(4, '(2a)')   'CONN ZN  CL  0.   2.6 ', trim(conn_radius_string)
   if(atom_Zn .and. atom_N)  WRITE(4, '(2a)')   'CONN ZN  N   0.   2.6 ', trim(conn_radius_string)
   if(atom_Zn .and. atom_O)  WRITE(4, '(2a)')   'CONN ZN  O   0.   2.6 ', trim(conn_radius_string)
   if(atom_Zn .and. atom_P)  WRITE(4, '(2a)')   'CONN ZR  P   0.   2.6 ', trim(conn_radius_string)
   if(atom_Zn .and. atom_S)  WRITE(4, '(2a)')   'CONN ZN  S   0.   2.6 ', trim(conn_radius_string)
   if(atom_Zn .and. atom_Se) WRITE(4, '(2a)')   'CONN ZN  SE  0.   2.6 ', trim(conn_radius_string)

   if(atom_Zr .and. atom_Cl) WRITE(4, '(2a)')   'CONN ZR  CL  0.   2.9 ', trim(conn_radius_string)
   if(atom_Zr .and. atom_N)  WRITE(4, '(2a)')   'CONN ZR  N   0.   2.9 ', trim(conn_radius_string)
   if(atom_Zr .and. atom_O)  WRITE(4, '(2a)')   'CONN ZR  O   0.   2.9 ', trim(conn_radius_string)
   if(atom_Zr .and. atom_P)  WRITE(4, '(2a)')   'CONN ZR  P   0.   2.6 ', trim(conn_radius_string)
   if(atom_Zr .and. atom_S)  WRITE(4, '(2a)')   'CONN ZR  S   0.   2.6 ', trim(conn_radius_string)
   if(atom_Zr .and. atom_Se) WRITE(4, '(2a)')   'CONN ZN  SE  0.   2.6 ', trim(conn_radius_string)


   !WRITE(4, '(a)') '!MOLECULE'
   !WRITE(4, '(a)') '!CONN Ni O 0. 2.3'
   if(create_FST_POLY) then
    if(atom_conn%type(1:1) /= '?') then
     WRITE(4, '(3a)') 'POLY ', trim(u_case(atom_conn%label)), ' COLOR 1 1 0 0.5 EDGES RADIUS 2. EDGECOL 1 1 0'
    else
     WRITE(4, '(a)')  'POLY ??? COLOR 1 1 0 0.5 EDGES RADIUS 2. EDGECOL 1 1 0'
    end if
   else
    WRITE(4, '(a)') '!POLY Ni1 COLOR 1 1 0 0.5 EDGES RADIUS 1. EDGECOL 1 1 0'
   end if

  close (UNIT=4)


  call write_info('')
  write(message_text, '(3a)') '   >> ', trim(created_FST_file), ' (for FP_studio) file has been created.'
  call write_info(trim(message_text))
  call write_info('')

  if(launch_FP_studio) then
   call system ('FP_studio '//trim(created_FST_file))
  end if


 RETURN
end subroutine create_FST_file

!------------------------------------------------------------------------
subroutine create_INS_file(input_file, extension)
 USE cryscalc_module
! USE CFML_IO_FORMATS, ONLY : Write_Shx_Template
 USE SHELX_module
 USE IO_module

 implicit none
  CHARACTER (LEN=*), INTENT(IN)        :: input_file, extension
  ! local variables
  INTEGER                              :: i, long, atom_order
  CHARACTER(LEN=256)                   :: created_INS_file

  if(debug_proc%level_2)  call write_debug_proc_level(2, "create_INS_file")

  !if(len_trim(input_file) >= 12) then
  ! if(input_file(1:12) == 'cryscalc.ins') then
  !  created_INS_file = 'cryscalc.ins'
  ! else
  !  i= INDEX(input_file, '.')
  !  WRITE(created_INS_file, '(a)') input_file(1:i-1)//'_'//TRIM(extension)//'.INS'
  ! end if
  !else
  ! i= INDEX(input_file, '.')
  ! WRITE(created_INS_file, '(a)') input_file(1:i-1)//'_'//TRIM(extension)//'.INS'
  !endif

  i = index(input_file, '.')
  WRITE(created_INS_file, '(a)') input_file(1:i-1)//'_'//TRIM(extension)//'.INS'
  if(len_trim(input_file) >=12) then
   if(input_file(1:12) == 'cryscalc.ins') created_INS_file = 'cryscalc.ins'
  end if

  !call Write_SHX_template(TRIM(created_INS_file), 2, TRIM(main_title), wavelength, INT(Z_unit), crystal_cell, SPG, Atm_list)
  close(unit=4)
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
   !do i=1, SPG%multip   ! whole set of sym. op.
   ! write(4, '(2a)') 'SYMM ', SPG%SymopSymb(i)
   !end do
   do i=2, SPG%NumOps    ! reduced set of sym. op.
    write(4, '(2a)') 'SYMM ', SPG%SymopSymb(i)
   end do

   IF(nb_atoms_type == 0) then
    WRITE(4, '(a)') 'SFAC ?'
    WRITE(4, '(a)') 'UNIT ?'
   ELSE
    IF(nb_atoms_type < 10) then
     WRITE(fmt_SFAC, '(a,i1,a)') "(a,", nb_atoms_type,"(1x,a6))"
     WRITE(fmt_UNIT, '(a,i1,a)') "(a,", nb_atoms_type,"(1x,F6.1))"
    ELSE
     WRITE(fmt_SFAC, '(a,i2,a)') "(a,", nb_atoms_type,"(1x,a6))"
     WRITE(fmt_UNIT, '(a,i2,a)') "(a,", nb_atoms_type,"(1x,F6.1))"
    ENDIF
    WRITE(4, fmt_SFAC) 'SFAC ', (SFAC%type(i)  ,i=1,nb_atoms_type)
    WRITE(4, fmt_UNIT) 'UNIT ', (SFAC%number(i),i=1,nb_atoms_type)
   endif

   write(4, "(a)")    'ACTA'
   write(4, "(a)")    'BOND $H'
   write(4, "(a)")    'WPDB -2'
   write(4, "(a)")    'L.S. 4'
   write(4, "(a)")    'WGHT 0.2'
   write(4, "(a)")    'FVAR 0.1'

             !---- Weight ----!

  do i=1,nb_atom
    !write(UNIT=4,'(a,1x,a4,5F10.5)')  atom_label(i), atom_typ (i) , new_atom_coord(1:3,i), 10.+atom_occ(i) ,  atom_Ueq(i)
    call get_atom_order(TRIM(atom_typ(i)), atom_order)
    !if(atom_adp_equiv(i) < 0.00001 .and. atom_Biso(i) > 0.00001) atom_adp_equiv(i) = atom_Biso(i)/(8*pi**2.)
    !if(atom_adp_equiv(i) < 0.00001) atom_adp_equiv(i) = 0.03
    !if(atom_order == 0) then
    ! write(4,'(a,1x,a,5F10.5)')  atom_label(i), '   ?', atom_coord(1:3,i), 10.+atom_occ(i) ,   atom_adp_equiv(i)
    !else
    ! write(4,'(a,1x,I4,5F10.5)')  atom_label(i), atom_order , atom_coord(1:3,i), 10.+atom_occ(i) ,   atom_adp_equiv(i)
    !endif

    long = len_trim(atom_typ(i))
    if (long == 1 .and. atom_typ(i)(1:1) == 'H' .and. no_H) cycle

    if(atom_adp_equiv(i) < 0.00001) atom_adp_equiv(i) = 8*pi**2. * 0.03
    if(atom_order == 0) then
     write(4,'(a,1x,a,5F10.5)')  atom_label(i), '   ?', atom_coord(1:3,i), 10.+atom_occ(i) ,  atom_adp_equiv(i)/(8*pi**2.)
    else
     write(4,'(a,1x,I4,5F10.5)')  atom_label(i), atom_order , atom_coord(1:3,i), 10.+atom_occ(i) ,   atom_adp_equiv(i)/(8*pi**2.)
    endif
  end do

  WRITE(4, '(a)') 'HKLF  4'
  WRITE(4, '(a)') 'END'


  close (UNIT=4)

  call write_info('')
  if(no_H) then
   write(message_text, '(3a)') '   >> ', trim(created_INS_file), ' file has been created (H atoms have been excluded).'
  else
   write(message_text, '(3a)') '   >> ', trim(created_INS_file), ' file has been created.'
  end if
  call write_info(trim(message_text))
  call write_info('')

  return
end subroutine create_INS_file

!------------------------------------------------------------------------
subroutine create_PCR_file(input_file, extension)
 USE cryscalc_module
 USE Pattern_profile_module

! USE CFML_IO_FORMATS, ONLY : Write_Shx_Template
 USE SHELX_module
 USE IO_module

 implicit none
  CHARACTER (LEN=*), INTENT(IN)        :: input_file, extension
  ! local variables
  INTEGER                              :: i, long
  CHARACTER(LEN=256)                   :: created_PCR_file

  if(debug_proc%level_2)  call write_debug_proc_level(2, "create_PCR_file")

  long = len_trim(input_file)
  if(long >= 12) then
   if(input_file(1:12) == 'cryscalc.pcr') then
    created_PCR_file = 'cryscalc.pcr'
   end if
  else
   i= INDEX(input_file, '.')
   WRITE(created_PCR_file, '(a)') input_file(1:i-1)//'_'//TRIM(extension)//'.PCR'
  endif

  if(X_rays) then   ! RX D8 CSM
   pattern = X_pattern
   PV = X_PV
  else
   pattern = N_pattern
   PV = N_PV
  end if


  OPEN(UNIT=4, FILE=TRIM(created_PCR_file))
   WRITE(4, '(a)')          'COMM PCR file created by CRYSCALC'
   WRITE(4, '(a)')          '!Job Npr Nph Nba Nex Nsc Nor Dum Iwg Ilo Ias Res Ste Nre Cry Uni Cor Opt Aut'
   write(4, '(i4,a)')       Pattern%job, '   5   1   0   0   0   0   1   0   0   0   0   0   0   0   0   0   1   1'

   write(4, '(a)')          '!'
   write(4, '(a)')          '!Ipr Ppl Ioc Mat Pcr Ls1 Ls2 Ls3 NLI Prf Ins Rpa Sym Hkl Fou Sho Ana'
   write(4, '(a)')          '   0   0   1   0   1   0   4   0   0   1   6   1   1   0   0   1   1'
   write(4, '(a)')          '!'
   write(4, '(a)')          '! Lambda1  Lambda2    Ratio    Bkpos    Wdt    Cthm     muR   AsyLim   Rpolarz  2nd-muR -> Patt# 1'
   write(4, '(2F9.6,a)')    pdp_simu%wave, pdp_simu%wave, '  1.00000   70.000  6.0000  0.7899  0.0000  160.00    0.0000  0.0000'
   write(4, '(a)')          '!'
   write(4, '(a)')          '!NCY  Eps  R_at  R_an  R_pr  R_gl     Thmin       Step       Thmax    PSD    Sent0'
   write(4, '(a,3F10.4,a)') '  1  0.10  1.00  1.00  1.00  1.00', Pattern%Xmin, Pattern%step, Pattern%Xmax, '   0.000   0.000'
   write(4, '(a)')          '!'
   write(4, '(a)')          '!'
   write(4, '(a)')          '   0    !Number of refined parameters'
   write(4, '(a)')          '!'
   write(4, '(a)')          '!  Zero    Code    SyCos    Code   SySin    Code  Lambda     Code MORE ->Patt# 1'
   write(4, '(a)')          ' 0.00000   00.0  0.00000    0.0  0.00000    0.0 0.000000    0.00   0'
   write(4, '(a)')          '!   Background coefficients/codes  for Pattern#  1  (Polynomial of 6th degree)'
   write(4, '(F12.3,a)')    Pattern%background, '    0.000       0.000       0.000       0.000       0.000'
   write(4, '(a)')          '       00.00       00.00       00.00       00.00       00.00       00.00'
   write(4, '(a)')          '!-------------------------------------------------------------------------------'
   write(4, '(a)')          '!  Data for PHASE number:   1  ==> Current R_Bragg for Pattern#  1:     0.00'
   write(4, '(a)')          '!-------------------------------------------------------------------------------'
   if(len_trim(main_title) /=0) then
    write(4, '(a)')         TRIM(main_title)
   else
    write(4, '(a)')          'Phase #1:'
   end if
   write(4, '(a)')          '!'
   write(4, '(a)')          '!Nat Dis Ang Pr1 Pr2 Pr3 Jbt Irf Isy Str Furth       ATZ    Nvk Npr More'
   write(4, '(I4,a)')       nb_atom, '   0   0 0.0 0.0 1.0   0   0   0   0   0          0.000   0   5   0'
   write(4, '(2a)')         trim(SPG%SPG_symb), '             <--Space group symbol'
   write(4, '(a)')          '!Atom   Typ       X        Y        Z     Biso       Occ     In Fin N_t Spc /Codes'

   do i=1, nb_atom
    write(4, '(2a,5(1x,F8.5), a)')  &
                            atom_label(i), atom_typ(i), atom_coord(1:3,i), atom_Biso(i), atom_occ(i), &
                            '   0   0   0    0'
    write(4, '(a)')         '                 0.00     0.00     0.00     0.00     0.00'
   end do

   write(4, '(a)')        '!-------> Profile Parameters for Pattern #  1'
   write(4, '(a)')        '!  Scale        Shape1      Bov      Str1      Str2      Str3   Strain-Model'
   write(4, '(F10.6,F10.5,a)')   Pattern%scale, PV%eta0, '   0.00000   0.00000   0.00000   0.00000       0'
   write(4, '(a)')         '     0.000     0.000     0.000     0.000     0.000     0.000'
   write(4, '(a)')         '!       U         V          W           X          Y        GauSiz   LorSiz Size-Model'
   write(4, '(4(1x,F10.6),a)')  PV%U, PV%V, PV%W, PV%eta1, '   0.000000   0.000000   0.000000    0'
   write(4, '(a)')         '      0.000      0.000      0.000      0.000      0.000      0.000      0.000'
   write(4, '(a)')         '!     a          b         c        alpha      beta       gamma'
   write(4, '(6(1x,F10.6))')  unit_cell%param(1:3), unit_cell%param(4:6)
   write(4, '(a)')         '    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000'
   write(4, '(a)')         '!  Pref1    Pref2      Asy1     Asy2     Asy3     Asy4      S_L      D_L'
   write(4, '(a)')         '0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000'
   write(4, '(a)')         '   0.00     0.00     0.00     0.00     0.00     0.00     0.00     0.00'

  close (UNIT=4)

  call write_info('')
  write(message_text, '(3a)') '   >> ', trim(created_PCR_file), ' file has been created.'
  call write_info(trim(message_text))
  call write_info('')

  return
end subroutine create_PCR_file


!------------------------------------------------------------------------
subroutine create_TIDY_file(input_file, extension)
 USE cryscalc_module
 USE IO_module
 USE macros_module, only : u_case

 implicit none
  CHARACTER (LEN=*), INTENT(IN)        :: input_file, extension
  ! local variables
  INTEGER                              :: i, i1
  CHARACTER(LEN=256)                   :: TIDY_file
  real, dimension(3)                   :: coord

  if(debug_proc%level_2)  call write_debug_proc_level(2, "create_TIDY_file")

  i1= INDEX(input_file, '.')
  if(i1 /=0) then
   WRITE(TIDY_file, '(a)') input_file(1:i1-1)//'_'//TRIM(extension)//'_tidy.dat'
  else
   WRITE(TIDY_file, '(a)') 'cryscalc_tidy.dat'
  end if
  OPEN(UNIT=4, FILE=TRIM(TIDY_file))

   write(4, '(1x,a)')          u_case(trim(SPG%SPG_symb))
   write(4, '(6F10.4)')        unit_cell%param(1:3), unit_cell%param(4:6)
   write(4, '(a)')             ''

   do i=1,nb_atom
    coord(1) = atom_coord(1,i)
    coord(2) = atom_coord(2,i)
    coord(3) = atom_coord(3,i)
    call inside_coord(coord)
    write(4,'(a4,4F10.5)')  atom_label(i), coord(1:3), atom_occ(i)
   end do
   write(4, '(a)') 'END'
   close (unit=4)

   call write_info('')
   write(message_text, '(3a)') '   >> ', trim(TIDY_file), ' file has been created.'
   call write_info(trim(message_text))
   call write_info('')

  return

end subroutine create_TIDY_file

!------------------------------------------------------------------------
subroutine create_SOLVE_input_files

 call create_SIR_file
 call create_SHELXS_file
 call create_SUPERFLIP_file

 return
end subroutine create_SOLVE_input_files

!------------------------------------------------------------------------
subroutine create_SIR_file()
 USE cryscalc_module, ONLY : SPG, molecule, unit_cell, debug_proc
 USE HKL_module,      ONLY : HKL_file
 USE IO_module
 implicit none
  INTEGER :: i1
  LOGICAL :: ok

 if(debug_proc%level_2)  call write_debug_proc_level(2, "create_SIR_file")

 call mandatory_keyword('CELL', 'SIR97', ok)
 if(.not. ok) return
 call mandatory_keyword('SPGR', 'SIR97', ok)
 if(.not. ok) return
 call mandatory_keyword('FILE', 'SIR97', ok)
 if(.not. ok) return
 call mandatory_keyword('CHEM', 'SIR97', ok)
 if(.not. ok) return


 i1 = INDEX(HKL_file%name, '.')

  OPEN(UNIT=51, FILE='cryscalc_SIR97.in')
   WRITE(51, '(a)' )       '%window'
   WRITE(51, '(2a)')       '%structure ',  HKL_file%name(1:i1-1)
   WRITE(51, '(a)' )       '%init'
   WRITE(51, '(4a)')       '%job ',  HKL_file%name(1:i1-1), '   in ', TRIM(SPG%SPG_symb)
   WRITE(51, '(a)' )       '%data'
   WRITE(51, '(a,6F12.5)') '      Cell         ',  unit_cell%param(1:6)
   WRITE(51, '(2a)')       '      Space        ',  TRIM(SPG%SPG_symb)
   WRITE(51, '(2a)')       '      Content      ',  molecule%content
   WRITE(51, '(2a)')       '      Reflections  ',  TRIM(HKL_file%HKL)
   WRITE(51, '(a)' )       '      Format (3i4,2F8.2)'
   WRITE(51, '(a)' )       '      Fosquare'
   WRITE(51, '(a)' )       '%normal'
   WRITE(51, '(a)' )       '%invariants'
   WRITE(51, '(a)' )       '%phase'
   WRITE(51, '(a)' )       '%fourier'
   WRITE(51, '(a)' )       '%menu'
   WRITE(51, '(2a)')       '      shelx  ',   HKL_file%name(1:i1-1)//'.res'
   WRITE(51, '(a)' )       '%end'
  close (51)

  call write_info('')
  call write_info('   >> cryscalc_SIR97.in input file for SIR97 has been created')
  call write_info('')



 RETURN
end subroutine create_SIR_file

!------------------------------------------------------------------------
subroutine create_SHELXS_file()
 USE cryscalc_module, ONLY : main_title, keyword_SIZE, SPG, molecule, unit_cell, wavelength, &
                             keyword_Zunit, Z_unit, &
                             crystal, nb_atoms_type, SFAC, Create_INS, debug_proc
 USE HKL_module,      ONLY : HKL_file
 USE SHELX_module

 USE IO_module
 implicit none
  INTEGER :: i, i1
  LOGICAL :: ok


 if(debug_proc%level_2)  call write_debug_proc_level(2, "create_SHELXS_file")

 call mandatory_keyword('CELL', 'SHELXS', ok)
 if(.not. ok) return
 call mandatory_keyword('SPGR', 'SHELXS', ok)
 if(.not. ok) return
 call mandatory_keyword('FILE', 'SHELXS', ok)
 if(.not. ok) return
 call mandatory_keyword('CHEM', 'SHELXS', ok)
 if(.not. ok) return

 i1 = INDEX(HKL_file%name, '.')


   OPEN(UNIT=51, FILE='cryscalc_SHELXS.ins')


    WRITE(51, '(2a)')                                   'TITL ', TRIM(main_title)
    if ((wavelength < 0.01)) then
     WRITE(51, '(a,F8.5, 3(1x,F8.4),3(1x,F8.3))')        'CELL ', 0.71073, unit_cell%param(1:3), unit_cell%param(4:6)
    else
     WRITE(51, '(a,F8.5, 3(1x,F8.4),3(1x,F8.3))')        'CELL ', wavelength, unit_cell%param(1:3), unit_cell%param(4:6)
    endif
    IF(keyword_ZUNIT)  then
     WRITE(51, '(a,F8.2,3(1x,F8.4),3(1x,F8.3))')   'ZERR ', Z_unit, unit_cell%param_ESD(1:3), &
                                                                    unit_cell%param_ESD(4:6)
    else
     WRITE(51, '(a,F8.2,3(1x,F8.4),3(1x,F8.3))')   'ZERR ', molecule%Z_unit, unit_cell%param_ESD(1:3), &
                                                                             unit_cell%param_ESD(4:6)
    end if

   call get_SHELX_nlatt()
   WRITE(51, '(a,I6)')   'LATT ', n_latt
   do i=2, SPG%NumOps    ! reduced set of sym. op.
    write(51, '(2a)') 'SYMM ', SPG%SymopSymb(i)
   end do

   IF(nb_atoms_type == 0) then
    WRITE(51, '(a)') 'SFAC ?'
    WRITE(51, '(a)') 'UNIT ?'
   ELSE
    IF(nb_atoms_type < 10) then
     WRITE(fmt_SFAC, '(a,i1,a)') "(a,", nb_atoms_type,"(1x,a6))"
     WRITE(fmt_UNIT, '(a,i1,a)') "(a,", nb_atoms_type,"(1x,F6.1))"
    ELSE
     WRITE(fmt_SFAC, '(a,i2,a)') "(a,", nb_atoms_type,"(1x,a6))"
     WRITE(fmt_UNIT, '(a,i2,a)') "(a,", nb_atoms_type,"(1x,F6.1))"
    ENDIF
    WRITE(51, fmt_SFAC) 'SFAC ', (SFAC%type(i)  ,i=1,nb_atoms_type)
    WRITE(51, fmt_UNIT) 'UNIT ', (SFAC%number(i),i=1,nb_atoms_type)
   endif

   WRITE(51, '(a,F8.0)') 'TEMP ', Create_INS%temperature
   IF(keyword_SIZE) then
    WRITE(51, '(a,3(1x,F6.3))')      'SIZE  ', crystal%size(1:3)
   endif
   WRITE(51, '(a)') 'PATT'
   WRITE(51, '(a)') 'HKLF  4'
   WRITE(51, '(a)') 'END'




  call write_info('')
  call write_info('   >> cryscalc_SHELXS.ins input file for SHELXS/T has been created')
  call write_info('')



 RETURN
end subroutine create_SHELXS_file

!------------------------------------------------------------------------
 subroutine create_SUPERFLIP_file
  use cryscalc_module, only : main_title, unit_cell, SPG, molecule, wavelength
  USE HKL_module,      ONLY : HKL_file
  use macros_module,   only : replace_car, replace_car1
  use IO_module,       only : write_info

  implicit none
   character(len=48)              :: spf_symop
   integer                        :: i
   logical                        :: ok


  call mandatory_keyword('CELL', 'SUPERFLIP', ok)
  if(.not. ok) return
  call mandatory_keyword('SPGR', 'SUPERFLIP', ok)
  if(.not. ok) return
  call mandatory_keyword('FILE', 'SUPERFLIP', ok)
  if(.not. ok) return
  call mandatory_keyword('CHEM', 'SUPERFLIP', ok)
  if(.not. ok) return


  OPEN(UNIT=51, FILE='cryscalc_SUPERFLIP.in')
      WRITE(51,'(a)') '# Superflip input file from CRYSCALc'
      WRITE(51,'(a)') ''
      WRITE(51,'(2a)') 'title ',trim(main_title)
      WRITE(51,'(a)') 'perform CF'
      WRITE(51,'(a)') 'terminal yes'
      WRITE(51,'(a)') ''
      WRITE(51,'(a)') '# keywords describing the file output'
      WRITE(51,'(a)') ''
      WRITE(51,'(a)') 'outputfile superflip.xplor'
      WRITE(51,'(a)') 'outputformat xplor'
      WRITE(51,'(a)') ''
      WRITE(51,'(a)') '# keywords describing crystallographic data'
      WRITE(51,'(a)') ''
      WRITE(51,'(a)') 'dimension 3'
      WRITE(51,'(a)') 'realdimension 3'
      WRITE(51,'(a)') 'voxel  AUTO'
      WRITE(51,'(a,6f10.4)')'cell ',unit_cell%param(1:6)
      WRITE(51,'(a,F10.5)') 'lambda ', wavelength
      WRITE(51,'(a)') ''
      WRITE(51,'(2a)') '# Space group  ',trim(SPG%SPG_symb)
      WRITE(51,'(a)') ''
      WRITE(51,'(a)') 'symmetry'
      do i=1, SPG%Multip
       WRITE(spf_symop, '(a)') trim(SPG%SymopSymb(i))
       spf_symop = replace_car(spf_symop, ",", "  ")
       spf_symop = replace_car1(spf_symop, "x", "x1")
       spf_symop = replace_car1(spf_symop, "y", "x2")
       spf_symop = replace_car1(spf_symop, "z", "x3")
       WRITE(51,'(a)') spf_symop
      end do
      WRITE(51,'(a)') 'endsymmetry'
      WRITE(51,'(a)') 'centers'
      WRITE(51,'(a)') ' 0.00000 0.00000 0.00000'
      WRITE(51,'(a)') 'endcenters'
      WRITE(51,'(a)') ''
      WRITE(51,'(a)') '# keywords influencing the algorithm'
      WRITE(51,'(a)') ''
      WRITE(51,'(a)') 'derivesymmetry yes'
      WRITE(51,'(a)') 'searchsymmetry average'
      WRITE(51,'(a)') 'delta AUTO'
      WRITE(51,'(a)') 'weakratio   0.00'
      WRITE(51,'(a)') 'biso        2.500'
      WRITE(51,'(a)') 'randomseed  AUTO'
      WRITE(51,'(a)') 'maxcycles 10000'
      WRITE(51,'(a)') 'normalize no'

      !if(len_trim(molecule%content) /=0) then
      ! WRITE(51,'(2a)') 'composition ', trim(molecule%content)
      if(len_trim(molecule%formula) /=0) then
       WRITE(51,'(2a)') 'composition ', trim(molecule%formula)
      else
       WRITE(51,'(a)') '#composition C62 H46 N14 (cell content !)'
      endif
      WRITE(51,'(a)') 'polish yes'
      WRITE(51,'(a)') 'coverage yes'
      WRITE(51,'(a)') 'convergencemode normal       0.80'
      WRITE(51,'(a)') 'expandedlog no'
      WRITE(51,'(a)') 'bestdensities 1'
      WRITE(51,'(a)') 'repeatmode  1'
      WRITE(51,'(a)') ''
      WRITE(51,'(a)') '# Reflection list, as intensities (F**2)'
      WRITE(51,'(a)') ''
      WRITE(51,'(a)') 'dataformat shelx'
      WRITE(51,'(2a)') 'fbegin  ', TRIM(HKL_file%HKL)
      WRITE(51,'(a)') 'endf'


  call write_info('')
  call write_info('   >> cryscalc_SUPERFLIP.in input file for SUPERFLIP has been created')
  call write_info('')

  return
 end subroutine create_SUPERFLIP_file

!------------------------------------------------------------------------
 subroutine mandatory_keyword(input_string, prog, ok)
  use cryscalc_module, only  : keyword_CELL, keyword_SPGR, keyword_CONT, keyword_SFAC_UNIT, keyword_CHEM, keyword_ZUNIT, &
                               keyword_FILE, SPG
  use IO_module,       only  : write_info
  implicit none
  character (len=*), intent(in)  :: input_string
  character (len=*), intent(in)  :: prog
  logical,           intent(out) :: ok

  ok = .false.

 select case(trim(input_string))
  case ('CELL')
   IF(.NOT. keyword_CELL) then
    call write_info('')
    call write_info('  >>> CELL keyword is mandatory to create an input file for '//trim(prog)//' <<<')
    call write_info('')
   else
    ok = .true.
   endif

  case ('SPGR')
   !IF(.NOT. keyword_SPGR) then
   if(SPG%numSPG ==0) then
    call write_info('')
    call write_info('  >>> SPGR keyword is mandatory to create an input file for '//trim(prog)//' <<<')
    call write_info('')
   else
    ok = .true.
   endif

  case ('FILE')
   IF(.NOT. keyword_FILE) then
    call write_info('')
    call write_info('  >>> FILE keyword is mandatory to create an input file for '//trim(prog)//' <<<')
    call write_info('')
   else
    ok = .true.
   endif

  case ('CHEM')
   IF(.NOT. keyword_CONT .and. .not. keyword_SFAC_UNIT ) then
    IF(.NOT. keyword_CHEM .or. (keyword_CHEM .and. .not. keyword_ZUNIT)) then
     call write_info('')
     call write_info('  >>> Atom content (keywords: CONT, SFAC/UNIT, CHEM/ZUNIT) is mandatory to create an input file ' &
                     // trim(prog) // ' <<<')
     call write_info('')
    else
     ok = .true.
    endif
   endif

 end select




  return
 end subroutine mandatory_keyword

!------------------------------------------------------------------------
subroutine get_SHELX_nlatt()
 USE cryscalc_module, ONLY : SPG, debug_proc
 use SHELX_module,    only : n_latt

 if(debug_proc%level_2)  call write_debug_proc_level(2, "get_SHELX_nlatt")


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


 return
END subroutine get_SHELX_nlatt

