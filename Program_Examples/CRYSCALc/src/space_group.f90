!     Last change:  TR   14 Feb 2007    2:20 pm
subroutine space_group_info
 use cryscalc_module,                ONLY : ON_SCREEN, space_group_symbol, symm_op_string,        &
                                            keyword_SYMM, symm_op_xyz, symm_op_mat, nb_symm_op,   &
                                            WRITE_SPG_info, Write_SPG_info_all, WRITE_SPG_EXTI,   &
                                            write_SPG_subgroups, message_text,                    &
                                            SPG, keyword_create_CIF,                              &
                                            keyword_CELL, unit_cell, CIF_cell_measurement,        &
                                            CIF_diffrn_reflns, input_PCR, known_cell_esd, keyword_read_INS, &
											keyword_ATOM_list, nb_atom, debug_proc
 
 use CFML_crystallographic_symmetry, ONLY : set_spacegroup, write_spacegroup,  Symmetry_Symbol,   &
                                            Wyckoff_Type, Wyck_Pos_Type,                          &
                                            get_T_SubGroups, space_group_type,                    &
                                            Lattice_trans, similar_transf_SG
                                            
 Use CFML_Reflections_Utilities,     ONLY : search_extinctions
 USE macros_module,                  ONLY : replace_car
 USE IO_module,                      ONLY : write_info


 implicit none
! TYPE (Wyckoff_Type)               :: Wyckoff
 CHARACTER (LEN=256)               :: read_line
 INTEGER                           :: i, i1, i2, i_error
 character (len=100), dimension(24):: texto
 CHARACTER (LEN=40)                :: input_line
 REAL                              :: esd
 CHARACTER (LEN=6)                 :: esd_string

 ! variables for sub_groups
  real, dimension (3,3)                 :: trans
  real, dimension (3  )                 :: orig
  type(space_group_type), dimension(48) :: Subgroup
  type(space_group_type)                :: SpGn
  integer                               :: nsg, j, ng, l
  logical                               :: enantio, chiral, polar
  CHARACTER (len=32)                    :: MOVE_string

  if(debug_proc%level_2)  call write_debug_proc_level(2, "SPACE_GROUP_INFO")

 !open (unit=3, file="sg.out", status="replace")
 call set_spacegroup(space_group_symbol, SPG)

 IF( SPG%NumSpg ==0)  then
  call write_info('')
  call write_info(' >>> Erroneous space group !!!')
  call write_info('')
  return
 endif

 IF(WRITE_SPG_info) then
  CLOSE(UNIT=3)
  open (unit=3, file="sg.out", status="replace")
  IF(write_SPG_info_all) then
   call write_spacegroup(SPG,  iunit=3, full=.true.)
  else
   call write_spacegroup(SPG,  iunit=3)
  endif
  close (UNIT=3)

!  call system('echo off')
!  call system('type sg.out')
!  close(unit=2)
!  call system('copy cryscalc.log + sg.out')
!  call system('echo on')

  if(ON_SCREEN) then
   OPEN(UNIT=3, FILE="sg.out")
   do
    READ(3, '(a)', IOSTAT=i_error) read_line
    IF(i_error < 0) exit
    call write_info(TRIM(read_line))
   END do
   CLOSE(UNIT=3)
  endif 
  
  CLOSE(UNIT=2)
  open (UNIT=2, FILE='cryscalc.log', ACTION="write",position = "append")
 else
  IF(.NOT. input_PCR .and. .not. keyword_read_INS) then    
   call write_info('')   
   WRITE(message_text, '(2a)')   '          > HM symbol : ', SPG%SPG_symb
   call write_info(TRIM(message_text))
   WRITE(message_text, '(a,i3)') '                 IT # : ', SPG%NumSpg
   call write_info(TRIM(message_text))
   if (SPG%centred == 2) then     ! groupe centro.
    WRITE(message_text, '(a,i3)') '              centric : yes (-1 at origin)'
   elseif (SPG%centred == 0) then ! groupe centro.
    WRITE(message_text, '(a,i3)') '              centric : yes (-1 not at origin)'
   elseif(SPG%centred == 1) then  ! groupe non centro
    WRITE(message_text, '(a,i3)') '              centric : no'
   endif
   call write_info(TRIM(message_text))

   if( SPG%centred == 1) then
    call test_chiral(SPG%NUmSPG, chiral)
    if(chiral) then
     WRITE(message_text, '(a,i3)') '               chiral : yes'
    else
     WRITE(message_text, '(a,i3)') '               chiral : no'
    endif
    call write_info(TRIM(message_text))
   endif	
   
   WRITE(message_text, '(a,i3)') ' general multiplicity : ', SPG%multip
   call write_info(TRIM(message_text))
   
   call write_info('')   
   call write_info('    Enter SG_INFO keyword for more details about selected space group.')
   call write_info('')
  endif
 END if

 ! liste le 'reduced set' des operateurs de symetrie
 IF(WRITE_SPG_info)   call write_symm_op_reduced_set()             ! space_group.F90

 IF(write_spg_exti) then
  call test_hkl_equiv_condition(SPG%NUmSPG)
  CLOSE(UNIT=3)
  open (unit=3, file="sg.out", status="replace")
  call search_extinctions(SPG, 3)
  close (UNIT=3)
!  call system('echo off')
!  call system('type sg.out')
!  close(unit=2)
!  call system('copy cryscalc.log + sg.out')
!  call system('echo on')


  if(ON_SCREEN) then
   OPEN(UNIT=3, FILE="sg.out")
   do
    READ(3, '(a)', IOSTAT=i_error) read_line
    IF(i_error < 0) exit
    call write_info(TRIM(read_line))
   enddo
  endif
  CLOSE(UNIT=3) 
  CLOSE(UNIT=2)
  open (UNIT=2, FILE='CRYSCALC.log', ACTION="write",position = "append")
 endif


 ! operateurs de symétrie du groupe
 do i=1,SPG%Multip
   texto(1)=" "
   call Symmetry_Symbol(SPG%SymopSymb(i),texto(1))

   input_line = trim(SPG%SymopSymb(i))
   !call replace_car(input_line, ',', ' ')
   input_line = replace_car(input_line, ',', ' ')
   i1=INDEX(input_line, ' ')
   symm_op_string(1,i) = input_line(1:i1-1)
   i2=INDEX(TRIM(input_line), ' ' , back=.TRUE.)
   symm_op_string(2,i) = input_line(i1+1:i2-1)
   symm_op_string(3,i) = input_line(i2+1:)

   keyword_SYMM   = .true.
   symm_op_xyz    = .true.
   symm_op_mat    = .false.

 end do
 nb_symm_op = SPG%Multip
 call get_op_mat  ! permet de recuperer les parties rot. et transl. des operateurs de symm.

 IF(keyword_create_CIF) then
  call write_CIF_file('UNIT_CELL_INFO')
  call write_CIF_file('SPACE_GROUP')
  !do i=2, SPG%Multip
  ! WRITE(CIF_string(1), '(a)')   SPG%SymopSymb(i)
  ! call write_CIF_file('SYM_OP', 1, CIF_string(1))
  !end do
  call write_CIF_file('')


  IF(keyword_CELL) then
   IF(unit_cell%volume < 0.1) call volume_calculation('out')
   !call write_CIF_file('UNIT_CELL_INFO')
   !IF(known_CELL_esd) then
   ! call write_CIF_file('CELL_PARAM_ESD')
   !else
   ! call write_CIF_file('CELL_PARAM')
   !endif

   ! CIF_Thmin, CIF_Thmax
    call write_CIF_file("CIF_THMIN,THMAX")
	
	IF(nb_atom /=0 .AND. ON_SCREEN)  call write_atom_list
  endif
 endif

 ! Wyckoff positions
 !call write_info('')
 !call write_info(' >> Wyckoff positions:')
 !write(message_text,*)'            .num_orbit: ', SPG%Wyckoff%num_orbit
 !call write_info(message_text))
 !call write_info(' ')
 !do i=1, SPG%Wyckoff%num_orbit
 ! write(message_text,*) '  site multiplicity: ', SPG%wyckoff%orbit(i)%multp
 ! call write_info(message_text))
 ! write(message_text,*) '  site norb        : ', SPG%wyckoff%orbit(i)%norb
 ! call write_info(message_text))
 ! do j=1,SPG%wyckoff%orbit(i)%norb,3
 ! write(message_text, *) '          ', SPG%wyckoff%orbit(i)%str_orbit(j:j+2)
 ! call write_info(message_text))
 ! END do
 !end do


 ! determination des regles d'extinctions du groupe d'espace
 !if (WRITE_SPG_exti) call search_extinctions(SPG)

 ! sub-groups
  if (write_SPG_subgroups) then
   ! trans = identity
   trans(:,:) = 0.
   trans(1,1) = 1.
   trans(2,2) = 1.
   trans(3,3) = 1.
   ! origine 
   orig(:) = 0.
   
   !> Construct the subgroup of SPG that is compatible
   call similar_transf_SG(trans,orig,SpG,SpGn)

   !> with the transformation matrix and change of origin give above
   call write_spacegroup(SPGn,full=.true.)

   !> Determine all subgroups of the new space group
   call get_T_SubGroups(SPGn,SubGroup,nsg)
   
   write(message_text, "(a)") " => LIST of Translationengleische Subgroups: "
   call write_info('')
   call write_info(trim(message_text))
   call write_info('')

   do i=1,nsg
    j=SPGn%Multip/SubGroup(i)%multip
    ng=SubGroup(i)%numops
    write(message_text,fmt="(4a,i2,30a)") " => ", SubGroup(i)%Spg_Symb, SubGroup(i)%hall,&
                                          " Index: [",j,"]   ->  { ", (trim(SubGroup(i)%SymopSymb(l))//" : ",l=1,ng-1),&
                                          trim(SubGroup(i)%SymopSymb(ng))," }    ", trim(SubGroup(i)%centre)
                                    
    call write_info(trim(message_text))

   end do   
  endif

 
  if(WRITE_SPG_info .AND. SPG%centred ==1 .and. on_screen) then
    call test_enantio(SPG%NumSPG, enantio)    
	call test_chiral(SPG%NUmSPG, chiral)
    call test_polar(SPG%NUmSPG, polar)	
	call write_info("")   
	if(chiral) then
	 call write_info(" => Chiral space group:          yes")
	else
	 call write_info(" => Chiral space group:          no")
	endif
	if(polar) then
	 call write_info(" => Polar space group:           yes")
	else
	 call write_info(" => Polar space group:           no")
	endif
    if(enantio) then
     call write_info(" => Enantiomorphic space group:  yes")
 	 call write_info(" => Inverse structure in SHELXL: MOVE 1 1 1 -1 and invert translation parts of the s.o.")
    else	
	 call get_MOVE_string(Move_string)
     call write_info(" => Enantiomorphic space group:  no")
	 call write_info(" => Inverse structure in SHELXL: "//trim(MOVE_string))
    endif 
   endif 	

 return
end subroutine space_group_info

!-----------------------------------------------------------------------------------
subroutine get_SITE_info(input_string)
 USE cryscalc_module 
 USE IO_module,                      ONLY : write_info
 USE CFML_crystallographic_symmetry, ONLY : set_spacegroup, Get_orbit, Get_Multip_Pos, &
                                            symmetry_symbol, Wyckoff_Type, Get_stabilizer, &
											Get_symSymb
 USE macros_module,                  ONLY : u_case
 USE CFML_Math_General,              ONLY : acosd, set_epsg, set_epsg_default

 !USE CFML_constants,                 ONLY : sp 
 USE CFML_GlobalDeps,                 ONLY : sp, cp

 implicit none
  character (len=*), intent(in)      :: input_string
 ! local variables
  INTEGER                            :: i, j, k
  REAL (kind=cp), DIMENSION(3)       :: r
  INTEGER                            :: site_multiplicity
  INTEGER                            :: Wyck_mult

  type(wyckoff_type)   :: wyckoff

  real (kind=cp) , dimension(3,192)  :: orbit
  CHARACTER (LEN=80)                 :: text80
  CHARACTER (LEN=40)                 :: Symm_Symb
  integer, dimension(192)            :: effsym  ! Pointer to effective symops
  integer, dimension(192)            :: ss_ptr  ! Pointer to stabilizer    >> oct. 2011
  !integer                            :: ss_ptr  
  integer                            :: order
  real(kind=cp),     dimension(3,48) :: atr
  integer, dimension(3,3)            :: Rsym
  real,    dimension(3)      :: tt

  
  LOGICAL, DIMENSION(nb_atom)        :: write_site_info_label
  
  ! pour determiner les constraintes sur les Bij 
  REAL , DIMENSION(6)                :: beta_ij 
  integer                            :: codini
  real(kind=cp), dimension(6)        :: codes  !codewords for positions

  if(debug_proc%level_2)  call write_debug_proc_level(2, "GET_SITE_INFO ("//trim(input_string)//")")

  call set_epsg(0.001_cp)  !In mathematical comparisons -> CFML_Math_General

  !call set_spacegroup(space_group_symbol, SPG)     ! <<<<<<<<<

  IF(.NOT. site_info_all_atoms) then
    do i=1,nb_atom
     write_site_info_label(i) = .false.

     do j=1, nb_atom_site_info
      IF(u_case(atom_label(i)) == u_case(site_info_label(j)) ) then
       write_site_info_label(i) = .true.
       cycle
      endif
     END do

    END do
  else
   write_site_info_label(1:nb_atom) = .true.
  endif


  do i=1, nb_atom
   IF(.NOT. write_site_info_label(i)) cycle

   text80 = ''
   r(1:3)= atom_coord(1:3,i)  ! coord. atomiques

   call write_info(' ')
   !WRITE(message_text,'(a,3F10.5)') '  ATOM POSITIONS: ',r(1:3)
   WRITE(message_text,'(a,i4, a6)') '  ATOM #: ',i, trim(atom_label(i))
   call write_info(trim(message_text))
   call write_info(' ')


 !  call Get_Orbit(r, SPG, Wyck_mult, orbit, effsym, ss_ptr )

   call Get_Orbit(r, SPG, wyckoff%num_orbit, orbit, effsym)
!   SITE_multiplicity = Get_Multip_Pos(r, SPG)

   ! oct. 2011
   call Get_stabilizer(r, SPG, order, ss_ptr, atr)

   !WRITE(message_text,'(a,I5)')     '     . site multiplicity: ', wyckoff%num_orbit
   !call write_info(trim(message_text))
   
   site_multiplicity = Get_Multip_Pos(r, SPG)
   WRITE(message_text,'(a,I5)')     '     . site multiplicity: ', site_multiplicity
   call write_info(trim(message_text))
   
   WRITE(message_text,'(a,F12.8)')  '     . site occupancy   : ', real(site_multiplicity)/SPG%multip
   call write_info(trim(message_text))
   
   
   call write_info('     . list of symmetry operators and symmetry elements of the site point group:')

   !do j=1, site_multiplicity
   do j=1, order
     Rsym = SPG%SymOp(ss_ptr(j))%Rot
     tt   = SPG%SymOp(ss_ptr(j))%tr + atr(:,j)
     call Get_SymSymb(Rsym,tt,Symm_Symb)
     call symmetry_symbol(Symm_Symb,text80)  
	 write(message_text,'(7x, a,i3,a,a,a)') ' Operator ',j,': ', Symm_Symb,trim(text80)
	 call write_info(trim(message_text))
   end do
 


   ! orbit
    call write_info(' ')
    call write_info('     --> Complete orbit of the atom: ')
    call write_info(' ')

   effsym(1) = 1
   do j=1,  wyckoff%num_orbit
    SPG%SymopSymb(effsym(1)) = 'x,y,z' 
	if(j <10) then
     write(message_text,"(10x,a6,a,i1,a,3f9.5,tr4,a)") trim(atom_label(i)),'_',j,":   ",orbit(:,j),trim(SPG%SymopSymb(effsym(j)))
    elseif(j < 100) then
     write(message_text,"(10x,a6,a,i2,a,3f9.5,tr4,a)") trim(atom_label(i)),'_',j,":  ", orbit(:,j),trim(SPG%SymopSymb(effsym(j)))
    else
     write(message_text,"(10x,a6,a,i3,a,3f9.5,tr4,a)") trim(atom_label(i)),'_',j,": ",  orbit(:,j),trim(SPG%SymopSymb(effsym(j)))
    endif
    call write_info(trim(message_text))  
   END do
   
   if(len_trim(input_string)==3 .and. input_string(1:3) == 'pcr') then
    call write_info('')
	call write_info('!Atom Typ       X        Y        Z     Biso       Occ     In Fin N_t Spc /Codes')
    call write_info('!    beta11   beta22   beta33   beta12   beta13   beta23  /Codes')
    do j=1,  wyckoff%num_orbit
     write(message_text, '(a6,a,i1,2a,3f9.5,3x,a)') trim(atom_label(i)),'_',j,"   ",trim(atom_typ(i)), orbit(:,j), &
	                                                '   0.300    1.000   0   0   0'      
	 call write_info(trim(message_text))
	 write(message_text, '(17x,a)')   '0.000    0.000    0.000      0.000    0.000' 
	 call write_info(trim(message_text))
	end do
   end if
   
   if(len_trim(input_string)==7 .and. input_string(1:7) == 'pcr_mag') then
    call write_info('')
	call write_info('!Atom Typ  Mag Vek    X      Y      Z       Biso   Occ      Rx      Ry      Rz')
    call write_info('!     Ix     Iy     Iz    beta11  beta22  beta33   MagPh')
    do j=1,  wyckoff%num_orbit
     write(message_text, '(a6,a,i1,3a,3f8.5,a)') trim(atom_label(i)),'_',j,"   ?",trim(atom_typ(i)), '?   1   0  ', orbit(:,j), &
	                                                ' 0.30000 1.00000   0.000   0.000   0.000'      
	 call write_info(trim(message_text))
	 write(message_text, '(29x,a)')   '0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00' 
	 call write_info(trim(message_text))
	 write(message_text, '(20x,a)')   '0.000   0.000   0.000   0.000   0.000   0.000 0.00000' 
	 call write_info(trim(message_text))
	 write(message_text, '(20x,a)')   ' 0.00    0.00    0.00    0.00    0.00    0.00    0.00' 
	 call write_info(trim(message_text))
	end do
   end if
   

   
   ! constrains on ADP
  
   beta_ij(1:6) = 0.01 
   codini = 0
   codes(:) = 1.
   call Get_Atombet_Ctr(r, beta_ij, SPG, codini, codes)

   call write_info(' ')
   call write_info('     --> Constraints on anisotropic ADP (11 22 33 12 13 23): ')
   call write_info(' ')
   write(message_text, '(8x,6(2x,F6.2))') codes(1:6)
   call write_info(trim(message_text))  



  end do
  
  call set_epsg_default()
 return

end subroutine get_site_info
!---------------------------------------------------------------------------

subroutine write_symm_op_reduced_set()
 USE IO_module,                      ONLY : write_info
 use CFML_crystallographic_symmetry, ONLY : set_spacegroup, searchop, write_sym
 USE cryscalc_module,                ONLY : ON_SCREEN, SPG, message_text, debug_proc
 USE CFML_symmetry_tables,           ONLY : intsymoh,x_oh, intsymd6h,x_d6h

 implicit none
  CHARACTER(LEN=3)                :: string_numor
  INTEGER                         :: i, j, i1, i2
  INTEGER, DIMENSION(192)         :: indx_op

  if(debug_proc%level_2)  call write_debug_proc_level(2, "WRItE_SYMM_OP_REDUCED_SET")
  
  !open (unit=3, file="sg.out", status="replace")
  !do i=1, SPG%multip
  ! WRITE(*,*)
  ! WRITE(*,*) i, SPG%Symop(i)%rot(1,1:3),  SPG%Symop(i)%tr(1)
  ! WRITE(*,*) i, SPG%Symop(i)%rot(2,1:3),  SPG%Symop(i)%tr(2)
  ! WRITE(*,*) i, SPG%Symop(i)%rot(3,1:3),  SPG%Symop(i)%tr(3)
  !
  ! call write_sym(3,SPG%Symop(i)%rot,  )
  !end do

 ! do i=1, SPG%multip
 !  call searchop(SPG%Symop(i)%rot(1:3,1:3), i1, i2, indx_op(i))
 !  call write_sym(3, indx_op(i), SPG%Symop(i)%rot(1:3,1:3), SPG%Symop(i)%tr(1:3), 0., .false.)
 ! END do

  if(ON_SCREEN) then
   call write_info('')
   call write_info('  -> Reduced set of symmetry operators:')
   call write_info('')
   call write_info('  No.  IT   Symmetry symbol    Rotation part     Associated Translation')
  endif

  IF(.not. SPG%hexa) then
   i1 = 1
   i2 = 24
   do i=1, SPG%NumOps
    call searchop(SPG%Symop(i)%rot(1:3,1:3), i1, i2, indx_op(i))
    j = ABS(indx_op(i))
    if (j<0) j=j+24
    if(ON_SCREEN) then
     WRITE(message_text,'(i4,a,i2,a,a14,a,a14,a,3f8.4,a)') i,': (',j,') ', intsymoh(j),' --> ',x_oh(j),   &
	                                                       ' + {', SPG%Symop(i)%tr(1:3),'}'
     call write_info(TRIM(message_text))
    endif 
   end do

  else
   i1 = 25
   i2 = 36
   do i=1, SPG%NumOps
    call searchop(SPG%Symop(i)%rot(1:3,1:3), i1, i2, indx_op(i))
    j=ABS(indx_op(i))-24
    IF(j < 0) j=j+12
    if(ON_SCREEN) then
     WRITE(message_text,'(i4,a,i2,a,a14,a,a14,a,3f8.4,a)') i,': (',j,') ', intsymd6h(j),' --> ',x_d6h(j),  &
	                                                       ' + {',SPG%Symop(i)%tr(1:3),'}'
     call write_info(TRIM(message_text))
    endif 
   end do
  endif
  if(ON_SCREEN) call write_info('')



 return
end subroutine write_symm_op_reduced_set

!!========================================================================================================
!! extracted from fp_codes.F90
!!
    !!
    !!----  Subroutine Get_Atombet_Ctr(x,betas,Spgr,codini,codes,ord,ss,att,Ipr)
    !!----     real(kind=cp), dimension(3),     intent(in    ) :: x      !Atom position (fractional coordinates)
    !!----     real(kind=cp), dimension(6),     intent(in out) :: betas  !Anisotropic temperature factors
    !!----     type(Space_Group_type),          intent(in    ) :: Spgr   !Space Group
    !!----     Integer,                         intent(in out) :: codini !Last attributed parameter
    !!----     real(kind=cp), dimension(6),     intent(in out) :: codes  !codewords for positions
    !!----     integer,               optional, intent(in    ) :: ord    !Order of the stabilizer
    !!----     integer, dimension(:), optional, intent(in    ) :: ss     !Pointer to SymmOp. of stabilizer
    !!----     integer,               optional, intent(in    ) :: Ipr    !Printing unit for debug
    !!----
    !!----  Subroutine to get the appropriate constraints in the refinement codes of
    !!----  anisotropic atomic displacement(thermal) parameters.
    !!----  New algorithm based in the Wigner theorem.
    !!----  The matrix Bet = Sum { R Beta RT} displays the symmetry constraints to be
    !!----  applied to the anisotropic temperature factors. The sum runs over all rotational
    !!----  symmetry operators of the stabilizer of the particular atom position in the given
    !!----  space group.
    !!----  There are a total of 29 kind of relations that may appear in the Bet matrix:
    !!----
    !!----     1    A A A 0   0   0  -> m-3m, -43m, 432, m-3,23, 3[111].2[001]
    !!----     2    A A C 0   0   0  -> 4/mmm, -42m, 4mm, 422, 4/m, -4,4, 4[001]
    !!----     3    A B A 0   0   0  -> 4[010]
    !!----     4    A B B 0   0   0  -> 4[100]
    !!----     5    A A A D   D   D  -> -3m, 3m, 32, -3, 3   3[111]
    !!----     6    A A A D  -D  -D  -> 3[11-1]
    !!----     7    A A A D  -D   D  -> 3[1-11]
    !!----     8    A A A D   D  -D  -> 3[-111]
    !!----     9    A A C A/2 0   0  -> 6/mmm, -6m2, 6mm, 622, 6/m, 6,-6,-3m, 32,-3, 3:  h 3[001]
    !!----    10    A B C 0   0   0  -> mmm, mm2, 222  2[001] 2[100]
    !!----    11    A A C D   0   0  -> 2[001], 2[110]    w
    !!----    12    A B A 0   E   0  -> 2[010], 2[101]
    !!----    13    A B B 0   0   F  -> 2[100], 2[011]
    !!----    14    A B C B/2 0   0  -> 2[001], 2[100]    h
    !!----    15    A B C A/2 0   0  -> 2[001], 2[010]    h
    !!----    16    A B C D   0   0  -> 2/m, m, 2: 2[001] w
    !!----    17    A B C 0   E   0  -> 2[010]
    !!----    18    A B C 0   0   F  -> 2[100]
    !!----    19    A A C D   E  -E  -> 2[110]            w
    !!----    20    A A C D   E   E  -> 2[1-10]           w
    !!----    21    A B A D   E  -D  -> 2[101]
    !!----    22    A B A D   E   D  -> 2[10-1]
    !!----    23    A B B D  -D   F  -> 2[011]
    !!----    24    A B B D   D   F  -> 2[01-1]
    !!----    25    A B C B/2 F/2 F  -> 2[100]            h
    !!----    26    A B C A/2 0   F  -> 2[210]            h
    !!----    27    A B C B/2 E   0  -> 2[120]            h
    !!----    28    A B C A/2 E   E/2-> 2[010]            h
    !!----    29    A B C D   E   F  -> 1, -1
    !!----
    !!----   Updated: 14 July 2011
    !!----

    !Subroutine Get_Atombet_Ctr(x,betas,Spgr,codini,codes,ord,ss,att,Ipr)	
	Subroutine Get_Atombet_Ctr(x,betas,Spgr,codini,codes)	
	  USE CFML_globaldeps,                only : cp
	  USE CFML_crystallographic_symmetry, only : Space_group_type, Get_stabilizer, Sym_b_relations
	  
	  
	  implicit none
       real(kind=cp), dimension(3),             intent(in    ) :: x
       real(kind=cp), dimension(6),             intent(in out) :: betas
       type(Space_Group_type),                  intent(in    ) :: Spgr
       Integer,                                 intent(in out) :: codini
       real(kind=cp), dimension(6),             intent(in out) :: codes
       !real(kind=cp), dimension(:,:), optional, intent(in)     :: att
       !integer,                       optional, intent(in    ) :: ord, Ipr
       !integer, dimension(:),         optional, intent(in    ) :: ss

       ! Local variables
       character (len=1), dimension(6)   :: cdd
       real(kind=cp),     dimension(6)   :: multip
       integer                           :: i,j,order
       real(kind=cp)                     :: suma
       integer,           dimension(48)  :: ss_ptr
       integer,           dimension(6)   :: codd
       integer,           dimension(3,3) :: Rsym
       real(kind=cp),     dimension(3,3) :: bet,bett,Rs
       real(kind=cp),     dimension(3)   :: tr
       real(kind=cp),     dimension(6)   :: cod,multi
       real(kind=cp),     dimension(3,48):: atr
       real(kind=cp),     parameter      :: epss=0.01_cp

       suma=0.0
       do j=1,6
          suma=suma+abs(codes(j))
          cod(j)=int(abs(codes(j))/10.0_cp)             !Input Parameter number with sign
          multi(j)=mod(codes(j),10.0_cp)                !Input Multipliers
          if(cod(j) < 1.0 .and. abs(multi(j)) > epss)  then
               codini=codini+1
               cod(j) = real(codini)
          end if
       end do
       if(suma < epss) return  !No refinement is required

       !if(present(ord) .and. present(ss) .and. present(att)) then
       !  order=ord
       !  ss_ptr(1:order) = ss(1:ord)
       !  atr(:,1:order)  = att(:,1:ord)
       !else
         call get_stabilizer(x,Spgr,order,ss_ptr,atr)
       !end if

       bet=reshape((/17.0, 7.0,3.0,  &
                     7.0,13.0,5.0,  &
                     3.0, 5.0,11.0/),(/3,3/))
       bett=bet
       if (order > 1 ) then
          do j=2,order
             Rsym=Spgr%SymOp(ss_ptr(j))%Rot
             Rs=real(Rsym)
             bett=bett+ matmul(Rs,matmul(bet,transpose(Rs)))
          end do
       end if
       Rsym=nint(1000.0*bett)
       codd=(/Rsym(1,1),Rsym(2,2),Rsym(3,3),Rsym(1,2),Rsym(1,3),Rsym(2,3)/)
       cdd=(/'a','b','c','d','e','f'/)
       multip=1.0
       !Search systematically all the possible constraints

       if(codd(1) == codd(2) .and. codd(1) == codd(3)) then ! a a a
         if(codd(4) == codd(5) .and. codd(4) == codd(6) ) then ! a a a d d d
           if(codd(4) == 0) then
             cdd=(/'a','a','a','0','0','0'/)     ! 1 A A A 0   0   0
             multip=(/1.0,1.0,1.0,0.0,0.0,0.0/)
             betas(4:6)=0.0
             betas(2:3)=betas(1)
             cod(2:3)=cod(1); cod(4:6)=0.0
           else
             cdd=(/'a','a','a','d','d','d'/)     ! 5 A A A D   D   D
             multip=(/1.0,1.0,1.0,1.0,1.0,1.0/)
             betas(5:6)=betas(4)
             betas(2:3)=betas(1)
             cod(2:3)=cod(1); cod(5:6)=cod(4)
           end if
         else if(codd(4) == -codd(5) .and. codd(4) == -codd(6) ) then !a a a d -d -d
           cdd=(/'a','a','a','d','d','d'/)       ! 6 A A A D  -D  -D
           multip=(/1.0,1.0,1.0,1.0,-1.0,-1.0/)
           betas(5:6)=-betas(4)
           betas(2:3)=betas(1)
           cod(2:3)=cod(1); cod(5:6)=cod(4)
         else if(codd(4) == -codd(5) .and. codd(4) ==  codd(6) ) then !a a a d -d  d
           cdd=(/'a','a','a','d','d','d'/)       ! 7 A A A D  -D   D
           multip=(/1.0,1.0,1.0,1.0,-1.0, 1.0/)
           betas(5)=-betas(4); betas(6)=betas(4)
           betas(2:3)=betas(1)
           cod(2:3)=cod(1); cod(5:6)= cod(4)
         else if(codd(4) ==  codd(5) .and. codd(4) == -codd(6) ) then !a a a d  d -d
           cdd=(/'a','a','a','d','d','d'/)       ! 8 A A A D   D  -D
           multip=(/1.0,1.0,1.0,1.0, 1.0,-1.0/)
           betas(6)=-betas(4); betas(5)=betas(4)
           betas(2:3)=betas(1)
           cod(2:3)=cod(1); cod(5:6)= cod(4)
         end if

       else if(codd(1) == codd(2)) then ! a a c
         if(codd(4) == codd(5) .and. codd(4) == codd(6) .and. codd(4) == 0) then ! a a c 0 0 0
             cdd=(/'a','a','c','0','0','0'/)     ! 2 A A C 0   0   0
             multip=(/1.0,1.0,1.0,0.0,0.0,0.0/)
             betas(4:6)=0.0
             betas(2)=betas(1)
             cod(2)=cod(1); cod(4:6)= 0.0
         else if(codd(5) == codd(6) .and. codd(5) == 0) then ! a a c x 0 0
             if(codd(4) == codd(1)/2) then
               cdd=(/'a','a','c','a','0','0'/)     ! 9 A A C A/2 0   0
               multip=(/1.0,1.0,1.0,0.5,0.0,0.0/)
               betas(5:6)=0.0; betas(4)=betas(1)*0.5
               betas(2)=betas(1)
               cod(2)=cod(1); cod(4)= cod(1); cod(5:6)=0.0
             else
               cdd=(/'a','a','c','d','0','0'/)     !11 A A C D   0   0
               multip=(/1.0,1.0,1.0,1.0,0.0,0.0/)
               betas(5:6)=0.0
               betas(2)=betas(1)
               cod(2)=cod(1); cod(5:6)=0.0
             end if
         else
             if(codd(5) == codd(6)) then  ! a a c d e e
               cdd=(/'a','a','c','d','e','e'/)     !20 A A C D   E   E
               multip=(/1.0,1.0,1.0,1.0,1.0,1.0/)
               betas(6)=betas(5)
               betas(2)=betas(1)
               cod(2)=cod(1); cod(6)=cod(5)
             else if(codd(5) == -codd(6)) then  ! a a c d e -e
               cdd=(/'a','a','c','d','e','e'/)     !19 A A C D   E  -E
               multip=(/1.0,1.0,1.0,1.0,1.0,-1.0/)
               betas(6)=-betas(5)
               betas(2)=betas(1)
               cod(2)=cod(1); cod(6)=cod(5)
             end if
         end if

       else if(codd(1) == codd(3)) then ! a b a
         if(codd(4) == codd(6)) then    ! a b a d x d
           if(codd(4) == 0) then  ! a b a 0 x 0
             if(codd(5) == 0) then ! a b a 0 0 0
               cdd=(/'a','b','a','0','0','0'/)     ! 3 A B A 0   0   0
               multip=(/1.0,1.0,1.0,0.0,0.0,0.0/)
               betas(4:6)=0.0
               betas(3)=betas(1)
               cod(3)=cod(1); cod(4:6)=0.0
             else                  ! a b a 0 e 0
               cdd=(/'a','b','a','0','e','0'/)     !12 A B A 0   E   0
               multip=(/1.0,1.0,1.0,0.0,1.0,0.0/)
               betas(4)=0.0;  betas(6)=0.0
               betas(3)=betas(1)
               cod(3)=cod(1); cod(4)=0.0;  cod(6)=0.0
             endif
           else  !! a b a d e d
             cdd=(/'a','b','a','d','e','d'/)       !22 A B A D   E   D
             multip=(/1.0,1.0,1.0,1.0,1.0,1.0/)
             betas(6)=betas(4)
             betas(3)=betas(1)
             cod(3)=cod(1); cod(6)=cod(4)
          end if

         else if(codd(4) == -codd(6)) then ! a b a d e -d
           cdd=(/'a','b','a','d','e','d'/)         !21 A B A D   E  -D
           multip=(/1.0,1.0,1.0,1.0,1.0,-1.0/)
           betas(6)=-betas(4)
           betas(3)=betas(1)
           cod(3)=cod(1); cod(6)=cod(4)
         end if

       else if(codd(2) == codd(3)) then ! a b b
         if(codd(4) == codd(5)) then    ! a b b d d x
           if(codd(4) == 0) then  ! a b b 0 0 x
             if(codd(6) == 0) then ! a b b 0 0 0
               cdd=(/'a','b','b','0','0','0'/)     ! 4 A B B 0   0   0
               multip=(/1.0,1.0,1.0,0.0,0.0,0.0/)
               betas(4:6)=0.0
               betas(3)=betas(2)
               cod(3)=cod(2); cod(4:6)=0.0
             else                  ! a b b 0 0 f
               cdd=(/'a','b','b','0','0','f'/)     !13 A B B 0   0   F
               multip=(/1.0,1.0,1.0,0.0,0.0,1.0/)
               betas(4:5)=0.0
               betas(3)=betas(2)
               cod(3)=cod(2); cod(4:5)=0.0
             endif
           else  !! a b b d d f
             cdd=(/'a','b','b','d','d','f'/)       !24 A B B D   D   F
             multip=(/1.0,1.0,1.0,1.0,1.0,1.0/)
             betas(5)=betas(4)
             betas(3)=betas(2)
             cod(3)=cod(2); cod(5)=cod(4)
           end if
         else if(codd(4) == -codd(5)) then ! a b b d -d e
           cdd=(/'a','b','b','d','d','f'/)         !23 A B B D  -D   F
           multip=(/1.0,1.0,1.0,1.0,-1.0,1.0/)
           betas(5)=-betas(4)
           betas(3)=betas(2)
           cod(3)=cod(2); cod(5)=cod(4)
         end if

       else !Now a /= b /= c

         if(codd(4) == codd(5) .and. codd(4) == 0) then ! a b c 0 0 x
           if(codd(6) == 0) then ! a b c 0 0 0
             cdd=(/'a','b','c','0','0','0'/)          !10 A B C 0   0   0
             multip=(/1.0,1.0,1.0,0.0,0.0,0.0/)
             betas(4:6)=0.0
             cod(4:6)=0.0
           else
             cdd=(/'a','b','c','0','0','f'/)          !18 A B C 0   0   F
             multip=(/1.0,1.0,1.0,0.0,0.0,1.0/)
             betas(4:5)=0.0
             cod(4:5)=0.0
           end  if
         else if(codd(5) == codd(6) .and. codd(5) == 0) then  ! a b c x 0 0
           if(codd(4) == codd(1)/2) then ! a b c a/2 0 0
             cdd=(/'a','b','c','a','0','0'/)          !15 A B C A/2 0   0
             multip=(/1.0,1.0,1.0,0.5,0.0,0.0/)
             betas(5:6)=0.0; betas(4)=betas(1)*0.5
             cod(4)=cod(1); cod(5:6)=0.0
           else if(codd(4) == codd(2)/2) then    !a b c b/2 0 0
             cdd=(/'a','b','c','b','0','0'/)          !14 A B C B/2 0   0
             multip=(/1.0,1.0,1.0,0.5,0.0,0.0/)
             betas(5:6)=0.0; betas(4)=betas(2)*0.5
             cod(4)=cod(2); cod(5:6)=0.0
           else
             cdd=(/'a','b','c','d','0','0'/)          !16 A B C D   0   0
             multip=(/1.0,1.0,1.0,1.0,0.0,0.0/)
             betas(5:6)=0.0
             cod(5:6)=0.0
           end  if
         else if(codd(4) == codd(6) .and. codd(4) == 0) then !a b c 0 e 0
           cdd=(/'a','b','c','0','e','0'/)            !17 A B C 0   E   0
           multip=(/1.0,1.0,1.0,0.0,1.0,0.0/)
           betas(4)=0.0; betas(6)=0.0
           cod(4)=0.0; cod(6)=0.0
         else if(codd(4) == codd(1)/2 .and. codd(5) == 0) then !a b c a/2 0 f
           cdd=(/'a','b','c','a','0','f'/)            !26 A B C A/2 0   F
           multip=(/1.0,1.0,1.0,0.5,0.0,1.0/)
           betas(4)=betas(1)*0.5; betas(5)=0.0
           cod(4)=cod(1); cod(5)=0.0
         else if(codd(4) == codd(2)/2 .and. codd(6) == 0) then !a b c b/2 e 0
           cdd=(/'a','b','c','b','e','0'/)            !27 A B C B/2 E   0
           multip=(/1.0,1.0,1.0,0.5,1.0,0.0/)
           betas(4)=betas(2)*0.5; betas(6)=0.0
           cod(4)=cod(2); cod(6)=0.0
         else if(codd(4) == codd(2)/2 .and. codd(5) == codd(6)/2) then !a b c b/2 f/2 f
           cdd=(/'a','b','c','b','f','f'/)            !25 A B C B/2 F/2 F
           multip=(/1.0,1.0,1.0,0.5,0.5,1.0/)
           betas(4)=betas(2)*0.5; betas(5)=betas(6)*0.5
           cod(4)=cod(2); cod(5)=cod(6)
         else if(codd(4) == codd(1)/2 .and. codd(6) == codd(5)/2) then !a b c a/2 e e/2
           cdd=(/'a','b','c','a','e','e'/)            !28 A B C A/2 E   E/2
           multip=(/1.0,1.0,1.0,0.5,1.0,0.5/)
           betas(4)=betas(1)*0.5; betas(6)=betas(5)*0.5
           cod(4)=cod(1); cod(6)=cod(5)
         else
           cdd=(/'a','b','c','d','e','f'/)            !29 A B C D   E   F
           multip=(/1.0,1.0,1.0,1.0,1.0,1.0/)
         end if
       end if

       do j=1,6
         if(abs(multi(j)) < epss .or. cdd(j) == '0' ) then
           codes(j) = 0.0_cp
         else
           codes(j) = sign(1.0_cp, multip(j))*(abs(cod(j))*10.0_cp + abs(multip(j)) )
         end if
       end do
       !if(present(Ipr)) then
       !  write(Ipr,'(a,6f10.4)')        '     Codes on Betas       : ',codes
       !  Write(Ipr,'(a,6(a,1x),6f7.3)') '     Codes and multipliers:  ',cdd,multip
       !  Write(Ipr,'(a)')               '     Beta_TOT matrix:  '
       !  do i=1,3
       !   Write(Ipr,'(a,3f12.4)')       '                      ',bett(i,:)
       !  end do
       !end if
       return
    End Subroutine Get_Atombet_Ctr

!==========================================================================================

