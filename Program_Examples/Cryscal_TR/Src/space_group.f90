!     Last change:  TR   14 Feb 2007    2:20 pm
subroutine space_group_info
 use cryscal_module,                 ONLY : ON_SCREEN, space_group_symbol, symm_op_string,        &
                                            keyword_SYMM, symm_op_xyz, symm_op_mat, nb_symm_op,   &
                                            WRITE_SPG_info, Write_SPG_info_all, WRITE_SPG_EXTI,   &
                                            write_SPG_subgroups, message_text,                    &
                                            SPG, keyword_create_CIF,                              &
                                            keyword_CELL, unit_cell, CIF_cell_measurement,        &
                                            CIF_diffrn_reflns, input_PCR, known_cell_esd, keyword_read_INS
 
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
!  call system('copy cryscal.log + sg.out')
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
  open (UNIT=2, FILE='cryscal.log', ACTION="write",position = "append")
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
  CLOSE(UNIT=3)
  open (unit=3, file="sg.out", status="replace")
  call search_extinctions(SPG, 3)
  close (UNIT=3)
!  call system('echo off')
!  call system('type sg.out')
!  close(unit=2)
!  call system('copy cryscal.log + sg.out')
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
  open (UNIT=2, FILE='CRYSCAL.log', ACTION="write",position = "append")
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
   IF(unit_cell%volume < 0.1) call volume_calculation('no_out')
   IF(known_CELL_esd) then
    call write_CIF_file('CELL_PARAM_ESD')
   else
    call write_CIF_file('CELL_PARAM')
   endif

   ! CIF_Thmin, CIF_Thmax
    call write_CIF_file("CIF_THMIN,THMAX")
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


 return
end subroutine space_group_info

!-----------------------------------------------------------------------------------
subroutine get_SITE_info()
 USE cryscal_module 
 USE IO_module,                      ONLY : write_info
 USE CFML_crystallographic_symmetry, ONLY : set_spacegroup, Get_orbit, Get_Multip_Pos, &
                                            symmetry_symbol, Wyckoff_Type
 USE macros_module,                  ONLY : u_case
 USE CFML_Math_General,              ONLY : acosd
 !USE CFML_constants,                 ONLY : sp 
 USE CFML_GlobalDeps,                 ONLY : sp

 implicit none
 ! local variables
  INTEGER                           :: i, j, k
  REAL , DIMENSION(3)               :: r
  INTEGER                           :: site_multiplicity
  INTEGER                           :: Wyck_mult

  type(wyckoff_type)   :: wyckoff

  real , dimension(3,192)           :: orbit
  CHARACTER (LEN=80)                :: text80
  integer, dimension(192)           :: effsym  ! Pointer to effective symops
  integer, dimension(192)           :: ss_ptr  ! Pointer to stabilizer

  LOGICAL, DIMENSION(nb_atom)       :: write_site_info_label
  
  ! pour determiner les constraintes sur les Bij 
  REAL , DIMENSION(6)               :: beta_ij 
  integer                           :: codini
  real(kind=sp), dimension(6)       :: codes  !codewords for positions

  


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

!   SITE_multiplicity = Get_Multip_Pos(r, SPG)
!   WRITE(message_text,'(a,I5)')     '     . site multiplicity: ', SITE_multiplicity
!   call write_info(trim(message_text))


 !  call Get_Orbit(r, SPG, Wyck_mult, orbit, effsym, ss_ptr )

   call Get_Orbit(r, SPG, wyckoff%num_orbit, orbit, effsym, ss_ptr )
!   SITE_multiplicity = Get_Multip_Pos(r, SPG)
   WRITE(message_text,'(a,I5)')     '     . site multiplicity: ', wyckoff%num_orbit
   call write_info(trim(message_text))


   call write_info('     . list of symmetry operators and symmetry elements of the site point group:')


   do j=1,  wyckoff%num_orbit
    if(ss_ptr(j) == 0) cycle
    call symmetry_symbol(SPG%SymOp(ss_ptr(j)),text80)
    write(message_text,'(7x, a,i3,a,a,a)') ' Operator ',j,': ', SPG%SymopSymb(ss_ptr(j)),trim(text80)
    call write_info(trim(message_text))
   END do


   ! orbit
    call write_info(' ')
    call write_info('     --> Complete orbit of the atom: ')
    call write_info(' ')

   do j=1,  wyckoff%num_orbit
    if(j <10) then
     write(message_text,"(10x,a6,a,i1,a,3f9.5,tr4,a)") trim(atom_label(i)),'_',j,":   ",orbit(:,j),trim(SPG%SymopSymb(effsym(j)))
    elseif(j < 100) then
     write(message_text,"(10x,a6,a,i2,a,3f9.5,tr4,a)") trim(atom_label(i)),'_',j,":  ", orbit(:,j),trim(SPG%SymopSymb(effsym(j)))
    else
     write(message_text,"(10x,a6,a,i3,a,3f9.5,tr4,a)") trim(atom_label(i)),'_',j,": ",  orbit(:,j),trim(SPG%SymopSymb(effsym(j)))
    endif
    call write_info(trim(message_text))
   END do
   
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

 return

end subroutine get_site_info
!---------------------------------------------------------------------------

subroutine write_symm_op_reduced_set()
 USE IO_module,                      ONLY : write_info
 use CFML_crystallographic_symmetry, ONLY : set_spacegroup, searchop, write_sym
 USE cryscal_module,                 ONLY : ON_SCREEN, SPG, message_text
 USE CFML_symmetry_tables,           ONLY : intsymoh,x_oh, intsymd6h,x_d6h

 implicit none
  CHARACTER(LEN=3)                :: string_numor
  INTEGER                         :: i, j, i1, i2
  INTEGER, DIMENSION(192)         :: indx_op


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
!!----  subroutine get_atombet_ctr(x,betas,Spgr,codini,codes,ord,ss,debug)
!!----     real(kind=sp), dimension(3),     intent(in    ) :: x      !Atom position (fractional coordinates)
!!----     real(kind=sp), dimension(6),     intent(in out) :: betas  !Anisotropic temperature factors
!!----     type(Space_Group_type),          intent(in    ) :: Spgr   !Space Group
!!----     Integer,                         intent(in out) :: codini !Last attributed parameter
!!----     real(kind=sp), dimension(6),     intent(in out) :: codes  !codewords for positions
!!----     integer,               optional, intent(in    ) :: ord    !Order of the stabilizer
!!----     integer, dimension(:), optional, intent(in    ) :: ss     !Pointer to SymmOp. of stabilizer
!!----     integer,               optional, intent(in    ) :: debug  !Debug variable
!!----
!!----  Subroutine to get the appropriate constraints in the refinement codes of
!!----  anisotropic atomic displacement(thermal) parameters.
!!----
!!----  Uses: Space_Group_type, get_stabilizer, sym_b_relations
!!----
!!----  Updated: 1 April 2002
!!----
!!
!!    subroutine get_atombet_ctr(x,betas,Spgr,codini,codes,ord,ss,debug)
    subroutine get_atombet_ctr(x,betas,Spgr,codini,codes)
     USE CFML_Math_General,              ONLY  : acosd
     !USE CFML_constants,                 ONLY      : sp 
     USE CFML_GlobalDeps,                 ONLY : sp
     USE CFML_crystallographic_symmetry, ONLY  : Space_group_type, get_stabilizer, sym_b_relations
    
     implicit none
       real(kind=sp), dimension(3),     intent(in    ) :: x
       real(kind=sp), dimension(6),     intent(in out) :: betas
       !real, dimension(3),     intent(in    ) :: x
       !real, dimension(6),     intent(in out) :: betas
       
       type(Space_Group_type),          intent(in    ) :: Spgr
       Integer,                         intent(in out) :: codini
       real(kind=sp), dimension(6),     intent(in out) :: codes
       !integer,               optional, intent(in    ) :: ord, debug
       !integer, dimension(:), optional, intent(in    ) :: ss
       ! Local variables
       character (len=1), dimension(6)   :: cdd
       real(kind=sp),     dimension(6)   :: multip
       integer,           dimension(6,48):: p_mul=0,p_cod=0
       integer :: i,j,order
       real(kind=sp)    :: suma
       integer,           dimension(48) :: ss_ptr
       character (len=1), dimension(6)  :: strc
       character (len=1)                :: car
       integer,           dimension(6)  :: codd
       real(kind=sp),     dimension(6)  :: cod,mul,multi
       real(kind=sp),     parameter     :: epss=0.01

       suma=0.0
       do j=1,6
          suma=suma+abs(codes(j))
          cod(j)=int(abs(codes(j))/10.0)             !Input Parameter number with sign
          multi(j)=mod(codes(j),10.0)                !Input Multipliers
          if(cod(j) < 1.0 .and. abs(multi(j)) > epss)  then
               codini=codini+1
               cod(j) = real(codini)
          end if
       end do
       if(suma < epss) return  !No refinement is required

       !if(present(ord) .and. present(ss)) then
       !  order=ord
       !  ss_ptr(1:order) = ss(1:ord)
       !else
         call get_stabilizer(x,Spgr,order,ss_ptr)
       !end if

       strc=(/'a','b','c','d','e','f'/)
       cdd = strc
!      multip=(/1.0,1.0,1.0,1.0,1.0,1.0/)
       multip=multi
       if (order > 1 ) then
          do j=2,order
             call sym_b_relations(Spgr%SymopSymb(ss_ptr(j)),codd,mul)
             p_mul(:,j)= nint(mul*10.0)
             p_cod(:,j)=codd
             do i=1,6
                if (abs(mul(i)) <= epss) then
                   cdd(i) = "0"
                   multip(i)=0.0
                else
                    if(cdd(i) /= "0") then
                     cdd(i) = strc(codd(i))
                     multip(i)=mul(i)
                    end if
                end if
             end do
             strc=cdd
          end do

          if(cdd(4) == cdd(5)) then
              i=0
              ss_ptr=0
              do j=2,order
                if(p_cod(4,j) /= p_cod(5,j)) cycle
                i=i+1
                ss_ptr(i)=10*p_mul(4,j)+ p_mul(5,j)
              end do
              if(i > 1) then
                do j=2,i
                  if(ss_ptr(1) /= ss_ptr(j)) then  !incompatible constraints => multip should be zero
                    multip(4)=0.0
                    multip(5)=0.0
                  end if
                end do
              end if
          end if

          if(cdd(5) == cdd(6)) then
              i=0
              ss_ptr=0
              do j=2,order
                if(p_cod(5,j) /= p_cod(6,j)) cycle
                i=i+1
                ss_ptr(i)=10*p_mul(5,j)+ p_mul(6,j)
              end do
              if(i > 1) then
                do j=2,i
                  if(ss_ptr(1) /= ss_ptr(j)) then  !incompatible constraints => multip should be zero
                    multip(5)=0.0
                    multip(6)=0.0
                  end if
                end do
              end if
          end if

          do j=1,5
           do i=j+1,6
             if( cdd(j) == cdd(i)) then
                if( abs(multip(j)) < epss .or. abs(multip(i)) < epss) then
                  cod(j) =0.0
                  cod(i) =0.0
                  cdd(j) = "0"
                  cdd(i) = "0"
                  multip(i) = 0.0
                  multip(j) = 0.0
                end if
             end if
           end do
          end do

          car=cdd(1)
          do j=2,6
             if (cdd(j) == "0") then
                cod(j) = 0.0
                betas(j)= 0.0
             end if
             if (cdd(j) == car) then
                cod(j) = cod(1)
                betas(j)= betas(1)*multip(j)
             end if
             if (cdd(j) == "a") then
                cod(j) = cod(1)
                betas(j)= betas(1)*multip(j)
             end if
             if (cdd(j) == "b") then
                cod(j) = cod(2)
                betas(j)= betas(2)*multip(j)
             end if
             if (cdd(j) == "c") then
                cod(j) = cod(3)
                betas(j)= betas(3)*multip(j)
             end if
             if (cdd(j) == "d") then
                cod(j) = cod(4)
                betas(j)= betas(4)*multip(j)
             end if
             if (cdd(j) == "e") then
                cod(j) = cod(5)
                betas(j)= betas(5)*multip(j)
             end if
             if (cdd(j) == "f") then
                cod(j) = cod(6)
                betas(j)= betas(6)*multip(j)
             end if
          end do
       end if

        do j=1,6
          if(abs(multi(j)) < epss .or. cdd(j) == '0' ) then
            codes(j) = 0.0
          else
            codes(j) = sign(1.0, multip(j))*(abs(cod(j))*10.0 + abs(multip(j)) )
          end if
        end do
        !if(present(debug)) then
        ! write(98,'(a,6f10.4)')        '     Codes on Betas       : ',codes
        ! WRITE(98,'(a,6(a,1x),6f7.3)') '     Codes and multipliers:  ',cdd,multip
        !end if
        return
end subroutine get_atombet_ctr

!==========================================================================================