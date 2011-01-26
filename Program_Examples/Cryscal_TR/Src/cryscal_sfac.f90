!     Last change:  TR   17 Jul 2007    4:50 pm

!--------------------------------------------------------------------------
subroutine Calcul_SFAC_hkl
 USE Cryscal_module
 USE IO_module
 USE CFML_Reflections_Utilities,     ONLY : HKL_absent,  HKL_equiv, Reflection_list_type, Get_MaxNumRef, HKL_uni
 USE CFML_Structure_factors
 USE CFML_Crystallographic_Symmetry, ONLY : Get_Multip_Pos, set_spacegroup
 USE CFML_Crystal_Metrics,           ONLY : Set_Crystal_Cell
 

  implicit none
  integer                             :: i, n
  REAL, DIMENSION(3)                  :: HQ
  REAL                                :: Q_hkl, d_hkl, stl_hkl, angle_2theta, Z
  REAL                                :: a_star, b_star, c_star
  REAL                                :: alfa_rad, beta_rad, gama_rad
  REAL                                :: cos_alfa_star, cos_beta_star, cos_gama_star

  ! pour le calcul de facteur de structure
   type (Reflection_List_Type)        :: reflex_list_hkl
   TYPE(atom_list_type)               :: Atm
   REAL                               :: sf2
   real,dimension(3)                  :: vet
   real                               :: sn
   integer                            :: Num_ref
   real                               :: stl_min, stl_max
  !-------------------------------------------------

  if(nb_atom ==0) then
   call write_info(' No input atoms !!')
   call write_info('')
   return
  endif
  
  call set_spacegroup(space_group_symbol, SPG)  
  if(SPG%NumSpg == 0) then
   call write_info(' No defined space group !!')
   call write_info('')
  endif 
  
  

  IF(unit_cell%volume < 0.1) call volume_calculation('out')
  call set_crystal_Cell(unit_cell%param(1:3), unit_cell%param(4:6), crystal_cell)

  
  a_star = unit_cell%param(2)*unit_cell%param(3) * sind(unit_cell%param(4)) / unit_cell%Volume
  b_star = unit_cell%param(3)*unit_cell%param(1) * sind(unit_cell%param(5)) / unit_cell%Volume
  c_star = unit_cell%param(1)*unit_cell%param(2) * sind(unit_cell%param(6)) / unit_cell%Volume

  alfa_rad = unit_cell%param(4) * pi/180.
  beta_rad = unit_cell%param(5) * pi/180.
  gama_rad = unit_cell%param(6) * pi/180.

  cos_alfa_star = (cos(beta_rad) * cos(gama_rad) - cos(alfa_rad) ) / (sin(beta_rad) * sin(gama_rad))
  cos_beta_star = (cos(gama_rad) * cos(alfa_rad) - cos(beta_rad) ) / (sin(gama_rad) * sin(alfa_rad))
  cos_gama_star = (cos(alfa_rad) * cos(beta_rad) - cos(gama_rad) ) / (sin(alfa_rad) * sin(beta_rad))

 call write_info('')
 if(beam_type(1:8) == 'neutrons') then
  call write_info('    . Structure factor calculation (neutrons)')
 else
  call write_info('    . Structure factor calculation (X-rays)')
 endif
 call write_info('      ----------------------------------------')
 call write_info('')
 if(keyword_wave) then
  write(message_text, '(5x,a,F10.5)') '. Wavelength : ',wavelength
  call write_info(trim(message_text))
 endif
 

 ! construction de l'objet Atm
   Atm%natoms = nb_atom
   if(allocated(Atm%atom))  deallocate(Atm%atom)
   allocate (Atm%atom(nb_atom))
 
   do i=1, Atm%natoms
    Atm%atom(i)%Lab      = atom_label(i)
    Atm%atom(i)%ChemSymb = atom_type(i)
    Atm%atom(i)%SFACsymb = atom_type(i)
    Atm%atom(i)%Biso     = atom_Biso(i)
    Atm%atom(i)%occ      = atom_occ(i)
    Atm%atom(i)%x(1:3)   = atom_coord(1:3, i) 
    vet=Atm%atom(i)%x
    Atm%atom(i)%Mult=Get_Multip_Pos(vet,SpG)
   end do
  !------------------------------------------- 


 !IF(keyword_WAVE) then
 ! call write_info('          h     k     l            d_hkl(A)   stl_hkl(A-1)     Q_hkl(A-1)   2theta(deg.)')
 ! call write_info(' ')
 !else
 ! call write_info('          h     k     l            d_hkl(A)   stl_hkl(A-1)     Q_hkl(A-1)')
 ! call write_info(' ')
 !endif

  n = 1
  do i = 1, nb_hkl_SFAC_calc
    
    HQ(:) = int(H(: ,i))
    
  !  verif. que la reflexion est presente dans le groupe d'espace
    IF(hkl_absent(HQ, SPG)) then
     call write_info('')
     write (message_text, '(a,3I3,a)') ' Reflexion (', int(HQ(1:3)),') is absent in the current space group !'
     call write_info(trim(message_text))
     call write_info('')
     cycle
    end if

   
   Q_hkl = HQ(1)**2 * a_star**2 + HQ(2)**2 * b_star**2 + HQ(3)**2 * c_star**2              &
        + 2*HQ(1)*HQ(2) * a_star*b_star*cos_gama_star                                      &
        + 2*HQ(2)*HQ(3) * b_star*c_star*cos_alfa_star                                      &
        + 2*HQ(3)*HQ(1) * c_star*a_star*cos_beta_star

   d_hkl = 1/sqrt(Q_hkl)
   stl_hkl = 1/(2*d_hkl)
   Q_hkl = 4*pi*stl_hkl

   IF(keyword_WAVE) then
    Z = wavelength / (2 * d_hkl)
    IF (Z**2 < 1.) then
     angle_2theta = 2 * ATAN(Z / (SQRT(1. - Z **2))) * 180. / PI
    else
     angle_2theta = 180.
    endif
    angle_2theta = angle_2theta + shift_2theta
   !
   ! write(message_text, '(5x,3(1x,F5.2),5x,4(5x,F10.4))') HQ(1:3), d_hkl, stl_hkl, Q_hkl, angle_2theta
   !  call write_info(TRIM(message_text))
   ! 
   !else
   ! write(message_text, '(5x,3(1x,F5.2),5x,3(5x,F10.4))') HQ(1:3), d_hkl, stl_hkl, Q_hkl
   ! call write_info(TRIM(message_text))
   endif
  
   stl_min = stl_hkl * 0.99
   stl_max = stl_hkl * 1.01
   !num_ref = get_MaxNumRef( STL_max, crystal_cell%cellVol, STL_min, mult=SPG%multip) 
   num_ref = get_MaxNumRef( STL_max, crystal_cell%cellVol, STL_min) 
   
!write(*,*) ' vol : ', crystal_cell%cellVol
!write(*,*) ' stl : ', stl_min, stl_max
!write(*,*) ' numref : ', num_ref
!write(*,*) ' mult   : ', SPG%multip




   call HKL_uni(crystal_cell, SPG, .true., STL_min, STL_max, "s", num_ref, reflex_list_HKL)

!   if(beam_type(1:8) == 'neutrons') then
!    call Init_Structure_Factors(reflex_list_HKL,Atm,Spg,mode="NUC")
!   else
!    call Init_Structure_Factors(reflex_list_HKL,Atm,Spg,mode="XRA")
!   endif
!   call Init_calc_hkl_StrFactors(Atm)  ! << indispensable à calc_hkl_StrFactor
!   call Structure_Factors(Atm,SpG,reflex_list_HKL)
!
!  
!   sn =  reflex_list_HKL%ref(1)%s**2
!   
! !  write(*,*) ' sn = ', sn
!   if(beam_type(1:8) == 'neutrons') then
!    call Calc_hkl_StrFactor("S","N", int(HQ), sn,  Atm, SPG, sf2) 
!   else
!    call Calc_hkl_StrFactor("S","X", int(HQ), sn,  Atm, SPG, sf2)
!   endif

   sn = reflex_list_HKL%ref(1)%s**2
   if(beam_type(1:8) == "neutrons") then
    call Init_Structure_Factors(reflex_list_HKL, Atm, Spg, mode = "NUC") 
    call Structure_Factors(Atm, SPG, reflex_list_HKL, mode="NUC")
    call Init_calc_hkl_StrFactors(Atm, mode ="NUC")     
    call Calc_hkl_StrFactor("S", "N", int(HQ), sn, Atm, Spg, sf2)
    
   else   ! X-rays
    call Init_Structure_Factors(reflex_list_HKL, Atm, Spg) 
    call Structure_Factors(Atm, SPG, reflex_list_HKL)
    call Init_calc_hkl_StrFactors(Atm)     
    call Calc_hkl_StrFactor("S", "X", int(HQ), sn, Atm, Spg, sf2)
   
   endif  
    
   if(i == 1) then
   call write_info('')
   if(keyword_WAVE) then
    call write_info('           H   K   L   Mult      2Theta     SinTh/Lda         dspc         |Fc|         Phase          F-Real        F-Imag       |Fc|^2')
   else   
    call write_info('           H   K   L   Mult   SinTh/Lda          dspc         |Fc|         Phase          F-Real        F-Imag       |Fc|^2')
   endif
   call write_info('')
   endif
    do n=1, reflex_list_HKL%Nref
	 if(reflex_list_HKL%ref(n)%h(1) /= HQ(1) .and. &
	    reflex_list_HKL%ref(n)%h(2) /= HQ(2) .and. &
            reflex_list_HKL%ref(n)%h(3) /= HQ(3) ) then
	  cycle	
     endif		
     if(keyword_WAVE) then
      write(message_text,fmt="(i5,3x, 3i4,i5,7f14.5,f14.3)")  &
		     i, reflex_list_HKL%ref(n)%h,     reflex_list_HKL%ref(n)%mult,  angle_2theta,        &
                reflex_list_HKL%ref(n)%S,     0.5/reflex_list_HKL%ref(n)%S,   &
		        reflex_list_HKL%ref(n)%Fc,    reflex_list_HKL%ref(n)%Phase,   &
                reflex_list_HKL%ref(n)%a,     reflex_list_HKL%ref(n)%b,       &
			reflex_list_HKL%ref(n)%Fc*reflex_list_HKL%ref(n)%Fc
     else
      write(message_text,fmt="(i5,3x, 3i4,i5,6f14.5,f14.3)")  &
		     i, reflex_list_HKL%ref(n)%h,     reflex_list_HKL%ref(n)%mult,    &
                reflex_list_HKL%ref(n)%S,     0.5/reflex_list_HKL%ref(n)%S,   &
		        reflex_list_HKL%ref(n)%Fc,    reflex_list_HKL%ref(n)%Phase,   &
                reflex_list_HKL%ref(n)%a,     reflex_list_HKL%ref(n)%b,       &
			reflex_list_HKL%ref(n)%Fc*reflex_list_HKL%ref(n)%Fc
     endif
     call write_info(trim(message_text)) 				
    end do
  end do
 
  call write_info('')
 
 return
end subroutine Calcul_SFAC_hkl


!----------------------------------------------------------

subroutine generate_HKL()
 USE CFML_Reflections_Utilities,  ONLY : Get_MaxNumRef, HKL_gen, reflect_type, Reflection_list_type, reflection_type, HKL_uni
 USE CFML_Structure_factors
 Use CFML_Crystallographic_Symmetry, only: Get_Multip_Pos, Write_SpaceGroup
 USE CFML_Crystal_Metrics,        ONLY : Set_Crystal_Cell, Write_Crystal_Cell
 !USE CFML_Atom_TypeDef,          ONLY : atom_list_type
 !use Structure_Factor_Module,    ONLY : Calc_StrFactor
 !USE CFML_Structure_Factors,      ONLY : Structure_factors
 USE CFML_Math_General,           ONLY : sort
 USE cryscal_module
 USE IO_module

 implicit none
  INTEGER                                       :: i, ord
  LOGICAL                                       :: Friedel
  INTEGER                                       :: Num_ref
  TYPE(Reflect_type),         ALLOCATABLE, DIMENSION(:) :: reflex_HKL
  !TYPE(Reflection_list_type), ALLOCATABLE, DIMENSION(:) :: reflex_list_HKL
  type (Reflection_List_Type) :: reflex_list_hkl

  TYPE(Reflection_type),      ALLOCATABLE, DIMENSION(:) :: ref

  INTEGER,            ALLOCATABLE, DIMENSION(:) :: ordered_array
  REAL                                          :: STL_min, STL_max
  REAL                                          :: d_hkl, Z, angle_2theta, angle_theta, Q_hkl

  TYPE(atom_list_type)                          :: Atm
  real                                          :: sn
  REAL                                          :: sf2
  REAL, ALLOCATABLE, DIMENSION(:)               :: sf_2 ! (structure factor)**2
  complex                                       :: fc
  real,dimension(3)                             :: vet



  Friedel = .true.


 IF(keyword_WAVE) then
  call write_info('')
  WRITE(message_text, '(a,F8.5,a)')'    . hkl list generation (l=',wavelength,' A)'
  call write_info(TRIM(message_text))
  call write_info('      ---------------------------------')
  call write_info('')
 else
  call write_info('')
  call write_info('    . hkl list generation')
  call write_info('      -------------------')
  call write_info('')
 endif


 IF(.NOT. keyword_SPGR) then
  call write_info('')
  call write_info(' SPGR parameter (space group) is a mandatory parameter for this calculation !')
  call write_info('')
  return
 endif
 IF(.NOT. keyword_CELL) then
  call write_info('')
  call write_info(' CELL parameters (in A) are mandatory parameters for this calculation !')
  call write_info('')
  return
 endif

 IF((HKL_2THETA .or. HKL_THETA) .and. .not. keyword_WAVE) then
  call write_info('')
  call write_info(' WAVE parameter (wavelength in A) is a mandatory parameter for this calculation !')
  call write_info('')
  return
 endif

 IF(HKL_STL) then
  STL_min = X_min
  STL_max = X_max
 elseif(HKL_Q) then
  STL_min = X_min/(4.*PI)
  STL_max = X_max/(4.*PI)
 ELSEIF(HKL_D) then
  STL_min = 0.5/X_max
  STL_max = 0.5/X_min
 ELSEIF(HKL_2THETA) then
  STL_min  = SIN(X_min  * PI / 360) / wavelength
  STL_max  = SIN(X_max  * PI / 360) / wavelength
 ELSEIF(HKL_THETA) then
  STL_min  = SIN(X_min  * PI / 180) / wavelength
  STL_max  = SIN(X_max  * PI / 180) / wavelength
 endif



 call set_crystal_Cell(unit_cell%param(1:3), unit_cell%param(4:6), crystal_cell)

 ! estimation du nombre de reflexions dans un domaine donne de stl
 ! estimation necessaire pour dimensionner les tableaux
 !Num_ref  =  Get_MaxNumRef(STL_max, Unit_cell%volume, STL_min, SPG%Multip )
 Num_ref  =  Get_MaxNumRef(STL_max, Unit_cell%volume, STL_min)
 
! write(*,*) ' Stl_min stl_max : ', stl_min, stl_max
! write(*,*) ' Volume  : ', unit_cell%volume
! write(*,*) ' Num_ref : ', Num_Ref
! pause
 IF(ALLOCATED(reflex_HKL))        DEALLOCATE(reflex_HKL)
 !IF(ALLOCATED(reflex_list_HKL))   DEALLOCATE(reflex_list_HKL)
 IF(ALLOCATED(ref))               DEALLOCATE(ref)
 ALLOCATE(reflex_HKL(Num_ref))
 !ALLOCATE(reflex_list_HKL(Num_ref))
 ALLOCATE(ref(Num_ref))

 IF(ALLOCATED(ordered_array)) DEALLOCATE(ordered_array)
 ALLOCATE(ordered_array(Num_ref))

 IF(SPG%centred /= 2) friedel=.false.

 call HKL_gen(crystal_cell, SPG, Friedel, STL_min, STL_max, Num_ref, reflex_HKL )

 call write_info(' ')
 if(HKL_STL) then
  write(message_text, '(a,F8.5,a,F8.5,a,I6)') '   Nb of reflections in the ',X_min, ' - ', X_max,  &
                                             ' SinTheta/lambda (A-1) range: ' , Num_ref

 elseif(HKL_Q) then
  write(message_text, '(a,F8.5,a,F8.5,a,I6)') '   Nb of reflections in the ',X_min, ' - ', X_max,  &
                                             ' Q (A-1) range: ' , Num_ref
 elseif(HKL_D) then
  write(message_text, '(a,F8.5,a,F8.5,a,I6)') '   Nb of reflections in the ',X_min, ' - ', X_max,  &
                                             ' d_hkl(A) range: ' , Num_ref

 elseif(HKL_2theta) then
  write(message_text, '(a,F7.2,a,F7.2,a,I6)') '   Nb of reflections in the ',X_min, ' - ', X_max,  &
                                             ' 2Theta (deg.) range: ' , Num_ref

 elseif(HKL_theta) then
  write(message_text, '(a,F7.2,a,F7.2,a,I6)') '   Nb of reflections in the ',X_min, ' - ', X_max,  &
                                             ' Theta (deg.) range: ' , Num_ref
 endif
 call write_info(trim(message_text))
 call write_info(' ')


 if(Num_ref ==0) return


 call write_info('')
 WRITE(message_text, '(a,2i4)') '  h range: ', MINVAL(reflex_HKL(1:Num_ref)%h(1)), MAXVAL(reflex_HKL(1:Num_ref)%h(1))
 call write_info(TRIM(message_text))
 WRITE(message_text, '(a,2i4)') '  k range: ', MINVAL(reflex_HKL(1:Num_ref)%h(2)), MAXVAL(reflex_HKL(1:Num_ref)%h(2))
 call write_info(TRIM(message_text))
 WRITE(message_text, '(a,2i4)') '  l range: ', MINVAL(reflex_HKL(1:Num_ref)%h(3)), MAXVAL(reflex_HKL(1:Num_ref)%h(3))
 call write_info(TRIM(message_text))
 call write_info('')

 if(.not. write_HKL) then
  DEALLOCATE (reflex_HKL)
  deallocate (ordered_array)
  return
 endif


! call sort(Num_ref,reflex_HKL%S, ordered_array )
 call sort(reflex_HKL%S, Num_ref, ordered_array )

 ! calcul du facteur de structure -----------------------------
  if (nb_atom /=0 .and. keyword_WAVE) then
   ! 
        
   ! construction de l'objet Atm
   Atm%natoms = nb_atom
   if(allocated(Atm%atom))  deallocate(Atm%atom)
   allocate (Atm%atom(nb_atom))
 
   do i=1, Atm%natoms
    Atm%atom(i)%Lab      = atom_label(i)
    Atm%atom(i)%ChemSymb = atom_type(i)
    Atm%atom(i)%SFACsymb = atom_type(i)
    Atm%atom(i)%Biso     = atom_Biso(i)
    Atm%atom(i)%occ      = atom_occ(i)
    Atm%atom(i)%x(1:3)   = atom_coord(1:3, i) 
    vet=Atm%atom(i)%x
    Atm%atom(i)%Mult=Get_Multip_Pos(vet,SpG)
   end do

   
 
  
   
   !call Write_Crystal_Cell(crystal_cell,1)
   !call Write_SpaceGroup(SpG,1)
   
   if (ALLOCATED(sf_2)) DEALLOCATE(sf_2)
   ALLOCATE(sf_2(nb_atom))
   
   ! subroutine HKL_uni(CELL_objet, SPG_objet, Friedel, stl1, stl2, "r": d_spacing, num_ref, reflex_list_HKL)
   !num_ref = get_MaxNumRef( STL_max, crystal_cell%cellVol, mult=SPG%multip)
   !write(*,*) 'STLmax : ', STL_max
   !write(*,*)' VOL: ', crystal_cell%cellVOL
   !write(*,*) ' num_ref = ', num_ref 
   call HKL_uni(crystal_cell, SPG, .true., STL_min, STL_max, "s", num_ref, reflex_list_HKL)

   ! chargement des facteurs de diffusion
   ! call Init_structure_factors (HKL_objet, Atoms_objet, SPG_objet, mode, lun)
   ! mode = XRA : Xrays
   ! mode = NUC : neutrons
   ! lun (optionel) : unite de sortie
   if(beam_type(1:8) == 'neutrons') then    
    call Init_Structure_Factors(reflex_list_HKL,Atm,Spg,mode="NUC")
    call Structure_Factors(Atm,SpG,reflex_list_HKL, mode="NUC")
   else    
    call Init_Structure_Factors(reflex_list_HKL,Atm,Spg,mode="XRA")
    call Structure_Factors(Atm,SpG,reflex_list_HKL)
   endif
   
   
   if(HKL_2theta) then 
    call write_info('           H   K   L   Mult        2Theta     SinTh/Lda          dspc          |Fc|         Phase        F-Real        F-Imag        |Fc|^2')
   else 
    call write_info('           H   K   L   Mult     SinTh/Lda          dspc          |Fc|         Phase        F-Real        F-Imag        |Fc|^2')
   endif
   call write_info('')
   
   do i=1,reflex_list_HKL%Nref
    if(HKL_2theta) then  !calcul du 2theta
     Z = wavelength * reflex_list_HKL%ref(i)%S
     IF (Z**2 < 1.) then
      angle_2theta = 2 * ATAN(Z / (SQRT(1. - Z **2))) * 180. / PI
     else
      angle_2theta = 180.
     endif
     angle_2theta = angle_2theta + shift_2theta
     write(message_text,fmt="(i5,3x, 3i4,2x,i5,7f14.5,f14.3)")  &
			         i, reflex_list_HKL%ref(i)%h,     reflex_list_HKL%ref(i)%mult,    angle_2theta, &
					    reflex_list_HKL%ref(i)%S,     0.5/reflex_list_HKL%ref(i)%S,   &
			            reflex_list_HKL%ref(i)%Fc,    reflex_list_HKL%ref(i)%Phase,   &
                        reflex_list_HKL%ref(i)%a,     reflex_list_HKL%ref(i)%b,       &
			            reflex_list_HKL%ref(i)%Fc*reflex_list_HKL%ref(i)%Fc

    else
     write(message_text,fmt="(i5,3x, 3i4,2x,i5,6f14.5,f14.3)")  &
			         i, reflex_list_HKL%ref(i)%h,     reflex_list_HKL%ref(i)%mult,    &
                        reflex_list_HKL%ref(i)%S,     0.5/reflex_list_HKL%ref(i)%S,   &
			            reflex_list_HKL%ref(i)%Fc,    reflex_list_HKL%ref(i)%Phase,   &
                        reflex_list_HKL%ref(i)%a,     reflex_list_HKL%ref(i)%b,       &
			            reflex_list_HKL%ref(i)%Fc*reflex_list_HKL%ref(i)%Fc
    endif			            
    call write_info(trim(message_text)) 				
   end do
   call write_info('')

   !do i=1, Num_ref
    !ord = ordered_array(i)
    !sn  = reflex_HKL(ord)%S**2     ! (sinTheta/lambda)**2	 
!    write(*,*) ord, sn
	!pause
    !sn = reflex_list_HKL%ref(i)%s**2
	!call Calc_StrFactor('S','X',i,sn,Atm,SPG,sf2,fc=fc)  ! 'S' pour single crystal 
    !sf_2(i) = sf2
    !WRITE(message_text,*) reflex_HKL(ord)%h(1), reflex_HKL(ord)%h(2), reflex_HKL(ord)%h(3), &
    !                         reflex_list_HKL%ref(ord)%Fc, sf_2(ord)
    ! call write_info(trim(message_text))

	!WRITE(message_text,*) reflex_list_HKL%ref(i)%h(1), reflex_list_HKL%ref(i)%h(2), &
	!                      reflex_list_HKL%ref(i)%h(3), &
    !                      reflex_list_HKL%ref(i)%Fc, sf_2(i)
							 
    !call write_info(trim(message_text))
	!call write_info('')
   !END do
   
      !   do i=1, hkl%nref
      !   sn=hkl%ref(i)%s * hkl%ref(i)%s
      !   call Calc_StrFactor("S","X",i,sn,A,Spg,sf2,fc=fc)
      !   write(unit=lun,fmt="(3i4,i5,5f12.5,i8,f12.5)") hkl%ref(i)%h, hkl%ref(i)%mult, &
      !                           hkl%ref(i)%S, hkl%ref(i)%Fc, hkl%ref(i)%Phase,   &
      !                           real(fc), aimag(fc), i, sqrt(sf2)
      ! end do

 !
 !  call structure_factors(Atm, SPG, reflex_list_HKL,'X-rays', wavelength)
 !
   return
  endif
 !------------------------------------------------------------

 IF(HKL_STL .or. HKL_Q .or. HKL_D) then
  if(.not. keyword_WAVE) then
   call write_info('   num    h   k   l     m           STL      Q(A-1)        d(A)')
   call write_info('')
  else
   call write_info('   num    h   k   l     m           STL      Q(A-1)        d(A)   2theta(deg).')
   call write_info('')
  endif


  do i=1,Num_ref
   ord = ordered_array(i)
   d_hkl = 0.5/reflex_HKL(ord)%S
   Q_hkl = 4.*pi*reflex_HKL(ord)%S
   !IF(keyword_STRUCT_factor) call Calc_StrFactor()
   if(keyword_WAVE) then
    call calcul_2theta(reflex_HKL(ord)%S, angle_2theta)
    WRITE(message_text,'(i6, 1x,3I4,I6,6x,F8.5,4x,F8.5,4x,F8.5,8x,F7.3)') i,reflex_HKL(ord)%h(1), reflex_HKL(ord)%h(2), reflex_HKL(ord)%h(3),  &
                                                                   reflex_HKL(ord)%mult, reflex_HKL(ord)%S,  Q_hkl, d_hkl, angle_2theta
   else
    WRITE(message_text,'(i6, 1x,3I4,I6,6x,F8.5,4x,F8.5,4x,F8.5)') i,reflex_HKL(i)%h(1), reflex_HKL(i)%h(2), reflex_HKL(i)%h(3),  reflex_HKL(i)%mult,  &
                                                                  reflex_HKL(i)%S, Q_hkl, d_hkl
   endif
   call write_info(TRIM(message_text))
  end do



 elseif(HKL_2THETA) then
  call write_info('   num    h   k   l     m           STL      Q(A-1)        d(A)   2theta(deg).')
  call write_info('')

  !call sort(Num_ref,reflex_HKL%S, ordered_array )
  call sort(reflex_HKL%S, Num_ref, ordered_array )
  do i=1, num_ref
   ord = ordered_array(i)
   d_hkl = 0.5/reflex_HKL(ord)%S
   Q_hkl = 4.*pi*reflex_HKL(ord)%S

   call calcul_2theta(reflex_HKL(ord)%S, angle_2theta)

   WRITE(message_text,'(i6, 1x,3I4,I6,6x,F8.5,4x,F8.5,4x,F8.5,8x,F7.3)') i, reflex_HKL(ord)%h(1), reflex_HKL(ord)%h(2), reflex_HKL(ord)%h(3),  &
                                                                   reflex_HKL(ord)%mult, reflex_HKL(ord)%S, Q_hkl, d_hkl, angle_2theta
   call write_info(TRIM(message_text))
  end do

 elseif(HKL_THETA) then

  call write_info('   num    h   k   l     m           STL      Q(A-1)        d(A)    theta(deg).')
  call write_info('')

  do i=1,Num_ref
   ord = ordered_array(i)
   d_hkl = 0.5/reflex_HKL(ord)%S
   Q_hkl = 4.*pi*reflex_HKL(ord)%S

   call calcul_2theta(reflex_HKL(ord)%S, angle_2theta)
   angle_theta = 0.5 * angle_2theta

   WRITE(message_text,'(i6, 1x,3I4,I6,6x,F8.5,4x,F8.5,4x,F8.5,8x,F7.3)') i,reflex_HKL(ord)%h(1), reflex_HKL(ord)%h(2), reflex_HKL(ord)%h(3),  &
                                                                  reflex_HKL(ord)%mult, reflex_HKL(ord)%S, Q_hkl, d_hkl, angle_theta
   call write_info(TRIM(message_text))
  end do


 endif

 DEALLOCATE (reflex_HKL)
 deallocate (ordered_array)

 return
end subroutine generate_HKL


