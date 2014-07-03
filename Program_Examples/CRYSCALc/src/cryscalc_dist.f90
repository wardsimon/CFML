!     Last change:  TR   17 Jul 2007    4:50 pm

subroutine  calcul_distances
  use cryscalc_module, ONLY      : nb_dist_calc, new_coord, atom1_dist, atom2_dist, SP_value, nb_atom, message_text, &
                                   keyword_DIST_X, dist_coef, keyword_DIST_plus, dist_plus,                          &
								   keyword_DHA, dist_AH, debug_proc
  USE IO_module

  implicit none
  integer                 :: i
  real                    :: distance
  REAL, DIMENSION(3)      :: atom_1_coord, atom_2_coord
  REAL, DIMENSION(3)      :: delta
  CHARACTER(LEN=16)       :: label_1, label_2
  LOGICAL                 :: ok

  if(debug_proc%level_2)  call write_debug_proc_level(2, "calcul_distances")
  
  call write_info('')
  call write_info('    . DISTANCES CALCULATIONS')
  call write_info('    ------------------------')
  call write_info('')

  IF(nb_atom == 0) then
   call write_info('')
   WRITE(message_text, '(a)') ' ... ATOM list not defined ...'
   call write_info(TRIM(message_text))
   call write_info('')
   return
  endif

  do i = 1, nb_dist_calc
   label_1 = atom1_dist(i)
   label_2 = atom2_dist(i)

   ok = .false.
   call get_label_atom_coord('dist', i,1, ok)
   IF(.NOT. ok) cycle
   atom_1_coord(1:3) = new_coord(1:3)

   ok = .false.
   call get_label_atom_coord('dist', i,2, ok)
   IF(.NOT. ok) cycle
   atom_2_coord(1:3) = new_coord(1:3)

   call distance_calculation(atom_1_coord(:), atom_2_coord(:), distance)

   WRITE(message_text,'(10x,5a,F10.5, a)') '    . d(', label_1(1:6), ' - ', label_2(1:6), ') = ',distance, ' A'
   call write_info(TRIM(message_text))

   IF(keyword_dist_X .or. keyword_dist_plus .or. keyword_DHA) then
    delta(1) =  atom_2_coord(1) - atom_1_coord(1)
    delta(2) =  atom_2_coord(2) - atom_1_coord(2)
    delta(3) =  atom_2_coord(3) - atom_1_coord(3)
    call write_info('')
    !WRITE(message_text,'(15x,5a, F10.5)')  '     Delta_x(', label_1(1:6), ' - ', label_2(1:6), ') = ', delta(1)
    !call write_info(TRIM(message_text))
    !WRITE(message_text,'(15x,5a, F10.5)')  '     Delta_y(', label_1(1:6), ' - ', label_2(1:6), ') = ', delta(2)
    !call write_info(TRIM(message_text))
    !WRITE(message_text,'(15x,5a, F10.5)')  '     Delta_z(', label_1(1:6), ' - ', label_2(1:6), ') = ', delta(3)
    !call write_info(TRIM(message_text))

	if(keyword_dist_plus) then
	 dist_coef = (distance + dist_plus) / distance
	 write(message_text,'(10x,3a,F10.5)')  '      > Distance from ', label_2(1:6), ': ', dist_plus
	 call write_info(trim(message_text))
	elseif(keyword_DHA) then
     dist_coef = (distance - dist_AH) / distance	
	 write(message_text,'(10x,3a,F10.5)')  '      > Distance from ', label_2(1:6), ': ', dist_AH
	 call write_info(trim(message_text))
	end if 
    WRITE(message_text,'(10x,a, 3F10.5)')  '      > Coordinates of the new atom:', atom_1_coord(1) + dist_coef * delta(1), &
                                                                                   atom_1_coord(2) + dist_coef * delta(2), &
                                                                                   atom_1_coord(3) + dist_coef * delta(3)
    call write_info(TRIM(message_text))
    call write_info('')
   endif
  end do


  return
end subroutine    calcul_distances


!------------------------------------------------------------------------------
subroutine distance_calculation(coord_1, coord_2, distance)
 use cryscalc_module, only : SP_value
 implicit none
  real, dimension(3), intent(in)     :: coord_1
  real, dimension(3), intent(in)     :: coord_2
  real,               intent(inout)  :: distance
  
  
  call scalar_product(coord_1(:), coord_2(:), coord_1(:), coord_2(:))
  distance = SQRT(SP_value)

 return
end subroutine distance_calculation 

!------------------------------------------------------------------------------
subroutine scalar_product(p1, p2, p3, p4)
 use cryscalc_module,         only      :  unit_cell , SP_value
 !USE CFML_Math_General,     ONLY      :  sp
 !USE CFML_constants,         ONLY      : sp 
 USE CFML_GlobalDeps,        ONLY : sp

 implicit none
 REAL(KIND=sp), DIMENSION(3), INTENT(IN)                 :: p1, p2, p3, p4
 REAL, parameter                                         :: pi=3.1415926535897932
 REAL,          DIMENSION(6)                             :: cell


 cell(1:6) = unit_cell%param(1:6)


 SP_value  =    ( (p2(1) - p1(1))  * (p4(1) - p3(1)) )  * cell(1)**2                                            &
             +  ( (p2(2) - p1(2))  * (p4(2) - p3(2)) )  * cell(2)**2                                            &
             +  ( (p2(3) - p1(3))  * (p4(3) - p3(3)) )  * cell(3)**2                                            &

             +  ( (p2(1) - p1(1))  * (p4(2) - p3(2))                                                            &
             +    (p4(1) - p3(1))  * (p2(2) - p1(2)) )  * cell(1) * cell(2) * COS(cell(6)* PI / 180.)           &

             +  ( (p2(1) - p1(1))  * (p4(3) - p3(3))                                                            &
             +    (p4(1) - p3(1))  * (p2(3) - p1(3)) )  * cell(1) * cell(3) * COS(cell(5)* PI / 180.)           &

             +  ( (p2(2) - p1(2))  * (p4(3) - p3(3))                                                            &
             +    (p4(2) - p3(2))  * (p2(3) - p1(3)) )  * cell(2) * cell(3) * COS(cell(4)* PI / 180.)


 return
end subroutine scalar_product

!------------------------------------------------------------------

 subroutine calcul_connectivity_atom
  use cryscalc_module,                ONLY : keyword_create_CIF, CIF_DIST, CONN_all, CONN_all_X, CONN_self, nb_atom, atom_CONN,  &
                                             CONN_species, CONN_out, CONN_excluded, calcul_BVS,                                  &
                                             CONN_dmin, CONN_dmax, CONN_ang, CONN_out_condensed,                                 &
											 SPG, crystal_cell, Atm_list, Atom_label, Atom_typ,                                  &
                                             Atom_Biso, Atom_occ, Atom_occ_perc, Atom_coord, Atom_mult,                          &
											 unit_cell, space_group_symbol,                                                      &
											 create_SHAPE_file, poly_vol_calc,                                                   &
											 message_text, debug_proc, keep_bond_str_out
  use CFML_geometry_calc,             ONLY : calc_dist_angle , calc_dist_angle_sigma, Coord_info
  use CFML_atom_TypeDef,              ONLY : allocate_Atom_list, deallocate_Atom_list, Atoms_cell_type, write_atom_list
  use CFML_BVS_energy_calc
  USE CFML_Crystal_Metrics,           ONLY : Set_Crystal_Cell, write_crystal_cell
  USE CFML_crystallographic_symmetry, ONLY : set_spacegroup, Get_Multip_Pos, write_spacegroup
  USE IO_module
  USE macros_module,                  ONLY : u_case
  

  implicit none
   integer                            :: i, i1, i2, nd, na, n_out
   integer                            :: n_atm, n_atm1, n_ang
   integer                            :: l1, l2, len_string
   integer                            :: i_error
   character(len=256)                 :: read_line, adjusted_line
   character(len=8)                   :: current_atom, atome_type, ligand
   logical                            :: atom_ok
   ! BVS
   type (Atoms_Conf_List_Type)        :: Ac   
   real                               :: ttol=20.0
   real                               :: eps
   integer                            :: charge_q
   integer                            :: n_bond, n_eff
   real                               :: r_eff, dist_AB
   real, dimension(100)               :: dist
   LOGICAL                            :: read_distances, read_angles
   LOGICAL                            :: search_angles
   
   
   read_distances = .false.
   read_angles    = .false.
   search_angles  = .false.


   if(create_SHAPE_file .or. poly_vol_calc) then
    if(CONN_all) then
     do i=1, nb_atom
      atom_CONN%label = atom_label(i)            
      call calcul_connect()      
     end do
    elseif(CONN_ALL_X) then
	 l1 =len_trim(CONN_species)
	 do i=1, nb_atom
	  l2 = len_trim(atom_typ(i))
	  if (l1 == l2 .and. atom_typ(i)(1:l2) == CONN_species(1:l1)) then
	   atom_CONN%label = atom_label(i)
	   CONN_out(i) = .true.
	   call calcul_connect()	   
	  end if
	 end do
	else
     call calcul_connect()   
    endif
	return
   end if	
   
   if(debug_proc%level_2)  call write_debug_proc_level(2, "calcul_connectivity_atom")
   
   eps = 0.0001
   
  !
  if(.not. CONN_all) then
   atom_ok = .false.
   do i=1, nb_atom
    if(CONN_all_X) then
	 l1 = len_trim(CONN_species)
	 l2 = len_trim(atom_typ(i))	 	 
	 if(l1 /= l2) cycle
	 if(u_case(atom_typ(i)(1:l2)) /= u_case(CONN_species(1:l1))) cycle
	 atom_CONN%label = atom_label(i)
	 conn_out(i) = .true.
	end if
    call Get_specie_from_type(atom_typ(i), atome_type)
	if(u_case(atom_CONN%label) == u_case(atom_label(i))) then    
     ATOM_CONN%type = atome_type
     atom_ok = .true.
     if(.not. CONN_all_X) exit
    end if
   end do
    
   if(.not. atom_ok) then
    call write_info('')
    call write_info(' Wrong atom for connectivity calculations !')
    call write_info('')
    return
   endif
  endif 
  
   
!  !--- write CIF ---------------------------------------------------------------------
!   if(keyword_create_CIF) then
!    if (allocated(CIF_DIST%text)) deallocate(CIF_DIST%text)
!    allocate(CIF_DIST%text(CIF_DIST%max_text)) !Maximum number of distances
!    CIF_DIST%text(:)(1:132)=" "
!    call write_CIF_file('DIST')
!   end if   
!  !-----------------------------------------------------------------------------------



!  --- new : jan. 2010 -----------------
!  idem Bond_str

   ! construction de l'objet Crystal_CELL
   if(debug_proc%level_3) call write_debug_proc_level(3, "construct SPG")
   !call Set_crystal_cell(unit_cell%param(1:3), unit_cell%param(4:6), crystal_cell)
   call create_CELL_object()
   
    ! construction de l'objet SPG
   call set_spacegroup(space_group_symbol, SPG)

   ! construction de l'objet Atm_list
   if(debug_proc%level_3) call write_debug_proc_level(3, "construct Atm_list")
   Atm_list%natoms = nb_atom
   if(allocated(Atm_list%atom))  deallocate(Atm_list%atom)
   allocate (Atm_list%atom(nb_atom))
 
   do i=1, Atm_list%natoms
    if(atom_mult(i) < eps) atom_mult(i) = Get_Multip_Pos(atom_coord(1:3,i), SPG)
   
    Atm_list%atom(i)%Lab      = atom_label(i)
    !Atm_list%atom(i)%ChemSymb = atom_typ(i)
	call get_specie_from_type(atom_typ(i), Atm_list%atom(i)%ChemSymb)
	Atm_list%atom(i)%Biso     = atom_Biso(i)
    !Atm_list%atom(i)%occ      = atom_occ(i)
	
	Atm_list%atom(i)%x(1:3)   = atom_coord(1:3, i) 
	Atm_list%atom(i)%x_std(1:3)   = 0.    ! <<< indispensable pour le calcul des distances !!!
	atm_list%atom(i)%mult     = atom_mult(i)
	Atm_list%atom(i)%occ      = real(Atm_list%atom(i)%Mult)/real(SPG%Multip)   !! << new
	
	charge_q = 0
	call Get_charge_from_specie(atom_typ(i), charge_q)
	Atm_list%atom(i)%charge  = real(charge_q)
	Atm_list%atom(i)%AtmInfo = ''
   end do
   
   close(unit=55)
   open(unit=55, file="bond_str.out")
    !call write_Crystal_cell(crystal_cell, 55)
	!call write_atom_list(Atm_list, 0, 55)
	  
    !call Calc_dist_angle(CONN_dmax, 0., crystal_cell, SPG, Atm_list, 55)
	if(CONN_ang) then
	 call Calc_Dist_Angle_Sigma(CONN_dmax, CONN_dmax, crystal_cell, SPG, Atm_list, 55)        
	else
	 call Calc_Dist_Angle_Sigma(CONN_dmax, 0., crystal_cell, SPG, Atm_list, 55)        
	endif 
	
    call Calc_PDB(Atm_list, CONN_dmax, 55)
	
	if(calcul_BVS) then	
	 if(debug_proc%level_3) call write_debug_proc_level(3, "calcul BVS (CFML)")
	 call Allocate_Atoms_Conf_List(Atm_list%natoms, Ac)
	 Ac%atom = Atm_list%atom	 
     
	 call Species_on_List(Ac, SPG%Multip, ttol)	 
	 if(err_conf) then
	  call write_info(' ')
	  write(message_text, '(a)') '   Species on list: '//trim(ERR_conf_Mess)
	  call write_info(trim(message_text))
	  call write_info('    ==> PROGRAM BOND_STR finished in error!')
	  call write_info('')
	  return
	 end if
	 	 
	 ! setting tables for B and D0
      call set_Table_d0_b(Ac)	 
	  if(err_conf) then
	   call write_info(' ')
	   write(message_text, '(a)') '   Species on list: '//trim(ERR_conf_Mess)
	   call write_info(trim(message_text))
	   call write_info('    ==> PROGRAM BOND_STR finished in error!')
	   call write_info('')
	   return
	  end if
	  
	  call Calc_BVS(Ac, 55)
	  
	  deallocate (Ac%atom)
	end if
	deallocate (Atm_list%atom )
	
   close(unit=55)
   
   call write_info('')
   call write_info('')
   call write_info('   >> Atomic connectivity calculation (JRC-JGP Bond_Str routine) : ')
   call write_info('')
   if(CONN_excluded) then
    call write_info('    !! '//trim(ATOM_CONN%excluded)// " atoms are excluded from listing !!")
	call write_info('')
   end if
   
   if(debug_proc%level_3) call write_debug_proc_level(3, "Atom connectivity")
   open(unit=55, file="bond_str.out")
   n_atm  = 0
   n_atm1 = 0
   n_ang  = 0
   do
    read(55, '(a)',iostat=i_error) read_line
    if(i_error /= 0) exit
    adjusted_line = adjustl(read_line)
    if(adjusted_line(1:10) == '----------') then
     !call write_info(trim(read_line))
     if(CONN_all) then    
      do
       read(55, '(a)',iostat=i_error) read_line
       if(i_error /= 0) exit   
       if(CONN_out_condensed) then
	    adjusted_line = adjustl(read_line)
		
		if(adjusted_line(1:24) == '{--- BONDS DISTRIBUTIONS' .and. .not. search_angles) then
		 rewind(unit=55)
		 search_angles = .true.
		 cycle
		end if 
		
		if(.not. search_angles) then
		 if(adjusted_line(1:19) == 'Distances less than') then
		  !call write_info('')	    
		  !call write_info('   - ' //trim(adjusted_line))
		  !call write_info('')
		  n_atm = n_atm + 1
		  if(n_atm == 1) then
		   call write_info('')
		   call write_info('    ..... Interatomic distances .....')
		   call write_info('')
		   call write_CIF_file('DIST')
		  end if 
		 end if
         !if(index(adjusted_line, '(') /=0 .and. index(adjusted_line,'{')==0) then
		 if(read_line(30:31) == '):') then
		  ! compare avec atome a exclure
		  if(CONN_excluded) then
		   ligand = read_line(25:29)
		   ligand = adjustl(ligand)
		   len_string = len_trim(ATOM_conn%excluded)
		   if(ligand(1:len_string) == ATOM_conn%excluded(1:len_string)) cycle
		  end if
		  
          call write_info('    '//read_line(16:40))    ! distance
		  if(keyword_create_CIF)  call write_dist_CIF(read_line)		
		 end if 
		
		else  ! seconde lecture pour les angles
		 if(adjusted_line(1:22) == '-  Angles around atom:') then
		  !call write_info('')	    
		  !call write_info(trim(read_line))
		  !call write_info('')
		  n_atm1 = n_atm1 + 1
		  if(n_atm1 == 1) then
		   call write_info('')
		   call write_info('    ..... Interatomic angles .....')
		   call write_info('')
		   call write_CIF_file('ANG')		   
		  end if 
		 end if
		 if(adjusted_line(1:5) == 'Atm-1') then
		  n_ang = 0
		  cycle 
		 end if 

		 if(read_line(5:5) == '(' .and. read_line(13:13) == '(' .and. read_line(21:21) == '(') then
		  n_ang = n_ang +1 
          if(n_ang == 1) then 
		   call write_info(trim(read_line))
	       if(keyword_create_CIF) call write_ang_CIF(read_line)   ! angle
		  end if 
		 end if
		end if
		cycle
       else	   
        call write_info(trim(read_line))
	   end if	
      end do
      exit     
      
     else  ! 
      ! DISTANCES
	  n_out = 0
      do
       read(55, '(a)', iostat=i_error) read_line
       if(i_error /=0) exit
	   
	   i1 = index(read_line, 'Distances less than')
	   i2 = index(read_line, 'to atom:')
       
       if(i1/=0 .and. i2/=0) then
	    n_out = n_out + 1
		
	    read_distances = .true.
        read(read_line(index(read_line, ':')+1:), *) current_atom
		
        if(u_case(ATOM_CONN%label) == u_case(current_atom) .or. (CONN_all_X .and. CONN_out(n_out))) then
		 call write_info(trim(read_line))
         read(55, '(a)', iostat=i_error) read_line
         if(i_error /=0) exit
         call write_info(trim(read_line))         
		 
         nd = 0
         do
          read(55, '(a)', iostat=i_error) read_line
          if(i_error /=0) exit
          adjusted_line = adjustl(read_line)
          if(adjusted_line(1:10) == '----------') exit
		  if(index(read_line, '(') /=0) then 
		   read(read_line(32:40),*, iostat=i_error) dist_AB
     	   if(i_error /= 0 .or. dist_AB < CONN_dmin) cycle
		  ! compare avec atome a exclure
	 	   if(CONN_excluded) then
		    ligand = read_line(25:29)
		    ligand = adjustl(ligand)
		    len_string = len_trim(ATOM_conn%excluded)
		    if(ligand(1:len_string) == ATOM_conn%excluded(1:len_string)) cycle
		   end if

		   if(CONN_out_condensed) then
		    call write_info(read_line(16:40))
			if(keyword_create_CIF)  call write_dist_CIF(read_line)	 
			cycle
		   end if
		  else
		   if(CONN_out_condensed) cycle
		  end if 
		  if(CONN_SELF) then
		   if(index(read_line, 'Orig. extr. p.equiv.')/=0) call write_info(trim(read_line))
		   if(read_line(17:21) == read_line(25:29))        call write_info(trim(read_line))
		  else
           call write_info(trim(read_line))
		  end if
		  if(keyword_create_CIF .and. read_distances) then
		   nd = nd + 1
           if(nd==1) call write_CIF_file('DIST')
		   if(keyword_create_CIF)  call write_dist_CIF(read_line)		  
		  end if 
         end do
        endif 
       endif 
	   
	   if(CONN_ang) then
	   i1 = index(read_line, 'Angles around atom')     
       if(i1/=0) then
	    read_angles = .true.
        read(read_line(index(read_line, ':')+1:), *) current_atom
        if(u_case(ATOM_CONN%label) == u_case(current_atom)) then
         call write_info(trim(read_line))
         read(55, '(a)', iostat=i_error) read_line       ! ligne '---------------------------'
         if(i_error /=0) exit
         call write_info(trim(read_line))         
         
		 if(CONN_out_condensed) then
		  call write_info('')
		  if(keyword_create_CIF) call write_CIF_file('ANG')    ! CIF angle ²_loop 
		 end if 
         nd = 0		 
         do
          read(55, '(a)', iostat=i_error) read_line
          if(i_error /=0) exit
          adjusted_line = adjustl(read_line)
		  if(CONN_out_condensed) then
		   if (len_trim(adjusted_line) == 0) cycle
		   if(adjusted_line(1:3) == 'Atm') then
		    na = 0
		    cycle
		   end if	
		   if(adjusted_line(1:1) == '(') then
		    na = na + 1
			if (na ==1) then
			 call write_info(trim(read_line))
			 if(keyword_create_CIF) call write_ang_CIF(read_line)
			end if 
			cycle
		   end if
		  end if
		  
          if(adjusted_line(1:10) == '----------') exit
		  if(index(read_line, 'd12 =') /=0) then 
		   read(read_line(44:50),*, iostat=i_error) dist_AB
		   if(i_error /= 0 .or. dist_AB < CONN_dmin) cycle
		  end if 
		  if(CONN_SELF) then
		   if(index(read_line, 'Atm-1   Atm-2   Atm-3')/=0) call write_info(trim(read_line))
		   if(read_line(6:10) == read_line(22:26))          call write_info(trim(read_line))
		  else
           call write_info(trim(read_line))
		  end if
		  if(keyword_create_CIF .and. read_angles) then 
		   nd = nd + 1
           if(nd==1) call write_CIF_file('ANG')
		   if(trim(adjustl(read_line(14:18))) == u_case(current_atom)) call write_ang_CIF(read_line)		  
		  end if 
         end do
        endif 
       endif 
       end if
	   
      end do 
     
     ! BONDS DISTRIBUTIONS        	 
      rewind(unit=55)
      do
       read(55, '(a)', iostat=i_error) read_line
       if(i_error /=0) exit
       i1 = index(read_line, '{--- BONDS DISTRIBUTIONS')
       if(i1 /=0) then
	    if(debug_proc%level_3) call write_debug_proc_level(3, "bond distributions")
		if(CONN_out_condensed) call write_info('')
	    call write_info(trim(read_line))
		call write_info('')
        backspace(unit=55)
        backspace(unit=55)
        if(.not. CONN_all) then
         do
          read(55, '(a)', iostat=i_error) read_line
          if(i_error /=0) exit
          i1 = index(read_line, 'Bond type:')
		  n_eff = 0
          if(i1 /=0) then
           read(read_line(index(read_line, ':')+1:), *) current_atom	
           read(read_line(18:20), *) ligand
		   ligand = adjustl(ligand)
		   len_string = len_trim(ATOM_conn%excluded)
		   if(ligand(1:len_string) == ATOM_conn%excluded(1:len_string)) cycle

		   
           if(u_case(ATOM_CONN%type) == u_case(current_atom)) then		    
		    !call write_info(trim(read_line))
            if(CONN_self) then
			 
			 if(read_line(13:14) == read_line(18:19)) then
   			  call write_info(read_line)		
              read(55, '(a)', iostat=i_error) read_line
			  if(i_error/=0) exit
			  call write_info(trim(read_line))			  
			  do 
			   read(55, '(a)', iostat=i_error) read_line
			   if(i_error /=0) exit
			   if(len_trim(read_line) == 0) then
			    call write_info('')
				exit
			   end if			  			   
               !call write_info(trim(read_line))		
               read(read_line, *, iostat=i_error) n_bond, dist_AB
			   if(dist_AB > CONN_dmin) then				 
			    n_eff = n_eff + 1
				dist(n_eff) = dist_AB
			    call write_info(trim(read_line))						
			   end if	
			  end do			  			 
			 else
			  cycle
			 end if 			
			else
			 call write_info(trim(read_line))
			 read(55, '(a)', iostat=i_error) read_line
			 if(i_error/=0) exit
			 call write_info(trim(read_line))
             do 
              read(55, '(a)', iostat=i_error) read_line
              if(i_error /=0)              exit
              if(len_trim(read_line) == 0) then
               call write_info('')
               exit
              endif
              read(read_line, *, iostat=i_error) n_bond, dist_AB
			  if(dist_AB > CONN_dmin) then
			   n_eff = n_eff + 1
			   dist(n_eff) = dist_AB				
			   call write_info(trim(read_line))
			  end if				  
             end do			 
			end if 
			if(n_eff /=0) then			 
			 r_eff = (n_eff/sum(dist(1:n_eff)**(-3.)))**(0.33333333)
			 write(message_text, '(a,I4)')     '  >> n_eff = ', n_eff 
             call write_info(trim(message_text))	
		     write(message_text, '(a,1F15.5)') '  >> r_eff = ', r_eff 
             call write_info(trim(message_text))	
			 write(message_text, '(a,1F15.5)') '  >> <r>   = ', sum(dist(1:n_eff))/n_eff
             call write_info(trim(message_text))				 
            else			 
			 write(message_text, '(a,F6.2,a,F6.2,a)') '  >> No atomic distance between ', CONN_dmin, ' and ', CONN_dmax, ' A !'
             call write_info(trim(message_text))	
			end if 
            call write_info('')
           endif           
		  endif
         end do 
		 
		 rewind(unit=55)
		 do 
		  read(55, '(a)', iostat=i_error) read_line
		  if(i_error /=0) exit
		  i1 = index(read_line, '{--- BOND-VALENCE')
		  if(i1 /=0) then
		   if(debug_proc%level_3) call write_debug_proc_level(3, "bond valence")
		   call write_info(trim(read_line))
		   call write_info('')
		   backspace(unit=55)
		   backspace(unit=55)
		   do
		    read(55, '(a)', iostat=i_error) read_line
			if(i_error/=0) exit
			!call write_info(trim(read_line))
			i1 = index(read_line, '=> Bond-valence and coordination of atom:')
			if(i1/=0) then
			 read(read_line(index(read_line, ':')+1:), *) current_atom
			 if(u_case(ATOM_CONN%label) == u_case(current_atom)) then
			  call write_info(trim(read_line))
			  do 
			   read(55, '(a)', iostat=i_error) read_line
			   if(i_error /=0) exit
			   call write_info(trim(read_line))
			   i1 = index(read_line, '{r2=<sij-<sij>>rms}')
			   if(i1 /=0) exit
			  end do
			  exit
			 else
			  cycle
			 end if
            
			end if	
			
		   end do
		   
		  end if
		  
		 end do
		 
        else ! calcul de la connectivite pour tous les atomes
         do 
          read(55, '(a)', iostat = i_error) read_line
          if(i_error /=0) exit
          call write_info(trim(read_line))
         end do 
        end if ! 
       endif 
      end do
      
     end if
    end if
    
    
   end do
   close(unit=55) 
   if(.not. keep_bond_str_out) call system("del bond_str.out")

  ! ---------------------------------------------------------------------
  
  ! --- old ----
  ! if(CONN_all) then
  !  do i=1, nb_atom
  !   atom_CONN%label = atom_label(i)      
  !   call calcul_connect()      
  !  end do
  ! else
  !  call calcul_connect()   
  ! endif
  ! ------------
  
  if (allocated(CIF_dist%text)) deallocate(CIF_dist%text)
  
 return 
 end subroutine calcul_connectivity_atom


!----------------------------------------------------------------------------------


subroutine calcul_connect()  
 USE IO_module
 USE cryscalc_module, ONLY : message_text, nb_atom, atom_label, atom_coord,  &
                             SPG, keyword_create_CIF,                        &
                             atom_CONN, CONN_dmax, CONN_excluded,            &
                             new_coord, SP_value, atom1_dist, atom2_dist, unit_cell, &
                             CIF_unit, keyword_create_CIF, CIF_DIST, debug_proc, &
							 tmp_unit, main_title, crystal_cell,                 &
							 create_shape_file, poly_vol_calc, cartesian_frame, message_text 
                             
 
 USE CFML_crystallographic_symmetry, only : Get_orbit, Wyckoff_Type, ApplySO
 USE CFML_GlobalDeps,                ONLY : cp, sp
 USE CFML_crystal_metrics,           only : cart_vector
 USE CFML_Math_3D,                   only : Polyhedron_Volume

  implicit none
   REAL (kind=cp), DIMENSION(3)      :: atom_1_coord, atom_2_coord
   real (kind=cp), dimension(3)      :: r
   real, dimension(3)                :: tn
   integer                           :: i, ieq, it1, it2, it3, i_tr, m
   integer                           :: id , len_string, n_blank
   integer                           :: itx, ity, itz 
   CHARACTER(LEN=16)                 :: label_1, label_2
   LOGICAL                           :: ok
   real                              :: distance
   type(wyckoff_type)                :: wyckoff
   
   integer                                 :: lk, nn
   real (kind=cp),    dimension(3)         :: x0,x1
   real (kind=cp), dimension(3,Spg%multip) :: uu
   character(len=16)                       :: transla

   real ,   dimension(3,192)               :: orbit
   CHARACTER (LEN=80)                      :: text80
   integer, dimension(192)                 :: effsym  ! Pointer to effective symops
   integer, dimension(192)                 :: ss_ptr  ! Pointer to stabilizer

   real (kind=sp),    dimension(3)     :: cart_vect
   real,              dimension(3, 30) :: poly_atom_coord
   real,              dimension(3, 30) :: poly_atom_coord_cart
   real,              dimension(3)     :: central_atom_coord_cart
   character (len=8), dimension(30)    :: poly_atom_label
   integer                             :: poly_atom_number
   character (len=256)                 :: string, shape_file
   real                                :: poly_vol    ! volume du polyedre
   character (len=128)                 :: fmt_

   poly_atom_number = 0
   
  if(debug_proc%level_2)  then
   if(create_SHAPE_file) then
    call write_debug_proc_level(2, "calcul_connect and create shape.dat file")
   else
    call write_debug_proc_level(2, "calcul_connect")
   end if	
  end if
  
  atom1_dist(1) = atom_CONN%label
  ok = .false.
  call get_label_atom_coord('dist', 1,1, ok)
  if(.not. ok) return
  atom_1_coord(1:3) = new_coord(1:3)
  x0(1:3) = atom_1_coord(1:3)
  
  label_1 = atom_CONN%label  
  call write_info('')
  cart_vect(1:3) = cart_vector('D',atom_1_coord(:), crystal_cell)
  WRITE(message_text, '(2a,5x,3F10.5,5x,a,3f10.5,a)') '  >> Connectivity around : ', trim(atom_CONN%label), atom_1_coord(1:3), &
                                                      '{cart. coord. = ', cart_vect(1:3) ,'}'
  call write_info(trim(message_text))
  WRITE(message_text, '(a,F6.2,a)')  '     d_max : ', CONN_dmax, ' A'
  call write_info(trim(message_text))
  call write_info('')
  
  !if(create_SHAPE_file) then
   poly_atom_number = poly_atom_number + 1
   poly_atom_coord(1:3, poly_atom_number) = atom_1_coord(1:3)
   poly_atom_label(poly_atom_number)      = trim(atom_CONN%label)
  !end if
  
  itx = int(unit_cell%param(1) / CONN_dmax) +  1.5
  ity = int(unit_cell%param(2) / CONN_dmax) +  1.5
  itz = int(unit_cell%param(3) / CONN_dmax) +  1.5
  
  !--- write CIF ---------------------------------------------------------------------
   if(keyword_create_CIF) then
    if (allocated(CIF_DIST%text)) deallocate(CIF_DIST%text)
    allocate(CIF_DIST%text(CIF_DIST%max_text)) !Maximum number of distances
    CIF_DIST%text(:)(1:132)=" "
    call write_CIF_file('DIST')
   end if   
  !-----------------------------------------------------------------------------------


  ! boucle sur les differents atomes presents 
  id = 0
  do i=1, nb_atom
   atom2_dist = atom_label(i)   
   label_2 = atom2_dist(i)
   ! compare avec atome a exclure
   if(CONN_excluded) then
	len_string = len_trim(ATOM_conn%excluded)
	if(label_2(1:len_string) == ATOM_conn%excluded(1:len_string)) cycle
   end if   
   
   ok = .false.
   call get_label_atom_coord('dist', i,2, ok)
   IF(.NOT. ok) cycle
   atom_2_coord(1:3) = new_coord(1:3)
   
   lk=1
   uu(:,lk) = x0(:)
   
   r(1:3) = atom_2_coord(1:3)
      
	  

 ! boucle sur les atomes equivalents par symétrie
   do ieq = 1, SPG%Multip
   
     new_coord(1:3) = ApplySO(SPG%symop(ieq), r(1:3))
   ! boucle sur les atomes des mailles voisines
    do it1=-itx, itx 
     do it2=-ity, ity
      do_it3:do it3=-itz, itz
 
       Tn(1) = real(it1)
       Tn(2) = real(it2)
       Tn(3) = real(it3)
       atom_2_coord(1:3) = new_coord(1:3) + Tn(1:3)
       do nn=1, lk
        if(sum(abs(uu(:,nn) - atom_2_coord(:))) < 0.01) cycle do_it3 
       end do
       call distance_calculation(atom_1_coord(:), atom_2_coord(:), distance)
       
     
       if(distance < 0.1 .or. distance > CONN_dmax) cycle do_it3       
       !if(distance > CONN_dmax) cycle
       
       lk=lk+1
       uu(:,lk) = atom_2_coord(1:3)
       id = id + 1
	   
	   cart_vect(1:3) = cart_vector('D',atom_2_coord(:), crystal_cell)
	   
	   len_string = len_trim(SPG%SymopSymb(ieq))
	   n_blank = 22 - len_string
	   write(fmt_, '(a, i2,a)')  '(2x,i3,5a,F10.5, a, 5x,3f10.5,2x,a,3I3,3a,', n_blank  , 'x,a,3F10.5,a )'
       WRITE(message_text,fmt=trim(fmt_))    &
        id,'. d(', label_1(1:4), ' - ', label_2(1:4), ') = ',distance, ' A', atom_2_coord(:),  &
        '  T=[', int(tn(1:3)), ']  (', trim(SPG%SymopSymb(ieq)), ')' ,'{cart. coord. = ',      &
		cart_vect(1:3), '}'
       !WRITE(message_text,'(2x,i3,5a,F10.5, a, 5x,3f10.5,2x,a,3I3,3a,5x,a,3F10.5,a )')    &
       ! id,'. d(', label_1(1:4), ' - ', label_2(1:4), ') = ',distance, ' A', atom_2_coord(:),  &
       ! '  T=[', int(tn(1:3)), ']  (', trim(SPG%SymopSymb(ieq)), ')' ,'{cart. coord. = ',      &
	!	cart_vect(1:3), '}'
       
       call write_info(TRIM(message_text))       
	   
	   !if(create_SHAPE_file) then
	    poly_atom_number = poly_atom_number + 1
        poly_atom_coord(1:3, poly_atom_number) = atom_2_coord(1:3)
		poly_atom_label(poly_atom_number)      = label_2(1:4)
	   !end if
	   
       if(keyword_create_CIF) then
        CIF_dist%n_text = CIF_dist%n_text + 1
        if(it1==0 .and. it2==0 .and. it3==0 .and. effsym(ieq) == 1) then
         write(CIF_dist%text(CIF_dist%n_text),'(a,2x,a,2x,F10.5,a)') label_1(1:4), label_2(1:4), distance, '  .  ?'
        else
         write(CIF_dist%text(CIF_dist%n_text),'(a,2x,a,2x,F10.5,a,i3,a,3I1,a)') label_1(1:4), label_2(1:4), distance,  &
                        " ", effsym(ieq), "_", int(tn(1:3)+5.), " ?"
        endif 
       end if
   
       end do do_it3  ! it3 
     end do ! it2
    end do ! it1
   end do ! ieq
  end do ! nb_atom
  
   
  
  if (id==0) then
   write(message_text, '(3a,F6.2,a)') "   >>> No atom connected to ", trim(atom_CONN%label), " until dmax = ", CONN_dmax, &
                                      " A. Try increase dmax !"      
   call write_info(trim(message_text))    
  end if

  
  if(keyword_create_CIF)  call write_CIF_file('DIST_values')  
  
  ! create shape.dat
  if(create_SHAPE_file) then
   call  create_CELL_object
   !call write_info('') 
   !call write_info('   Cartesian frame type : '//cartesian_frame%type//' ' //trim(cartesian_frame%string))
   !call write_info('') 
 
   write(shape_file, '(3a)') 'shape_',trim(poly_atom_label(1)),'.dat'
   open(unit = tmp_unit, file=trim(shape_file))
   write(unit = tmp_unit, fmt='(a)')    '! Input file for SHAPE, created by CRYSCALC (TR/CDIFX-ISCR Rennes)'
   write(unit = tmp_unit, fmt='(2a)') '$ ', trim(main_title)
   write(unit = tmp_unit, fmt='(a)')  '! Ligands    central atom (line 1 in the atoms list)' 
   write(string, fmt='(2(1x,i3))') poly_atom_number-1, 1
   write(unit = tmp_unit, fmt='(a)') adjustl(string)
   write(unit = tmp_unit, fmt='(a)')  '! Polyedron will be compared to any kind of polyedra'
   write(unit = tmp_unit, fmt='(a)') '1 2 3 4 5 6 7 8 9 10 11 12 13'
   write(unit = tmp_unit, fmt='(a)')  '! Label for the polyedron'
   write(unit = tmp_unit, fmt='(2a)') 'POLY_',trim(poly_atom_label(1))
   write(unit = tmp_unit, fmt='(a)')  '! List of cartesian atomic coordinates' 
  end if
  
  if(create_SHAPE_file) call write_info('')
  do i= 1, poly_atom_number
   cart_vect(1:3) = cart_vector('D', poly_atom_coord(1:3, i), crystal_cell)
   if(create_SHAPE_file) then
    write(unit = tmp_unit, fmt='(a4,2x,3F12.6)')     poly_atom_label(i) , cart_vect(1:3) 	
   end if	
   poly_atom_coord_cart(1:3, i) = cart_vect(1:3)	
  end do	

  if(poly_vol_calc) then
  ! calcul du volume du polyedre
   central_atom_coord_cart(1:3) = poly_atom_coord_cart(1:3, 1)
   poly_vol = Polyhedron_volume(poly_atom_number-1,  poly_atom_coord_cart(1:3, 2:poly_atom_number), central_atom_coord_cart)
   call write_info('')
   write(message_text, '(5x,a,F15.5,a)') ' >> Polyedron volume = ', poly_vol, ' A3'
   call write_info(trim(message_text))
   call write_info('')
  end if 
   

  close(unit=tmp_unit)

  if(create_SHAPE_file) then
   call write_info('')
   call write_info('   > '//trim(shape_file)//' for SHAPE has been created ('//   &
                         'cartesian frame type: ' //trim(cartesian_frame%string)//').')
   call write_info('')   
  end if 
 !end if  

  if (allocated(CIF_dist%text)) deallocate(CIF_dist%text)
  
 return 
end subroutine calcul_connect

!----------------------------------------------------------------------------------------------
! distribution des distances
! extract from Bond_Str (JRC, JGP)
   Subroutine Calc_PDB(A,Dmax,lun)
    
   USE CFML_Atom_TypeDef,              ONLY  : Atom_list_Type
   USE CFML_Geometry_Calc,             ONLY  : Coord_Info
   USE CFML_atom_TypeDef,              ONLY  : allocate_Atom_list, deallocate_Atom_list
   USE CFML_Math_General,              ONLY  : sort
   USE MACROS_module,                  ONLY  : U_case
   USE cryscalc_module,                ONLY  : debug_proc

   
    implicit none
      !---- Arguments ----!
      type (Atom_list_Type), intent(in)         :: A
      real,                  intent(in)         :: Dmax
      integer,               intent(in)         :: lun

      !---- Local variables ----!
      integer, parameter                        :: Max_NSpecies=20
      character(len=2), dimension(Max_NSpecies) :: species
      character(len=2)                          :: car

      integer                                   :: n_spec
      integer, dimension(:,:),allocatable       :: bond_spec
      integer                                   :: i,j,n,n_c,n1,n2

      integer, parameter                        :: Max_dis=8000
      integer                                   :: n_dis,n_ini
      integer, dimension(Max_Dis)               :: indx
      real,dimension(Max_Dis)                   :: disbond
      real                                      :: dis_ini

      type (Atom_list_Type)                     :: Ac

	  
	  if(debug_proc%level_3)  call write_debug_proc_level(3, "calc_PDB")
	  
      !---- Init control ----!
      if (A%natoms <=0) return

      write(unit=lun,fmt="(/,a)")          "  ---------------------------------------------"
      write(unit=lun,fmt="(a,f6.3,a)")     "  {--- BONDS DISTRIBUTIONS (UP TO ",dmax, " ) ---}"
      write(unit=lun,fmt="(a)")            "  ---------------------------------------------"

      !---- Make a copy of A ----!
      call Allocate_Atom_List(A%natoms,Ac)
      Ac%Atom=A%atom

      !---- Calculate the number of different species in the List ----!
      n_spec=0
      species=" "
      at1:do n=1,Ac%natoms
         if (n_spec > 0) then
            car=u_case(Ac%Atom(n)%ChemSymb)
            do i=1,n_spec
               if (car == species(i)) cycle at1
            end do
         end if
         if (n_spec == Max_NSpecies) then
            write(unit=lun,fmt="(/,a)")  "Overflow the number of different species the program can use"
            call Deallocate_atom_list(Ac)
            return
         end if
         n_spec=n_spec+1
         species(n_spec)=u_case(Ac%Atom(n)%ChemSymb)
      end do at1

      !---- Assign the index in the Atom list respect to species ----!
      do n=1,Ac%natoms
         do i=1,n_spec
            if (u_case(Ac%Atom(n)%ChemSymb) == species(i)) then
               Ac%atom(n)%ind(1)=i
               exit
            end if
         end do
      end do

      !---- Creating Table for bonds between species ----!
      allocate(bond_spec(n_spec,n_spec))
      bond_spec=0

      do n=1,coord_info%natoms
         n1=Ac%atom(n)%ind(1)
         n_c=coord_info%coord_num(n)
         do i=1,n_c
            !j=coord_info%n_cooatm(n,i)
			j=coord_info%n_cooatm(i,n)   ! 05.2014 : avec nouvelle library CFML
            n2=Ac%atom(j)%ind(1)
			bond_spec(n1,n2)=bond_spec(n1,n2)+1
            bond_spec(n2,n1)=bond_spec(n2,n1)+1
         end do
      end do

      !---- Calculate the Distributions Information ----!
      do n1=1,n_spec
         do n2=1,n_spec
            if (bond_spec(n1,n2) == 0) cycle
            !if (n1 > n2) cycle   <<< TR nov. 2010

            n_dis=0
            indx=0
            disbond=0.0
            do n=1,coord_info%natoms
               if (Ac%atom(n)%ind(1) /= n1) cycle
               n_c=coord_info%coord_num(n)
               do i=1,n_c
                  !j=coord_info%n_cooatm(n,i)
				  j=coord_info%n_cooatm(i,n)   ! 05.2014 : avec nouvelle library CFML
                  if (Ac%atom(j)%ind(1) /= n2) cycle
                  if (n_dis+1 > Max_Dis) then
                     write(unit=lun,fmt="(a,i6/)")  "Overflow the number of distances the program can handle: ",max_dis
                     deallocate(bond_spec)
                     call Deallocate_atom_list(Ac)
                     return
                  end if
                  n_dis=n_dis+1
                  !disbond(n_dis)=coord_info%dist(n,i)
				  disbond(n_dis)=coord_info%dist(i,n)          ! 05.2014 : avec nouvelle library CFML
               end do
            end do

            call sort(disbond,n_dis,indx)

            write(unit=lun,fmt='(/,a,i4)') " Bond type: "//species(n1)//" - "//species(n2)//"      Num: ",n_dis
            write(unit=lun,fmt='(a)')   "  Num. Bonds         Distance"

            n_ini=1
            dis_ini=disbond(indx(n_ini))
            !do i=2,n_dis
			do i=n_ini, n_dis   ! modif. TR 19.04.2011
               if (abs(disbond(indx(i))-dis_ini) <= 0.001) then
                  if (i < n_dis) cycle
                  write(unit=lun,fmt='(3x,i3,12x,f12.4)') i-n_ini+1,disbond(indx(n_ini))
               else
                  write(unit=lun,fmt='(3x,i3,12x,f12.4)') i-n_ini,disbond(indx(n_ini))
                  n_ini=i
                  dis_ini=disbond(indx(i))
                  if (n_ini == n_dis) write(unit=lun,fmt='(3x,i3,12x,f12.4)') 1,dis_ini
               end if
            end do

         end do
      end do

      deallocate(bond_spec)
      call Deallocate_atom_list(Ac)

      return
   End Subroutine Calc_PDB

   
!--------------------------------------------------------
 subroutine Get_Charge_from_specie(label, q)
  implicit none
   character (len=*), intent(in) :: label
   integer,           intent(out)   :: q
   integer                          :: ip, im, ier
   
   
   
   ip = index(label, '+')
   select case(ip)
    case(0,1)  ! no sign +
	 im = index(label, '-')
	 select case(im)
	  case (0,1)
	   q = 0
	  
	  case (2) ! element with a single character (ex: F-1)
	   read(unit=label(3:), fmt='(i1)', iostat = ier) q
	    if (ier/=0) q = 0				
	
      case (3) ! element in the form: F1- or Br-1
	   read(unit=label(2:2), fmt='(i1)', iostat =  ier) q
	   if(ier /=0) then
	    read(unit=label(4:4), fmt='(i1)', iostat =  ier) q
		if(ier/=0) q = 0			   
	   end if
	   
	  case  (4)  ! element in the form: Br1-
        read(unit=label(3:3), fmt='(i1)', iostat =  ier) q	
        if(ier/=0) q = 0				
		
	 end select
	 q = -q  ! anions
	
    case(2) ! element with a single character (ex: C+4)
      read(unit=label(3:), fmt='(i1)', iostat = ier) q
	    if (ier/=0) q = 0		
		
	case(3) ! element in the form: C4+ or Fe+3
      read(unit=label(2:2), fmt='(i1)', iostat = ier) q
	   if(ier /=0) then
	    read(unit=label(4:4), fmt='(i1)', iostat =  ier) q
		if(ier/=0) q = 0	
	   end if
	   
	case  (4)  ! element in the form: Fe3+
        read(unit=label(3:3), fmt='(i1)', iostat =  ier) q	
        if(ier/=0) q = 0			
		
   end select
   
   
   
 
  return
 end subroutine Get_Charge_from_specie
 !--------------------------------------------------------
 subroutine Get_specie_from_type(label, atome_type)
  implicit none
   character (len=*), intent(in)  :: label
   character (len=*), intent(out) :: atome_type
   integer                        :: ip, im, i, long
   
   atome_type = label
   long = len_trim(label)
   
   ip = index(label, '+')
   if(ip /=0) then
    if(ip == long) then
     atome_type = label(1:long-2)
    else	 
	 atome_type = label(1:ip-1)  ! cation
	end if 
   else	
    im = index(label, '-')
	if(im /=0) then
	 if(im == long) then
	  atome_type = label(1:long-2)
	 else 
	  atome_type = label(1:im-1)
	 end if 
	endif 
   endif 
    
  return
 end subroutine Get_specie_from_type
!--------------------------------------------------------

subroutine write_dist_CIF(read_line)
 use cryscalc_module, only : CIF_unit
!    
  implicit none
  character (len=*), intent(in) :: read_line
  character (len=256)           :: new_line
  integer                       :: i1, i2
  character (len=8)             :: label_1, label_2
  real                          :: distance
  integer                       :: tx, ty, tz
  integer                       :: sym_op_nb

  ! label 1
  i1 = index(read_line, '(')  
  i2 = index(read_line, ')')
  if(i1 /=0 .and. i2/=0 .and. i2>i1) then
   read(read_line(i1+1:i2-1),*) label_1
   label_1 = adjustl(label_1)
  else
   return
  endif
  new_line = read_line(i2+1:)
 
  ! label 2
  i1 = index(new_line, '(')  
  i2 = index(new_line, ')')
  if(i1 /=0 .and. i2/=0 .and. i2>i1) then
   read(new_line(i1+1:i2-1),*) label_2
   label_2 = adjustl(label_2)
  else
   return
  endif
  
  new_line=new_line(i2+2:) 
  ! distance  
  read(new_line, *) distance
  
  ! tx, ty, tz
  i1 = index(new_line, '(')
  i2 = index(new_line, ',')
  if(i1 /=0 .and. i2 /=0 .and. i2 > i1) then
   read(new_line(i1+1: i2-1), *) tx
  else
   return
  endif
  new_line = new_line(i2+1:)  
  read(new_line, *) ty
  
  i1 = index(new_line, ",")
  i2 = index(new_line, ')')
  read(new_line(i1+1:i2-1), *) tz
  
  i1 =index(new_line, ')', back=.true.)
  if(i1/=0) then
   new_line = new_line(i1+1:)
   new_line = adjustl(new_line)
   if(new_line(1:5) == 'x,y,z') then
    sym_op_nb = 1
   else
    call get_sym_op_nb(new_line, sym_op_nb)
   endif
  else
   return
  end if 

  
  if(tx==0 .and. ty==0 .and. tz==0 .and. sym_op_nb == 1) then
   write(CIF_unit,'(a,2x,a,2x,F10.5,5x,a)') label_1(1:4), label_2(1:4), distance, '  . ?'
  else
   write(CIF_unit,'(a,2x,a,2x,F10.5,a,i3,a,3I1,a)') label_1(1:4), label_2(1:4), distance,  &
                        " ", sym_op_nb, "_", int(tx+5.), int(ty+5.), int(tz+5.), " ?"

  endif
  
 		
 return		
end subroutine write_dist_CIF		
!--------------------------------------------------------

subroutine get_sym_op_nb(input_string, sym_op_nb)
 use cryscalc_module, only : SPG
 implicit none
  character (len=*), intent(in)    :: input_string
  integer,           intent(out)   :: sym_op_nb
  integer                          :: i, long_1, long_2

 
 long_1 = len_trim(input_string)
 
 do i = 1, SPG%Multip
  long_2 = len_trim(SPG%SymopSymb(i))
  if (long_1 == long_2 .and. input_string(1:long_1) == SPG%SymopSymb(i) (1:long_2)) then
   sym_op_nb = i
   exit
  else
  end if  
 end do
 
 return
end subroutine get_sym_op_nb
!--------------------------------------------------------


subroutine write_ang_CIF(input_line)
 use cryscalc_module, only : CIF_unit

 implicit none
 character (len=*), intent(in) :: input_line
 character (len=256)           :: new_line
 integer                       :: i1, i2, i3
 character (len=8)             :: label_1, label_2, label_3
 character (len=16)            :: ang_string

 new_line = input_line
 
 ! label 1
 i1 = index(new_line, '(')
 i2 = index(new_line, ')')
 if(i1 /=0 .and. i2/=0 .and. i2>i1) then
  label_1 = new_line(i1+1:i2-1)
  label_1 = adjustl(label_1)
   new_line = new_line(i2+1:)
 else
  return
 endif
 
  ! label 2
 i1 = index(new_line, '(')
 i2 = index(new_line, ')')
 if(i1 /=0 .and. i2/=0 .and. i2>i1) then
  label_2 = new_line(i1+1:i2-1)
  label_2 = adjustl(label_2)
  new_line = new_line(i2+1:)
 else
  return
 endif

   ! label 3
 i1 = index(new_line, '(')
 i2 = index(new_line, ')')
 if(i1 /=0 .and. i2/=0 .and. i2>i1) then
  label_3 = new_line(i1+1:i2-1)
  label_3 = adjustl(label_3)
  new_line = new_line(i2+1:)
 else
  return
 endif
 
 ! angle
 i1 = index(new_line, ':')
 if(i1==0) return
 ang_string = new_line(i1+1:)
 ang_string = adjustl(ang_string)
  
 
 write(CIF_unit,'(a,2x,a,2x,a,2x,a, 2x,a)') label_1(1:4), label_2(1:4), label_3(1:4), ang_string(1:12), '.       .       ?'

 return
end subroutine write_ang_CIF


!--------------------------------------------------------
   