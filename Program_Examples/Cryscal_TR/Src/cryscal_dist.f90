!     Last change:  TR   17 Jul 2007    4:50 pm

subroutine  calcul_distances
  use cryscal_module, ONLY       : nb_dist_calc, new_coord, atom1_dist, atom2_dist, SP_value, nb_atom, message_text, &
                                   keyword_DIST_, dist_coef
  USE IO_module

  implicit none
  integer                 :: i
  real                    :: distance
  REAL, DIMENSION(3)      :: atom_1_coord, atom_2_coord
  REAL, DIMENSION(3)      :: delta
  CHARACTER(LEN=16)       :: label_1, label_2
  LOGICAL                 :: ok


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

   IF(keyword_dist_) then
    delta(1) =  atom_2_coord(1) - atom_1_coord(1)
    delta(2) =  atom_2_coord(2) - atom_1_coord(2)
    delta(3) =  atom_2_coord(3) - atom_1_coord(3)
    call write_info('')
    WRITE(message_text,'(15x,5a, F10.5)')  '     Delta_x(', label_1(1:6), ' - ', label_2(1:6), ') = ', delta(1)
    call write_info(TRIM(message_text))
    WRITE(message_text,'(15x,5a, F10.5)')  '     Delta_y(', label_1(1:6), ' - ', label_2(1:6), ') = ', delta(2)
    call write_info(TRIM(message_text))
    WRITE(message_text,'(15x,5a, F10.5)')  '     Delta_z(', label_1(1:6), ' - ', label_2(1:6), ') = ', delta(3)
    call write_info(TRIM(message_text))

    WRITE(message_text,'(10x,a, 3F10.5)')  '     . Coordinates of the new atom:', atom_1_coord(1) + dist_coef * delta(1), &
                                                                                  atom_1_coord(2) + dist_coef * delta(2), &
                                                                                  atom_1_coord(3) + dist_coef * delta(3)
    call write_info(TRIM(message_text))

   endif
  end do






  return
end subroutine    calcul_distances


!------------------------------------------------------------------------------
subroutine distance_calculation(coord_1, coord_2, distance)
 use cryscal_module, only : SP_value
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
 use cryscal_module,         only      :  unit_cell , SP_value
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
  use cryscal_module,                 ONLY : keyword_create_CIF, CIF_DIST, CONN_all, nb_atom, atom_CONN,  &
                                             CONN_dmax, SPG, crystal_cell, Atm_list, Atom_label, Atom_type,  &
                                             Atom_Biso, Atom_occ, Atom_occ_perc, Atom_coord, unit_cell, space_group_symbol
  use CFML_geometry_calc,             ONLY : calc_dist_angle 
  use CFML_atom_TypeDef,              ONLY : allocate_Atom_list, deallocate_Atom_list
  use CFML_BVS_energy_calc
  USE CFML_Crystal_Metrics,           ONLY : Set_Crystal_Cell
  USE CFML_crystallographic_symmetry, ONLY : set_spacegroup
  USE IO_module
  USE macros_module,                  ONLY : u_case


  implicit none
   integer                            :: i, i1, nd
   integer                            :: i_error
   character(len=256)                 :: read_line, adjusted_line
   character(len=6)                   :: current_atom
   logical                            :: atom_ok
   
   
  !
  if(.not. CONN_all) then
   atom_ok = .false.
   do i=1, nb_atom
    if(u_case(atom_CONN%label) == u_case(atom_label(i))) then    
     ATOM_CONN%type = atom_type(i)
     atom_ok = .true.
     exit
    end if
   end do
   if(.not. atom_ok) then
    call write_info('')
    call write_info(' Wrong atom for connectivity calculations !')
    call write_info('')
    return
   endif
  endif 
   
   
  !--- write CIF ---------------------------------------------------------------------
   if(keyword_create_CIF) then
    if (allocated(CIF_DIST%text)) deallocate(CIF_DIST%text)
    allocate(CIF_DIST%text(CIF_DIST%max_text)) !Maximum number of distances
    CIF_DIST%text(:)(1:132)=" "
    call write_CIF_file('DIST')
   end if   
  !-----------------------------------------------------------------------------------



!  --- new : jan. 2010 -----------------
!  idem Bond_str

   ! construction de l'objet Crystal_CELL
  call Set_crystal_cell(unit_cell%param(1:3), unit_cell%param(4:6), crystal_cell)
   
   ! construction de l'objet SPG
  call set_spacegroup(space_group_symbol, SPG)

   ! construction de l'objet Atm_list
   Atm_list%natoms = nb_atom
   if(allocated(Atm_list%atom))  deallocate(Atm_list%atom)
   allocate (Atm_list%atom(nb_atom))
 
   do i=1, Atm_list%natoms
    Atm_list%atom(i)%Lab      = atom_label(i)
    Atm_list%atom(i)%ChemSymb = atom_type(i)
    Atm_list%atom(i)%Biso     = atom_Biso(i)
    Atm_list%atom(i)%occ      = atom_occ(i)
    Atm_list%atom(i)%x(1:3)   = atom_coord(1:3, i) 
   end do
   
   close(unit=55)
   open(unit=55, file="bond_str.out")
    call Calc_dist_angle(CONN_dmax, 0., crystal_cell, SPG, Atm_list, 55)
    call Calc_PDB(Atm_list, CONN_dmax, 55)
    deallocate (Atm_list%atom )
   close(unit=55)
   
   
   call write_info('')
   call write_info('')
   call write_info('   >> Atomic connectivity calculation (JRC-JGP Bond_Str routine) : ')
   call write_info('')
   open(unit=55, file="bond_str.out")
   do
    read(55, '(a)',iostat=i_error) read_line
    if(i_error /= 0) exit
    adjusted_line = adjustl(read_line)
    if(adjusted_line(1:10) == '----------') then
     call write_info(trim(read_line))
     if(CONN_all) then
      do
       read(55, '(a)',iostat=i_error) read_line
       if(i_error /= 0) exit      
       call write_info(trim(read_line))
      end do
      exit     
      
     else  ! 
      ! DISTANCES
      do
       read(55, '(a)', iostat=i_error) read_line
       if(i_error /=0) exit
       i = index(read_line, 'to atom:')
       
       if(i/=0) then
        read(read_line(index(read_line, ':')+1:), *) current_atom
        if(u_case(ATOM_CONN%label) == u_case(current_atom)) then
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
          call write_info(trim(read_line))
         end do
        endif 
       endif 
      end do 
     
     ! BONDS DISTRIBUTIONS                
      rewind(unit=55)
      do
       read(55, '(a)', iostat=i_error) read_line
       if(i_error /=0) exit
       i1 = index(read_line, '{--- BONDS DISTRIBUTIONS')
       if(i1 /=0) then
        backspace(unit=55)
        backspace(unit=55)
        if(.not. CONN_all) then
         do
          read(55, '(a)', iostat=i_error) read_line
          if(i_error /=0) exit
          i1 = index(read_line, 'Bond type:')
          if(i1 /=0) then
           read(read_line(index(read_line, ':')+1:), *) current_atom
           if(u_case(ATOM_CONN%type) == u_case(current_atom)) then
            call write_info(read_line)
            do 
             read(55, '(a)', iostat=i_error) read_line
             if(i_error /=0)              exit
             if(len_trim(read_line) == 0) then
              call write_info('')
              exit
             endif
             call write_info(trim(read_line))
            end do
           endif
          endif
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
   call system("del bond_str.out")

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
 USE cryscal_module, ONLY : message_text, nb_atom, atom_label, atom_coord,  &
                            SPG, keyword_create_CIF,                        &
                            atom_CONN, CONN_dmax,                           &
                            new_coord, SP_value, atom1_dist, atom2_dist, unit_cell, &
                            CIF_unit, keyword_create_CIF, CIF_DIST 
                             
 
 USE CFML_crystallographic_symmetry, only : Get_orbit, Wyckoff_Type, ApplySO
 USE CFML_GlobalDeps,        ONLY : cp


  implicit none
   REAL (kind=cp), DIMENSION(3)      :: atom_1_coord, atom_2_coord
   real (kind=cp), dimension(3)      :: r
   real, dimension(3)       :: tn
   integer                 :: i, ieq, it1, it2, it3, i_tr, m
   integer                 :: id 
   integer                 :: itx, ity, itz 
   CHARACTER(LEN=16)       :: label_1, label_2
   LOGICAL                 :: ok
   real                    :: distance
   type(wyckoff_type)      :: wyckoff
   
   integer                 :: lk, nn
   real (kind=cp),    dimension(3)         :: x0,x1
   real (kind=cp), dimension(3,Spg%multip) :: uu
   character(len=16)                 :: transla

   real , dimension(3,192)           :: orbit
   CHARACTER (LEN=80)                :: text80
   integer, dimension(192)           :: effsym  ! Pointer to effective symops
   integer, dimension(192)           :: ss_ptr  ! Pointer to stabilizer



  atom1_dist(1) = atom_CONN%label
  ok = .false.
  call get_label_atom_coord('dist', 1,1, ok)
  if(.not. ok) return
  atom_1_coord(1:3) = new_coord(1:3)
  x0(1:3) = atom_1_coord(1:3)
  
  label_1 = atom_CONN%label  
  call write_info('')
  WRITE(message_text, '(2a,5x,3F10.5)') '  >> Connectivity around : ', trim(atom_CONN%label), atom_1_coord(1:3)
  call write_info(trim(message_text))
  WRITE(message_text, '(a,F6.2,a)')  '     d_max : ', CONN_dmax, ' A'
  call write_info(trim(message_text))
  call write_info('')

  itx = int(unit_cell%param(1) / CONN_dmax) +  1
  ity = int(unit_cell%param(2) / CONN_dmax) +  1
  itz = int(unit_cell%param(3) / CONN_dmax) +  1


  ! boucle sur les differents atomes presents 
  id = 0
  do i=1, nb_atom
   atom2_dist = atom_label(i)
   label_2 = atom2_dist(i)
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
       WRITE(message_text,'(2x,i3,5a,F10.5, a, 5x,3f10.5,2x,a,3I3,3a )')    &
        id,' . d(', label_1(1:4), ' - ', label_2(1:4), ') = ',distance, ' A', atom_2_coord(:),  &
        '  T=[', int(tn(1:3)), ']  (', trim(SPG%SymopSymb(ieq)), ')'                  
       
       call write_info(TRIM(message_text))       

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
            j=coord_info%n_cooatm(n,i)
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
                  j=coord_info%n_cooatm(n,i)
                  if (Ac%atom(j)%ind(1) /= n2) cycle
                  if (n_dis+1 > Max_Dis) then
                     write(unit=lun,fmt="(a,i6/)")  "Overflow the number of distances the program can handle: ",max_dis
                     deallocate(bond_spec)
                     call Deallocate_atom_list(Ac)
                     return
                  end if
                  n_dis=n_dis+1
                  disbond(n_dis)=coord_info%dist(n,i)
               end do
            end do

            call sort(disbond,n_dis,indx)

            write(unit=lun,fmt='(/,a,i4)') " Bond type: "//species(n1)//" - "//species(n2)//"      Num: ",n_dis
            write(unit=lun,fmt='(a)')   "  Num. Bonds         Distance"

            n_ini=1
            dis_ini=disbond(indx(n_ini))
            do i=2,n_dis
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
