

subroutine search_SPGR(input_string)
 USE IO_module,                      ONLY : write_info
 USE cryscalc_module,                ONLY : keyword_file, crystal_system, unit_cell, Get_SPGR, check_cent, check_all_sg, &
                                            check_only_C, only_monoclinic_angle, space_group_symbol, SPG, search_mono_criteria, &
                                            debug_proc
 USE HKL_module
 USE CFML_symmetry_tables,           ONLY : spgr_info, Set_Spgr_Info, Remove_Spgr_Info
 USE CFML_crystallographic_symmetry, ONLY : Space_Group_Type, set_spacegroup
 USE CFML_Reflections_Utilities,     ONLY : HKL_Absent
 USE CFML_Math_General,              ONLY : Sort
 USE CFML_string_utilities,          ONLY : l_case
 USE Macros_module,                  ONLY : replace_car



 implicit none
  character (len=*),   intent(inout)            :: input_string
  integer, parameter                            :: Num_Spgr_Info = 612
  integer, parameter                            :: nb_Monoc  =  15
  integer, parameter                            :: nb_Ortor  = 163
  integer, parameter                            :: nb_Hexag  = 527
  integer, parameter                            :: nb_cubic  = 554
  integer, parameter                            :: nb_tetra  = 410
  integer, parameter                            :: nb_Trigo  = 495
  integer                                       :: i, i1, i2, j, m, n, i_get
  integer                                       :: long
  integer                                       :: n_good, n_absent, numg
  real                                          :: F2_threshold
  integer,         dimension(1000)              :: num_group, num_abs
  integer(kind=4), dimension(1000)              :: ptr
  character(len=12)                             :: hms,system_n, HM_s
  character(len=15)                             :: hall
  character(len=256)                            :: message
  character(len=1)                              :: centric_text
  INTEGER, DIMENSION(3)                         :: ref_H
  character (len=6), dimension(3)               :: hms_str
  logical                                       :: angle_mono_alfa, angle_mono_beta, angle_mono_gamma

  !logical                                       :: check_cent=.true.
  logical                                       :: absent
  type(Space_Group_Type)                        :: Spacegroup
  real, parameter                               :: eps = 0.001


  if(debug_proc%level_2)  call write_debug_proc_level(2, "search_SPGR")

  IF(.NOT. keyword_file) then
  call write_info('')
  call write_info('  !! FILE keyword mandatory for the search group procedure !!')
  call write_info('')
  return
 endif




! cell parameters
! symmetry
! hkl file

      !---- Define SpaceGroup ----!
      call Set_Spgr_Info()


      !---- Cryst. System ----!
      Select Case(crystal_system)
         case("MONO")
            i1=nb_monoc
            i2=nb_ortor-1
            system_n="Monoclinic"
         case("ORTHO")
            i1=nb_ortor
            i2=nb_tetra-1
            system_n="Orthorhombic"
         case("TETRA")
            i1=nb_tetra
            i2=nb_trigo-1
            system_n="Tetragonal"
         case("RHOMB", "TRIG")
            i1=nb_trigo
            i2=nb_hexag-1
            system_n="Rhombohedral"
         case("HEXA")
            i1=nb_hexag
            i2=nb_cubic-1
            system_n="Hexagonal"
         case("CUB")
            i1=nb_cubic
            i2=Num_Spgr_Info
            system_n="Cubic"
         case("TRICL", "TRIC")
            call write_info('')
            call write_info('  > triclinic case: only P1 and P-1 space groups available !')
            call write_info('')
            return
         case default
            call write_info('')
            call write_info('  > Unknown crystal system. CHECK_GROUP can not be executed !')
            call write_info('')
            return
      End Select


    ! ------ octobre 2016 : recherche de l'angle monoclinique ------------
     if(only_monoclinic_angle) then
      angle_mono_alfa  = .false.
      angle_mono_beta  = .false.
      angle_mono_gamma = .false.
      long = len_trim(system_n)
      if(long == 10) then
       if(system_n(1:10) == 'Monoclinic') then
        if(abs(unit_cell%param(4) - 90.) > search_mono_criteria(1) .and. &
           abs(unit_cell%param(5) - 90.) < eps .and.  abs(unit_cell%param(6) - 90.) < eps) angle_mono_alfa  = .true.
        if(abs(unit_cell%param(5) - 90.) > search_mono_criteria(1) .and. &
           abs(unit_cell%param(4) - 90.) < eps .and.  abs(unit_cell%param(6) - 90.) < eps) angle_mono_beta  = .true.
        if(abs(unit_cell%param(6) - 90.) > search_mono_criteria(1) .and. &
           abs(unit_cell%param(4) - 90.) < eps .and.  abs(unit_cell%param(5) - 90.) < eps) angle_mono_gamma  = .true.
       end if
      end if
     end if
    ! --------------------------------------------------------------------


     if(input_string(1:3) == "out") then
      call write_info('')
      call write_info( "   ----------------------------------------------------------------------- ")
      call write_info( "     PROGRAM CHECK_GROUP: attempt to select the possible space groups from ")
      call write_info( "                          a single crystal list of structure factors       ")
      call write_info( "   ----------------------------------------------------------------------- ")
      call write_info( "     Author: J.Rodriguez-Carvajal (version 0.01, based on CrysFML)         ")
      call write_info( "   ----------------------------------------------------------------------- ")
      call write_info( "    " )
     end if

      n_good   = 0
      HKL_good = 0
      do i=1,n_ref
        if (n_sig*sig_F2(i) < F2(i)) then
          HKL_good(i)=1
          n_good=n_good+1
        end if
      end do

      if (n_good < 5) then
       if(input_string(1:3) == "out") then
        call write_info('Too few GOOD reflections! Humm..., re-measure your crystal or another one!')
       end if
       return
      endif


      F2_threshold = threshold*MAXVAL(F2(1:n_ref))
      m=0

do_group: do i=i1,i2
          if ( m /= 0 ) then !Check if the number of the group is already selected
            do j=1,m
               if ( spgr_info(num_group(j))%N == spgr_info(i)%N ) cycle do_group
            end do
         end if
         hms=adjustl(spgr_info(i)%HM)
         hall=spgr_info(i)%hall
         if ( hms(1:1) /= "P" .and. .not. check_cent ) cycle do_group ! Skip centred groups
         call set_spacegroup(hall,Spacegroup,Force_Hall="y")

         n_absent=0
         do j=1,n_ref
           ref_H(1) = INT(h(j))
           ref_H(2) = INT(k(j))
           ref_H(3) = INT(l(j))
           absent=Hkl_Absent(ref_H(:), Spacegroup)
           if (absent .and. (HKL_good(j) == 0 .or. F2(j) < F2_threshold) ) n_absent=n_absent+1
           if (absent .and. F2(j) > F2_threshold .and. HKL_good(j) == 1 ) cycle do_group !Group not allowed
         end do

         ! Passing here means that all reflections are allowed in the group -> Possible group!
         !if(.not. check_cent .and. hms(1:1)/='P')   cycle do_group

         if(.not. check_all_sg) then
          if(check_cent   .and. .not. hms(1:1)/='P' .and. .not. unit_cell%H_M(1:1)/='?') cycle do_group
          if(check_only_C .and. .not. hms(1:1)/='P' .and. .not. unit_cell%H_M(1:1)/='?') cycle do_group
         end if

         m=m+1
         num_group(m)= i
         num_abs(m)  = n_absent


         !!if(unit_cell%H_M(1:1) /='?') then
         ! if(check_all_sg .or. hms(1:1) == unit_cell%H_M(1:1)) then
         !  m=m+1
         !  num_group(m)= i
         !  num_abs(m)  = n_absent
         ! !endif
         !else
         ! m=m+1
         ! num_group(m)= i
         ! num_abs(m)  = n_absent
         !endif
      end do  do_group

    if(input_string(1:3) == "out") then
     call write_info( "   ")
     write(unit=message,fmt="( a,i5,a,F4.1,a)")      " => Number of good reflections  : ",n_good, '  (F2/sig > ',n_sig,')'
     call write_info(trim(message))
     write(unit=message,fmt="( a,f12.4)")            "    Maximum intensity           : ",MAXVAL(F2(:))
     call write_info(trim(message))
     write(unit=message,fmt="( a,f12.4,a,F5.2,a)")   "    Minimum (for observed)      : ",F2_threshold, '  (',threshold,'*MAX(F2))'
     call write_info(trim(message))
     write(unit=message,fmt='(2a)')                  "    Crystal system tested       : ", system_n
     call write_info(trim(message))
     write(unit=message,fmt="( a,i3   )")            "    Number of Space Group tested: ",i2-i1+1
     call write_info(trim(message))
     if(.not. check_cent) then
      write(unit=message,fmt="( a,i3   )")            "    Search restriction          : only Primitive space groups"
     elseif(check_only_C) then
      write(unit=message,fmt="( a,i3   )")            "    Search restriction          : only centered space groups"
     else
      write(unit=message,fmt="( a,i3   )")            "    Search restriction          : no"
     end if
     call write_info(trim(message))
     call write_info( "   ")
     call write_info( "   ")

     if(m ==0) then
      write(unit=message,fmt="(a)") " => NO POSSIBLE SPACE GROUPS ! Strange !!"
      call write_info(trim(message))
      return
     else
      write(unit=message,fmt="(a,i3,a)") " => LIST OF POSSIBLE SPACE GROUPS, a total of ",m," groups are possible"
      call write_info(trim(message))
     endif
     call write_info( "   ")
     call write_info( "     ------------------------------------------------------------------------------------")
     call write_info( "     Number(IT)      Hermann-Mauguin Symbol     Hall Symbol       Absences     Centric   ")
     call write_info( "     ------------------------------------------------------------------------------------")
     call write_info( "   ")
    end if ! fin de la condition 'if input_string(1:3) == "out"

     call sort(num_abs,m,ptr)

     n = 0
     i_get = 0
     do i=m,1,-1
      j=num_group(ptr(i))
      hms=adjustl(spgr_info(j)%HM)
      hms(2:)=l_case(hms(2:))
      hall=spgr_info(j)%hall
      numg=spgr_info(j)%N

      ! ---------------------- oct. 2016 ----------------------------
      ! cas d'un systeme monoclinic :
      !  . exclusion d'un groupe si l'axe monoclinique est different
      !    de l'angle monoclinique (metrique)
      if(only_monoclinic_angle) then
      long = len_trim(system_n)
      if(long == 10) then
       if(system_n(1:10) == 'Monoclinic') then
        HM_S = replace_car(hms, "/", "_")
        read(HM_S(2:), *) hms_str(1), hms_str(2), hms_str(3)
        if(angle_mono_alfa      .and. .not. angle_mono_beta .and. .not. angle_mono_gamma) then
         if(hms_str(1)(1:1) == "1" .and. hms_str(2)(1:1) == "1") cycle
         if(hms_str(1)(1:1) == "1" .and. hms_str(3)(1:1) == "1") cycle
        elseif(angle_mono_beta  .and. .not. angle_mono_alfa .and. .not. angle_mono_gamma) then
         if(hms_str(1)(1:1) == "1" .and. hms_str(2)(1:1) == "1") cycle
         if(hms_str(2)(1:1) == "1" .and. hms_str(3)(1:1) == "1") cycle
        elseif(angle_mono_gamma .and. .not. angle_mono_alfa .and. .not. angle_mono_beta) then
         if(hms_str(1)(1:1) == "1" .and. hms_str(3)(1:1) == "1") cycle
         if(hms_str(2)(1:1) == "1" .and. hms_str(3)(1:1) == "1") cycle
        end if
       end if
      end if
      end if
      ! --------------------------------------------------------------

      !if(unit_cell%H_M(1:1) /= "?") then
      if(check_cent .and. unit_cell%H_M(1:1) /= "?") then
       if(check_only_C .and. .not. hms(1:1)/='P') cycle
       if (check_all_sg .or. hms(1:1) == unit_cell%H_M(1:1)) then  ! uniquement les groupes du même reseau de Bravais
        n = n + 1
        if (n == 1) i_get = j
        !write(unit=message,fmt="(i3,a1,i6,4a,i8)") m-i+1, '.', numg,"                   ",hms,"          ",hall,num_abs(ptr(i))
        call set_spacegroup(hms,Spacegroup)
        centric_text = "Y"
        if(Spacegroup%centred == 1) centric_text = "N"
        write(unit=message,fmt="(i3,a1,i6,4a,i8,5x,a)") n, '.', numg,"                   ",hms,"          ", hall, &
                                                        num_abs(ptr(i)), centric_text
        if(input_string(1:3) == "out") call write_info(trim(message))
       ! -------- nov. 2016 ------------
       elseif(check_only_C .and. hms(1:1)/='P') then
        n = n + 1
        if (n == 1) i_get = j
        call set_spacegroup(hms,Spacegroup)
        centric_text = "Y"
        if(Spacegroup%centred == 1) centric_text = "N"
        write(unit=message,fmt="(i3,a1,i6,4a,i8,5x,a)") n, '.', numg,"                   ",hms,"          ", hall, &
                                                        num_abs(ptr(i)), centric_text
        if(input_string(1:3) == "out") call write_info(trim(message))
       ! --------------------------------
       endif
      else
       if(.not. check_cent    .and.       hms(1:1)/='P') cycle
       if(      check_only_C  .and. .not. hms(1:1)/='P') cycle
       n = n + 1
       if (n == 1) i_get = j
       !write(unit=message,fmt="(i3,a1,i6,4a,i8)")   m-i+1, '.', numg,"                   ",hms,"          ",hall,num_abs(ptr(i))
       call set_spacegroup(hms,Spacegroup)
       centric_text = "Y"
       if(Spacegroup%centred == 1) centric_text = "N"
       write(unit=message,fmt="(i3,a1,i6,4a,i8,5x,a)") n, '.', numg,"                   ",hms,"          ", hall, &
                                                       num_abs(ptr(i)), centric_text
       if(input_string(1:3) == "out") call write_info(trim(message))
      endif
     end do


     if(get_SPGR) then
      !space_group_symbol=adjustl(spgr_info(num_group(ptr(m)))%HM)
      if(i_get /= 0) then
       space_group_symbol = adjustl(spgr_info(i_get)%HM)
       space_group_symbol(2:) = l_case(space_group_symbol(2:))
       if (input_string(1:3) == "out") call write_info( "")
       write(unit=message,fmt="(a,a)") "   => MOST PROBABLE SPACE GROUP: ", trim(Space_group_symbol)
       call write_info(trim(message))
       call write_info( "")
       call Set_spacegroup(space_group_symbol, SPG)
      else
       space_group_symbol = '?'
       call Set_spacegroup(space_group_symbol, SPG)
       call write_info( "   ")
       write(unit=message,fmt="(a)") "   => NO SPACE GROUP FOUND !!"
       call write_info(trim(message))
      end if
     endif

 return
end subroutine search_SPGR


