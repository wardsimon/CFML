

subroutine search_SPGR
 USE IO_module,                      ONLY : write_info
 USE cryscal_module,                 ONLY : crystal_system
 USE HKL_module
 USE CFML_symmetry_tables,           ONLY : spgr_info, Set_Spgr_Info, Remove_Spgr_Info
 USE CFML_crystallographic_symmetry, ONLY : Space_Group_Type, set_spacegroup
 USE CFML_Reflections_Utilities,     ONLY : HKL_Absent
 USE CFML_Math_General,              ONLY : Sort 
 USE CFML_string_utilities,          ONLY : l_case



 implicit none
  integer, parameter                            :: Num_Spgr_Info = 612
  integer, parameter                            :: nb_Monoc  =  15
  integer, parameter                            :: nb_Ortor  = 163
  integer, parameter                            :: nb_Hexag  = 527
  integer, parameter                            :: nb_cubic  = 554
  integer, parameter                            :: nb_tetra  = 410
  integer, parameter                            :: nb_Trigo  = 495
  integer                                       :: i, i1, i2, j, m
  integer                                       :: n_good, n_absent, numg
  real                                          :: F2_threshold
  integer,         dimension(1000)              :: num_group, num_abs
  integer(kind=4), dimension(1000)              :: merit, ptr
  character(len=12)                             :: hms,system_n
  character(len=15)                             :: hall
  character(len=256)                            :: message
  INTEGER, DIMENSION(3)                         :: ref_H


  logical                                       :: check_cent=.true.
  logical                                       :: absent
  type(Space_Group_Type)                        :: Spacegroup



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
      
      call write_info('')
      call write_info( "   ----------------------------------------------------------------------- ")
      call write_info( "     PROGRAM CHECK_GROUP: attempt to select the possible space groups from ")
      call write_info( "                          a single crystal list of structure factors       ")
      call write_info( "   ----------------------------------------------------------------------- ")
      call write_info( "     Author: J.Rodriguez-Carvajal (version 0.01, based on CrysFML)         ")
      call write_info( "   ----------------------------------------------------------------------- ")
      call write_info( "    " )
  
      n_good   = 0
      HKL_good = 0
      do i=1,n_ref 
        if (n_sig*sig_F2(i) < F2(i)) then
          HKL_good(i)=1
          n_good=n_good+1
        end if
      end do
      
      if (n_good < 5) then
       call write_info('Too few GOOD reflections! Humm..., re-measure your crystal!')
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
         m=m+1
         num_group(m)= i
         num_abs(m)  = n_absent
      end do  do_group

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
     call write_info( "   ")

     call write_info( "   ")
     write(unit=message,fmt="(a,i3,a)") " => LIST OF POSSIBLE SPACE GROUPS, a total of ",m," groups are possible"
     call write_info(trim(message))
     call write_info( "   ")
     call write_info( "     ---------------------------------------------------------------------")
     call write_info( "     Number(IT)      Hermann-Mauguin Symbol     Hall Symbol       Absences")
     call write_info( "     ---------------------------------------------------------------------")
     call write_info( "   ")

     call sort(num_abs,m,ptr)
     do i=m,1,-1
        j=num_group(ptr(i))
        hms=adjustl(spgr_info(j)%HM)
        hms(2:)=l_case(hms(2:))
        hall=spgr_info(j)%hall
        numg=spgr_info(j)%N
        write(unit=message,fmt="(i10,4a,i8)")  numg,"                   ",hms,"          ",hall,num_abs(ptr(i))
        call write_info(trim(message))
     end do

      
 return
end subroutine search_SPGR