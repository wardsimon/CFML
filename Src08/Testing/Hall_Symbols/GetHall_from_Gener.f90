!!----
!!---- Program: Recover Magnetic_Hall symbols
!!----
!!---- JGP May2019               
!!
Program Magnetic_Hall
   !---- Use Modules ----!
   Use CFML_Globaldeps
   Use CFML_gSpaceGroups

   !---- Variables ----!
   implicit none
   
   character(len=50)                            :: str_BNS, str_Hall, str
   character(len=60), dimension(:), allocatable :: gen
   integer                                      :: n, ngen
   integer, dimension(3)                        :: shift
   
   !> Magnetic Symmetry Symbols (Proposed)
   call set_Shubnikov_info()
   
   call system("cls")
   print*,'-----------------------------------------------------'
   print*,' Testing Magnetic Symmetry Groups using Hall symbols'
   print*,'-----------------------------------------------------'
      
   !> Main
   do n=1,1651   !19,66,103
      select case (n)
         case (1)
            print*,'---- Triclinic Space groups ----'
         case (8)  
            print*,' '
            print*,'---- Monoclinic Space groups ----'
         case (99)   
            print*,' '
            print*,'---- Orthorrhombic Space groups ----'
         case (661)
            print*,' '
            print*,'---- Tetragonal Space groups ----'  
         case (1231)  
            print*,' '
            print*,'---- Trigonal Space groups ----' 
         case (1339)     
            print*,' '
            print*,'---- Hexagonal Space groups ----'
         case (1503)     
            print*,' '
            print*,'---- Cubic Space groups ----'   
      end select
      
      
      !> Info from Database
      str_BNS =adjustl(shubnikov_info(n)%BNS)
      str_Hall=adjustl(shubnikov_info(n)%MHall)         
      
      !> Magnetic Hall symbol
      ngen=0
      gen=" "
      call Get_Generators_from_Hall(str_Hall, ngen, Gen)
      
      !> Test routine
      select case (n)
         case (1259:1261)
            shift=[0,0,1]
            str=Get_Hall_from_Generators(Ngen, Gen, shift)
         case (1262) 
            shift=[0,0,-1]  
            str=Get_Hall_from_Generators(Ngen, Gen, shift)         
         case (1267:1269)   
            shift=[0,0,-1]
            str=Get_Hall_from_Generators(Ngen, Gen, shift)
         case (1270) 
            shift=[0,0,1]  
            str=Get_Hall_from_Generators(Ngen, Gen, shift)
         case (1386:1390)   
            shift=[0,0,-1]  
            str=Get_Hall_from_Generators(Ngen, Gen, shift)
         case (1391:1400) 
            shift=[0,0,1]
            str=Get_Hall_from_Generators(Ngen, Gen, shift)
         case (1401:1407) 
            shift=[0,0,-1]  
            str=Get_Hall_from_Generators(Ngen, Gen, shift)
         case (1408:1409) 
            shift=[0,0,1]
            str=Get_Hall_from_Generators(Ngen, Gen, shift)      
         case default
            str=Get_Hall_from_Generators(Ngen, Gen)
      end select      
      
      if (trim(pack_string(Str)) /= trim(pack_string(str_Hall))) then
         write(unit=*,fmt='(i6,t18,a,t35,a,t60,a,t85,a)') n, trim(str_BNS), trim(str_Hall), trim(str),  'ERROR'
      else
         write(unit=*,fmt='(i6,t18,a,t35,a,t60,a,t85,a)') n, trim(str_BNS), trim(str_Hall), trim(str),  'OK'
      end if   
      
   end do   

End Program Magnetic_Hall
