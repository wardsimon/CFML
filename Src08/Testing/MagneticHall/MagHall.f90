!!----
!!---- Program: Test Magnetic_Hall
!!----
!!---- JGP May2019               
!!
Program Magnetic_Hall
   !---- Use Modules ----!
   Use CFML_Extra

   !---- Variables ----!
   implicit none
   
   character(len=50)                            :: str_BNS, str_Hall
   character(len=60), dimension(:), allocatable :: gen
   integer                                      :: n, k
   integer                                      :: ik,ngen
   real(kind=cp), dimension(3)                  :: shift
   
   type(spg_type)                               :: SpG
    
   integer                                      :: i,j,nt

   !> Read Magnetic Symmetry Symbols (Proposed)
   call read_magnetic_symbols()
   
   !> Magnetic database from FullProf
   call Read_Magnetic_Data()
   if (err_CFML%IErr /= 0) then
      print*, trim(err_cfml%msg)
      stop
   end if   
   
   call system("cls")
   print*,'-----------------------------------------------------'
   print*,' Testing Magnetic Symmetry Groups using Hall symbols'
   print*,'-----------------------------------------------------'
      
   !> Main
   nt=0
   do n=1,1651
      !> Info from Database
      str_BNS=adjustl(spacegroup_label_bns(n))
         
      ik=0
      do k=1,1651
         if (trim(str_BNS) /= trim(Mag_symb(k)%BNS)) cycle
         ik=k
         exit
      end do
      if (ik == 0) then
         print*,'    ---> Problem with this BNS Symbol: '//trim(str_BNS)
         stop
      end if     
      
      hexa=.false.
      k=index(str_BNS,'6')
      if (k > 0) hexa=.true.
      
      !> Magnetic Hall symbol
      str_Hall=trim(Mag_symb(ik)%HallM)
      ngen=0
      gen=" "
      call Get_Generators_from_Hall(str_Hall, ngen, Gen, Shift)
      
      !> Constructor
      call Init_SpaceG(SpG)
      call Group_Constructor(gen,SpG)
      if (Err_CFML%Ierr /= 0) then
         print*,'    --->'//trim(Err_CFML%Msg)
         cycle
      end if   

      !> Identify group
      call Identify_Group(Spg)
      
      str_BNS=pack_string(str_BNS)
      nt=nt+1
      
      if (trim(spg%shu_symb) /= trim(str_BNS)) then
         write(unit=*,fmt='(i6,i6,t18,a,t40,a,t60,a)') nt, ik,str_BNS, spg%shu_symb, 'ERROR'
      else
         write(unit=*,fmt='(i6,i6,t18,a,t40,a,t60,a)') nt, ik, str_BNS, spg%shu_symb, 'OK'
      end if   
      
   end do   
             
End Program Magnetic_Hall
