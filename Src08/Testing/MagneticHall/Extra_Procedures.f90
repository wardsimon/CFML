!!----
!!----
!!----
Module CFML_Extra
   !---- Use Modules ----!
   Use CFML_Maths, only: determ
   Use CFML_Rational
   Use CFML_Symmetry_Tables
   Use CFML_gSpaceGroups
   Use CFML_Strings
       
   !---- Variables ----!
   character(len=132), dimension(:), allocatable :: ffile
   
   Type :: MAG_SYMB_TYPE
      character(len=40) :: BNS
      character(len=40) :: OG
      character(len=40) :: STD
      character(len=40) :: HallM
      character(len=15) :: ID_BNS
      character(len=15) :: ID_OG
   End Type MAG_SYMB_TYPE
   type(Mag_Symb_Type), dimension(1651) :: Mag_Symb
   
   Contains
   
   !!----
   !!----
   !!----
   !!---- 15/05/19
   !!
   Subroutine Read_Magnetic_Symbols()
      !---- Local Variables ----!
      integer                          :: i, n, nlines, ic, iv
      
      
      !> Reading File
      nlines=Number_Lines('MSGsymbols_proposal.txt')
      if (nlines ==0) then
         print*,'>>> 0 lines readed in MSGsymbols_proposal.txt'
         return
      end if 
      if (nlines /= 1652) then
         print*,">>> The number of lines readed was not correct!"
         return
      end if   
      
      if (allocated(ffile)) deallocate(ffile)
      allocate(ffile(nlines))
      ffile=' '
      call reading_lines('MSGsymbols_proposal.txt',nlines,ffile) 
      
      do i=2,nlines
         read(unit=ffile(i)(61:70), fmt='(i4)') n
         if (n <= 0 .or. n > 1651) then
            print*,">>> N. IT was wrong!"
            exit
         end if   
         Mag_symb(n)%BNS=adjustl(ffile(i)(16:30))    
         Mag_symb(n)%OG =adjustl(ffile(i)(46:60))
         Mag_symb(n)%STD=adjustl(ffile(i)(76:105))   
         Mag_symb(n)%HallM=adjustl(ffile(i)(106:))
         Mag_symb(n)%ID_BNS=adjustl(ffile(i)(1:15))
         Mag_symb(n)%ID_OG=adjustl(ffile(i)(31:45))
      end do   
      
      deallocate(ffile)
   End Subroutine Read_Magnetic_Symbols
   
End Module CFML_Extra