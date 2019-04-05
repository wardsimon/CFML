Program Testing01
    !---- Use modules ----!
    Use CFML_GlobalDeps
    
    implicit none
    
    !---- Local Variables ----!
    character(len=132) :: linea
    
    !> Header
    print*,"-------------------------------"
    print*,"----    CrysFML Testing    ----"
    print*,"---- JGP-JRC        Test01 ----"
    print*,"-------------------------------"
    
    !> Operating System
    print*, " Operative System: ",ops_name

    !> Compiler
    print*, " Fortran Compiler: ",compiler
    
    !> Precision
    if (cp == sp) then
       print*, ' Default precision for CrysFML is set to: SIMPLE'
    else
       print*, ' Default precision for CrysFML is set to: DOUBLE'
    end if        
    print*," " 
    !> Numeric
    print*, " Maximum value control: ",v_huge
    print*, " Minimum value control: ",v_tiny
    print*, " Epsilon value control: ",v_epsi
    
    
    !> DirInfo
    !do
    !   write(unit=*,fmt='(a)',advance="no") "Directory Path?:"
    !   read(unit=*,fmt='(a)') linea
    !   if (len_trim(linea)<=0) exit
    !   if (Directory_Exists(trim(linea))) then
    !      print*,'This directory exists!'
    !   else
    !      print*,"Please, check the path because I don't found this directory!"
    !   end if  
    !end do
    

End Program Testing01
