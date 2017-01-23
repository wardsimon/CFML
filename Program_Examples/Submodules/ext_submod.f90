 Program extract_submodules
   use CFML_GlobalDeps
   use CFML_String_Utilities, only: l_case,u_case
   use CFML_IO_Formats, only: file_list_type, File_To_FileList, Err_Form, Err_Form_Mess
   implicit none
   character(len=*),dimension(3),parameter :: remov=["del   ","rm -f ","rm -f "]
   character(len=256) :: line, Cmdline, mod_name,fileins,fileout,proc_name,aux_name
   character(len=60), dimension(:),   allocatable :: fun_names
   character(len=60), dimension(:),   allocatable :: sub_names
   character(len=60), dimension(:),   allocatable :: int_names
   integer,           dimension(:),   allocatable :: point_to_endfun,point_to_endsub
   integer,           dimension(:),   allocatable :: point_to_finterface
   integer,           dimension(:),   allocatable :: point_to_sinterface
   integer,           dimension(:,:), allocatable :: point_to_fun,point_to_sub
   integer,           dimension(:,:), allocatable :: point_to_subtxt
   integer,           dimension(:,:), allocatable :: point_to_funtxt
   type(file_list_type) :: file_list

   integer :: i,j,k,nfun,nsub,i_contains,n_interf
   logical :: existe
   !
   call Get_Command(Command=Cmdline)
   call Get_Command_Argument(1,fileins)
   inquire(file=trim(fileins),exist=existe)
   write(*,"(a)")' ---------------------------'
   write(*,"(a)")' PROGRAM: Extract_Submodules'
   write(*,"(a)")' ---------------------------'
   if (.not. existe) then
      write(*,"(a)")' Input File not found!'
      call Stop_Program()
   end if
   write(*,"(a)") " Converting "//trim(fileins)//"  to file_list"
   call File_To_FileList(fileins, file_list)
   call Set_Pointers()
   call Write_Main_Module()
   call Create_Submodules()

   Contains

    Subroutine Set_Pointers()
      integer :: n,i_subtxt,i_sub,i_endsub,mss,mse,mff,mfe
      integer :: mf,ms,i_funtxt,i_fun,i_endfun,ipar,i_int,nlong
      logical :: new_text, done_cont, done_mod
      ms=0; mf=0; done_cont=.false.; done_mod=.false.
      do n=1,file_list%nlines
        line=l_case(adjustl(file_list%line(n)))
        if(line(1:6) == "module" .and. .not. done_mod) then
           nlong=len_trim(file_list%line(n))
           i_fun=index(file_list%line(n)(1:nlong)," ",back=.true.)
           mod_name=file_list%line(n)(i_fun+1:)
           done_mod=.true.
        end if
        if(line(1:1) == "!") cycle
        if(line(1:8) == "contains" .and. .not. done_cont) then
          i_contains=n
          done_cont=.true.
          cycle
        end if
        if(line(1:8) == "function") then
          mf=mf+1; cycle
        end if
        if(line(1:4) == "pure" .and. index(line,"function") /= 0) then
           mf=mf+1; cycle
        end if
        if(line(1:9) == "elemental" .and. index(line,"function") /= 0) then
          mf=mf+1;cycle
        end if
        if(line(1:9) == "recursive" .and. index(line,"function") /= 0) then
          mf=mf+1;cycle
        end if
        if(line(1:10) == "subroutine") then
          ms=ms+1; cycle
        end if
        if(line(1:4) == "pure" .and. index(line,"subroutine") /= 0) then
          ms=ms+1; cycle
        end if
        if(line(1:9) == "recursive" .and. index(line,"subroutine") /= 0) then
          ms=ms+1; cycle
        end if
        if(line(1:9) == "elemental" .and. index(line,"subroutine") /= 0) ms=ms+1
      end do

      n_interf=0
      do n=1,i_contains-1
        line=l_case(adjustl(file_list%line(n)))
        if(line(1:1) == "!") cycle
        if(line(1:9) == "interface") then
          n_interf=n_interf+1
        end if
      end do

      allocate(fun_names(mf),sub_names(ms),point_to_sub(2,ms),    &
               point_to_endsub(ms),point_to_subtxt(2,ms), &
               point_to_fun(2,mf),   point_to_endfun(mf), &
               point_to_funtxt(2,mf),int_names(n_interf), &
               point_to_finterface(mf),point_to_sinterface(ms))

      point_to_sub=0; point_to_endsub=0; point_to_subtxt=0
      point_to_fun=0; point_to_endfun=0; point_to_funtxt=0
      point_to_finterface=0; point_to_sinterface=0
      i_int=0
      do n=1,i_contains-1
        line=l_case(adjustl(file_list%line(n)))
        if(line(1:1) == "!") cycle
        if(line(1:9) == "interface") then
          i_int=i_int+1
          line=adjustl(file_list%line(n))
          i=index(line," ")
          int_names(i_int)=adjustl(line(i:))
          !write(*,"(i4,a)") i_int,"  "//trim(int_names(i_int))
        end if
      end do

      mss=0; mse=0; mff=0; mfe=0
      new_text=.true.


      do n=i_contains+1,file_list%nlines
        line=l_case(adjustl(file_list%line(n)))
        !write(*,"(a,6i3)") " mss,mse,mff,mfe",mss,mse,mff,mfe
        if(line(1:3) == "end") then

          i_endfun=index(line(1:14)," function")
          if(i_endfun /= 0) then
            new_text=.true.
            mfe=mfe+1
            if(mfe > mf) cycle
            point_to_endfun(mfe)=n
            cycle
          end if
          i_endsub=index(line(1:15)," subroutine")
          if(i_endsub /= 0) then
            new_text=.true.
            mse=mse+1
            if(mse > ms) cycle
            point_to_endsub(mse)=n
            cycle
          end if

        else

          if(line(1:8) == "function" .or.  &
             (line(1:4) == "pure" .and. index(line,"function") /= 0) .or. &
             (line(1:9) == "recursive" .and. index(line,"function") /= 0) .or. &
             (line(1:9) == "elemental" .and. index(line,"function") /= 0)) then
            new_text=.true.
            mff=mff+1;   if(mff > mf) cycle
            i_fun=index(file_list%line(n),"function")
            if(i_fun == 0) i_fun=index(file_list%line(n),"Function")
            ipar=index(file_list%line(n),"(")
            fun_names(mff) = adjustl(file_list%line(n)(i_fun+8:ipar-1))
            point_to_fun(1,mff)=n; cycle
          end if
          if(line(1:10) == "subroutine" .or. &
             (line(1:4) == "pure" .and. index(line,"subroutine") /= 0) .or. &
             (line(1:9) == "recursive" .and. index(line,"subroutine") /= 0) .or. &
             (line(1:9) == "elemental" .and. index(line,"subroutine") /= 0)) then
            new_text=.true.
            mss=mss+1; if(mss > ms) cycle
            i_sub=index(file_list%line(n),"subroutine")
            if(i_sub == 0) i_sub=index(file_list%line(n),"Subroutine")
            ipar=index(file_list%line(n),"(")
            sub_names(mss) = adjustl(file_list%line(n)(i_sub+10:ipar-1))
            point_to_sub(1,mss)=n
            cycle
          end if
        end if
      end do

      !Treatment of the text lines
      !Functions
      do i=1,mf
        n=point_to_fun(1,i)
        point_to_funtxt(2,i)=n-1
        proc_name=l_case(fun_names(i))
        doj1: do j=n-1,1,-1
          line=l_case(adjustl(file_list%line(j)))
          nlong=len_trim(line)
          if(line(1:4) == "!!--" .or. nlong == 0) cycle
          do k=j+1,j+80 !no more than 80 comment lines
            line=l_case(adjustl(file_list%line(k)))
            if(line(1:4) /= "!!--") cycle
            if(index(line," function ") /= 0 .and. index(line,trim(proc_name)) /= 0) then
              point_to_funtxt(1,i)=k
              exit doj1
            end if
          end do
        end do doj1
        if(point_to_funtxt(1,i) == 0) point_to_funtxt(1,i)= point_to_funtxt(2,i)
      end do
      !Subroutines
      do i=1,ms
        n=point_to_sub(1,i)
        point_to_subtxt(2,i)=n-1
        proc_name=l_case(sub_names(i))
        doj2:do j=n-1,1,-1
          line=l_case(adjustl(file_list%line(j)))
          nlong=len_trim(line)
          if(line(1:4) == "!!--" .or. nlong == 0) cycle
          do k=j+1,j+80 !no more than 80 comment lines
            line=l_case(adjustl(file_list%line(k)))
            if(line(1:4) /= "!!--") cycle
            if(index(line," subroutine ") /= 0 .and. index(line,trim(proc_name)) /= 0) then
              point_to_subtxt(1,i)=k
              exit doj2
            end if
          end do
        end do doj2
        if(point_to_subtxt(1,i) == 0) point_to_subtxt(1,i)= point_to_subtxt(2,i)
      end do
      !Complete the pointers to subroutines and functions for
      !lately create the interfaces.
      !Functions
      do i=1,mf
        n=point_to_fun(1,i)
        line=l_case(adjustl(file_list%line(n)))
        i_fun=index(line,"result")
        i_fun=index(line(i_fun+6:),"(")
        proc_name=line(i_fun+1:)
        i_fun=index(proc_name,")")
        proc_name=proc_name(1:i_fun-1)
        do j=n+1,n+80 !no more that 80 lines of interface
          line=l_case(adjustl(file_list%line(j)))
          nlong=len_trim(line)
          if(index(line,"intent") /= 0 .or. line(1:1) == "!") cycle
          if(index(line,trim(proc_name)) == 0) then !The result variable is the end of the interface
            point_to_fun(2,i)=j               !It does not work if the result clause is not given
          else if(nlong == 0 .or. line(1:1) == "!") then
            point_to_fun(2,i)=j
          else
            point_to_fun(2,i)=j-1 !?
          end if
          exit
        end do
      end do
      !Subroutines
      do i=1,ms
        n=point_to_sub(1,i)
        do j=n+1,n+80 !no more that 80 lines of interface
          line=l_case(adjustl(file_list%line(j)))
          if(line(1:1) == "!") cycle
          if(index(line,"intent") /= 0) then
            point_to_sub(2,i)=j !The last line with "intent" gives the pointer
            cycle
          end if
          exit
        end do
      end do
      nsub=mss; nfun=mff

      write(*,"(a,i3)") " Number of overload.inter.: ",n_interf
      write(*,"(a,i3)") " Number of functions      : ",nfun
      write(*,"(a,i3)") " Number of end functions  : ",mfe
      write(*,"(a,i3)") " Number of subroutines    : ",nsub
      write(*,"(a,i3)") " Number of end subroutines: ",mse
      write(*,*) "  FUNCTIONS"
      do i=1,mff
        write(*,"(a)") trim(file_list%line(point_to_fun(1,i)))//"  => "//trim(file_list%line(point_to_endfun(i)))//&
        "  => "//trim(file_list%line(point_to_funtxt(1,i)))//"  => "//trim(fun_names(i))
      end do
      write(*,*) "  SUBROUTINES"
      do i=1,mss
        write(*,"(a)") trim(file_list%line(point_to_sub(1,i)))//"  => "//trim(file_list%line(point_to_endsub(i)))// &
        "  => "//trim(file_list%line(point_to_subtxt(1,i)))//"  => "//trim(sub_names(i))
      end do
      !Associate each function and subroutine to one of the overloaded procedure
      !The association is done by the principal name of the overloaded procedure
      do i=1,n_interf
         proc_name=l_case(int_names(i))
         !write(*,*) trim(proc_name)
         do j=1,nfun
          aux_name=l_case(fun_names(j))
          k=index(aux_name,trim(proc_name))
          !write(*,*) " Testing:    "//trim(aux_name)//" with "//trim(proc_name), k
          if(k /= 0) then
            point_to_finterface(j)=i
            !write(*,"(a,i4,2i8)") "  OK =>   "//trim(fun_names(j)),point_to_finterface(j),point_to_fun(1,j),point_to_endfun(j)
          end if
         end do
         do j=1,nsub
          aux_name=l_case(sub_names(j))
          k=index(aux_name,trim(proc_name))
          !write(*,*) " Testing:    "//trim(aux_name)//" with "//trim(proc_name), k
          if(k /= 0) then
            point_to_sinterface(j)=i
            !write(*,"(a,i4,2i8)") "  OK =>   "//trim(sub_names(j)),point_to_sinterface(j),point_to_sub(1,j),point_to_endsub(j)
          end if
         end do
      end do
    End Subroutine Set_Pointers

    Subroutine Write_Main_Module()
      integer :: i_mod, n
      open(newunit=i_mod,file=trim(mod_name)//".f08",status="replace",action="write")
      write(*,"(a)") "  Writing the main module: "//trim(mod_name)//".f08"
      do i=1,i_contains-1
        write(unit=i_mod,fmt="(a)") trim(file_list%line(i))
      end do
      !Creation of interfaces of individual procedures
      write(unit=i_mod,fmt="(a)")
      write(unit=*,fmt="(a)")  "  Writing Function Interfaces: "
      write(unit=i_mod,fmt="(a)")  "   ! FUNCTION INTERFACES: "
      do i=1,nfun
        write(unit=i_mod,fmt="(a)")
        !Write text part
        j=point_to_funtxt(1,i); k=point_to_funtxt(2,i)
        do n=j,k
          write(unit=i_mod,fmt="(a)") trim(file_list%line(n))
        end do
        write(unit=i_mod,fmt="(a)")  "   Interface "//trim(fun_names(i))
        j=point_to_fun(1,i); k=point_to_fun(2,i)
        do n=j,k
          write(unit=i_mod,fmt="(a)")"    Module "//adjustl(trim(file_list%line(n)))
        end do
        n=point_to_endfun(i)
        write(unit=i_mod,fmt="(a)") trim(file_list%line(n))
        write(unit=i_mod,fmt="(a)") "   End Interface "//trim(fun_names(i))
      end do
      write(unit=i_mod,fmt="(a)")
      write(unit=*,fmt="(a)")  "  Writing Subroutine Interfaces: "
      write(unit=i_mod,fmt="(a)")  "   ! SUBROUTINE INTERFACES: "
      do i=1,nsub
        write(unit=i_mod,fmt="(a)")
        !Write text part
        j=point_to_subtxt(1,i); k=point_to_subtxt(2,i)
        do n=j,k
          write(unit=i_mod,fmt="(a)") trim(file_list%line(n))
        end do
        write(unit=i_mod,fmt="(a)") "   Interface "//trim(sub_names(i))
        j=point_to_sub(1,i); k=point_to_sub(2,i)
        do n=j,k
          write(unit=i_mod,fmt="(a)") "    Module "//adjustl(trim(file_list%line(n)))
        end do
        n=point_to_endsub(i)
        write(unit=i_mod,fmt="(a)") trim(file_list%line(n))
        write(unit=i_mod,fmt="(a)") "   End Interface "//trim(sub_names(i))
      end do

      write(unit=i_mod,fmt="(a)")
      do i=file_list%nlines, file_list%nlines-10,-1
        line=l_case(adjustl(file_list%line(i)))
        if(line(1:3) == "end") then
          write(unit=i_mod,fmt="(a)") trim(file_list%line(i))
          exit
        end if
      end do

      close(unit=i_mod)
    End Subroutine Write_Main_Module

    Subroutine Create_Submodules()
      character(len=256) :: directory_name, submod_name
      character(len=1)   :: key
      integer :: i_sub,i_sum
      logical :: filled
      !Create first the overloaded procedures in a single
      !submodule. Check if the subdirectory ./Module_Name exist
      !if not  abort the program
      if(directory_exists(trim(mod_name))) then
        directory_name=trim(mod_name)//OPS_SEP
      else
        !This really does not work because "system" is not run synchronously
        !Just running the program a second time (the directory is really created
        !before running ext_submod.
        call system("mkdir "//trim(mod_name))
        call system("dir *.")
        write(unit=*,fmt="(a)") "  The directory: "//trim(mod_name)//" has been created!"
        write(unit=*,fmt="(/,a)") " => Press <enter> to continue "
        read(unit=*,fmt="(a)") key
      end if
      write(unit=*,fmt="(a)") " Writing Submodules in the directory: "//trim(mod_name)
      i_sum=0
      do i=1,n_interf
         filled=.false.
         open(newunit=i_sub,file=trim(directory_name)//trim(int_names(i))//".f90",status="replace",action="write")
         write(unit=i_sub,fmt="(a)") "  Submodule ("//trim(mod_name)//") "//trim(int_names(i))
         i_sum=i_sum+1
         write(unit=*,fmt="(i5,a)") i_sum,"  Submodule ("//trim(mod_name)//") "//trim(int_names(i))
         write(unit=i_sub,fmt="(a)") "   Contains"

         do j=1,nfun
            !write(*,"(i5,a,i5)") j, "  "//trim(fun_names(j)), point_to_finterface(j)
            if(point_to_finterface(j) /= i) cycle
            k=point_to_fun(1,j)
            write(unit=i_sub,fmt="(a)") "   "
            write(unit=i_sub,fmt="(a)") "    Module "//adjustl(trim(file_list%line(k)))
            do k=point_to_fun(1,j)+1, point_to_endfun(j)
              write(unit=i_sub,fmt="(a)") trim(file_list%line(k))
            end do
            filled=.true.
         end do
         do j=1,nsub
            !write(*,"(i5,a,i5)") j, "  "//trim(sub_names(j)), point_to_sinterface(j)
            if(point_to_sinterface(j) /= i) cycle
            k=point_to_sub(1,j)
            write(unit=i_sub,fmt="(a)") "   "
            write(unit=i_sub,fmt="(a)") "    Module "//adjustl(trim(file_list%line(k)))
            do k=point_to_sub(1,j)+1, point_to_endsub(j)
              write(unit=i_sub,fmt="(a)") trim(file_list%line(k))
            end do
            filled=.true.
         end do
         write(unit=i_sub,fmt="(a)") "   "
         write(unit=i_sub,fmt="(a)") "  End Submodule "//trim(int_names(i))
         flush(unit=i_sub)
         close(unit=i_sub)
         if(.not. filled) then
          !call execute_command_line(command[,wait,exitstat,cmdstat,cmdmsg])
          call system(remov(OPS)//trim(directory_name)//trim(int_names(i))//".f90")
         end if
      end do
      !Now create the submodules corresponding to independent functions and subroutines
      !Functions
      do j=1,nfun
         if(point_to_finterface(j) /= 0) cycle
         open(newunit=i_sub,file=trim(directory_name)//trim(fun_names(j))//".f90",status="replace",action="write")
         write(unit=i_sub,fmt="(a)") "  Submodule ("//trim(mod_name)//") "//trim(fun_names(j))
         i_sum=i_sum+1
         write(unit=*,fmt="(i5,a)") i_sum,"  Submodule ("//trim(mod_name)//") "//trim(fun_names(j))
         write(unit=i_sub,fmt="(a)") "   Contains"
         k=point_to_fun(1,j)
         write(unit=i_sub,fmt="(a)") "   "
         write(unit=i_sub,fmt="(a)") "    Module "//adjustl(trim(file_list%line(k)))
         do k=point_to_fun(1,j)+1, point_to_endfun(j)
           write(unit=i_sub,fmt="(a)") trim(file_list%line(k))
         end do
         write(unit=i_sub,fmt="(a)") "   "
         write(unit=i_sub,fmt="(a)") "  End Submodule "//trim(fun_names(j))
         flush(unit=i_sub)
         close(unit=i_sub)
      end do
      !subroutines
      do j=1,nsub
         if(point_to_sinterface(j) /= 0) cycle
         open(newunit=i_sub,file=trim(directory_name)//trim(sub_names(j))//".f90",status="replace",action="write")
         i_sum=i_sum+1
         write(unit=*,fmt="(i5,a)") i_sum,"  Submodule ("//trim(mod_name)//") "//trim(sub_names(j))
         write(unit=i_sub,fmt="(a)") "  Submodule ("//trim(mod_name)//") "//trim(sub_names(j))
         write(unit=i_sub,fmt="(a)") "   Contains"
         k=point_to_sub(1,j)
         write(unit=i_sub,fmt="(a)") "   "
         write(unit=i_sub,fmt="(a)") "    Module "//adjustl(trim(file_list%line(k)))
         do k=point_to_sub(1,j)+1, point_to_endsub(j)
           write(unit=i_sub,fmt="(a)") trim(file_list%line(k))
         end do
         write(unit=i_sub,fmt="(a)") "   "
         write(unit=i_sub,fmt="(a)") "  End Submodule "//trim(sub_names(j))
         flush(unit=i_sub)
         close(unit=i_sub)
      end do
    End Subroutine Create_Submodules

    Subroutine Stop_Program()
      character(len=1) :: key
      write(unit=*,fmt="(/,a)") " => Press <enter> to finish "
      read(unit=*,fmt="(a)") key
      stop
    End Subroutine Stop_Program

 End Program extract_submodules