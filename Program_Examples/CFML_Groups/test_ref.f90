program test_reflections

  use CFML_Rational_Arithmetic
  use CFML_ssg_datafile
  use CFML_GlobalDeps
  use CFML_Crystal_Metrics
  use CFML_SuperSpaceGroups
  Use CFML_String_Utilities, only : Pack_String,u_case
  use CFML_IO_Formats,       only : File_To_FileList,File_List_Type
  implicit none

  character(len=80)             :: line,fmto,fileout,message,cfl_file,filcod,database_path
  type(sReflection_List_Type)   :: ref
  type(Crystal_Cell_Type)       :: cell
  type(SuperSpaceGroup_Type)    :: SSG
  type(File_List_Type)          :: File_list
  type(kvect_info_type)         :: kinfo
  real(kind=cp), dimension(3)   :: abc,albega,h
  !real(kind=cp), dimension(3,6) :: kv
  real(kind=cp), dimension(3,6) :: kcomp
  real(kind=cp)                 :: sintlmax=0.5,tini,tfin
  integer                       :: nk,i,Dd,i_out,m,j,narg
  !integer,       dimension(6)   :: nharm=1
  !real(kind=cp), dimension(6)   :: sintl
  logical :: ok,arggiven,powder,mag=.true.

  fileout="Test_Reflections.out"
  database_path="c:\CrysFML\Program_Examples\Testing_SuperSpace"
  call Read_SSG(trim(database_path),ok,message)
  if(.not. ok) then
    write(*,"(a)") "   !!! "//message//" !!!"
    stop
  end if
  narg=COMMAND_ARGUMENT_COUNT()
  powder=.false.
  call cpu_time(tini)
  if(narg > 0) then
      call GET_COMMAND_ARGUMENT(1,filcod)
      i=index(filcod,'.cfl',back=.true.)
      if( i /= 0) then
        cfl_file=filcod
        filcod=filcod(1:i-1)
      else
        cfl_file=trim(filcod)//".cfl"
      end if
      fileout=trim(filcod)//"_ref.out"
      call File_To_FileList(cfl_file,File_list)
      call Readn_Set_SuperSpace_Group(database_path,File_list%line,cell,SSG,kinfo)

      if(Err_ssg) then
        write(unit=*,fmt="(a)") trim(Err_ssg_mess)
        stop
      end if
      arggiven=.true.
      nk=SSG%nk
      Dd=SSG%nk+3
      do i=1,File_list%nlines
        line=u_case(adjustl(File_list%line(i)))
        if(line(1:6) == "POWDER") then
          powder=.true.
          exit
        end if
      end do
      call Write_SSG(SSG,full=.true.,kinfo=kinfo)
      call k_SSG_compatible(SSG,kcomp)
      do i=1,6
        write(*,"(i4,a,3f10.5,a)") i," [",kcomp(:,i)," ]"
      end do
  end if

  open(newunit=i_out,file=fileout,status="replace",action="write")
  write(i_out,"(/,a)") " ======================================="
  write(i_out,"(a)")   " Program: Test_Reflections in SuperSpace"
  write(i_out,"(a,/)") " ======================================="


  do

    if(.not. arggiven) then
       write(*,"(a)",advance="no") " => Enter cell parameters (<cr>=exit): "
       read(*,"(a)") line
       if(len_trim(line)==0) exit
       read(line,*) abc,albega
       call set_crystal_cell(abc,albega,cell)
       write(*,"(//,a)",advance="no") " => Enter the number of the SSG: "
       read(*,*) m
       if(m <= 0) exit
       if(m > 16697) then
         write(*,"(a)") " => There are only 16697 superspace groups in the database! "
         cycle
       end if
       call Set_SSG_Reading_Database(database_path,m,SSG,ok,message,"x1x2x3")
       if(.not. ok) then
         write(*,"(a)") "   !!! "//message//" !!!"
         stop
       end if
       call Write_SSG(SSG,full=.true.,kinfo=kinfo)
       Dd=size(SSG%Op(1)%Mat,dim=1)-1
       nk=Dd-3
       write(*,"(/,a,i3,/)") " => Number of Propagation Vectors: ",nk
       if(nk /= 0) then
        call Allocate_kvect_info(nk,kinfo)
        !Dd=Dd+nk
        do i=1,nk
          write(*,"(a,i2,a)",advance="no") " => Enter propagation vector # ",i," number of harmonics and maximum SinTheta/Lambda: "
          read(*,*) kinfo%kv(:,i), kinfo%nharm(i),kinfo%sintlim(i)
        end do
       end if
    end if
    call Write_SSG(SSG,iunit=i_out,full=.true.,kinfo=kinfo)
    fmto="(i8,tr4, i4,f12.4,tr4,2i2,tr4,3f10.5)"
    if(nk > 0) then
      !do i=1,nk
      ! write(i_out,"(a,i2,a,3f10.5,a,i2,a,f10.5)") "     Propagation vector # ",i," : [",kinfo%kv(:,i)," ], Number of harmonics:  ",kinfo%nharm(i), "   Max_SinT/L: ",kinfo%sintlim(i)
      !end do
      if(powder) then
        call  Generate_Reflections(Cell,sintlmax,Ref%nref,Ref%ref,mag,kinfo,order="yes",SSG=SSG, powder=powder)
      else
        call  Generate_Reflections(Cell,sintlmax,Ref,mag,kinfo,order="yes",SSG=SSG)
      end if
    else
      if(powder) then
        call  Generate_Reflections(Cell,sintlmax,Ref%nref,Ref%ref,mag,SSG=SSG,order="yes",powder=powder)
      else
        call  Generate_Reflections(Cell,sintlmax,Ref,mag,order="yes",SSG=SSG)
      end if
    end if
    write(fmto(9:9),"(i1)") Dd
    write(i_out,"(/,a,i8)")" => Total number of reflections: ", Ref%nref
    write(i_out,"(/,a)")   "       List of Reflections  "
    write(i_out,"(a,/)")   "       ===================  "
    do i=1,Ref%nref
      h=Ref%ref(i)%h(1:3)
      do j=1,nk
        h=h+Ref%ref(i)%h(3+j)*kinfo%kv(:,j)
      end do
      !write(*,fmto)     i, Ref%ref(i)%h,Ref%ref(i)%s,Ref%ref(i)%mult,Ref%ref(i)%imag,h
      write(i_out,fmto) i, Ref%ref(i)%h,Ref%ref(i)%s,Ref%ref(i)%mult,Ref%ref(i)%imag,h
    end do
    if(arggiven) exit
  end do
  call cpu_time(tfin)
  write(*,"(a,f10.5,a)") " CPU-time: ",tfin-tini," seconds"
end program test_reflections