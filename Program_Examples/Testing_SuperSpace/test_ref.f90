program test_reflections

  use CFML_Rational_Arithmetic
  use CFML_ssg_datafile
  use CFML_GlobalDeps
  use CFML_Crystal_Metrics
  use CFML_SuperSpaceGroups
  Use CFML_String_Utilities, only : Pack_String
  implicit none

  character(len=80)             :: line,fmto,fileout,message
  type(sReflection_List_Type)   :: ref
  type(Crystal_Cell_Type)       :: cell
  type(SuperSpaceGroup_Type)    :: SSG
  real(kind=cp), dimension(3)   :: abc,albega
  real(kind=cp), dimension(3,6) :: kv
  real(kind=cp)                 :: sintlmax=0.5
  integer                       :: nk,i,Dd,i_out,m
  integer,       dimension(6)   :: nharm=1
  real(kind=cp), dimension(6)   :: sintl
  logical :: ok

  fileout="Test_Reflections.out"
  call Read_SSG(" ",ok,message)
  if(.not. ok) then
    write(*,"(a)") "   !!! "//message//" !!!"
    stop
  end if

  open(newunit=i_out,file=fileout,status="replace",action="write")
  write(i_out,"(/,a)") " ======================================="
  write(i_out,"(a)")   " Program: Test_Reflections in SuperSpace"
  write(i_out,"(a,/)") " ======================================="


  do
  	write(*,"(a)",advance="no") " => Enter cell parameters (<cr>=exit): "
  	read(*,"(a)") line
  	if(len_trim(line)==0) exit
  	read(line,*) abc,albega
  	call set_crystal_cell(abc,albega,cell)
  	call write_crystal_cell(cell,i_out)

    write(*,"(//,a)",advance="no") " => Enter the number of the SSG: "
    read(*,*) m
    if(m <= 0) exit
    if(m > 16697) then
      write(*,"(a)") " => There are only 16697 superspace groups in the database! "
      cycle
    end if
    call Set_SSG_Reading_Database(m,SSG,ok,message)
    if(.not. ok) then
      write(*,"(a)") "   !!! "//message//" !!!"
      stop
    end if

    call Write_SSG(SSG,full=.true.)
    call Write_SSG(SSG,iunit=i_out,full=.true.)
    Dd=size(SSG%SymOp(1)%Mat,dim=1)-1
    nk=Dd-3
  	write(*,"(/,a,i3,/)") " => Number of Propagation Vectors: ",nk
  	write(i_out,"(/,a,i3,/)") " => Number of Propagation Vectors: ",nk
  	!     123456789012345678
  	fmto="(i8,tr4, i4,f12.4,tr4,2i2)"
  	if(nk /= 0) then
  		!Dd=Dd+nk
  		do i=1,nk
  			write(*,"(a,i2,a)",advance="no") " => Enter propagation vector # ",i," number of harmonics and maximum SinTheta/Lambda: "
  			read(*,*) kv(:,i), nharm(i),sintl(i)
  			write(i_out,"(a,i2,a,3f10.5,a,i2,a,f10.5)") "     Propagation vector # ",i," : [",kv(:,i)," ], Number of harmonics:  ",nharm(i), "   Max_SinT/L: ",sintl(i)
  		end do
  		call  Gen_SReflections(Cell,sintlmax,Ref%nref,Ref%ref,nk,nharm,kv,sintl,order="yes")
      !     Gen_SReflections(Cell,sintlmax,Num_Ref,Reflex,nk,nharm,kv,maxsinl,order,SSG,powder)
      !call  Gen_SReflections(Cell,sintlmax,Ref%nref,Ref%ref,nk,nharm,kv,sintl,SSG=SSG,powder="powder")
  	else
  	  call  Gen_SReflections(Cell,sintlmax,Ref%nref,Ref%ref,order="yes")
  	end if
    write(fmto(9:9),"(i1)") Dd
    write(i_out,"(/,a,i8)")" => Total number of reflections: ", Ref%nref
    write(i_out,"(/,a)")   "       List of Reflections  "
    write(i_out,"(a,/)")   "       ===================  "
  	do i=1,Ref%nref
  		write(*,fmto)     i, Ref%ref(i)%h,Ref%ref(i)%s,Ref%ref(i)%mult,Ref%ref(i)%imag
  		write(i_out,fmto) i, Ref%ref(i)%h,Ref%ref(i)%s,Ref%ref(i)%mult,Ref%ref(i)%imag
  	end do

  end do

end program test_reflections