program test_reflections

  use CFML_GlobalDeps
  use CFML_Crystal_Metrics
  use CFML_SuperSpaceGroups
  Use CFML_String_Utilities, only : Pack_String
  Use CFML_Math_General,     only : invert_matrix
  implicit none

  character(len=80)             :: line,fmto,fileout
  type(sReflection_List_Type)   :: ref
  type(Crystal_Cell_Type)       :: cell
  real(kind=cp), dimension(3)   :: abc,albega
  real(kind=cp), dimension(3,6) :: kv
  real(kind=cp)                 :: sintlmax=0.5
  integer                       :: nk,i,Dd,i_out
  integer,       dimension(6)   :: nharm=1
  real(kind=cp), dimension(6)   :: sintl

  fileout="Test_Reflections.out"
  open(newunit=i_out,file=fileout,status="replace",action="write")
  write(i_out,"(/,a)") " ======================================="
  write(i_out,"(a)")   " Program: Test_Reflections in SuperSpace"
  write(i_out,"(a,/)") " ======================================="
  do
  	write(*,"(a)",advance="no") " => Enter cell parameters: "
  	read(*,"(a)") line
  	if(len_trim(line)==0) exit
  	read(line,*) abc,albega
  	call set_crystal_cell(abc,albega,cell)
  	call write_crystal_cell(cell,i_out)

  	write(*,"(a)",advance="no") " => Enter the number of propagation vectors (<= 8): "
  	read(*,"(a)") line
  	if(len_trim(line)==0) exit
  	read(line,*) nk
  	write(i_out,"(/,a,i3)") " => Number of Propagation Vectors: ",nk
  	Dd=3
  	fmto="( i4,f12.4)"
  	if(nk /= 0) then
  		Dd=Dd+nk
  		do i=1,nk
  			write(*,"(a,i2,a)",advance="no") " => Enter propagation vector # ",i," number of harmonics and maximum SinTheta/Lambda: "
  			read(*,*) kv(:,i), nharm(i),sintl(i)
  			write(i_out,"(a,i2,a,3f10.5,a,i2,a,f10.5)") "     Propagation vector # ",i," : [",kv(:,i)," ], Number of harmonics:  ",nharm(i), "   Max_SinT/L: ",sintl(i)
  		end do
  		call  Gen_SReflections(Cell,sintlmax,Ref%nref,Ref%ref,nk,nharm,kv,sintl)
  	else
  		call  Gen_SReflections(Cell,sintlmax,Ref%nref,Ref%ref)
  	end if
    write(fmto(2:2),"(i1)") Dd
    write(i_out,"(/,a,i8)")" => Total number of reflections: ", Ref%nref 
    write(i_out,"(/,a)")   "       List of Reflections  "
    write(i_out,"(a,/)")   "       ===================  "
  	do i=1,Ref%nref
  		write(i_out,fmto) Ref%ref(i)%h,Ref%ref(i)%s
  	end do
  	
  end do

end program test_reflections