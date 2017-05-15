program test_reflections

  use CFML_GlobalDeps
  use CFML_Crystal_Metrics
  use CFML_SuperSpaceGroups
  Use CFML_String_Utilities, only : Pack_String
  Use CFML_Math_General,     only : invert_matrix
  implicit none

  character(len=80)             :: line,fmto
  type(sReflection_List_Type)   :: ref
  type(Crystal_Cell_Type)       ::cell
  real(kind=cp), dimension(3)   :: abc,albega
  real(kind=cp), dimension(3,6) :: kv
  real(kind=cp)                 :: sintlmax=0.5
  integer                       :: nk,i,Dd
  integer, dimension(6)         :: nharm=2

  
  do
  	write(*,"(a)",advance="no") " => Enter cell parameters: "
  	read(*,"(a)") line
  	if(len_trim(line)==0) exit
  	read(line,*) abc,albega
  	call set_crystal_cell(abc,albega,cell)

  	write(*,"(a)",advance="no") " => Enter the number of propagation vectors (<= 8): "
  	read(*,"(a)") line
  	if(len_trim(line)==0) exit
  	read(line,*) nk
  	Dd=3
  	fmto="( i4,f12.4)"
  	if(nk /= 0) then
  		Dd=Dd+nk
  		do i=1,nk
  			write(*,"(a)",advance="no") " => Enter propagation vector # ",i," : "
  			read(*,*) kv(:,i)
  		end do
  		call  Gen_SReflections(Cell,sintlmax,Ref%nref,Ref%ref,nk,nharm,kv)
  	else
  		call  Gen_SReflections(Cell,sintlmax,Ref%nref,Ref%ref)
  	end if
    write(fmto(2:2),"(i1)") Dp
  	do i=1,Ref%nref
  		write(*,fmto) Ref%ref(i)%h,Ref%ref(i)%s
  	end do
  	
  end do

end program test_reflections