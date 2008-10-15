module WCFGL_glbond 
!------------------------------------------------------------------------------
  
  use WCFGL_geometry
  use WCFGL_metrix
  use WCFGL_glatom, only : gl_atom
  
  implicit none 
  
  type ::  gl_bond
    real                          :: color1(4), color2(4) 
    real, dimension(:,:), pointer :: xf       
  end type gl_bond
  
  real, parameter, private                       :: eps=0.01
  
  contains 
!------------------------------------------------------------------------------
  subroutine new_gl_bonds(atom1,atom2,minvalue,maxvalue,bond_out,bond_symbol,color)
    type(gl_atom), intent(in), dimension(:)              :: atom1, atom2
    real         , intent(in)                            :: minvalue
    real         , intent(in)                            :: maxvalue
    type(gl_bond), intent(out),dimension(:), allocatable :: bond_out
    logical      , intent(in)                            :: bond_symbol
    real         , intent(in), optional                  :: color(4)
    ! Local variables 
    integer :: h1,k1,l1,i,j, count, n1, n2, bondcount, n2max, neqmax
    real    :: blen
    type(gl_bond), dimension(:), allocatable :: bondbuffer
    real, dimension(:,:), allocatable :: buffer 
    
    bondcount=0
    
    allocate(bondbuffer(size(atom1)*size(atom2)))
    
    do n1=1, size(atom1)
      if (.not.(bond_symbol)) then
        if (atom1(1)%label==atom2(1)%label) then
          n2max=n1
        else
          n2max=size(atom2)
        end if
      else
        if (atom1(1)%symbol==atom2(1)%symbol) then
          n2max=n1
        else
          n2max=size(atom2)
        end if
      end if
      do n2=1,n2max
        count=0
        allocate(buffer(6,27*size(atom1(n1)%xf_eq)*size(atom2(n2)%xf_eq)))
        do i=1, size(atom1(n1)%xf_eq, DIM=2)
          do h1=-1,1 
            do k1=-1,1
              do l1=-1,1
                if (.norm.(atom1(n1)%xf-atom2(n2)%xf)<=eps) then 
                  neqmax=i
                else
                  neqmax=size(atom2(n2)%xf_eq, DIM=2) 
                end if
                do j=1,neqmax        
                  blen=.norm.(f2c(atom1(n1)%xf_eq(:,i)-atom2(n2)%xf_eq(:,j)-real((/h1,k1,l1/))))
                  if ((blen<=maxvalue).and.(blen>=minvalue).and.(blen>eps)) then
                    count=count+1
                    buffer(1:3,count)=atom1(n1)%xf_eq(:,i)
                    buffer(4:6,count)=atom2(n2)%xf_eq(:,j)+real((/h1,k1,l1/))
                  end if 
                end do
              end do 
            end do
          end do
        end do 
        
        if (count>0) then
          bondcount=bondcount+1
          allocate(bondbuffer(bondcount)%xf(6,count))
          bondbuffer(bondcount)%xf(:,:)=buffer(:,1:count)
        else 
          bondcount=bondcount+1
          allocate(bondbuffer(bondcount)%xf(6,1))
          bondbuffer(bondcount)%xf(:,:)=0.0
        end if 
          
        if (.not.(present(color))) then 
          bondbuffer(bondcount)%color1   =atom1(n1)%color
          bondbuffer(bondcount)%color2   =atom2(n2)%color
        else 
          bondbuffer(bondcount)%color1   =color
          bondbuffer(bondcount)%color2   =color
        end if

        if (allocated(buffer)) deallocate(buffer) 
      end do
    end do  
    
    allocate(bond_out(bondcount))
      bond_out(:)=bondbuffer(1:bondcount)
    
    if (allocated(bondbuffer)) deallocate(bondbuffer)
     
  return
    
  end subroutine new_gl_bonds
!------------------------------------------------------------------------------
end module WCFGL_glbond