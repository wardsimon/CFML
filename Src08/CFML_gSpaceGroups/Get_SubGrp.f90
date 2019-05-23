SubModule (CFML_gSpaceGroups) SPG_022
   Contains
   
   !!----
   !!---- GET_SUBGROUPS
   !!----
   !!----   This subroutine provides all subgroups of a generic space group
   !!----   It is supposed that the operators have been already ordered, so
   !!----   that if the group is centrosymmetric the position of the corresponding
   !!----   operator is Numops+1, moreover if it is centred the symmetry operators
   !!----   are orderd as:
   !!----   [1...Numops] {-1|t}*[1...Numops] {1|t1}*[1...Numops] {-1|t}*[1 ... Numops]}
   !!----   {1|t2}*[1...Numops] {-1|t}*[1 ... Numops]} ....
   !!----   where t1,t2, ... are the centring translations
   !!----   This ordering facilitates the calculation of the major part of subgroups
   !!----   (excluding those for which less centring translations are considered)
   !!----   If indexg is present only the subgroups of index=indexg are output
   !!----
   !!---- 20/04/2019
   !!
   Module Subroutine Get_SubGroups(SpG, SubG, nsg, indexg, point)
      !---- Arguments ----!
      type(Spg_Type),                  intent( in) :: SpG
      type(Spg_Type),dimension(:),     intent(out) :: SubG
      integer,                          intent(out) :: nsg
      integer,                 optional,intent(in)  :: indexg
      logical, dimension(:,:), optional,intent(out) :: point
      
      !--- Local variables ---!
      integer  :: i,L,j,k,d, nc,maxg,ng,kb, nla, nop,ns_1,ns_2,n_nc_group 
      logical  :: newg 
      character (len=40), dimension(:),allocatable :: gen
      character (len=40), dimension(30)            :: gen_lat
      character (len=40)                           :: gen_cent

      type(Symm_Oper_Type)                         :: Op_cent
      type(Symm_Oper_Type), dimension(30)          :: Op_lat
      type(Spg_Type), dimension(:), allocatable   :: sG
      integer, dimension(size(SubG))               :: index_sg, ind

      !> Init
      maxg=size(SubG)
      allocate(gen(SpG%multip))
      d=SpG%d
      
      !> Construct first the generators of centring translations
      ng=0; nc=0
      nop=SpG%numops !number of symmetry operators excluding lattice centrings & centre of symmetry
      if (SpG%centred /= 1) then
         nop=nop*2 !!number of symmetry operators excluding lattice centrings
         nc=SpG%Numops+1  !Position of the centre of symmetry if it exist
         gen_cent=SpG%Symb_Op(nc)
         call Allocate_Symm_Op(SpG%d,Op_cent)
         Op_cent=SpG%Op(nc)
      end if
      nla=0
      if (SpG%num_lat > 0) then
         do i=1,SpG%num_lat
            ng=ng+1
            gen_lat(ng)= SpG%Symb_Op(1+nop*i)
            Op_lat(ng)= SpG%Op(1+nop*i)
         end do
         nla=ng
      end if
      
      !> First work with the Numops operators to determine the subgroups, the other subgroups
      !> will be obtained adding progressively the rest of generators (centre of symmetry and
      !> lattice centrings.
      L=0
      
      !> Determine first the groups with only one rotational generator
      ng=1
      do i=2,SpG%numops
         gen(1) = SpG%Symb_Op(i)
         L=L+1
         if (L > maxg) then
            nsg=maxg
            return
         end if
         newg=.true.
         call Group_Constructor(gen(1:ng),SubG(L))
         if (Err_CFML%Ierr == 0) then
            do k=1,L-1
               if (SubG(L) == SubG(k)) then
                  newg=.false.
                  exit
               end if
            end do
            if (.not. newg) L=L-1
         else
            L=L-1
         end if
      end do
      ns_1=L
      
      !> Determine now the groups with two rotational generators
      if (SpG%numops > 2) then
         ng=2
         do i=2,SpG%numops-1
            gen(1) = SpG%Symb_Op(i)
            do j=i+1,SpG%numops
               gen(2)=SpG%Symb_Op(j)
               L=L+1
               if (L > maxg) then
                  nsg=maxg
                  return
               end if
               newg=.true.
               call Group_Constructor(gen(1:ng),SubG(L))
               if (Err_CFML%Ierr == 0) then
                  do k=1,L-1
                     if (SubG(L) == SubG(k)) then
                        newg=.false.
                        exit
                     end if
                  end do
                  if (.not. newg) L=L-1
               else
                  L=L-1
               end if
            end do
         end do
         ns_2=L-ns_1
      end if
      nsg=L
      n_nc_group=L

      !> Determine now the new groups adding a centre of symmetry (without lattice centring) if it exists
      if (SpG%centred /= 1) then !This doubles the number of groups
         do i=1,n_nc_group
            L=L+1
            !> List of generators for the new group
            kb= SpG%numops+1
            SubG(L)%generators_list=SubG(i)%generators_list//";"//trim(SpG%Symb_Op(kb))
            newg=.true.
            call Group_Constructor(SubG(L)%generators_list,SubG(L))
            if (Err_CFML%Ierr == 0) then
               do k=1,L-1
                  if (SubG(L) == SubG(k)) then
                     newg=.false.
                     exit
                  end if
               end do
               if (.not. newg) L=L-1
         
            else
               L=L-1
            end if
         end do
      end if
      nsg=L
      n_nc_group=L

      !> Determine now the rest of groups adding the lattice translations if they exist in the
      !> original space group
      if (SpG%num_lat > 0) then
         do j=1,nla
            do i=1,n_nc_group
               L=L+1
               SubG(L)%generators_list=SubG(i)%generators_list//";"//trim(gen_lat(j))
               newg=.true.
               call Group_Constructor(SubG(L)%generators_list,SubG(L))
               if (Err_CFML%Ierr == 0) then
                  do k=1,L-1
                     if (SubG(L) == SubG(k)) then
                        newg=.false.
                        exit
                     end if
                  end do
                  if (.not. newg) L=L-1
               else
                  L=L-1
               end if
            end do
            n_nc_group=n_nc_group+L
         end do
      end if

      nsg=L
      if (present(point)) then
         point=.false.
         do j=1,nsg
            L=1
            do i=1,SpG%multip
               do k=L,SubG(j)%multip
                  if (SubG(j)%Symb_Op(k) == SpG%Symb_Op(i)) then
                     point(i,j) = .true.
                     L=k+1
                     exit
                  end if
               end do
            end do
         end do
      end if
      
      if (present(indexg)) then
         k=0
         do L=1,nsg
            index_sg(L)=Spg%multip/SubG(L)%multip
            if (index_sg(L) == indexg) then
               k=k+1
               ind(k)=L
            end if
         end do
        
         nsg=k
         if (nsg /= 0) then
            allocate(sG(nsg))
            do k=1,nsg
               sG(k)=SubG(ind(k))
            end do
          
            do L=1,nsg
               SubG(L)=sG(L)
            end do
         end if
      end if
   End Subroutine Get_SubGroups
   
End SubModule SPG_022  