!!----
!!----
!!----
Submodule (CFML_gSpaceGroups) SPG_023
   Contains
   !!----
   !!---- GET_SUBGROUPS_SUBGEN
   !!----
   !!---- 20/04/19
   !!
   Module Subroutine Get_SubGroups_Subgen(SpG, SubG, nsg, indexg)
      !---- Arguments ----!
      type(Spg_Type),                   intent( in) :: SpG
      type(Spg_Type),dimension(:),      intent(out) :: SubG
      integer,                          intent(out) :: nsg
      integer,                 optional,intent(in)  :: indexg
      
      !--- Local variables ---!
      integer  :: i,L,j,k,d, nc, mp,maxg,ngen,nla,n,nop,idx,ng       
      character (len=40), dimension(:),allocatable :: gen,gen_min
      character (len=40), dimension(30)            :: gen_lat
      character (len=256),dimension(:),allocatable :: list_gen
      character (len=40)                           :: gen_cent, gen_aux

      type(Symm_Oper_Type)                         :: Op_cent, Op_aux
      type(Symm_Oper_Type), dimension(30)          :: Op_lat


      !> Test if generators are available
      !> construct a procedure for selecting the minimal set of generators
      if (len_trim(SpG%generators_list) ==  0) then 
         Err_CFML%Ierr = 1
         Err_CFML%Msg  = "Get_Subgroups_Subgen@SPACEG: A list of generators of the parent group is needed!"
         return
      end if
       
      maxg=size(SubG)
      allocate(gen(SpG%multip))
      d=SpG%d
      ng=0; nc=0
      nop=SpG%numops !number of symmetry operators excluding lattice centrings & centre of symmetry
      if (SpG%centred /= 1) then
         nop=nop*2        ! number of symmetry operators excluding lattice centrings
         nc=SpG%Numops+1  ! Position of the centre of symmetry if it exist
         gen_cent=SpG%Symb_Op(nc)
         call Allocate_Symm_Op(SpG%d,Op_cent)
         Op_cent=SpG%Op(nc)  !Operator corresponding to the centre of symmetry
      end if
      
      nla=0
      if (SpG%num_lat > 0) then
         do i=1,SpG%num_lat
            ng=ng+1
            gen_lat(ng)= SpG%Symb_Op(1+nop*i)
            Op_lat(ng)= SpG%Op(1+nop*i)  !Operators corresponding to the lattice centring vectors
         end do
         nla=ng
      end if
       
      !> Formation of the list of possible generators extracted from list of generators of 
      !> the input Group
      call Get_Gener_From_Str(SpG%generators_list,d,ngen,gen_min)
      
      !>Purge the list of operators eliminating centre of symmetry and Lattice Translations
      n=0
      ng=ngen
      if (allocated(gen)) deallocate(gen)
      if (SpG%centred /= 1) ng=ngen*2
      allocate(gen(ng))
      gen=" "
      do i=1,ngen
         gen_aux=gen_min(i)
         Op_aux=Get_Op_from_Symb(gen_aux)
         if (is_inversion_centre(Op_aux)) cycle
         if (is_lattice_centring(Op_aux)) cycle
         n=n+1
         gen(n) = gen_aux
         if (SpG%centred /= 1) then
            n=n+1
            Op_aux=Op_aux*Op_cent
            gen(n)=Get_Symb_from_Op(Op_aux)
         end if
      end do
      ngen=n
      mp=2**(ngen+2)
      if (allocated(list_gen)) deallocate(list_gen)
      allocate(list_gen(mp))

      L=0
      if(ngen >= 3) then
         do i=1,ngen-2
            do j=i+1,ngen-1
               do k=j+1,ngen
                  L=L+1
                  list_gen(L)=trim(gen(i))//";"//trim(gen(j))//";"//trim(gen(k))
               end do
            end do
         end do
      end if
      
      if (ngen >= 2) then
         do i=1,ngen-1
            do j=i+1,ngen
               L=L+1
               list_gen(L)=trim(gen(i))//";"//trim(gen(j))
            end do
         end do
      end if
      
      if (ngen >= 1) then
         do i=1,ngen
            L=L+1
            list_gen(L)=trim(gen(i))
         end do
      end if
      mp=L
      
      if (SpG%num_lat > 0)  then
         do j=1,SpG%num_lat
            do i=1,mp
               L=L+1
               list_gen(L)=trim(list_gen(i))//";"//trim(gen_lat(j))
            end do
         end do
      end if
      nsg=L

      !> Now generate the subgroups
      n=0
      do L=1,nsg
         n=n+1
         call Group_Constructor(List_gen(L),SubG(n))
         if (present(indexg)) then
            idx=SpG%multip/SubG(n)%multip
            if (idx /= indexg) then
               n=max(0,n-1)
               cycle
            end if
         end if
         do i=n-1,1,-1
            if (SubG(n) == SubG(i)) then
               n=max(0,n-1)
               exit
            end if
         end do
      end do
      nsg=n
   End Subroutine Get_SubGroups_Subgen
    
End Submodule SPG_023   