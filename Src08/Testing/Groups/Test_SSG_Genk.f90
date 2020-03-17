  Module getting_ssg
    use CFML_GlobalDeps
    use CFML_gSpaceGroups
    use CFML_Rational
    use CFML_Propagation_Vectors
    use CFML_Maths,             only : Equal_Matrix
    use CFML_Symmetry_Tables,   only : Get_HM_Compact_HM
    use CFML_Strings,           only : String_Fraction_2Dig

    contains

    Subroutine Set_Intersection_Groups(SpGs,SpG)
      Type (SpG_Type),dimension(:), intent(in)   :: SpGs
      Type (SpG_Type),              intent(out)  :: SpG
      !--- Local Variables ---!
      integer :: i,j,k,ng,ipos,n
      integer,dimension(1) :: iip
      character(len=80),dimension(192) :: gen
      logical,dimension(size(SpGs(:))) :: estak

      !write(*,"(/a/)") " => Entering subroutine Set_Intersection_Groups "

      n=size(SpGs(:))
      if(n == 1) then
         SpG=SpGs(1)
         return
      end if
      iip=minloc(SpGs(:)%multip)
      ipos=iip(1)
      ng=1
      gen(1)="x,y,z,1"
      estak=.false.
      estak(ipos)=.true.

      !write(*,"(3(a,i3))") " => Number of space groups: ",n, "  Position: ",ipos,"  Multiplicity:",SpGs(ipos)%multip

      do_ext:do i=2,SpGs(ipos)%multip
        ng=ng+1
        gen(ng)=SpGs(ipos)%Symb_Op(i)
        do j=1,n
           if(j == ipos) cycle
           estak(j)=.false.
           do k=2,SpGs(j)%multip
             !write(*,"(2i4,a,a)") k,ng, "   "//trim(gen(ng))//"   "//trim(SpGs(j)%Symb_Op(k))
             if(trim(SpGs(j)%Symb_Op(k)) == trim(gen(ng))) then
               estak(j)=.true.
               exit
             end if
           end do
        end do
        !write(*,*) estak
        k=count(estak(1:n))
        if(k /= n) then
          !write(*,"(a,i3)") "  Operator: "//trim(gen(ng))//"  Discarded  ng=",ng-1
          ng=ng-1
          cycle
        end if
        !Passing here means that the operator is common to all space groups
        if(ng > 192) exit
      end do do_ext
      call Set_SpaceGroup(" ",SpG,ng,gen)

    End Subroutine Set_Intersection_groups

    ! This function returns the maximum order
    ! (defined as the exponent to become the Identity)
    ! of a group of rotational operators,optionally the
    ! order of each element is also provided

    Subroutine get_Order(SpG, maxorder,g_order)
      class(SpG_type),                           intent(in)    :: SpG
      integer,                                   intent(out)   :: maxorder
      integer,dimension(:),optional,allocatable, intent(out)   :: g_order

      integer, dimension(3,3)  :: ID, CM, TM
      integer                  :: cnt, i, j

      if(present(g_order)) then
        allocate(g_order(SpG%multip))
        g_order(1)=1
        g_order(2:)=0
      end if
      TM=SpG%Op(1)%Mat(1:3,1:3) ! Current Matrix
      CM=SpG%Op(1)%Mat(1:3,1:3) ! Initiallized to the Identity
      ID=SpG%Op(1)%Mat(1:3,1:3) ! Identity matrix is always the first of the SG

      maxorder=0

      do_ext:do i=2, SpG%Multip
         CM=SpG%op(i)%Mat(1:3,1:3)  ! Current Matrix
         TM=SpG%op(1)%Mat(1:3,1:3)  ! Test Matrix (initiallized to the Id)
         cnt=0
         do j=2,6 !order is never greater than 6
             cnt=cnt+1
             TM=matmul(TM,CM)
             if(Equal_Matrix(TM,ID,3)) then
                if(present(g_order) ) g_order(i)=cnt
                if(cnt > maxorder) maxorder=cnt
                cycle do_ext
             end if
        end do
      end do do_ext

    End Subroutine get_Order

    !!----
    !!---- Subroutine Get_Generators_From_SpGSymbol(SpG,gen,point_op,ngen)
    !!----    class(SpG_Type),                intent (in) :: SpG
    !!----    Character(len=*), dimension(:), intent(out) :: gen
    !!----    integer,          dimension(:), intent(out) :: point_op
    !!----    integer,                        intent(out) :: ngen
    !!----
    !!----    This subroutine provides the generators of the space group that
    !!----    are explicitly written in the Hermann-Mauguin symbol of the space group.
    !!----    The generators of the lattice are ignored. The generators gen(i),are written
    !!----    in the Jones faithful representation. There are ngen generators and the
    !!----    integer vector "point_op" contains the index of the corresponding operator in
    !!----    the list of the total SpG%Multip operators
    !!----
    !!----   Update: February - 2017
    !!----
    Subroutine Get_Generators_From_SpGSymbol(SpG,gen,point_op,ngen)
      class(SpG_Type),                intent (in) :: SpG
      Character(len=*), dimension(:), intent(out) :: gen
      integer,          dimension(:), intent(out) :: point_op
      integer,                        intent(out) :: ngen
      !--- Local variables ---!
      integer :: i,j,n,m
      integer,       dimension(3,3)  :: s
      real(kind=cp), dimension(3)    :: t
      character(len=20),dimension(5) :: sgen,msgen
      logical,          dimension(10):: done
      character(len=6)               :: symbg
      character(len=20)              :: spg_symb
      character(len=60)              :: op_symb
      character(len=2), dimension(11):: screw=(/"21","31","32","41","42","43","61","62","63","64","65"/)
      character(len=1), dimension(11):: rm_screw=(/"2","3","3","4","4","4","6","6","6","6","6"/)

      point_op=0; done=.false.
      if(len_trim(SpG%SPG_Symb) == 0) then
        spg_symb= Get_HM_Compact_HM(SpG%BNS_symb)
      else
        spg_symb=adjustl(SpG%SPG_Symb(2:))
      end if
      write(*,"(a)") " => Symbol of the space group: "//trim(spg_symb)
      i=index(spg_symb,":")
      if( i /= 0) spg_symb=spg_symb(1:i-1)
      j=index(spg_symb," ")
      n=1
      sgen(n)=spg_symb(1:j-1)
      i=index(sgen(n),"/")
      if(i /= 0) then
        sgen(n+1)=sgen(n)(i+1:)
        sgen(n)=sgen(n)(1:i-1)
        n=n+1
      end if
      n=n+1
      spg_symb=spg_symb(j+1:)
      if(len_trim(spg_symb) /= 0) then
        j=index(spg_symb," ")
        sgen(n)= spg_symb(1:j-1)
        i=index(sgen(n),"/")
        if(i /= 0) then
          sgen(n+1)=sgen(n)(i+1:)
          sgen(n)=sgen(n)(1:i-1)
          n=n+1
        end if
        spg_symb = spg_symb(j+1:)
        if(len_trim(spg_symb) /= 0) then
          n=n+1
          sgen(n)=spg_symb
          i=index(sgen(n),"/")
          if(i /= 0) then
            sgen(n+1)=sgen(n)(i+1:)
            sgen(n)=sgen(n)(1:i-1)
            n=n+1
          end if
        end if
      else
        n=n-1
      end if
      m=0
      do i=1,n
        if(sgen(i)(1:1) == "1") cycle
        m=m+1
        msgen(m)=sgen(i)
      end do
      do i=1,m
        sgen(i)=msgen(i)
      end do

      ngen = m
      !Remove the second number for screw axes
      do j=1,ngen
        do i=1,11
          if(trim(sgen(j)) == screw(i)) then
            sgen(j)= rm_screw(i)
            exit
          end if
        end do
      end do

      do i=2,SpG%Multip
         s=SpG%Op(i)%Mat(1:3,1:3)
         t=SpG%Op(i)%Mat(1:3,4)
         op_symb = Symmetry_Symbol(s,t)
         j=index(op_symb," ")
         symbg=op_symb(1:j-1)
         do j=1,ngen
           if(index(symbg,trim(sgen(j))) /= 0  .and. .not. done(j)) then
             point_op(j)=i
             done(j)=.true.
             exit
           end if
         end do
         if(all(done(1:ngen))) exit
      end do
      do j=1,ngen
        i=point_op(j)
        gen(j)=SpG%Symb_Op(i)
      end do

    End Subroutine Get_Generators_From_SpGSymbol

    Subroutine Get_Symbol_SSG_from_Operators()
      !Translations along extra coordinates allowed for the different order of operators
      !  t =  0  1/2  1/3  -1/3   1/4  -1/4   1/6   -1/6
      !Symb   0   s    t   -t      q    -q     h     -h
    End Subroutine Get_Symbol_SSG_from_Operators

  End Module getting_ssg

    Program Test_ssg_Genk
      Use CFML_GlobalDeps
      Use CFML_gSpaceGroups
      Use CFML_gSpaceGroups
      Use CFML_Propagation_Vectors, only: K_Star, Write_Group_K, Set_Gk, Group_k_Type
      Use CFML_SuperSpace_Database
      Use CFML_IOForm, only : Read_CFL_SpG,Read_CFL_Cell,Read_kinfo
      use CFML_Rational
      use getting_ssg

      implicit none

      integer                         :: nkv,i,maxorder,ngen !,j,k,m,multip
      character(len=50)               :: str !,forma
      real(kind=cp), dimension(3,12)  :: kv
      type(Group_k_Type),dimension(12):: Gk
      type(SpG_Type)                  :: SpG,intSpG
      type(SpG_Type),dimension(12)    :: Grpk
      logical                         :: ext=.true.
      integer,           dimension(:),allocatable :: order,point_op
      character(len=40), dimension(:),allocatable :: gen

      do
        write(*,"(a)",advance="no") " => Enter the number (or the symbol) of a space group: "
        read(*,"(a)") str
        if(len_trim(str) == 0) exit
        write(*,"(a)",advance="no") " => Enter the number of propagation vectors (1 integer): "
        read(*,*) nkv
        do i=1,nkv
           write(*,"(a,i2,a)",advance="no") " => Enter the propagation vector # ",i,": "
           read(*,*) kv(:,i)
        end do

        call Set_SpaceGroup(str,SpG)
        if(Err_CFML%Ierr /= 0) then
          write(unit=*,fmt="(a)") " => "//trim(Err_CFML%Msg)
          cycle
        end if
        call Write_SpaceGroup_Info(SpG)

        do i=1,nkv
          call K_Star(kv(:,i),SpG,Gk(i),ext)
          call Set_Gk(Gk(i),Grpk(i),ext)
          write(*,"(//a,3f10.5/)") " => GROUP OF THE PROPAGATION VECTOR: ",kv(:,i)
          call Write_Group_K(Gk(i))
          write(*,"(//a/)") " => FULL PROPAGATION VECTOR SPACE GROUP INCLUDING -k (Extended Little Group): "
          call Write_SpaceGroup_Info(Grpk(i))
        end do
        call Set_Intersection_Groups(Grpk(1:nkv),intSpG)
        write(*,"(//a)")  "  INTERSECTION SPACE GROUP FOR GETTING POSSIBLE SUPERSPACE GROUPS"
        call Write_SpaceGroup_Info(intSpG)
        call get_Order(intSpG, maxorder, order)
        if(allocated(point_op)) deallocate(point_op)
        allocate(point_op(intSpG%multip))
        if(allocated(gen)) deallocate(gen)
        allocate(gen(intSpG%multip))
        call Get_Generators_From_SpGSymbol(intSpG,gen,point_op,ngen)
        write(*,"(//a,i2)")  "  VISIBLE GENERATORS OF SPACE GROUP: ",ngen
        do i=1,ngen
          write(*,"(i4,a,t70,i3)") i,"  "//gen(i)//"  "//trim(intSpG%Symb_Op(point_op(i))),point_op(i)
        end do
        write(*,"(//a,i2)")  "  MAXIMUM ORDER OF INTERSECTION SPACE GROUP: ", maxorder
        do i=1,intSpG%multip
          write(*,"(i4,a,t50,i2)") i,"  "//trim(intSpG%Symb_Op(i)),order(i)
        end do
        !Determine now the possible superspace groups supposing that the propagation vectors are purely magnetic
        !That is magnetic superspace groups of type 4. Just adding the operation 1'(0 0 0 1/2) and adding in the
        !addiditonal dimensions a phase of different types depending on the order of operators
        !Take up to three generators of the Extended Little Group
        !Determine the order
      end do

    End Program Test_ssg_Genk