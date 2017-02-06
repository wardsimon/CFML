    Module ssg_datafile
      ! read data for (3+d)-dimensional superspace groups (d=1,2,3)
      implicit none

      public :: Read_SSG
      logical, public            :: database_read=.false.
      integer, parameter, public ::m_cen=16, m_ncl=322, m_ngs=16697, m_ops=48, &
                                   m_dim=6, m_cond=50, m_qv=3 !D=3+d (d=3)
      !
      integer                     :: nclasses=0         ! number of Bravais classes
      integer,   dimension(m_ncl) :: iclass_nmod=0      ! for each Bravais class: number of modulation q vectors
      integer,   dimension(m_ncl) :: iclass_number=0    ! class number
      integer,   dimension(m_ncl) :: iclass_spacegroup=0! basic space group of lattice
      integer,   dimension(m_ncl) :: iclass_nstars=0    ! number of different stars of q
      integer, dimension(3,m_ncl) :: iclass_nmodstar=0  ! number of modulation q vectors for each star

      character(len=5),  dimension(m_ncl)          :: class_nlabel=" " ! class number label: 1.1, 1.2, 1.3, etc.
      character(len=50), dimension(m_ncl)          :: class_label=" "  !   class label
      integer,           dimension(3,3,m_qv,m_ncl) :: iclass_qvec=0    !   q vectors
      integer,           dimension(m_ncl)          :: iclass_ncentering=0 ! number of centering translations
      integer,  dimension(m_dim+1,m_cen,m_ncl)     :: iclass_centering=0  ! centering translations: D=3+d integers followed by a common denominator
      integer                             :: ngroups=0           ! number of superspace groups
      integer,           dimension(m_ngs) :: igroup_number=0    ! for each superspace group,  group number
      integer,           dimension(m_ngs) :: igroup_class=0     ! Bravais class
      integer,           dimension(m_ngs) :: igroup_spacegroup=0! Basic space group
      character(len=13), dimension(m_ngs) :: group_nlabel=" " !   group number label: 1.1.1.1, 2,1,1,1, etc.
      character(len=60), dimension(m_ngs) :: group_label=" "  !   group label
      integer,           dimension(m_ngs) :: igroup_nops      !   number of operators
      !   (d+1)x(d+1) augmented matrix for each operator in supercentered setting
      !   common denominator in element (d+1,d+1)
      integer, dimension(m_dim+1,m_dim+1,m_ops,m_ngs) :: igroup_ops=0
      integer,                       dimension(m_ngs) :: igroup_nconditions=0     ! number of reflection conditions
      integer,    dimension(m_dim,m_dim,m_cond,m_ngs) :: igroup_condition1=0  ! matrix representation of righthand side
      integer,        dimension(m_dim+1,m_cond,m_ngs) :: igroup_condition2=0  ! vector representation of lefthand side
      integer,  dimension(m_ngs) :: pos_group !position in the file of the groups
      integer,  dimension(m_ncl) :: pos_class !position in the file of the Bravais classes
      integer,  dimension(m_ncl+m_ngs) :: pos_all !position in the file of all

      contains

      Subroutine Read_SSG(ok,mess)
        logical,          intent(out) :: ok
        character(len=*), intent(out) :: mess
        !
        integer :: i,j,k,m,n,imax,nmod,iclass
        integer :: i_db, ier, line_pos,L
        character(len=512) :: ssg_file
        character(len=4)   :: line
         !imax=0
         if(database_read) return
         i_db=1; ier=0
         ok=.true.
         mess=" "
         ! open data file
         ssg_file='ssg_datafile.txt'
         open(unit=i_db,file=ssg_file,status='old',action='read',position='rewind',iostat=ier)
         if(ier /= 0) then
           ok=.false.
           mess= 'Error opening the database file: '//trim(ssg_file)
           return
         end if
         L=0
         do i=1,2526
           read(i_db,"(a)") line
           if(line(1:1) == '"') then
             L=L+1
             pos_class(L)= i-1
           end if
         end do
         L=0
         do i=2527,300000
           read(i_db,"(a)",iostat=ier) line
           if (ier /= 0) exit
           if(line(1:1) == '"') then
             L=L+1
             pos_group(L)= i-1
           end if
         end do
         rewind(i_db)
         ! skip heading
         read(i_db,*)
         read(i_db,*)
         read(i_db,*)
         ! read number of Bravais classes
         read(i_db,*) nclasses
         ! read each Bravais class
         do m=1,nclasses
           read(i_db,*)n,iclass_nmod(m),iclass_number(m), iclass_spacegroup(m),iclass_nstars(m), &
                     (iclass_nmodstar(i,m),i=1,iclass_nstars(m))
           nmod=iclass_nmod(m)
           if(n /= m) then
             ok=.false.
             write(mess,"(a,i3)") 'Error in ssg_datafile @reading Bravais class #: ',m
             return
           end if

           read(i_db,*)class_nlabel(m),class_label(m)
           line_pos=line_pos+1
           read(i_db,*)(((iclass_qvec(i,j,k,m),i=1,3),j=1,3),k=1,nmod)
           read(i_db,*)iclass_ncentering(m)
           read(i_db,*)((iclass_centering(i,j,m),i=1,nmod+4),j=1,iclass_ncentering(m))
         end do

         ! read number of superspace groups
         read(i_db,*)ngroups
         ! read each superspace group

         do m=1,ngroups
           !write(6,'(i5)')m
           read(i_db,*)n,igroup_number(m),igroup_class(m),igroup_spacegroup(m)
           if(n /= m)then
             ok=.false.
             write(mess,"(a,i3)") 'Error in ssg_datafile @reading group#: ',m
             return
           end if
           iclass=igroup_class(m)
           nmod=iclass_nmod(iclass)
           read(i_db,*)group_nlabel(m),group_label(m)
           read(i_db,*)igroup_nops(m)
           read(i_db,*)(((igroup_ops(i,j,k,m),i=1,nmod+4),j=1,nmod+4), k=1,igroup_nops(m))
           read(i_db,*)igroup_nconditions(m)
           if(igroup_nconditions(m) > 0)then
               line_pos=line_pos+igroup_nconditions(m)-1
               read(i_db,*)(((igroup_condition1(i,j,k,m),i=1,nmod+3),j=1,nmod+3),(igroup_condition2(j,k,m),j=1,nmod+4), &
                          k=1,igroup_nconditions(m))
           end if
         end do
         database_read=.true.
      End Subroutine Read_SSG

      Subroutine Read_single_SSG(num,ok,Mess)
        integer,          intent(in)  :: num
        logical,          intent(out) :: ok
        character(len=*), intent(out) :: mess
        !
        integer :: i,j,k,n,m,i_pos,n_skip,nmod,i_db,ier,iclass
        character(len=512) ssg_file,pos_file
         ok=.true.
         mess=" "
         ! open data file
         ssg_file='ssg_datafile.txt'
         open(newunit=i_db,file=ssg_file,status='old',action='read',position='rewind',iostat=ier)
         if(ier /= 0) then
           ok=.false.
           mess= 'Error opening the database file: '//trim(ssg_file)
           return
         end if
         pos_file='class+group_pos.txt'
         open(newunit=i_pos,file=pos_file,status='old',action='read',position='rewind',iostat=ier)
         if(ier /= 0) then
           ok=.false.
           mess= 'Error opening the database file: '//trim(pos_file)
           return
         end if
         read(unit=i_pos,fmt=*) !skip class line
         read(unit=i_pos,fmt=*) pos_class
         read(unit=i_pos,fmt=*) !skip group line
         read(unit=i_pos,fmt=*) pos_group
         close(unit=i_pos)
         read(i_db,*)
         read(i_db,*)
         read(i_db,*)
         ! read number of Bravais classes
         read(i_db,*) nclasses
         ! read each Bravais class
         do m=1,nclasses
           read(i_db,*)n,iclass_nmod(m),iclass_number(m), iclass_spacegroup(m),iclass_nstars(m), &
                     (iclass_nmodstar(i,m),i=1,iclass_nstars(m))
           nmod=iclass_nmod(m)
           read(i_db,*)class_nlabel(m),class_label(m)
           read(i_db,*)(((iclass_qvec(i,j,k,m),i=1,3),j=1,3),k=1,nmod)
           read(i_db,*)iclass_ncentering(m)
           read(i_db,*)((iclass_centering(i,j,m),i=1,nmod+4),j=1,iclass_ncentering(m))
         end do
         rewind(i_db)
         !write(*,"(10i8)") pos_group
         n_skip=pos_group(num)-1
         !write(*,"(a,i12)") "Skipping ",n_skip
         do i=1,n_skip
           read(unit=i_db,fmt=*)
         end do
         m=num
         read(i_db,*)n,igroup_number(m),igroup_class(m),igroup_spacegroup(m)
         if(n /= m)then
           ok=.false.
           write(mess,"(a,2i5)") 'Error in ssg_datafile @reading group#: ',m,n
           return
         end if
         iclass=igroup_class(m)
         nmod=iclass_nmod(iclass)
         read(i_db,*)group_nlabel(m),group_label(m)
         read(i_db,*)igroup_nops(m)
         read(i_db,*)(((igroup_ops(i,j,k,m),i=1,nmod+4),j=1,nmod+4), k=1,igroup_nops(m))
      End Subroutine Read_single_SSG

    End Module ssg_datafile


    Module CFML_SuperSpaceGroups
      use ssg_datafile
      use CFML_Rational_Arithmetic

      Implicit None

      public :: operator (*)

      interface operator (*)
        module procedure multiply_ssg_symop
      end interface

      !!---- Type, public   :: SSym_Oper_Type
      !!----
      !!---- Rational matrix of special type dimension (3+d+1,3+d+1)
      !!---- The matrix of the superspace operator is extended with
      !!---- a column containing the translation in supespace plus
      !!---- a zero 3+d+1 row and + or -1 for time reversal at position
      !!---- (3+d+1,3+d+1). In order to limit the operators to the
      !!---- factor group the translations, a modulo-1 is applied
      !!---- in the multiplication of two operators.
      Type, public   :: SSym_Oper_Type
        Type(rational), allocatable, dimension(:,:) :: Mat
      End Type SSym_Oper_Type

      Type, public        :: SuperSpaceGroup_Type
        logical           :: standard_setting=.true.  !true or false
        Character(len=60) :: SSG_symbol=" "
        Character(len=60) :: SSG_Bravais=" "
        Character(len=13) :: SSG_nlabel=" "
        Character(len=20) :: Parent_spg=" "
        Character(len=80) :: trn_from_parent=" "
        Character(len=80) :: trn_to_standard=" "
        character(len=80) :: Centre=" "    ! Alphanumeric information about the center of symmetry
        integer           :: d=1           !(d=1,2,3, ...) number of q-vectors
        integer           :: Parent_num=0  ! Number of the parent Group
        integer           :: Bravais_num=0 ! Number of the Bravais class
        integer           :: Num_Lat=0     ! Number of lattice points in a cell
        integer           :: Num_aLat=0    ! Number of anti-lattice points in a cell
        integer           :: Centred=1     ! Centric or Acentric [ =0 Centric(-1 no at origin),=1 Acentric,=2 Centric(-1 at origin)]
        integer           :: NumOps=0      ! Number of reduced set of S.O. (removing lattice centring and anticentrings and centre of symmetry)
        integer           :: Multip=0      ! General multiplicity
        real,   allocatable,dimension(:,:)  :: Latt_trans    ! Lattice translations (3+d,Num_lat)
        real,   allocatable,dimension(:,:)  :: aLatt_trans   ! Lattice anti-translations (3+d,Num_alat)
        real,   allocatable,dimension(:)    :: Centre_coord  ! Fractional coordinates of the inversion centre (3+d)
        real,   allocatable,dimension(:,:)  :: kv            !k-vectors (3,d)
        integer,allocatable,dimension(:)    :: time_rev      !Multip
        integer,allocatable,dimension(:,:,:):: Om            !Operator matrices (3+d+1,3+d+1,Multip) common denominator at (4+d,4+d)
        character(len=80),allocatable,    dimension(:)  :: SymopSymb  ! Alphanumeric Symbols for SYMM
        type(SSym_Oper_Type),allocatable, dimension(:)  :: SymOp      ! Crystallographic symmetry operators

      End Type SuperSpaceGroup_Type


    Contains

      Subroutine Allocate_SSG_SymmOps(d,multip,SymOp)
        integer,                                        intent(in)      :: d,multip
        type(SSym_Oper_Type),allocatable, dimension(:), intent(in out)  :: SymOp
        integer :: i,Di
        if(allocated(SymOp)) deallocate(SymOp)
        allocate(SymOp(multip))
        Di=3+d+1
        do i=1,multip
          allocate(SymOp(i)%Mat(Di,Di))
          SymOp(i)%Mat=0//1
        end do
      End Subroutine Allocate_SSG_SymmOps

      Subroutine Set_SSG_Reading_Database(num,ssg,ok,Mess)
        integer,                    intent(in)  :: num
        type(SuperSpaceGroup_Type), intent(out) :: ssg
        Logical,                    intent(out) :: ok
        character(len=*),           intent(out) :: Mess
        !
        integer :: i,j,d,nmod,D
        if(.not. database_read)  then
         call Read_single_SSG(num,ok,Mess)
         read(i_db,*)group_nlabel(m),group_label(m)
         read(i_db,*)igroup_nops(m)
         read(i_db,*)(((igroup_ops(i,j,k,m),i=1,nmod+4),j=1,nmod+4), k=1,igroup_nops(m))
         read(i_db,*)iclass_ncentering(m)
         read(i_db,*)((iclass_centering(i,j,m),i=1,nmod+4),j=1,iclass_ncentering(m))
        end if
        nmod=iclass_nmod(iclass)
        D=3+nmod+1
        ssg%standard_setting=.true.
        iclass=igroup_class(num)

        ssg%SSG_symbol=group_label(num)
        ssg%SSG_Bravais=class_label(iclass)
        ssg%SSG_nlabel=group_nlabel(num)
        ssg%Parent_spg=" "
        ssg%trn_from_parent="a,b,c;0,0,0"
        ssg%trn_to_standard="a,b,c;0,0,0"
        ssg%Centre=" "             ! Alphanumeric information about the center of symmetry
        ssg%d=iclass_nmod(iclass)  !(d=1,2,3, ...) number of q-vectors
        ssg%Parent_num=0           ! Number of the parent Group
        ssg%Bravais_num=0          ! Number of the Bravais class
        ssg%Num_Lat=iclass_ncentering(num)     ! Number of lattice points in a cell
        ssg%Num_aLat=0                ! Number of anti-lattice points in a cell
        ssg%Centred=1                 ! Centric or Acentric [ =0 Centric(-1 no at origin),=1 Acentric,=2 Centric(-1 at origin)]
        ssg%NumOps=igroup_nops(num)        ! Number of reduced set of S.O. (removing lattice centring and anticentrings and centre of symmetry)
        ssg%Multip==igroup_nops(num)* iclass_ncentering(num)     ! General multiplicity
        Allocate(ssg%Latt_trans(D-1,ssg%Num_Lat)
        nmod=iclass_nmod(iclass)
        do j=1,iclass_ncentering(num)
           ssg%Latt_trans(1:nmod+3,j)= real(iclass_centering(1:nmod+3,j,num))/real(iclass_centering(nmod+4,j,num))
        end do
        call Allocate_SSG_SymmOps(d,multip,SymOp)
        do i=1,ssg%NumOps
           SymOp(i)%Mat= igroup_ops(1:nmod+4,1:nmod+4,i,num)//igroup_ops(nmod+4,nmod+4,i,num)
        end do

      End Subroutine Set_SSG_Reading_Database

      function multiply_ssg_symop(Op1,Op2) result (Op3)
        type(SSym_Oper_Type), intent(in) :: Op1,Op2
        type(SSym_Oper_Type)             :: Op3
        integer :: n,d
        n=size(Op1%Mat,dim=1)
        d=n-1
        Op3%Mat=matmul(Op1%Mat,Op2%Mat)
        Op3%Mat(1:d,n)=mod(Op3%Mat(1:d,n),1)
      end function multiply_ssg_symop

    End Module CFML_SuperSpaceGroups

    Program read_ssg_datafile
      use ssg_datafile
      implicit none

      integer :: iclass,nmod,i,j,k,m
      character(len=20)  :: forma
      character(len=130) :: message
      logical :: ok

      !call Read_SSG(ok,message)
      !if(.not. ok) then
      !  write(*,"(a)") "   !!! "//message//" !!!"
      !  stop
      !end if
      !!
      !open(unit=1,file="class+group_pos.txt",status="replace",action="write")
      !write(1,"(a,i6)") "Bravais Classes positions in database, Number of classes ",m_ncl
      !write(1,"(10i7)") pos_class
      !write(1,"(a,i6)") "Space Group positions in database, Number of groups ",m_ngs
      !write(1,"(10i7)") pos_group
      !close(unit=1)
      !stop
      do
        write(*,"(a)",advance="no") " => Enter the number of the SSG: "
        read(*,*) m
        if(m == 0) exit
        Call Read_single_SSG(m,ok,Message)
        !Call Read_SSG(ok,Message)
        !if(.not. ok) then
        !  write(*,"(a)") "   !!! "//trim(message)//" !!!"
        !  stop
        !end if
        write(*,"(3(a,i5))") " Group number:",igroup_number(m), " Bravais class:",igroup_class(m), "  Basic Group #:",igroup_spacegroup(m)
        iclass=igroup_class(m)
        nmod=iclass_nmod(iclass)
        write(*,"(a,tr4,a)")  group_nlabel(m), group_label(m)
        write(*,"(a,i3)") " Number of operators:", igroup_nops(m)
        do k=1,igroup_nops(m)
          write(*,"(a,i3)") " Operator #",k
          forma="(  i4)"
          write(forma(2:3),"(i2)") nmod+4
          do i=1,nmod+4
            write(*,forma) igroup_ops(i,1:nmod+4,k,m)
            !write(*,*)(((igroup_ops(i,j,k,m),i=1,nmod+4),j=1,nmod+4), k=1,igroup_nops(m))
          end do
        end do
        write(*,"(a,i3)") " Reflection conditions: ",igroup_nconditions(m)
        if(igroup_nconditions(m) > 0)then
         do k=1,igroup_nconditions(m)
           write(*,*)((igroup_condition1(i,j,k,m),i=1,nmod+3),j=1,nmod+3),(igroup_condition2(j,k,m),j=1,nmod+4)
         end do
        end if
      end do

    End Program read_ssg_datafile
