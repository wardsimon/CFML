!!----
!!---- Program: Super Space Group Info (SSPG_INFO)
!!----          Made after the Space Group Info program in  CFML.
!!----          Restricted to Magnetic SSGs in particular.
!!----
!!---- Author: R. Martinez
!!---- Revision: 
!!

!-----------------------------------------------------------------

    Module ssg_datafile
      ! read data for (3+d)-dimensional superspace groups (d=1,2,3)
      implicit none

      public :: Read_SSG
      logical, public            :: database_read=.false.
      integer, parameter, public :: m_cen=16, m_ncl=322, m_ngs=16697, m_ops=48, &
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
        integer :: i_db, ier,L
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
      private
      public :: Allocate_SG_SymmOps, Allocate_SSG_SymmOps, Set_SSG_Reading_Database, Write_SSG, Gen_Group

      public :: operator (*)
      interface operator (*)
        module procedure multiply_ssg_symop
      end interface
      private :: multiply_ssg_symop

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
        logical                                          :: standard_setting=.true.  !true or false
        Character(len=60)                                :: SSG_symbol=" "
        Character(len=60)                                :: SSG_Bravais=" "
        Character(len=13)                                :: SSG_nlabel=" "
        Character(len=20)                                :: Parent_spg=" "
        Character(len=80)                                :: trn_from_parent=" "
        Character(len=80)                                :: trn_to_standard=" "
        character(len=80)                                :: Centre="Acentric" ! Alphanumeric information about the center of symmetry
        Character(len=1)                                 :: SPG_Lat="P"   ! Symbol of the lattice
        integer                                          :: d=1           !(d=1,2,3, ...) number of q-vectors
        integer                                          :: Parent_num=0  ! Number of the parent Group
        integer                                          :: Bravais_num=0 ! Number of the Bravais class
        integer                                          :: Num_Lat=0     ! Number of lattice points in a cell
        integer                                          :: Num_aLat=0    ! Number of anti-lattice points in a cell
        integer                                          :: Centred=1     ! Centric or Acentric [ =0 Centric(-1 no at origin),=1 Acentric,=2 Centric(-1 at origin)]
        integer                                          :: NumOps=0      ! Number of reduced set of S.O. (removing lattice centring and anticentrings and centre of symmetry)
        integer                                          :: Multip=0      ! General multiplicity
        real,                allocatable,dimension(:,:)  :: kv            ! k-vectors (3,d)
        integer,             allocatable,dimension(:)    :: time_rev      ! Time Reversal
        type(rational),      allocatable,dimension(:,:)  :: Latt_trans    ! Lattice translations (3+d,Num_lat)
        type(rational),      allocatable,dimension(:,:)  :: aLatt_trans   ! Lattice anti-translations (3+d,Num_alat)
        type(rational),      allocatable,dimension(:)    :: Centre_coord  ! Fractional coordinates of the inversion centre (3+d)
        type(rational),      allocatable,dimension(:,:,:):: Om            ! Operator matrices (3+d+1,3+d+1,Multip) common denominator at (4+d,4+d)
        type(SSym_Oper_Type),allocatable, dimension(:)   :: SymOp         ! Crystallographic symmetry operators
        character(len=80),   allocatable,dimension(:)    :: SymopSymb     ! Alphanumeric Symbols for SYMM
      End Type SuperSpaceGroup_Type

      logical, public :: Err_ssg
      character(len=80), public :: Err_ssg_mess

    Contains

      ! need this to convert SG operators into rationals
      ! to decouple both problems: aperiodicity and superspace
      Subroutine Allocate_SG_SymmOps(d,multip,SymOp) ! but d=0 always in this case
        integer,                                        intent(in)      :: d,multip
        type(SSym_Oper_Type),allocatable, dimension(:), intent(in out)  :: SymOp
        integer :: i,Dp

        if(allocated(SymOp)) deallocate(SymOp)
        allocate(SymOp(multip))
        Dp=3+d+1 ! 
        do i=1,multip
          allocate(SymOp(i)%Mat(Dp,Dp))
          SymOp(i)%Mat=0//1
        end do
      End Subroutine Allocate_SG_SymmOps


      Subroutine Allocate_SSG_SymmOps(d,multip,SymOp)
        integer,                                        intent(in)      :: d,multip
        type(SSym_Oper_Type),allocatable, dimension(:), intent(in out)  :: SymOp
        integer :: i,Dp

        if(allocated(SymOp)) deallocate(SymOp)
        allocate(SymOp(multip))
        Dp=3+d+1
        do i=1,multip
          allocate(SymOp(i)%Mat(Dp,Dp))
          SymOp(i)%Mat=0//1
        end do
      End Subroutine Allocate_SSG_SymmOps

      Subroutine Set_SSG_Reading_Database(num,ssg,ok,Mess)
        integer,                    intent(in)  :: num
        type(SuperSpaceGroup_Type), intent(out) :: ssg
        Logical,                    intent(out) :: ok
        character(len=*),           intent(out) :: Mess
        !
        integer :: i,j,nmod,Dp,D,iclass
        type(rational), dimension(:,:), allocatable :: Inv
        character(len=15) :: forma
        logical :: inv_found

        if(.not. database_read)  then
           call Read_single_SSG(num,ok,Mess)
           if(.not. ok) return
        end if
        iclass=igroup_class(num)
        nmod=iclass_nmod(iclass)
        D=3+nmod
        Dp=D+1
        ssg%standard_setting=.true.
        ssg%SSG_symbol=group_label(num)
        ssg%SSG_Bravais=class_label(iclass)
        ssg%SSG_nlabel=group_nlabel(num)
        i=index(ssg%SSG_symbol,"(")
        ssg%Parent_spg=ssg%SSG_symbol(1:i-1)
        ssg%trn_from_parent="a,b,c;0,0,0"
        ssg%trn_to_standard="a,b,c;0,0,0"
        ssg%Centre="Acentric"                  ! Alphanumeric information about the center of symmetry
        ssg%d=nmod                             !(d=1,2,3, ...) number of q-vectors
        ssg%Parent_num=igroup_spacegroup(num)  ! Number of the parent Group
        ssg%Bravais_num=igroup_class(num)      ! Number of the Bravais class
        ssg%Num_Lat=iclass_ncentering(iclass)  ! Number of lattice points in a cell
        if( ssg%Num_Lat > 1 .and. ssg%SSG_symbol(1:1) == "P" ) then
          ssg%SPG_Lat="Z"
        else
          ssg%SPG_Lat=ssg%SSG_symbol(1:1)
        end if


        ssg%Num_aLat=0                         ! Number of anti-lattice points in a cell
        ssg%Centred= 1                         ! Centric or Acentric [ =0 Centric(-1 no at origin),=1 Acentric,=2 Centric(-1 at origin)]
        ssg%NumOps=igroup_nops(num)            ! Number of reduced set of S.O. (removing lattice centring and anticentrings and centre of symmetry)
        ssg%Multip=igroup_nops(num)* max(iclass_ncentering(iclass),1)     ! General multiplicity

        Allocate(ssg%Latt_trans(D,ssg%Num_Lat))
        Allocate(ssg%aLatt_trans(D,ssg%Num_aLat))
        Allocate(ssg%kv(3,nmod))
        Allocate(ssg%time_rev(ssg%Multip))
        Allocate(ssg%Centre_coord(D))
        call Allocate_SSG_SymmOps(ssg%d,ssg%Multip,ssg%SymOp)

        ssg%kv=0.0; ssg%Latt_trans=0//1; ssg%aLatt_trans=0//1; ssg%time_rev=1; ssg%Centre_coord=0//1
        nmod=iclass_nmod(iclass)
        !forma="(a,i3, i4)"
        !write(forma(7:7),"(i1)") Dp
        do j=1,iclass_ncentering(iclass)
           !write(*,forma) " Vector #",j, iclass_centering(1:Dp,j,iclass)
           ssg%Latt_trans(1:D,j)= rational_simplify(iclass_centering(1:D,j,iclass)//iclass_centering(Dp,j,iclass))
        end do
        do i=1,ssg%NumOps
           ssg%SymOp(i)%Mat= rational_simplify(igroup_ops(1:Dp,1:Dp,i,num)//igroup_ops(Dp,Dp,i,num))
        end do

        !Look for a centre of symmetry
        inv_found=.false.
        allocate(Inv(D,D))
        Inv=0//1
        do i=1,D
          Inv(i,i) = -1//1
        end do
        do i=1,ssg%NumOps
          if(equal_rational_matrix(Inv,ssg%SymOp(i)%Mat(1:D,1:D))) then
            if(sum(ssg%SymOp(i)%Mat(1:D,Dp)) == 0//1) then
              ssg%Centred=2
              ssg%Centre_coord=0//1
              ssg%Centre="Centrosymmetric with centre at origin"
            else
              ssg%Centred=0
              ssg%Centre_coord=ssg%SymOp(i)%Mat(1:D,Dp)/2
              write(unit=ssg%Centre,fmt="(a)") "Centrosymmetric with centre at : "//print_rational(ssg%Centre_coord)
            end if
            exit
          end if
        end do

      End Subroutine Set_SSG_Reading_Database

      function multiply_ssg_symop(Op1,Op2) result (Op3)
        type(SSym_Oper_Type), intent(in) :: Op1,Op2
        type(SSym_Oper_Type)             :: Op3
        integer :: n,d,i
        n=size(Op1%Mat,dim=1)
        d=n-1
        Op3%Mat=matmul(Op1%Mat,Op2%Mat)
        Op3%Mat(1:d,n)=mod(Op3%Mat(1:d,n),1)
        do i=1,d
          if(Op3%Mat(i,n) < 0//1) Op3%Mat(i,n) = Op3%Mat(i,n) + 1
        end do
      end function multiply_ssg_symop

      Subroutine Write_SSG(SpaceGroup,iunit,full)
        type(SuperSpaceGroup_Type), intent(in) :: SpaceGroup
        integer, optional,          intent(in) :: iunit
        logical, optional,          intent(in) :: full
        !---- Local variables ----!
        integer :: lun,i,j,k,D,Dp,nlines
        character(len=40),dimension(:),  allocatable :: vector
        character(len=40),dimension(:,:),allocatable :: matrix
        character(len=15) :: forma
        integer,  parameter                      :: max_lines=192
        character (len=100), dimension(max_lines):: texto
        logical                                  :: print_latt

        lun=6
        print_latt=.false.
        if (present(iunit)) lun=iunit
        if (present(full))  print_latt=.true.

        D=3+SpaceGroup%d
        Dp=D+1
        allocate(vector(D),matrix(Dp,Dp))
        write(unit=lun,fmt="(/,a,/)")  " "
        write(unit=lun,fmt="(/,/,a)")          "        Information on SuperSpace Group: "
        write(unit=lun,fmt="(a,/ )")           "        ------------------------------- "
        write(unit=lun,fmt="(a,a)")           " =>   Number of Space group: ", trim(SpaceGroup%SSG_nlabel)
        write(unit=lun,fmt="(a,a)")           " => SuperSpace Group Symbol: ", trim(SpaceGroup%SSG_symbol)
        write(unit=lun,fmt="(a,a)")           " =>    Bravais Class Symbol: ", trim(SpaceGroup%SSG_Bravais)
        write(unit=lun,fmt="(a,a)")           " =>      Parent Space Group: ", trim(SpaceGroup%Parent_spg)
        if(.not. SpaceGroup%standard_setting) then
          write(unit=lun,fmt="(a,a)")         " =>     Transf. from Parent: ", trim(SpaceGroup%trn_from_parent)
          write(unit=lun,fmt="(a,a)")         " =>     Transf. to Standard: ", trim(SpaceGroup%trn_to_standard)
        end if
        write(unit=lun,fmt="(a,a)")           " =>          Centrosymmetry: ", trim(SpaceGroup%Centre)
        write(unit=lun,fmt="(a,a)")           " =>         Bravais Lattice: ", trim(SpaceGroup%SPG_Lat)

        write(unit=lun,fmt="(a,i3)")          " => Number of  Parent Group: ", SpaceGroup%Parent_num
        write(unit=lun,fmt="(a,i3)")          " => Number of Bravais Class: ", SpaceGroup%Bravais_num
        write(unit=lun,fmt="(a,i3)")          " => Number of q-vectors (d): ", SpaceGroup%d
        write(unit=lun,fmt="(a,i3)")          " =>   # of Centring Vectors: ", max(SpaceGroup%Num_Lat-1,0)
        write(unit=lun,fmt="(a,i3)")          " => # Anti-Centring Vectors: ", SpaceGroup%Num_aLat

        !write(unit=lun,fmt="(a,a)")           " =>          Crystal System: ", trim(SpaceGroup%CrystalSys)
        !write(unit=lun,fmt="(a,a)")           " =>              Laue Class: ", trim(SpaceGroup%Laue)
        !write(unit=lun,fmt="(a,a)")           " =>             Point Group: ", trim(SpaceGroup%Pg)


        write(unit=lun,fmt="(a,i3)")          " =>  Reduced Number of S.O.: ", SpaceGroup%NumOps
        write(unit=lun,fmt="(a,i3)")          " =>    General multiplicity: ", SpaceGroup%Multip
        !write(unit=lun,fmt="(a,i4)")          " =>  Generators (exc. -1&L): ", SpaceGroup%num_gen
        if (SpaceGroup%centred == 0) then
           write(unit=lun,fmt="(a,a)")        " =>               Centre at: ", print_rational(SpaceGroup%Centre_coord)
        end if
        if (print_latt .and. max(SpaceGroup%Num_lat-1,0) > 0) then
           texto(:) (1:100) = " "
           write(unit=lun,fmt="(a,i3)")       " =>        Centring vectors: ",max(SpaceGroup%Num_lat-1,0)
           nlines=1
           forma="(a,i2,a,   a)"
           write(unit=forma(9:11),fmt="(i3)") D
           do i=2,SpaceGroup%Num_lat
              vector=print_rational(SpaceGroup%Latt_trans(:,i))
              if (mod(i-1,2) == 0) then
                 write(unit=texto(nlines)(51:100),fmt=forma) &
                                            " => Latt(",i-1,"): ",(trim(vector(j))//" ",j=1,D)
                 nlines=nlines+1
              else
                 write(unit=texto(nlines)( 1:50),fmt=forma)  &
                                            " => Latt(",i-1,"): ",(trim(vector(j))//" ",j=1,D)
              end if
           end do
           do i=1,nlines
              write(unit=lun,fmt="(a)") texto(i)
           end do
        end if
        !Writing of the rational operator matrices
        forma="( a8)"
        write(forma(2:2),"(i1)") Dp
        do i=1,SpaceGroup%NumOps
          matrix=print_rational(SpaceGroup%SymOp(i)%Mat)
          write(unit=lun,fmt="(a,i3)") "  Rational Operator #",i
          do j=1,Dp
             write(unit=lun,fmt=forma) (trim( Matrix(j,k))//" ",k=1,Dp)
          end do
        end do

      End Subroutine Write_SSG

      !This subroutine assumes that Op contains the identity as the first operator, followed
      !by few generators. The value of ngen icludes also the identity
      Subroutine Gen_Group(ngen,Op,multip)
        integer,                            intent(in)     :: ngen
        type(SSym_Oper_Type), dimension(:), intent(in out) :: Op
        integer,                            intent(out)    :: multip
        !--- Local variables ---!
        integer :: i,j,k,n, nt,max_op
        type(SSym_Oper_Type) :: Opt
        logical, dimension(size(Op),size(Op)) :: done

        max_op=size(Op)
        done=.false.
        done(1,:) = .true.
        done(:,1) = .true.
        nt=ngen

        do_ext:do
          n=nt
          do i=1,n
            do_j:do j=1,n
              if(done(i,j)) cycle
              Opt=Op(i)*Op(j)
              do k=1,nt
                if(equal_rational_matrix(Opt%Mat,Op(k)%Mat)) cycle do_j
              end do
              done(i,j)=.true.
              nt=nt+1
              if(nt > max_op) then
                exit do_ext
              end if
              Op(nt)=Opt
            end do do_j
          end do
          if ( n == nt) exit do_ext
        end do do_ext
        multip=nt
      End Subroutine Gen_Group

    End Module CFML_SuperSpaceGroups



!-----------------------------------------------------------------
!---> Main
!-----------------------------------------------------------------



Program SSPG_Info

   !---- Use Modules ----!
   
   use CFML_Crystallographic_Symmetry, only: Space_Group_Type, set_SpaceGroup, &
                                             Write_SpaceGroup, Symmetry_Symbol, &
                                             Get_Generators_From_SpGSymbol,     &
                                             err_symm, err_symm_mess
   use CFML_GlobalDeps,                only: cp
   use CFML_Propagation_Vectors,       only: Group_k_Type, Set_Gk, K_Star, &
                                             Write_Group_k, k_EQUIV
   use CFML_String_Utilities,          only: Get_Fraction_1Dig
   use CFML_Math_General,              only: Determinant, Equal_Matrix, & 
                                             Equal_Vector, Modulo_Lat, Zbelong
   use CFML_SuperSpaceGroups
   use ssg_datafile
   use CFML_Rational_Arithmetic

   !---- Variables ----!
   implicit none

   !>
   !> final result must be written as a SuperSpaceGroup_Type
   !> use notation similar to SG write()
   !> type(SuperSpaceGroup_Type),dimension(:),allocatable :: SSG_out 
   !> use write_ssg()
   !>


   !---variables for  STEP 1: 
   !                  Prompts the user to type a SG and a  k-vector
   !                  Computes the SG's Small Group (or Group of k)
   !
   character(len=20)             :: spgr
   type(Space_Group_type)        :: grp_espacial, SPGk_in, SPGk_out
   type(Group_k_Type)            :: Gk_out
   character(len=1)              :: default_example ! len=1 ('y' or 'n')  
   real(kind=cp), dimension(3)   :: vec_k
   integer                       :: size_Gk
   !
   !                  For the retrieval of the SG generators
   integer,           dimension(10) :: point_op
   character(len=50), dimension(10) :: gen
   character(len=70)                :: op_symb
   integer :: ngen



   !---variables for  STEP 2: 
   !                  Compute the SG's Small Group Extended for the k-vector.
   !                  Cast the result to a vector of type(SSym_Oper_Type)
   !
   integer :: size_Gkmink
   integer :: i,j,k




   !---variables for  STEP 3: 
   !                  Compute the reciprocal vectors {H(R)}
   !                  that fulfill conditions 1 and  2, namely the small group 
   !                  (G_{k}) and extended small group (G_{k,-k}), respectively
   !
   real(kind=cp), dimension(3)                :: vec_k_input
   real(kind=cp), dimension(3)                :: vec_k_prime
   real(kind=cp), dimension(3)                :: vec_k_condition1
   real(kind=cp), dimension(3)                :: vec_k_condition2
   real(kind=cp), dimension(:,:), allocatable :: vecs_H_R
   Type (Space_Group_Type)                    :: SPGk_input



   !---variables for  STEP 4:
   !                  Compute the possible values t4 and store in vector
   !
   integer                             :: m,n, max_order, size_vec_all_orders
   real                                :: t4
   real, dimension(:), allocatable     :: vec_T4_possible_values
   real, dimension(4,2)                :: vec_T4_possible_combinations 
   integer, dimension(:), allocatable  :: vec_all_orders



   !---variables for  STEP 5:
   !                  Compound the (3+d+1) matrices as a type(SSym_Oper_Type)
   !                  leaving just t4(=44) indetermined. 
   !                  Note: include the extra operator 'closeness of SSG'
   !                  Note: this step has been skipped and merged with 6
   !
   type(SSym_Oper_Type), dimension(:), allocatable  :: setOp 
   type(SSym_Oper_Type), dimension(:), allocatable  :: setOpGen 
   type(SSym_Oper_Type), dimension(:), allocatable  :: setOpGenTest 
   character(len=80),    dimension(5,5)             :: cMat1, cMat2



   !---variables for  STEP 6:
   !                  this step has been skipped and merged with 7 
   ! 
   integer i_size, f_size, ss_size, sg_size, gg_size
   integer, dimension(:), allocatable    :: gg_ID    ! indexes of elements of a group
   integer, dimension(:), allocatable    :: sset_ID  ! indexes of a subset of a group


   !---variables for  STEP 7:
   !                  work with just the generators of the group plus 1'
   !
   integer ii, nGenerators, nmax, ncmb, ngen_ssg
   type (rational),   dimension(5,5)               :: MAT_IDENTITY 
   character(len=40), dimension(:,:), allocatable  :: matrix   
   character(len=15)                               :: forma
   Type (SuperSpaceGroup_Type)                     :: SSG_out



   !---- Procedure ----!
   do
      WRITE(*,*) " "
      WRITE(*,*) " "
      WRITE(*,*) " ========================================================================"
      WRITE(*,*) " "
      WRITE(*,*) " STEP 1: find the SMALL GROUP"
      WRITE(*,*) " "

      ! Overriding with a default input example
      !
      !write(unit=*,fmt="(a)") " => Enter a space group: "
      !write(unit=*,fmt="(a)",advance="no") " => Space Group (HM/Hall symbol or number): "
      !read(unit=*,fmt="(a)") spgr
      !overriding: default example
      !spgr="P n m a"
      !WRITE(*,*) "your input:", spgr

      !if (len_trim(spgr)==0) exit
      !write(unit=*,fmt="(a)",advance="no") " => Enter a a k-vector: "
      !read(unit=*,fmt=*) vec_k
      !overriding: default example
      !vec_k=(/ 0.321, 0.0, 0.5 /)

      spgr=""
      vec_k=0
      write(unit=*,fmt="(a)") " => Would you like the default Space Group and h-vector"
      write(unit=*,fmt="(a)") " => P n m a"
      write(unit=*,fmt="(a)") " => 0.321 0 0.5"
      write(unit=*,fmt="(a)") " => y/n?"
      read(unit=*,fmt=*) default_example
      if(default_example=="y") then
          write(unit=*,fmt="(a)") "Using default example"
          spgr="P n m a"
          vec_k=(/ 0.321, 0.0, 0.5 /)
      else
          !> New prompt
          write(unit=*,fmt="(a)") " => Enter a space group: "
          write(unit=*,fmt="(a)",advance="no") " => Space Group (HM/Hall symbol or number): "
          read(unit=*,fmt="(a)") spgr
          if (len_trim(spgr)==0) exit
          write(unit=*,fmt="(a)",advance="no") " => Enter a a k-vector: "
          read(unit=*,fmt=*) vec_k
      end if

      !> Setting and writing the SG Information
      !> Casting the SG variable to decouple this part of the code from the following
      WRITE(*,*) " setting and writing grp_espacial:"
      call Set_SpaceGroup(spgr, grp_espacial)
      call Write_SpaceGroup(grp_espacial, full=.true.)
      size_Gk=grp_espacial%multip  
      SPGk_in=grp_espacial ! casting


      !> Get the generators of the SG from its SG symbol
      call Get_Generators_From_SpGSymbol(grp_espacial,gen,point_op,ngen)
      write(unit=*,fmt="(/,a,i3,a)") " => Generators of the Space Group: # ",grp_espacial%NumSpg,"  "//trim(grp_espacial%SPG_Symb)
      do i=1,ngen
        j=point_op(i)
        call Symmetry_Symbol(grp_espacial%SymopSymb(j),op_symb)
        write(unit=*,fmt="(a,i3,a)") "    Generator #",i, "  "//gen(i)(1:30)//"Symbol: "//trim(op_symb), "    SYMM(",j,")"        
      end do
      write(unit=*,fmt="(a)") "   "

      

      
      WRITE(*,*) " ========================================================================"
      WRITE(*,*) " "
      WRITE(*,*) " STEP 2: find the SMALL GROUP EXTENDED "
      WRITE(*,*) " "

      !> Computing the K-Star (i.e. the small group or G_{k})
      !> NOTE: ext=.true. in K_Star() provides the EG (extended group, G_{k,-k})
      !> The output will be loaded into SPGk_out:
      WRITE(*,*) " "
      call K_Star(vec_k, SPGk_in, Gk_out, ext=.true.)  ! ext=.true. gives G_{k,-k} 
      call Write_Group_k(Gk_out)
      call Set_Gk(Gk_out, SPGk_out, ext=.true.)
      call Write_SpaceGroup(SPGk_out, full=.true.) ! SPGk_out is now the EG
      size_Gkmink=SPGk_out%multip

 


      WRITE(*,*) " "
      WRITE(*,*) " ========================================================================"
      WRITE(*,*) " "
      WRITE(*,*) " STEP 3: find the H(R) vectors and R_S matrices "
      WRITE(*,*) " "


      !> allocate H(R) vector with the number of operators (SPGk_out%multip)
      WRITE(*,*) "number of  R_S matrices = ", SPGk_out%multip
      if(allocated(vecs_H_R)) deallocate(vecs_H_R) 
      allocate(vecs_H_R(4,SPGk_out%multip)) 
      vecs_H_R=-10

      do j=1, SPGk_out%multip+0
          vec_k_prime=matmul(vec_k, SPGk_out%SymOp(j)%Rot) 
          vec_k_condition1=vec_k_prime-vec_k ! k1=k'-k  <=>  k' = kR = +k + H(R)   
          vec_k_condition2=vec_k_prime+vec_k ! k2=k'+k  <=>  k' = kR = -k + H(R)

          ! if k1=H(R) belongs to reciprocal lattice, then ...
          if(        k_EQUIV(vec_k_prime, vec_k,  SPGk_out%SPG_lat)    ) then  
               vecs_H_R(1:3,j)=vec_k_condition1
               vecs_H_R(4:4,j)=1 
          else if(   k_EQUIV(vec_k_prime, -vec_k, SPGk_out%SPG_lat)    ) then
               vecs_H_R(1:3, j)=vec_k_condition2
               vecs_H_R(4:4, j)=-1     
          else
               WRITE(*,*) "Error computing the H(R)!"
          end if
      end do

      WRITE(*,*) ""

      !> printing the H(R_j) vectors plus the 4th component RI
      do j=1, SPGk_out%multip
      WRITE(*,*) "H(R)|RI=", vecs_H_R(1:4,j)
      end do


      WRITE(*,*) " "
      WRITE(*,*) " ========================================================================"
      WRITE(*,*) " "
      WRITE(*,*) " STEP 4: find the t4(=m/n) values"
      WRITE(*,*) " "
      !NOTE: n = maximum order of the set of matrices
      !NOTE: m_values  = (/0,1,2,3,4,6/) ! not needed
      !NOTE: there's no need to compute the possibilities
      !      just loop 1st loop from  m=0,...,max_order-1
      !      and  loop 2nd loop from  n=1,...,max_order 
      !
      !NOTE: example with rational type
      max_order              = GET_MAX_ORDER(SPGk_out) ! maximum order possible order is 6
      vec_T4_possible_values = GET_POSSIBLE_t4_VALUES(max_order) 
      !WRITE(*,*) " m values = (/0,1,2,3,4,6/)"
      !WRITE(*,*) " n values = ", max_order
      WRITE(*,*) " possible t4 values  = ", vec_T4_possible_values
      size_vec_all_orders=SPGk_out%multip
      allocate(vec_all_orders(size_vec_all_orders))
      vec_all_orders=GET_GROUP_ORDERS(SPGk_out) 
      WRITE(*,*) " vector of SG orders = ", vec_all_orders ! orders to reach to Identity



!> TO DO: automate the cartesian product (with repeat=number of values) 
!>        of possible found t4 values
!>        product('AB', repeat=2) AA AB BA BB
!>
!>        Explicit manual filling:
!>
vec_T4_possible_combinations(1,:)=(/0.,   0./)
vec_T4_possible_combinations(2,:)=(/0.,   0.5/)
vec_T4_possible_combinations(3,:)=(/0.5,  0./)
vec_T4_possible_combinations(4,:)=(/0.5,  0.5/)



WRITE(*,*) " "
WRITE(*,*) " ========================================================================"
WRITE(*,*) " "
WRITE(*,*) " STEP 5:"
WRITE(*,*) " "
WRITE(*,*) " done"

 
WRITE(*,*) " "
WRITE(*,*) " ========================================================================"
WRITE(*,*) " "
WRITE(*,*) " STEP 6:"
WRITE(*,*) " "
WRITE(*,*) " done"




WRITE(*,*) " "
WRITE(*,*) " ========================================================================"
WRITE(*,*) " "
WRITE(*,*) " STEP 7:"
WRITE(*,*) " "

!> Making the supermatrices. Leaving t4 position free (as 44).
!> Casting the Tr and H_R to setOp (of type SSG, with rational matrices).
!> Number of allocated operators  = 2*SPGk_out%multip

!> The known set of generators for this Pnma example are the 6,7, and 8 (nma) 
!> but we exclude, for unknown reasosn, the n

!> The generator set for the SSG is:
!> The super-identity
!> The generator set of the SG
!> The closure matrix for the SSG

call Allocate_SSG_SymmOps(1, 2*SPGk_out%multip, setOpGen) ! double number is a general rule?

MAT_IDENTITY=0//1
do i=1, 5 ! automate superdimension
MAT_IDENTITY(i,i)=1//1
end do

ngen_ssg=0        ! number of generatos of the SSG is num of gen of SG plus 2 (Id plus Closure)
ngen_ssg=ngen+2-1 ! removing the n operator!  automate this
do k=1, ngen_ssg 
 if(k==1) then ! the Id for the SSG
  setOpGen(k)%Mat(:,:)=MAT_IDENTITY
 else if(k==ngen_ssg) then ! the closure for the SSG
  setOpGen(k)%Mat(:,:)    =0//1
  setOpGen(k)%Mat(1:3,1:3)=SPGk_out%SymOp(1)%Rot(:,:)//1 ! ID is the 1st (sumbatrix Rot 3*3)
  setOpGen(k)%Mat(4,4)    =1
  setOpGen(k)%Mat(5,5)    =-1 ! def. of 1' must have inv
  setOpGen(k)%Mat(4,5)    =0.5
 else
  j=point_op(k) ! j is here the number of SG generator from the SG generator pointer
  !
  if(j==1) then
    write(*,*) "excluding the trivial generator n, which is the first one (j=1)" ! automate
  endif
  setOpGen(k)%Mat(:,:)    =0//1
  setOpGen(k)%Mat(1:3,1:3)=SPGk_out%SymOp(j)%Rot(:,:)//1 ! ID sumbatrix Rot 3*3
  setOpGen(k)%Mat(4,1:4)  =vecs_H_R(1:4,j)
  setOpGen(k)%Mat(1:3,5)  =SPGk_out%SymOp(j)%Tr(1:3) ! sic, at the fifth column
  setOpGen(k)%Mat(4,5)    =44//1 ! location of t4
  setOpGen(k)%Mat(5,5)    =1//1  ! sic, just to complete
 end if
enddo


WRITE(*,*) " "
WRITE(*,*) " Checking SSG Generators---------------------------------------------"
do k=1, ngen_ssg !size(setOpGen)
WRITE(*,*) " "
WRITE(*,*) " setOpGen", k
cMat1=print_rational(setOpGen(k)%Mat(:,:))
 do i=1,5
   write(*,"(5a16)") (trim( cMat1(i,j))//"   ",j=1,5)
 end do
enddo
WRITE(*,*) " End Checking SSG Generators------------------------------------------"
WRITE(*,*) " "


! We have the possible t4 values in a vector (cartesian products after removing integers except zero).
! For each possible combination, and in the Pnma case, operators m-a (= 7-8) must be filled with possible t4 values
! starts the cycle of combinations for:
!
! For the case of Pnma, the generator set is n, m, a. But n is trivial. For the other two generators:
!
!               m   a              
! Pnma1'    t4: 0   0              (substitute the t4 values    0-0 in operators m-a)
! Pnma1'    t4: 0.5 0              (substitute the t4 values  0.5-0 in operators m-a)
! Pnma1'    t4: 0   0.5            (...)
! Pnma1'    t4: 0.5 0.5            (...)
!

!--->
!--->



! they are really cartesian products, not 'combinations'
WRITE(*,*) ""
WRITE(*,*) "-----------c vector: ", size(vec_T4_possible_combinations, dim=1)
WRITE(*,*) "-----------c vector: ", vec_T4_possible_combinations(1,:)
WRITE(*,*) "-----------c vector: ", vec_T4_possible_combinations(2,:)
WRITE(*,*) "-----------c vector: ", vec_T4_possible_combinations(3,:)
WRITE(*,*) "-----------c vector: ", vec_T4_possible_combinations(4,:)
WRITE(*,*) ""


!> Computing all the compatible SSGs 
!> their total number being the number of combinations
!>
do ii=1, size(vec_T4_possible_combinations, dim=1) 
  WRITE(*,*) "================================================== Starting: t4 Combination Num.--->", ii
  WRITE(*,*) "=== Combining numbers: ",vec_T4_possible_combinations(ii,:)
  !>
  !> Substitute the t4 value in the corresponding generator
  !> TO DO: a generalization of this step 
  !> 
  setOpGen(2)%Mat(4,5)=vec_T4_possible_combinations(ii,1)
  setOpGen(3)%Mat(4,5)=vec_T4_possible_combinations(ii,2)
  !>
  WRITE(*,*) "=== Number of SSG generators before expansion =", ngen_ssg
  k=0 ! to get the number of elements after expansion
  !>
  Call Gen_Group(ngen_ssg, setOpGen, k)
  Write(*,*) "=== Number of SSG elements after expansion    = ", k 
  if(allocated(matrix)) deallocate(Matrix)
  allocate(Matrix(5,5)) !automate
  forma="( a8)"
  write(forma(2:2),"(i1)") 5
  do i=1,k !all elements
     matrix=print_rational(setOpGen(i)%Mat)
     write(unit=*,fmt="(a,i3)") "  Rational Operator #",i
     do j=1,5 !Dp, automate
        write(unit=*,fmt=forma) (trim( Matrix(j,k))//" ",k=1,5)!Dp automate
     end do
  end do
  WRITE(*,*) "==="
  WRITE(*,*) "================================================== End of:   t4 Combination Num.--->",ii
  WRITE(*,*) ""
end do


!> writing final results in SSG info format
!>
WRITE(*,*) ""
WRITE(*,*) ""
SSG_out%SymOp=setOpGen
call Write_SSG(SSG_out) 
WRITE(*,*) ""
WRITE(*,*) ""


end do ! prompting loop



!--- FUNCTIONS--------------------------------------------

contains

function GET_MAX_ORDER(SPG_in) result(maxorder)
! this function returns the maximum order 
! (defined as the exponent to become the Identity) 
! of a group of rotational operators

type(Space_Group_type), intent(in)    :: SPG_in
!logical,   optional,    intent(in)    :: full
integer                               :: maxorder

real(kind=cp), dimension(3,3)     :: ID, CM, TM 
integer                           :: cnt, i, j
!integer, dimension(SPG_in%multip) :: vec_all_orders

TM=SPG_in%Symop(1)%rot ! Current Matrix
CM=SPG_in%Symop(1)%rot ! Initiallized to the Identity
ID=SPG_in%Symop(1)%rot ! Identity matrix is always the first of the SG

maxorder=0

out:do i=1, SPG_in%multip
   CM=SPG_in%Symop(i)%rot ! Current Matrix
   TM=SPG_in%Symop(1)%rot ! Test Matrix (initiallized to the Id)
   cnt=0 
   in:do j=2,6 !order is never greater than 6
       cnt=cnt+1
       TM=matmul(TM,CM)
       if(Equal_Matrix(TM,ID,3)) then
          !cnt is the order of the i-operator; load out
          !vec_all_orders(i)=cnt
          if(cnt > maxorder) then
             maxorder=cnt
          end if
          cycle out
       else
          continue
       end if
  
  end do in
end do out

!if (present(full))  then
!  return vec_all_orders
!else
  return
!end if
end function GET_MAX_ORDER

!-----------------------------------------------------

function GET_POSSIBLE_t4_VALUES(maximum_order) result(vec_possible_t4_values)
! this function returns a vector with the possible t4 values
! corresponding to a maximum order (matrix potentiation)
!
integer, intent(in)             :: maximum_order
real, dimension(:), allocatable :: vec_possible_t4_values
!
integer :: init_size, inst_size, max_order, num, den, final_size
real    :: t4
real, dimension(:), allocatable  :: vectorA
real, dimension(:), allocatable  :: vectorB


init_size=0  ! start with a zero-sized vector
inst_size=init_size

max_order=maximum_order

!> using vectorB as copy buffer
allocate(vectorA(init_size))
allocate(vectorB(init_size))
vectorA=-11
vectorB=-11


do num=0, max_order-1
   do den=1, max_order
      t4=real(num)/real(den)
      !WRITE(*,*) "a possible t4 is: ", num, " / ",den, " = ", t4
      
      if( Zbelong( t4 ) .and. (t4/=0)) then
         !WRITE(*,*) "        exclude ", t4 
         !cycle
         continue
      else
         !WRITE(*,*) "        accept ", t4
         if ( ANY( vectorA==t4 ) ) then
            !WRITE(*,*) "        is present: do nothing "
            !cycle
            continue
         else
            !WRITE(*,*) "        is impresent: increase vec size and add "
             
            !> BEFORE: copy A -> B
            if(allocated(vectorB)) deallocate(vectorB)
            allocate(vectorB(inst_size))
            vectorB=vectorA
            vectorA(inst_size)=t4 
            
            !> DURING: increase size, reallocate A
            inst_size=inst_size+1
            if(allocated(vectorA)) deallocate(vectorA)
            allocate(vectorA(inst_size))
            
            !> AFTER: copy B -> A, add new to A 
            vectorA(1:inst_size-1)=vectorB ! sic, since different sizes
            vectorA(inst_size)=t4
            
            !WRITE(*,*) "initial vector------------------------------------------------------- a     ", vectorA
            !WRITE(*,*) "initial vector------------------------------------------------------- b     ", vectorB
            
         end if
      end if
      !        
   end do ! num
end do ! den

! casting the results
final_size=size(vectorA)
if(allocated(vec_possible_t4_values)) deallocate(vec_possible_t4_values)
allocate(vec_possible_t4_values(final_size))
vec_possible_t4_values=vectorA

return

end function GET_POSSIBLE_t4_VALUES

!-----------------------------------------------------
!-----------------------------------------------------
function GET_GROUP_ORDERS(SPG_in) result(vec_all_orders)
! this function returns the maximum order 
! (defined as the exponent to become the Identity) 
! of a group of rotational operators

type(Space_Group_type), intent(in)    :: SPG_in
!logical,   optional,    intent(in)    :: full
integer, dimension(SPG_in%multip+1)     :: vec_all_orders

real(kind=cp), dimension(3,3)     :: ID, CM, TM 
integer                           :: cnt, i, j, maxorder

TM=SPG_in%Symop(1)%rot ! Current Matrix
CM=SPG_in%Symop(1)%rot ! Initiallized to the Identity
ID=SPG_in%Symop(1)%rot ! Identity matrix is always the first of the SG

maxorder=0

out:do i=1, SPG_in%multip
   CM=SPG_in%Symop(i)%rot ! Current Matrix
   TM=SPG_in%Symop(1)%rot ! Test Matrix (initz. to Id)
   cnt=0 
   in:do j=2,6 !order is never greater than 6
       cnt=cnt+1
       TM=matmul(TM,CM)
       if(Equal_Matrix(TM,ID,3)) then
          !cnt is the order of the i-operator; load out
          vec_all_orders(i)=cnt
          if(cnt > maxorder) then
             maxorder=cnt
          end if
          cycle out
       else       
          continue
       end if
  
  end do in
end do out

return

end function GET_GROUP_ORDERS



!----------------------------------------------------------------------------

End Program SSPG_Info
