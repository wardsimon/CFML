    Module ssg_datafile
      ! read data for (3+d)-dimensional superspace groups (d=1,2,3)
      implicit none

      public :: Read_SSG, Read_single_SSG
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
      use CFML_String_Utilities, only: pack_string
      use CFML_Rational_Arithmetic

      Implicit None
      private
      public :: Allocate_SSG_SymmOps, Set_SSG_Reading_Database, Write_SSG, Gen_Group, &
                Get_SSymSymb_from_Mat

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
        type(rational),      allocatable,dimension(:,:,:):: Om            !Operator matrices (3+d+1,3+d+1,Multip) common denominator at (4+d,4+d)
        type(SSym_Oper_Type),allocatable, dimension(:)   :: SymOp         ! Crystallographic symmetry operators
        character(len=80),   allocatable,dimension(:)    :: SymopSymb     ! Alphanumeric Symbols for SYMM
      End Type SuperSpaceGroup_Type

      logical, public :: Err_ssg
      character(len=80), public :: Err_ssg_mess

    Contains

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

    Subroutine Get_SSymSymb_from_Mat(Mat,Symb,x1x2x3_type)
       !---- Arguments ----!
       type(rational),dimension(:,:), intent( in) :: Mat
       character (len=*),             intent(out) :: symb
       logical, optional,             intent( in) :: x1x2x3_type

       !---- Local Variables ----!
       character(len=*),dimension(10),parameter :: xyz=(/"x","y","z","t","u","v","w","p","q","r"/)
       character(len=*),dimension(10),parameter :: x1x2x3=(/"x1 ","x2 ","x3 ","x4 ","x5 ","x6 ","x7 ","x8 ","x9 ","x10"/)
       character(len=3),dimension(10) :: x_typ
       character(len= 15)               :: car
       character(len= 40),dimension(10) :: sym
       integer                          :: i,j,Dp,d
       Dp=size(Mat,dim=1)
       d=Dp-1
       x_typ=xyz
       if(present(x1x2x3_type)) then
         if(x1x2x3_type) x_typ=x1x2x3
       end if
       !---- Main ----!
       symb=" "
       do i=1,d
          sym(i)=" "
          do j=1,d
             if(Mat(i,j) == 1) then
                sym(i) = trim(sym(i))//"+"//trim(x_typ(j))
             else if(Mat(i,j) == -1) then
                sym(i) =  trim(sym(i))//"-"//trim(x_typ(j))
             else if(Mat(i,j) /= 0) then
               car=" "
               write(unit=car,fmt="(i3,a)") int(Mat(i,j)),trim(x_typ(j))
               if(Mat(i,j) > 0) car="+"//trim(car)
               sym(i)=trim(sym(i))//pack_string(car)
             end if
          end do
          !Write here the translational part for each component
          if (Mat(i,Dp) /= 0) then
             car=adjustl(print_rational(Mat(i,Dp)))
             if(car(1:1) == "-") then
                sym(i)=trim(sym(i))//trim(car)
             else
                sym(i)=trim(sym(i))//"+"//trim(car)
             end if
          end if
          sym(i)=adjustl(sym(i))
          if(sym(i)(1:1) == "+")  then
            sym(i)(1:1) = " "
            sym(i)=adjustl(sym(i))
          end if
          sym(i)=pack_string(sym(i))
       end do
       symb=sym(1)
       do i=2,d
         symb=trim(symb)//","//trim(sym(i))
       end do
       car=print_rational(Mat(Dp,Dp))
       symb=trim(symb)//","//trim(car)
    End Subroutine Get_SSymSymb_from_Mat

    Subroutine Write_SSG(SpaceGroup,iunit,full,x1x2x3_typ)
      type(SuperSpaceGroup_Type), intent(in) :: SpaceGroup
      integer, optional,          intent(in) :: iunit
      logical, optional,          intent(in) :: full
      logical, optional,          intent(in) :: x1x2x3_typ
      !---- Local variables ----!
      integer :: lun,i,j,k,D,Dp,nlines
      character(len=40),dimension(:),  allocatable :: vector
      character(len=40),dimension(:,:),allocatable :: matrix
      character(len=15) :: forma
      integer,  parameter                      :: max_lines=192
      character (len=100), dimension(max_lines):: texto
      logical                                  :: print_latt,xyz_typ

      lun=6
      print_latt=.false.
      if (present(iunit)) lun=iunit
      if (present(full))  print_latt=.true.
      xyz_typ=.false.
      if(present(x1x2x3_typ)) then
        if(x1x2x3_typ) xyz_typ=.true.
      end if

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
        call Get_SSymSymb_from_Mat(SpaceGroup%SymOp(i)%Mat,texto(1),xyz_typ)
        write(unit=lun,fmt="(a,i3,a)") "  Operator # ",i,"  "//trim(texto(1))
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

    Program read_ssg_datafile
      use ssg_datafile
      use CFML_SuperSpaceGroups
      use CFML_Rational_Arithmetic

      implicit none

      integer :: iclass,nmod,i,j,k,m,multip,Dp
      character(len=20)  :: forma
      character(len=130) :: message
      logical :: ok
      type(SuperSpaceGroup_Type) :: SSpaceGroup
      character(len=40),dimension(:,:),allocatable :: matrix
      character(len=60) :: Operator_Symbol

      call Read_SSG(ok,message)
      if(.not. ok) then
        write(*,"(a)") "   !!! "//message//" !!!"
        stop
      end if
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
        if(m <= 0) exit
        if(m > 16697) then
          write(*,"(a)") " => There are only 16697 superspace groups in the database! "
          cycle
        end if
        !Call Read_single_SSG(m,ok,Message)
        !Call Read_SSG(ok,Message)
        !if(.not. ok) then
        !  write(*,"(a)") "   !!! "//trim(message)//" !!!"
        !  stop
        !end if
        !write(*,"(3(a,i5))") " Group number:",igroup_number(m), " Bravais class:",igroup_class(m), "  Basic Group #:",igroup_spacegroup(m)
        !iclass=igroup_class(m)
        !nmod=iclass_nmod(iclass)
        !write(*,"(a,tr4,a)")  group_nlabel(m), group_label(m)
        !write(*,"(a,i3)") " Number of operators:", igroup_nops(m)
        !do k=1,igroup_nops(m)
        !  write(*,"(a,i3)") " Operator #",k
        !  forma="(  i4)"
        !  write(forma(2:3),"(i2)") nmod+4
        !  do i=1,nmod+4
        !    write(*,forma) igroup_ops(i,1:nmod+4,k,m)
        !    !write(*,*)(((igroup_ops(i,j,k,m),i=1,nmod+4),j=1,nmod+4), k=1,igroup_nops(m))
        !  end do
        !end do
        call Set_SSG_Reading_Database(m,SSpaceGroup,ok,message)
        if(.not. ok) then
          write(*,"(a)") "   !!! "//message//" !!!"
          stop
        end if
        call Write_SSG(SSpaceGroup,full=.true.)
        Dp=size(SSpaceGroup%SymOp(1)%Mat,dim=1)
        if(SSpaceGroup%Centred == 2) then
          SSpaceGroup%SymOp(4)%Mat = -SSpaceGroup%SymOp(1)%Mat
          SSpaceGroup%SymOp(4)%Mat(Dp,Dp) = 1
        end if
        Call Gen_Group(4,SSpaceGroup%SymOp,multip)
        !Writing of the rational operator matrices
        Write(*,*) " Number of generated operators: ",multip
        if(allocated(matrix)) deallocate(Matrix)
        allocate(Matrix(Dp,Dp))
        forma="( a8)"
        write(forma(2:2),"(i1)") Dp
        do i=1,multip
          call Get_SSymSymb_from_Mat(SSpaceGroup%SymOp(i)%Mat,Operator_Symbol,.true.)
          write(unit=*,fmt="(a,i3,a)") "  Operator # ",i,"  "//trim(Operator_Symbol)

          !matrix=print_rational(SSpaceGroup%SymOp(i)%Mat)
          !write(unit=*,fmt="(a,i3)") "  Rational Operator #",i
          !do j=1,Dp
          !   write(unit=*,fmt=forma) (trim( Matrix(j,k))//" ",k=1,Dp)
          !end do
        end do
        !write(*,"(a,i3)") " Reflection conditions: ",igroup_nconditions(m)
        !if(igroup_nconditions(m) > 0)then
        ! do k=1,igroup_nconditions(m)
        !   write(*,*)((igroup_condition1(i,j,k,m),i=1,nmod+3),j=1,nmod+3),(igroup_condition2(j,k,m),j=1,nmod+4)
        ! end do
        !end if
      end do

    End Program read_ssg_datafile
