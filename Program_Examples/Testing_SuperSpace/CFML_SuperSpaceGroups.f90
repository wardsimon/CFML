  Module CFML_SuperSpaceGroups
    use CFML_GlobalDeps,       only: sp,dp,cp,tpi
    use CFML_String_Utilities, only: pack_string, Get_Separator_Pos
    use CFML_Math_General,     only: sort, trace, iminloc, SVDcmp
    use CFML_Crystal_Metrics,  only: Crystal_Cell_Type
    use CFML_Crystallographic_Symmetry, only: Space_Group_Type, Set_SpaceGroup,Get_Generators_From_SpGSymbol, &
                                              set_Intersection_SPG, Write_SpaceGroup
    use CFML_Propagation_Vectors, only: K_Star, Set_Gk, Group_k_Type
    use CFML_ssg_datafile
    use CFML_Rational_Arithmetic_test
    use Matrix_Mod

    Implicit None
    private
    public :: Allocate_SSG_SymmOps, Set_SSG_Reading_Database, Write_SSG, Gen_Group, &
              Get_SSymSymb_from_Mat,Get_Mat_From_SSymSymb, Gen_SReflections, Set_SSGs_from_Gkk, &
              Gen_SSGroup


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
      integer                                          :: MagType       ! 1: No time inversion present, 2: Paramagnetic, 3: 1' is not present, 4: 1' is present
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
      !type(rational),      allocatable,dimension(:,:,:):: Om            !Operator matrices (3+d+1,3+d+1,Multip) common denominator at (4+d,4+d)
      type(SSym_Oper_Type),allocatable, dimension(:)   :: SymOp         ! Crystallographic symmetry operators
      character(len=80),   allocatable,dimension(:)    :: SymopSymb     ! Alphanumeric Symbols for SYMM
    End Type SuperSpaceGroup_Type

    !!----
    !!---- TYPE: sReflect_Type
    !!----
    !!---- Type, Public :: sReflect_Type
    !!----    integer,dimension(:),allocatable :: H       ! H
    !!----    integer                          :: Mult=0  ! mutiplicity
    !!----    real(kind=cp)                    :: S=0.0   ! Sin(Theta)/lambda=1/2d
    !!----    integer                          :: imag=0  !=0 nuclear reflection, 1=magnetic, 2=both
    !!---- End Type sReflect_Type
    !!----
    !!----

    Type, Public :: sReflect_Type
       integer,dimension(:),allocatable :: H       ! H
       integer                          :: Mult=0  ! mutiplicity
       real(kind=cp)                    :: S=0.0   ! Sin(Theta)/lambda=1/2d
       integer                          :: imag=0  !=0 nuclear reflection, 1=magnetic, 2=both
    End Type sReflect_Type

    Type, Public, extends(sReflect_Type) :: sReflection_Type
       real(kind=cp)        :: Fo    ! Observed Structure Factor
       real(kind=cp)        :: Fc    ! Calculated Structure Factor
       real(kind=cp)        :: SFo   ! Sigma of  Fo
       real(kind=cp)        :: Phase ! Phase in degrees
       real(kind=cp)        :: A     ! real part of the Structure Factor
       real(kind=cp)        :: B     ! Imaginary part of the Structure Factor
    End Type sReflection_Type

    Type, Public, extends(sReflection_Type) :: gReflection_Type
       real(kind=cp)                    :: mIvo          ! Observed modulus of the Magnetic Interaction vector
       real(kind=cp)                    :: sigma_mIvo    ! Sigma of observed modulus of the Magnetic Interaction vector
       complex(kind=cp),dimension(3)    :: msF           ! Magnetic Structure Factor
       complex(kind=cp),dimension(3)    :: mIv           ! Magnetic Interaction vector
    End Type gReflection_Type

    !!----
    !!---- TYPE :: SREFLECTION_LIST_TYPE
    !!--..
    !!---- Type, public :: sReflection_List_Type
    !!----    integer                                        :: NRef ! Number of Reflections
    !!----    type(sReflect_Type),allocatable,dimension(:) :: Ref  ! Reflection List
    !!---- End Type sReflection_List_Type
    !!----
    !!---- Update: February - 2005
    !!
    Type, public :: sReflection_List_Type
       integer                                          :: NRef  ! Number of Reflections
       class(sReflect_Type),allocatable, dimension(:)   :: Ref   ! Reflection List
    End Type sReflection_List_Type


    logical,            public :: Err_ssg
    character(len=180), public :: Err_ssg_mess
    real(kind=cp), parameter, private :: eps_ref  = 0.0002_cp
    integer, private, parameter :: max_mult=1024

  Contains

    !!---- function multiply_ssg_symop(Op1,Op2) result (Op3)
    !!----
    !!---- The arguments are two SSym_Oper_Type operators and the
    !!---- result is another SSym_Oper_Type operator. This implements
    !!---- the operator (*) between SSym_Oper_Type objects making the
    !!---- rational multiplications of the D+1 matrices and taking
    !!---- only positive translations modulo-1.
    !!----
    !!----   Created : February 2017 (JRC)
    function multiply_ssg_symop(Op1,Op2) result (Op3)
      type(SSym_Oper_Type), intent(in) :: Op1,Op2
      type(SSym_Oper_Type)             :: Op3
      integer :: n,d,i
      n=size(Op1%Mat,dim=1)
      allocate(Op3%Mat(n,n))  !needed for f95
      d=n-1
      Op3%Mat=matmul(Op1%Mat,Op2%Mat) !automatic allocation in f2003
      !write(*,*) Op3%Mat(1:d,n)
      Op3%Mat(1:d,n)=mod(Op3%Mat(1:d,n),1_ik)
      !call rational_modulo_lat(Op3%Mat(1:d,n))
      ! do i=1,d
      ! 	  write(*,*) Op3%Mat(i,n)
      !  	Op3%Mat(i,n)=mod(Op3%Mat(i,n),1_ik)
      ! end do
       do i=1,d
       	 do
            if(Op3%Mat(i,n) < 0_ik//1_ik) then
            	Op3%Mat(i,n) = Op3%Mat(i,n) + 1_ik
            else
            	exit
            end if
         end do
      end do
    end function multiply_ssg_symop

    Pure Subroutine reduced_translation(Mat)
    	type(rational), dimension(:,:), intent( in out) :: Mat
    	integer :: d,i,n
    	n=size(Mat,dim=1)
    	d=n-1
      do i=1,d
      	 do
           if(Mat(i,n) < 0_ik) then
           	Mat(i,n) = Mat(i,n) + 1_ik
           else
           	exit
           end if
         end do
      	 do
           if(Mat(i,n) > 1_ik) then
           	Mat(i,n) = Mat(i,n) - 1_ik
           else
           	exit
           end if
        end do
      end do
    End Subroutine reduced_translation


    Subroutine Allocate_SSG_SymmOps(d,multip,SymOp)
      integer,                                        intent(in)      :: d,multip
      type(SSym_Oper_Type),allocatable, dimension(:), intent(in out)  :: SymOp
      integer :: i,Dd

      if(allocated(SymOp)) deallocate(SymOp)
      allocate(SymOp(multip))
      Dd=3+d+1
      do i=1,multip
        allocate(SymOp(i)%Mat(Dd,Dd))
        SymOp(i)%Mat=0_ik//1_ik
      end do
    End Subroutine Allocate_SSG_SymmOps

    !This subroutine assumes that Op contains the identity as the first operator, followed
    !by few non-equal generators. The value of ngen icludes also the identity
    Subroutine Gen_Group(ngen,Op,multip,table)
      integer,                                        intent(in)     :: ngen
      type(SSym_Oper_Type), dimension(:),             intent(in out) :: Op
      integer,                                        intent(out)    :: multip
      integer, dimension(:,:), allocatable, optional, intent(out)    :: table
      !--- Local variables ---!
      integer :: i,j,k,n, nt,max_op
      type(SSym_Oper_Type) :: Opt
      logical, dimension(size(Op),size(Op)) :: done
      integer, dimension(size(Op),size(Op)) :: tb
      max_op=size(Op)
      done=.false.
      done(1,:) = .true.
      done(:,1) = .true.
      tb(1,:) = [(i,i=1,max_op)]
      tb(:,1) = [(i,i=1,max_op)]
      nt=ngen
      Err_ssg=.false.
      Err_ssg_mess=" "

      do_ext:do
        n=nt
        do i=1,n
          do_j:do j=1,n
            if(done(i,j)) cycle
            Opt=Op(i)*Op(j)
            do k=1,nt
              if(equal_rational_matrix(Opt%Mat,Op(k)%Mat)) then
                tb(i,j)=k
                done(i,j)=.true.
                cycle do_j
              end if
            end do
            done(i,j)=.true.
            nt=nt+1
            if(nt > max_op) then
              nt=nt-1
              exit do_ext
            end if
            tb(i,j)=nt
            Op(nt)=Opt
          end do do_j
        end do
        if ( n == nt) exit do_ext
      end do do_ext
      if(any(done(1:nt,1:nt) .eqv. .false. ) ) then
      	Err_ssg=.true.
      	Err_ssg_mess="Table of SSG operators not exhausted! Increase the expected order of the group!"
      end if
      if(nt == max_op) then
        Err_ssg=.true.
      	Err_ssg_mess= "Max_Operations reached! "//trim(Err_ssg_mess)
      end if
      
      multip=nt
      if(present(table)) then
        allocate(Table(multip,multip))
        Table=tb
      end if
    End Subroutine Gen_Group

    Subroutine Gen_SSGroup(ngen,gen,SSG,x1x2x3_type,table)
      integer,                             intent(in)  :: ngen
      type(SSym_Oper_Type), dimension(:),  intent(in)  :: gen
      type(SuperSpaceGroup_Type),          intent(out) :: SSG
      character(len=*), optional,          intent(in)  :: x1x2x3_type
      integer, dimension(:,:), allocatable, optional, intent(out) :: table
      !--- Local variables ---!
      integer :: i,j,k,n, nt, Dd, Dex, ngeff, nlat,nalat, d
      type(SSym_Oper_Type), dimension(:), allocatable :: Op
      type(SSym_Oper_Type) :: Opt
      type(Rational),dimension(size(gen(1)%Mat,dim=1),size(gen(1)%Mat,dim=2)) :: identity
      integer, dimension(max_mult) :: ind_lat,ind_alat
      logical :: esta

      nt=ngen
      Err_ssg=.false.
      Err_ssg_mess=" "
      SSG%standard_setting=.false.
      Dex=size(gen(1)%Mat,dim=1) 
      Dd=Dex-1
      d=Dd-3
      SSG%d= Dex-4
      call Allocate_SSG_SymmOps(d,max_mult,Op)
      call Identity_Matrix(Dex,identity)
      Op(1)%Mat=identity
      j=0
      do i=1,ngen
      	if(equal_rational_matrix(gen(1)%Mat,identity)) then
      		j=i
      		exit
      	end if
      end do
      k=1
      do i=1,ngen
      	if(i == j) cycle
      	k=k+1
      	Op(k)=gen(i)
      end do
      ngeff=k
      nt=ngeff
      
      if(present(table)) then
        call Gen_Group(ngeff,Op,nt,table)
      else
        call Gen_Group(ngeff,Op,nt)
      end if
      
      if(Err_ssg) then
      	write(unit=*,fmt="(a)") " => ERROR !  "//trim(Err_ssg_mess)
      	return
      end if

      if(allocated(SSG%SymOp)) deallocate(SSG%SymOp)
      allocate(SSG%SymOp(nt))
      if(allocated(SSG%SymOpSymb)) deallocate(SSG%SymOpSymb)
      allocate(SSG%SymOpSymb(nt))
      if(allocated(SSG%Centre_coord))  deallocate(SSG%Centre_coord)
      allocate(SSG%Centre_coord(Dd))
      if(allocated(SSG%time_rev))  deallocate(SSG%time_rev)
      allocate(SSG%time_rev(nt))
      SSG%SymOp(:)=Op(1:nt)
      SSG%multip=nt
      do i=1,nt
      	call Get_SSymSymb_from_Mat(SSG%SymOp(i)%Mat,SSG%SymOpSymb(i),"x1x2x3")
      end do
      
      !Search for an inversion centre
      j=0
      do i=2,nt
      	if(equal_rational_matrix(SSG%SymOp(i)%Mat(1:Dd,1:Dd),-identity(1:Dd,1:Dd))) then
      		j=i
      		exit
      	end if
      end do
      SSG%Centre_coord=0_ik//1_ik
      if( j /= 0) then
      	k=0
      	do i=1,Dd
      		if(SSG%SymOp(i)%Mat(Dex,i) /= 0_ik//1_ik) then
      			k=i
      			exit
      		end if
      	end do
      	if(k == 0) then
      		SSG%Centre="Centric with centre at origin"
      		SSG%Centred=2
      	else
      		SSG%Centre="Centric with centre NOT at origin"
      		SSG%Centred=0
      		SSG%Centre_coord(:)=SSG%SymOp(k)%Mat(Dex,1:Dd)/2_ik
      	end if
      else
      	SSG%Centre="Acentric"
      	SSG%Centred=1
      end if

      nlat=0
      nalat=0
      do i=2,nt
      	if(equal_rational_matrix(SSG%SymOp(i)%Mat(1:Dd,1:Dd),identity(1:Dd,1:Dd))) then
      	  if(SSG%SymOp(i)%Mat(Dex,Dex) == 1_ik) then
      	  	 nlat=nlat+1
      	  	 ind_lat(nlat)=i
      	  else
      	  	 nalat=nalat+1
      	  	 ind_alat(nalat)=i
      	  end if
      	end if
      end do
      SSG%Num_Lat=nlat+1
      SSG%Num_aLat=nalat
      
      !Determine the type of magnetic group and select centring translations and anti-translations
      !if(any()) <- implement any for rational arrays      
      do i=1,nt
      	SSG%time_rev(i)=SSG%SymOp(i)%Mat(Dex,Dex)
      end do
      if(.not. any(SSG%time_rev == -1)) SSG%MagType=1
      
      if(nalat /= 0) then
      	 SSG%MagType=4
      	 allocate(SSG%aLatt_trans(Dd,nalat))
      	 SSG%aLatt_trans=0_ik//1_ik
      	 do i=1,nalat
      	 	  j=ind_alat(i)
      	 	  SSG%aLatt_trans(:,i)=SSG%SymOp(j)%Mat(1:Dd,Dex)
         end do
      else if(SSG%MagType /= 1 ) then
      	 SSG%MagType=3
      end if
      
      if(nlat /=  0) then
      	 SSG%SPG_Lat="X"
      	 allocate(SSG%Latt_trans(Dd,nlat+1))
      	 SSG%Latt_trans=0_ik//1_ik
      	 do i=1,nlat
      	 	  j=ind_lat(i)
      	 	  SSG%Latt_trans(:,i+1)=SSG%SymOp(j)%Mat(1:Dd,Dex)
         end do
      end if
            
    End Subroutine Gen_SSGroup

    Subroutine Set_SSGs_from_Gkk(SpG,nk,kv)!,ssg,nss)
      type(Space_Group_Type),                                intent(in)  :: SpG
      integer,                                               intent(in)  :: nk
      real(kind=cp),dimension(:,:),                          intent(in)  :: kv
      !type(SuperSpaceGroup_Type), dimension(:), allocatable, intent(out) :: ssg
      !integer,                                               intent(out) :: nss
      !--- Local variables ---!
      character(len=132) :: line
      integer :: i,j,k,Dd
      type(Space_Group_Type),dimension(nk) :: Gkks !extended Little Groups
      type(Space_Group_Type)               :: Gkk  !extended Little Group (Intersection of Litte Groups for all nk)
      type(SuperSpaceGroup_Type)           :: trial_ssg
      Type (Group_k_Type)                  :: Gk
      real(kind=cp), dimension(3+nk)       :: tr
      integer,       dimension(3+nk,3+nk)  :: Mat

      !Initializing variables
      Err_ssg=.false.
      Err_ssg_mess=" "
      Dd=3+nk+1 !Dimension of the extended matrices of ssg

      !Determine the extended little Groups
      do i=1,nk
         call K_Star(kv(:,i),SpG,Gk,.true.)
         call Set_Gk(Gk,Gkks(i),.true.)
         call Write_Spacegroup(Gkks(i),Full=.true.)
      end do
      if(nk > 1) then !Determine the intersection of the space groups
        call set_Intersection_SPGt(Gkks,Gkk)
      else
        Gkk=Gkks(1)
      end if
      call Write_Spacegroup(Gkk,Full=.true.)
      !Derive possible superspace groups from Gkk
      !First: colorless groups of type 1
      !Second: black and white groups of type 3
      !Third: black and white groups of type 4

    End Subroutine Set_SSGs_from_Gkk


    Subroutine Set_Intersection_SPGt(SpGs,SpG)
      Type (Space_Group_Type),dimension(:), intent(in)   :: SpGs
      Type (Space_Group_Type),              intent(out)  :: SpG
      !--- Local Variables ---!
      integer :: i,j,k,ng,ipos,n
      character(len=40),dimension(192) :: gen
      logical,dimension(size(SpGs(:))) :: estak

      write(*,"(/a/)") " => Entering subroutine Set_Intersection_SPGt "
      ipos=iminloc(SpGs(:)%multip)
      ng=1
      gen(1)="x,y,z"
      n=size(SpGs(:))
      estak=.false.
      estak(ipos)=.true.

      write(*,"(3(a,i3))") " => Number of space groups: ",n, "  Position: ",ipos,"  Multiplicity:",SpGs(ipos)%multip

      do_ext:do i=2,SpGs(ipos)%multip
        ng=ng+1
        gen(ng)=SpGs(ipos)%SymopSymb(i)

        do j=1,n
           if(j == ipos) cycle
           estak(j)=.false.
           do k=2,SpGs(j)%multip
           	 write(*,"(2i4,a,a)") k,ng, "   "//trim(gen(ng))//"   "//trim(SpGs(j)%SymopSymb(k))
             if(trim(SpGs(j)%SymopSymb(k)) == trim(gen(ng))) then
             	 estak(j)=.true.
             	 exit
             end if
           end do
        end do
        write(*,*) estak
        k=count(estak(1:n))
        if(k /= n) then
          write(*,"(a,i3)") "  Operator: "//trim(gen(ng))//"  Discarded  ng=",ng-1
          ng=ng-1
          cycle
        end if
        !Passing here means that the operator is common to all space groups
        if(ng > 192) exit
      end do do_ext
      call Set_SpaceGroup(" ",SpG,gen,ng,Mode="GEN")

    End Subroutine Set_Intersection_SPGt



    Subroutine Set_SSG_Reading_Database(num,ssg,ok,Mess,x1x2x3_type)
      integer,                    intent(in)  :: num
      type(SuperSpaceGroup_Type), intent(out) :: ssg
      Logical,                    intent(out) :: ok
      character(len=*),           intent(out) :: Mess
      character(len=*),optional,  intent(in)  :: x1x2x3_type
      !
      integer :: i,j,nmod,Dd,D,iclass,m
      type(rational), dimension(:,:), allocatable :: Inv
      type(SSym_Oper_Type)                        :: transla
      character(len=15) :: forma,xyz_typ
      logical :: inv_found

      if(.not. ssg_database_allocated)  then
         call Read_single_SSG(num,ok,Mess)
         if(.not. ok) return
      end if
      iclass=igroup_class(num)
      nmod=iclass_nmod(iclass)
      D=3+nmod
      Dd=D+1
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


      ssg%MagType= 1                         ! No time-reversal is associated with the symmetry operators
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
      Allocate(ssg%SymopSymb(ssg%Multip))

      ssg%kv=0.0; ssg%Latt_trans=0_ik//1_ik; ssg%aLatt_trans=0_ik//1_ik; ssg%time_rev=1_ik; ssg%Centre_coord=0_ik//1_ik
      nmod=iclass_nmod(iclass)
      !forma="(a,i3, i4)"
      !write(forma(7:7),"(i1)") Dd
      do j=1,ssg%Num_Lat
         !write(*,forma) " Vector #",j, iclass_centering(1:Dd,j,iclass)
         ssg%Latt_trans(1:D,j)= rational_simplify(iclass_centering(1:D,j,iclass)//iclass_centering(Dd,j,iclass))
      end do
      do i=1,ssg%NumOps
         ssg%SymOp(i)%Mat= rational_simplify(igroup_ops(1:Dd,1:Dd,i,num)//igroup_ops(Dd,Dd,i,num))
      end do

      !Look for a centre of symmetry
      inv_found=.false.
      allocate(Inv(D,D),transla%Mat(Dd,Dd))
      Inv=0//1; transla%Mat=0//1
      do i=1,D
        Inv(i,i) = -1//1
      end do
      do i=1,Dd
        transla%Mat(i,i) = 1//1
      end do

      do i=1,ssg%NumOps
        if(equal_rational_matrix(Inv,ssg%SymOp(i)%Mat(1:D,1:D))) then
          if(sum(ssg%SymOp(i)%Mat(1:D,Dd)) == 0//1) then
            ssg%Centred=2
            ssg%Centre_coord=0//1
            ssg%Centre="Centrosymmetric with centre at origin"
          else
            ssg%Centred=0
            ssg%Centre_coord=ssg%SymOp(i)%Mat(1:D,Dd)/2_ik
            write(unit=ssg%Centre,fmt="(a)") "Centrosymmetric with centre at : "//print_rational(ssg%Centre_coord)
          end if
          exit
        end if
      end do
      !Extend the symmetry operator to the whole set of lattice centring and
      !set the symmetry operators symbols
      m=ssg%NumOps
      do i=2,ssg%Num_Lat
        transla%Mat(1:D,Dd)=ssg%Latt_trans(1:D,i)
        do j=1,ssg%NumOps
           m=m+1
           ssg%SymOp(m)=ssg%SymOp(j)*transla
        end do
      end do
      if(m /= ssg%Multip) then
        Err_ssg=.true.
        Err_ssg_mess="Error extending the symmetry operators for a centred cell"
      end if
      !Get the symmetry symbols
      xyz_typ="xyz"
      if(present(x1x2x3_type)) xyz_typ=x1x2x3_type
      do i=1,ssg%Multip
        call Get_SSymSymb_from_Mat(ssg%SymOp(i)%Mat,ssg%SymOpSymb(i),xyz_typ)
      end do

    End Subroutine Set_SSG_Reading_Database


    Subroutine Get_Mat_From_SSymSymb(Symb,Mat)
      character(len=*),                intent(in)  :: Symb
      type(rational),dimension(:,:),   intent(out) :: Mat
      !---- local variables ----!
      integer :: i,j,k,Dd, d, np,ns, n,m,inv,num,den,ind,ier
      character(len=*),dimension(10),parameter :: xyz=(/"x","y","z","t","u","v","w","p","q","r"/)
      character(len=*),dimension(10),parameter :: x1x2x3=(/"x1 ","x2 ","x3 ","x4 ","x5 ","x6 ","x7 ","x8 ","x9 ","x10"/)
      character(len=*),dimension(10),parameter :: abc=(/"a","b","c","d","e","f","g","h","i","j"/)
      character(len=3),dimension(10)           :: x_typ
      character(len=len(Symb)), dimension(size(Mat,dim=1))  :: split
      character(len=len(Symb))                              :: string,pSymb,translation,transf
      character(len=10),        dimension(size(Mat,dim=1))  :: subst
      integer,                  dimension(size(Mat,dim=1)-1):: pos,pn
      logical                                               :: abc_transf

      err_ssg=.false.
      err_ssg_mess=" "
      Dd=size(Mat,dim=1)
      d=Dd-1
      abc_transf=.false.
      x_typ=xyz
      i=index(Symb,"x2")
      if(i /= 0) x_typ=x1x2x3
      i=index(Symb,"b")
      if(i /= 0) then
        x_typ=abc
        abc_transf=.true.
      end if
      Mat=0//1
      if(abc_transf) then
        i=index(Symb,";")
        if(i == 0) return !check for errors
        translation=pack_string(Symb(i+1:)//",0")
        pSymb=pack_string(Symb(1:i-1)//",1")
        call Get_Separator_Pos(translation,",",pos,np)
        if(np /= d-1) then
          err_ssg=.true.
          err_ssg_mess="Error in the origin coordinates"
          return
        end if
        j=1
        do i=1,d
          split(i)=translation(j:pos(i)-1)
          j=pos(i)+1
        end do

        do i=1,d
            k=index(split(i),"/")
            if(k /= 0) then
              read(unit=split(i)(1:k-1),fmt=*) num
              read(unit=split(i)(k+1:),fmt=*) den
            else
              den=1
              read(unit=split(i),fmt=*) num
            end if
            Mat(i,Dd)= num//den
        end do
        split=" "

      else

        pSymb=pack_string(Symb)

      end if

      pos=0
      call Get_Separator_Pos(pSymb,",",pos,np)
      if(np /= d) then
      	if(np == d-1) then
      		pSymb=trim(pSymb)//",1"
      		np=np+1
      		pos(np)=len_trim(pSymb)-1
        else
          err_ssg=.true.
          err_ssg_mess="Error in the symbol of the operator"
          return
        end if
      end if

      read(unit=pSymb(pos(np)+1:),fmt=*,iostat=ier) inv
      if(ier == 0) then
        Mat(Dd,Dd)=inv//1
      else
        Mat(Dd,Dd)=1//1
      end if
      j=1
      do i=1,d
        split(i)=pSymb(j:pos(i)-1)
        j=pos(i)+1
        !write(*,"(i3,a)")  i, "  "//trim(split(i))
      end do
      !Analysis of each item
                      !   |     |  |  |  <- pos(i)
      !items of the form  -x1+3/2x2-x3+x4+1/2
                      !   |  |     |  |  |  <- pn(m)
                      !   -  +3/2  -  +  +1/2
      do i=1,d
        string=split(i)
        pn=0; m=0
        do j=1,len_trim(string)
          if(string(j:j) == "+" .or. string(j:j) == "-" ) then
            m=m+1
            pn(m)=j
          end if
        end do
        if(m /= 0) then
          if(pn(1) == 1) then
            n=0
          else
            subst(1)=string(1:pn(1)-1)
            n=1
          end if
          do k=1,m-1
            n=n+1
            subst(n)=string(pn(k):pn(k+1)-1)
          end do
          if(pn(m) < len_trim(string)) then
            n=n+1
            subst(n)=string(pn(m):)
          end if
        else
          n=1
          subst(1)=string
        end if

        !write(*,"(11a)") " ---->  "//trim(string),("  "//trim(subst(k)),k=1,n)

        ! Now we have n substrings of the form +/-p/qXn  or +/-p/q or +/-pXn/q or Xn or -Xn
        ! Look for each substring the item x_typ and replace it by 1 or nothing
        ! m is reusable now
        do j=1,n
          k=0
          do m=1,d
             k=index(subst(j),trim(x_typ(m)))
             if(k /= 0) then
               ind=m
               exit
             end if
          end do
          if(k == 0) then !pure translation
            k=index(subst(j),"/")
            if(k /= 0) then
              read(unit=subst(j)(1:k-1),fmt=*) num
              read(unit=subst(j)(k+1:),fmt=*) den
            else
              den=1
              read(unit=subst(j),fmt=*) num
            end if
            Mat(i,Dd)= num//den
          else  !Component m of the row_vector
            !suppress the symbol
            subst(j)(k:k+len_trim(x_typ(ind))-1)=" "
            if( k == 1 ) then
              subst(j)(1:1)="1"
            else if(subst(j)(k-1:k-1) == "+" .or. subst(j)(k-1:k-1) == "-") then
              subst(j)(k:k)="1"
            else
              subst(j)=pack_string(subst(j))
            end if
            !Now read the integer or the rational
            k=index(subst(j),"/")
            if( k /= 0) then
              read(unit=subst(j)(1:k-1),fmt=*) num
              read(unit=subst(j)(k+1:),fmt=*) den
            else
              den=1
              read(unit=subst(j),fmt=*) num
            end if
            Mat(i,ind)=num//den
          end if
        end do
      end do

    End Subroutine Get_Mat_From_SSymSymb


    Subroutine Get_SSymSymb_from_Mat(Mat,Symb,x1x2x3_type)
       !---- Arguments ----!
       type(rational),dimension(:,:), intent( in) :: Mat
       character (len=*),             intent(out) :: symb
       character(len=*), optional,    intent( in) :: x1x2x3_type

       !---- Local Variables ----!
       character(len=*),dimension(10),parameter :: xyz=(/"x","y","z","t","u","v","w","p","q","r"/)
       character(len=*),dimension(10),parameter :: x1x2x3=(/"x1 ","x2 ","x3 ","x4 ","x5 ","x6 ","x7 ","x8 ","x9 ","x10"/)
       character(len=*),dimension(10),parameter :: abc=(/"a","b","c","d","e","f","g","h","i","j"/)
       character(len=3),dimension(10) :: x_typ
       character(len= 15)               :: car
       character(len= 40)               :: translation
       character(len= 40),dimension(10) :: sym
       integer                          :: i,j,Dd,d,k
       logical                          :: abc_type
       Dd=size(Mat,dim=1)
       d=Dd-1
       x_typ=xyz
       abc_type=.false.
       if(present(x1x2x3_type)) then
         Select Case (trim(x1x2x3_type))
          Case("xyz")
            x_typ=xyz
          Case("x1x2x3")
            x_typ=x1x2x3
          Case("abc")
            x_typ=abc
            abc_type=.true.
          Case Default
          	x_typ=xyz
         End Select
       end if
       !---- Main ----!
       symb=" "
       translation=" "
       do i=1,d
          sym(i)=" "
          do j=1,d
             if(Mat(i,j) == 1_ik) then
                sym(i) = trim(sym(i))//"+"//trim(x_typ(j))
             else if(Mat(i,j) == -1_ik) then
                sym(i) =  trim(sym(i))//"-"//trim(x_typ(j))
             else if(Mat(i,j) /= 0_ik) then
               car=adjustl(print_rational(Mat(i,j)))
               k=index(car,"/")
               if(k /= 0) then
                 if(car(1:1) == "1") then
                   car=trim(x_typ(j))//car(k:)
                 else if(car(1:2) == "-1") then
                   car="-"//trim(x_typ(j))//car(k:)
                 else
                   car=car(1:k-1)//trim(x_typ(j))//car(k:)
                 end if
               else
                 car=trim(car)//trim(x_typ(j))
               end if
               !write(unit=car,fmt="(i3,a)") int(Mat(i,j)),trim(x_typ(j))
               if(Mat(i,j) > 0_ik) car="+"//trim(car)
               sym(i)=trim(sym(i))//pack_string(car)
             end if
          end do
          !Write here the translational part for each component
          if (Mat(i,Dd) /= 0_ik) then
            car=adjustl(print_rational(Mat(i,Dd)))
            if(abc_type) then
              translation=trim(translation)//","//trim(car)
            else
               if(car(1:1) == "-") then
                  sym(i)=trim(sym(i))//trim(car)
               else
                  sym(i)=trim(sym(i))//"+"//trim(car)
               end if
            end if
          else
            if(abc_type) translation=trim(translation)//",0"
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
       if(abc_type)then
          symb=trim(symb)//";"//trim(translation(2:))
       else
          car=print_rational(Mat(Dd,Dd))
          symb=trim(symb)//","//trim(car)
       end if
    End Subroutine Get_SSymSymb_from_Mat

    Subroutine Write_SSG(SpaceGroup,iunit,full,x1x2x3_typ)
      type(SuperSpaceGroup_Type), intent(in) :: SpaceGroup
      integer, optional,          intent(in) :: iunit
      logical, optional,          intent(in) :: full
      logical, optional,          intent(in) :: x1x2x3_typ
      !---- Local variables ----!
      integer :: lun,i,j,k,D,Dd,nlines
      character(len=40),dimension(:),  allocatable :: vector
      character(len=40),dimension(:,:),allocatable :: matrix
      character(len=15)                            :: forma,xyz_typ
      integer,  parameter                          :: max_lines=192
      character (len=100), dimension(max_lines)    :: texto
      logical                                      :: print_latt

      lun=6
      print_latt=.false.
      if (present(iunit)) lun=iunit
      if (present(full))  print_latt=.true.
      xyz_typ="xyz"
      if(present(x1x2x3_typ)) then
        xyz_typ="x1x2x3"
      end if

      D=3+SpaceGroup%d
      Dd=D+1
      allocate(vector(D),matrix(Dd,Dd))
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
      write(unit=lun,fmt="(a,i3)")           " =>           Magnetic Type: ", SpaceGroup%MagType
      write(unit=lun,fmt="(a,a)")           " =>          Centrosymmetry: ", trim(SpaceGroup%Centre)
      write(unit=lun,fmt="(a,a)")           " =>         Bravais Lattice: ", "  "//trim(SpaceGroup%SPG_Lat)

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
      write(forma(2:2),"(i1)") Dd
      nlines=SpaceGroup%NumOps
      if(print_latt) nlines=SpaceGroup%Multip
      if(present(x1x2x3_typ)) then
         do i=1,nlines
           call Get_SSymSymb_from_Mat(SpaceGroup%SymOp(i)%Mat,texto(1),xyz_typ)
           write(unit=lun,fmt="(a,i3,a)") "  Operator # ",i,"  "//trim(texto(1))
         end do
      else
         do i=1,nlines
           write(unit=lun,fmt="(a,i3,a)") "  Operator # ",i,"  "//trim(SpaceGroup%SymOpSymb(i))
         end do
      end if
    End Subroutine Write_SSG

    Function H_S(hkl,Cell,nk,kv) result(s)
    	integer, dimension(:),                  intent(in) :: hkl
     	type(Crystal_Cell_Type),                intent(in) :: Cell
    	integer, optional,                      intent(in) :: nk
    	real(kind=cp),dimension(:,:), optional, intent(in) :: kv
    	real(kind=cp)   :: s
    	!--- Local variables ---!
    	real(kind=cp), dimension(3) :: h
    	integer :: i

    	h=hkl(1:3)
    	if(present(nk) .and. present(kv)) then
    	   do i=1,nk
    	    	h=h+hkl(3+i)*kv(:,i)
    	   end do
    	end if
      s= 0.5*sqrt( h(1)*h(1)*Cell%GR(1,1) +     h(2)*h(2)*Cell%GR(2,2) + &
                   h(3)*h(3)*Cell%GR(3,3) + 2.0*h(1)*h(2)*Cell%GR(1,2) + &
               2.0*h(1)*h(3)*Cell%GR(1,3) + 2.0*h(2)*h(3)*Cell%GR(2,3) )
    End Function H_S

    Function H_Equal(H,K) Result (Info)
       !---- Arguments ----!
       integer, dimension(:), intent(in) :: h,k
       logical                           :: info
       integer :: i

       info=.true.
       do i=1,size(h)
         if (h(i) /= k(i)) then
          info=.false.
          return
         end if
       end do

    End Function H_Equal


    Function H_Equiv(H,K,SSG,Friedel) Result (Info)
       !---- Arguments ----!
       integer, dimension(:),        intent(in)  :: h,k
       Type (SuperSpaceGroup_Type),  intent (in) :: SSG
       logical, optional,        intent(in)      :: Friedel
       logical                                   :: info

       !---- Local Variables ----!
       integer                                :: i, nops,Dd
       integer, dimension(size(h))            :: hh
       Integer,dimension(size(h),size(h))     :: Mat

       info=.false.
       nops= SSG%numops
       if(SSG%centred /= 1) nops=nops*2

       Dd=size(h)
       do i=1,nops
          Mat=SSG%SymOp(i)%Mat(1:Dd,1:Dd)
          hh = matmul(h,Mat)
          if (H_equal(k,hh)) then
             info=.true.
             exit
          end if
          if (present(Friedel)) then
             if (Friedel) then
                if (h_equal(k,-hh)) then
                   info=.true.
                   exit
                end if
             end if
          end if
       end do

       return
    End Function H_Equiv

    Function H_Lat_Absent(h,Latt,n) result(info)
       integer,        dimension(:),  intent (in) :: h
       type(rational), dimension(:,:),intent (in) :: Latt
       integer,                       intent (in) :: n
       logical                                    :: info

       !---- Local Variables ----!
       integer                          :: k,i
       real(kind=cp),dimension(size(h)) :: Lat
       real(kind=cp)                    :: r1,r2

       info=.false.
       do i=2,n
          Lat=Latt(:,i)
          r1=dot_product(Lat,real(h))
          r2=nint(r1)
          k=nint(2.0*r1)
          if(mod(k,2) /= 0) info=.true.
          exit
       end do
       return
    End Function H_Lat_Absent

    Function H_Absent_SSG(H,SSG) Result(Info)
       !---- Arguments ----!
       integer, dimension(:),            intent (in) :: h
       Type (SuperSpaceGroup_Type),      intent (in) :: SSG
       logical                                       :: info

       !---- Local Variables ----!
       integer, dimension(size(h))        :: k
       integer                            :: i,Dd
       real(kind=cp)                      :: r1,r2
       Integer,dimension(size(h),size(h)) :: Mat
       real(kind=cp),dimension(size(h))   :: tr

       info=.false.
       Dd=size(h)
       do i=1,SSG%multip
          Mat=SSG%SymOp(i)%Mat(1:Dd,1:Dd)
          k = matmul(h,Mat)
          if (h_equal(h,k)) then
             tr=SSG%SymOp(i)%Mat(Dd+1,1:Dd)
             r1=dot_product(tr,real(h))
             r2=nint(r1)
             if (abs(r1-r2) > eps_ref) then
                info=.true.
                exit
             end if
          end if
       end do
    End Function H_Absent_SSG

    Function mH_Absent_SSG(H,SSG) Result(Info)
       !---- Arguments ----!
       integer, dimension(:),       intent (in) :: h
       Type (SuperSpaceGroup_Type), intent (in) :: SSG
       logical                                  :: info

       !---- Local Variables ----!
       integer,       dimension(size(h))        :: k
       integer,       dimension(size(h),size(h)):: Mat
       real(kind=cp), dimension(size(h))        :: tr
       integer                                  :: i,n_id,Dd
       real(kind=cp)                            :: r1 !,r2

       info=.false.
       Dd=size(h)
       Select Case(SSG%MagType)
         Case(1,3)
           n_id=0
           do i=1,SSG%Multip
              Mat=SSG%SymOp(i)%Mat(1:Dd,1:Dd)
              tr= SSG%SymOp(i)%Mat(Dd+1,1:Dd)
              k = matmul(h,Mat)
              if (h_equal(h,k)) then
                 r1=cos(tpi*dot_product(Tr,real(h)))
                 n_id=n_id+Trace(Mat*nint(r1))
              end if
           end do
           if(n_id == 0) info=.true.
         Case(2)
           info=.true.
         Case(4)
           info=.true.
           do i=1,SSG%Num_aLat
              tr=SSG%aLatt_trans(:,i)
              r1=cos(tpi*dot_product(tr,real(h)))
              if(nint(r1) == -1) info=.false.
           end do
       End Select
    End Function mH_Absent_SSG
    
    !!----
    !!---- Function  H_Mult(H, Spacegroup,Friedel)
    !!----    integer, dimension(:),    intent(in) :: h
    !!----    Type (SSpace_Group_Type), intent(in) :: SpaceGroup
    !!----    Logical,                  intent(in) :: Friedel
    !!----
    !!----    Calculate the multiplicity of the reflection
    !!----
    !!----    Created: June - 2017
    !!
    Function H_Mult(H,Spacegroup,Friedel) Result(N)
       !---- Arguments ----!
       integer, dimension(:),       intent (in)  :: h
       Type (SuperSpaceGroup_Type), intent (in)  :: SpaceGroup
       Logical,                     intent (in)  :: Friedel
       integer                                   :: N

       !---- Local Variables ----!
       logical                                      :: esta
       integer, dimension(size(h))                  :: k
       integer, dimension(size(h),size(h))          :: Mat
       integer                                      :: i,j,ng,Dd
       integer, dimension(size(h),SpaceGroup%numops):: klist

       ng=SpaceGroup%numops
       n=1
       Dd=size(h)

       !if NG = 0 (strange case), skip it, fix by Petr
       if (ng > 1) then
           klist(:,1)=h(:)

           do i=2,ng
           	  Mat=SpaceGroup%SymOp(i)%Mat(1:Dd,1:Dd)
              k = matmul(h,Mat)
              esta=.false.
              do j=1,n
                 if (h_equal(k,klist(:,j)) .or. h_equal(-k,klist(:,j))) then
                    esta=.true.
                    exit
                 end if
              end do
              if (esta) cycle
              n=n+1
              klist(:,n) = k
           end do
       end if
       if (Friedel .or. SpaceGroup%centred == 2) then
           n=n*2
       end if

       return
    End Function H_Mult    

    !!----
    !!---- Subroutine  Gen_SReflections(Cell,sintlmax,Num_Ref,Reflex,nk,nharm,kv,maxsinl,SSG,powder)
    !!----   type (Crystal_Cell_Type),                        intent(in)     :: Cell
    !!----   real(kind=cp),                                   intent(in)     :: sintlmax
    !!----   integer,                                         intent(out)    :: num_ref
    !!----   class (sReflect_Type), dimension(:), allocatable,intent(out)    :: reflex
    !!----   integer,                       optional,         intent(in)     :: nk
    !!----   integer,       dimension(:),   optional,         intent(in)     :: nharm
    !!----   real(kind=cp), dimension(:,:), optional,         intent(in)     :: kv
    !!----   real(kind=cp), dimension(:),   optional,         intent(in)     :: maxsinl
    !!----   type (SuperSpaceGroup_Type) ,  optional,         intent(in)     :: SSG
    !!----   character(len=*),              optional,         intent(in)     :: powder
    !!----
    !!----    Calculate unique reflections between two values of
    !!----    sin_theta/lambda.  The output is not ordered.
    !!----
    !!---- Created: March - 2016
    !!

    Subroutine Gen_SReflections(Cell,sintlmax,Num_Ref,Reflex,nk,nharm,kv,maxsinl,SSG,powder)
       !---- Arguments ----!
       type (Crystal_Cell_Type),                        intent(in)     :: Cell
       real(kind=cp),                                   intent(in)     :: sintlmax
       integer,                                         intent(out)    :: num_ref
       class (sReflect_Type), dimension(:), allocatable,intent(out)    :: reflex
       integer,                       optional,         intent(in)     :: nk
       integer,       dimension(:),   optional,         intent(in)     :: nharm
       real(kind=cp), dimension(:,:), optional,         intent(in)     :: kv
       real(kind=cp), dimension(:),   optional,         intent(in)     :: maxsinl
       type (SuperSpaceGroup_Type) ,  optional,         intent(in)     :: SSG
       character(len=*),              optional,         intent(in)     :: powder

       !---- Local variables ----!
       real(kind=cp)         :: sval !,vmin,vmax
       integer               :: h,k,l,hmin,kmin,lmin,hmax,kmax,lmax, maxref,i,j,indp,indj, &
                                maxpos, mp, iprev,Dd, nf, ia
       integer,      dimension(:),   allocatable :: hh,kk,nulo
       integer,      dimension(:,:), allocatable :: hkl,hklm
       integer,      dimension(:),   allocatable :: indx,ini,fin,itreat
       real(kind=cp),dimension(:),   allocatable :: max_s
       real(kind=cp),dimension(:),   allocatable :: sv,sm
       logical :: kvect

       Dd=3
       kvect=present(nk) .and. present(nharm) .and. present(kv)
       if(kvect) Dd=3+nk             !total dimension of the reciprocal space
       hmax=nint(Cell%cell(1)*2.0*sintlmax+1.0)
       kmax=nint(Cell%cell(2)*2.0*sintlmax+1.0)
       lmax=nint(Cell%cell(3)*2.0*sintlmax+1.0)
       hmin=-hmax; kmin=-kmax; lmin= -lmax
       maxref= (2*hmax+1)*(2*kmax+1)*(2*lmax+1)
       if(kvect) then
          do k=1,nk
            	 maxref=maxref*2*nharm(k)
          end do
          allocate(max_s(nk))
          if(present(maxsinl)) then
            max_s=maxsinl
          else
            max_s=sintlmax
          end if
       end if
       allocate(hkl(Dd,maxref),indx(maxref),sv(maxref),hh(Dd),kk(Dd),nulo(Dd))
       nulo=0

       num_ref=0
       !Generation of fundamental reflections
       ext_do: do h=hmin,hmax
          do k=kmin,kmax
             do l=lmin,lmax
                hh=0
                hh(1:3)=[h,k,l]
                sval=H_s(hh,Cell)
                if (sval > sintlmax) cycle
                if(present(SSG)) then
                   if (H_Lat_Absent(hh,SSG%Latt_trans,SSG%Num_Lat)) cycle
                end if
                num_ref=num_ref+1
                if(num_ref > maxref) then
                   num_ref=maxref
                   exit ext_do
                end if
                sv(num_ref)=sval
                hkl(:,num_ref)=hh
             end do
          end do
       end do ext_do

       !Generation of satellites
       if(kvect) then
       	 k=0
       	 do_ex: do
           k = k + 1
           if(k > nk) exit do_ex
            	nf=num_ref
      	      do i=1,nf
       	      	 hh=hkl(:,i)
       	      	 do ia=-1,1,2
       	   	 	     do j=1,nharm(k)
       	   	 	     	 hh(3+k)=ia*j
                     sval=H_s(hh,Cell,nk,kv)
                     if (sval > max_s(k)) cycle
                     if(present(SSG)) then
                        if (H_Lat_Absent(hh,SSG%Latt_trans,SSG%Num_Lat)) cycle
                     end if
                     num_ref=num_ref+1
                     if(num_ref > maxref) then
                        num_ref=maxref
                        exit do_ex
                     end if
                     sv(num_ref)=sval
                     hkl(:,num_ref)=hh
                   end do  !j
       	   	 	   end do !ia
       	      end do  !i
         end do do_ex
       end if

       call sort(sv,num_ref,indx)

       allocate(hklm(Dd,num_ref-1),sm(num_ref-1),ini(num_ref-1),fin(num_ref-1),itreat(num_ref-1))
       do i=2,num_ref  !Eliminate the reflection 0 0 0 0 0 ...
         j=indx(i)
         hklm(:,i-1)=hkl(:,j)
         sm(i-1)=sv(j)
       end do
       num_ref=num_ref-1
       deallocate(hkl,sv,indx)
       if(present(SSG) .and. present(powder)) then
         itreat=0; ini=0; fin=0
         indp=0
         do i=1,num_ref              !Loop over all reflections
           !write(*,"(i6,3i5,i8)") i, hklm(:,i),itreat(i)
           if(itreat(i) == 0) then   !If not yet treated do the following
             hh(:)=hklm(:,i)
             indp=indp+1  !update the number of independent reflections
             itreat(i)=i  !Make this reflection treated
             ini(indp)=i  !put pointers for initial and final equivalent reflections
             fin(indp)=i
             do j=i+1,num_ref  !look for equivalent reflections to the current (i) in the list
                 if(abs(sm(i)-sm(j)) > 0.000001) exit
                 kk=hklm(:,j)
                 if(h_equiv(hh,kk,SSG)) then     ! if  hh eqv kk
                   itreat(j) = i                 ! add kk to the list equivalent to i
                   fin(indp)=j
                 end if
             end do
           end if !itreat
         end do

         !Selection of the most convenient independent reflections
         allocate(hkl(Dd,indp),sv(indp),indx(indp))
         indx=2 !nuclear and magnetic contribution by default
         do i=1,indp
           maxpos=0
           indj=ini(i)
           iprev=itreat(indj)
           do j=ini(i),fin(i)
             if(iprev /= itreat(j)) cycle
             hh=hklm(:,j)
             mp=count(hh > 0)
             if(mp > maxpos) then
               indj=j
               maxpos=mp
             end if
           end do !j
           hkl(:,i)=hklm(:,indj)
           if(hkl(1,i) < 0) hkl(:,i)=-hkl(:,i)
           sv(i)=sm(indj)
         end do
         !Now apply systematic absences other than lattice type
         num_ref=0
         do i=1,indp
           hh=hkl(:,i)
           if(H_Absent_SSG(hh,SSG)) then
             if(mH_Absent_SSG(hh,SSG)) then
               cycle
             else
               indx(i)=1   !pure magnetic
             end if
           else
             if(mH_Absent_SSG(hh,SSG)) indx(i)=0  !pure nuclear
           end if
           num_ref=num_ref+1
           hklm(:,num_ref)=hh
           sm(num_ref) = sv(i)
         end do
       end if  !SSG and Powder
       !Final assignments
       if(allocated(reflex)) deallocate(reflex)
       allocate(reflex(num_ref))
       do i=1,num_ref
         reflex(i)%h    = hklm(:,i)
         reflex(i)%s    = sm(i)
         reflex(i)%mult = h_mult(hh,SSG,.true.)
         reflex(i)%imag = indx(i)
       end do
       return
    End Subroutine Gen_SReflections


  End Module CFML_SuperSpaceGroups
