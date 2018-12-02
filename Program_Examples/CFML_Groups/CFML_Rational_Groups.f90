  Module CFML_Rational_Groups
    Use CFML_GlobalDeps,   only : cp
    use CFML_Math_General, only : sort, Set_Epsg
    use CFML_Rational_Arithmetic
    use CFML_String_Utilities, only : pack_string, Get_Separator_Pos

    implicit none
    Private
    Public  :: Allocate_Operator, Allocate_Operators, Init_Group, Get_Operators_From_String
    Public  :: Group_Constructor, Get_Group_From_Generators, is_Lattice_vec, Get_dimension, &
               Get_Mat_From_Symb_Op, Get_Symb_Op_from_Mat, Reorder_Operators,  &
               print_group,Get_SubGroups,Allocate_Group,Get_Multiplication_Table, &
               Get_subgroups_from_Table

    integer, private :: maxnum_op=2048
    logical,           public :: Err_group
    character(len=256),public :: Err_group_mess

   !-------Type declarations

    Type, public :: Symm_Oper_Type
       integer       :: time_inv=1
       integer       :: dt=1  !determinant of the submatrix (1:3,1:3), it should be 1 or -1
       type(rational), allocatable, dimension(:,:) :: Mat
    End Type Symm_Oper_Type

    Type, public :: Group_type
      integer :: Multip
      integer :: d  !Dimension of operator matrices (for 2D, d=3, for 3D, d=4, etc.)
      type(Symm_Oper_Type),dimension(:),allocatable:: Op
      Character(len=80),   dimension(:),allocatable:: Symb_Op !":" doesn't work properly in gfortran
    End Type Group_type

    Type, extends(Group_type), public :: Spg_Type
      integer :: Numops
      integer :: centred !=0 Centric(-1 no at origin), =1 Acentric, =2 Centric(-1 at origin)
      integer :: mag_type
      integer :: num_lat
      integer :: num_alat
      type(rational),dimension(:),   allocatable :: centre_coord
      type(rational),dimension(:,:), allocatable :: Lat_tr
      type(rational),dimension(:,:), allocatable :: aLat_tr
    End Type Spg_Type

   !------ Interfaces declarations (overloaded procedures)
    public :: operator (*)
    interface operator (*)
      module procedure multiply_Symm_Oper
    end interface
    private :: multiply_Symm_Oper

    public :: operator (==)
    interface operator (==)
      module procedure equal_Symm_Oper
      module procedure equal_Group
    end interface
    private :: equal_group,equal_Symm_Oper

    interface Group_Constructor
      module procedure Group_Constructor_gen
      module procedure Group_Constructor_string
    end interface Group_Constructor
    Private :: Group_Constructor_gen,Group_Constructor_string

    !Private symbols for parsing operators in both senses: Symbol <-> Operator
    character(len=*),dimension(10),parameter :: xyz=(/"x","y","z","t","u","v","w","p","q","r"/)
    character(len=*),dimension(10),parameter :: x1x2x3=(/"x1 ","x2 ","x3 ","x4 ","x5 ","x6 ","x7 ","x8 ","x9 ","x10"/)
    character(len=*),dimension(10),parameter :: abc=(/"a","b","c","d","e","f","g","h","i","j"/)

    Contains

    Pure recursive Function rdet(a) result(acc)
      type(rational), dimension(:,:), intent(in) :: a
      type(rational) :: acc
      !Local variables
      type(rational), dimension(size(a,dim=1)-1, size(a,dim=1)-1) :: b
      type(rational) :: sgn
      integer :: i, n
      n=size(a,dim=1)
      if (n == 1) then
        acc = a(1,1)
      else
        acc = 0_ik//1_ik
        sgn = 1_ik/1_ik
        do i=1,n
          b(:, :(i-1)) = a(2:, :i-1)
          b(:, i:) = a(2:, i+1:)
          acc = acc + sgn * a(1, i) * rdet(b)
          sgn = sgn * (-1_ik/1_ik)
        end do
      end if
    End Function rdet

    subroutine lu(a,p)
      !   in situ decomposition, corresponds to LAPACK's dgebtrf
      real(8), intent(in out) :: a(:,:)
      integer, intent(out  )  :: p(:)
      integer                 :: n,i,j,k,kmax
      n = size(a,1)
      p = [ ( i, i=1,n ) ]
      do k = 1,n-1
          kmax = maxloc(abs(a(p(k:),k)),1) + k-1
          if (kmax /= k ) p([k, kmax]) = p([kmax, k])
          a(p(k+1:),k) = a(p(k+1:),k) / a(p(k),k)
          forall (j=k+1:n) a(p(k+1:),j) = a(p(k+1:),j) - a(p(k+1:),k) * a(p(k),j)
      end do
    end subroutine

    Pure function multiply_Symm_Oper(Op1,Op2) result (Op3)
      type(Symm_Oper_Type), intent(in) :: Op1,Op2
      type(Symm_Oper_Type)             :: Op3
      integer :: n,d,i
      n=size(Op1%Mat,dim=1)
      allocate(Op3%Mat(n,n))
      d=n-1
      Op3%Mat=matmul(Op1%Mat,Op2%Mat) !automatic allocation in f2003
      Op3%Mat(1:d,n)=mod(Op3%Mat(1:d,n),1_ik)
       do i=1,d
       	 do
            if(Op3%Mat(i,n) < 0_ik//1_ik) then
            	Op3%Mat(i,n) = Op3%Mat(i,n) + 1_ik
            else
            	exit
            end if
         end do
      end do
      Op3%time_inv=Op1%time_inv*Op2%time_inv
      Op3%dt=Op1%dt*Op2%dt
    End Function multiply_Symm_Oper

    Pure function equal_Symm_Oper(Op1,Op2) result (info)
      type(Symm_Oper_Type), intent(in) :: Op1,Op2
      logical                          :: info
      info=.false.
      if(Op1%time_inv == Op2%time_inv) then
        if(equal_rational_matrix(Op1%Mat,Op2%Mat)) info=.true.
      end if
    end function equal_Symm_Oper

    Pure Function equal_Group(Gr1,Gr2) result (info)
      type(Spg_Type), intent(in) :: Gr1,Gr2
      logical                    :: info
      integer :: i,j
      info=.false.
      if(Gr1%multip /= Gr2%multip) return
      do i=2,Gr1%multip
        do j=2,Gr2%multip
           if(Gr1%Op(i) == Gr2%Op(j)) then
             info=.true.
             exit
           end if
        end do
        if(.not. info) return
      end do
    End Function equal_Group

    Pure Function is_Lattice_vec(V,Ltr,nlat) Result(Lattice_Transl)
       !---- Argument ----!
       type(rational), dimension(:),   intent( in) :: v
       type(rational), dimension(:,:), intent( in) :: Ltr
       integer,                        intent( in) :: nlat
       logical                                     :: Lattice_Transl

       !---- Local variables ----!
       type(rational), dimension(size(v)) :: vec
       integer                            :: i

       Lattice_Transl=.false.

       if (IsInteger(v)) then       ! if v is an integral vector =>  v is a lattice vector
          Lattice_Transl=.true.
       else                       ! if not look for lattice type
          do i=1,nlat
            vec=Ltr(:,i)-v
            if (IsInteger(vec)) then
              Lattice_Transl=.true.
              return
            end if
          end do
       end if
    End Function is_Lattice_vec

    Function Symbol_Operator(Op) result(symb)
      type(Symm_Oper_Type) :: Op
      Character(len=80)    :: symb
      call Get_Symb_Op_from_Mat(Op%Mat,Symb,"xyz",Op%time_inv)
    End Function Symbol_Operator

    Function Get_dimension(Symbol) result(d)
       character(len=*), intent (in) :: Symbol
       integer                       :: d
       integer, dimension(10) :: Pos
       integer                :: np,i,j
       call Get_Separator_Pos(symbol,",",Pos,np)
       !Verify if time reversal symbol is provided
       read(unit=symbol(Pos(np)+1:),fmt=*,iostat=i) j !reading an integer supposed to be invt
       if( i == 0 ) then
         d=np+1     !time_inv provided  (remmeber for 2D d=3, for 3D d=4, etc.)
       else
         d=np+2     !time_inv not-provided
       end if
    End Function Get_dimension

    Subroutine Init_Group(maxop,epsg)
      integer,       optional, intent(in) :: maxop
      real(kind=cp), optional, intent(in) :: epsg
      Err_group=.false.
      Err_group_mess=" "
      if(present(maxop)) maxnum_op=maxop
      if(present(epsg)) then
        call Set_Epsg(epsg)
      else
        call Set_Epsg(0.001)
      end if
    End Subroutine Init_Group

   Subroutine Allocate_Operator(d,Op)
       integer,                       intent(in)     :: d
       type(Symm_Oper_Type), intent(in out) :: Op
       integer :: i
       if(allocated(Op%Mat)) deallocate(Op%Mat)
       allocate(Op%Mat(d,d))
       !Inititalize to identity matrix
       Op%Mat=0_ik//1_ik
       do i=1,d
         Op%Mat(i,i)=1_ik
       end do
       Op%time_inv=1
       Op%dt=1
    End Subroutine Allocate_Operator

    !!---- Subroutine Allocate_Operators(d,multip,Op)
    !!----    integer,              intent(in)     :: d,multip
    !!----    type(Symm_Oper_Type), dimension(:), allocatable, intent(in out) :: Op
    !!----
    !!----  Multip is the expected maximum number of operators
    !!----
    Subroutine Allocate_Operators(d,multip,Op)
       integer,                                         intent(in)     :: d,multip
       type(Symm_Oper_Type), dimension(:), allocatable, intent(in out) :: Op
       integer :: i
       if(allocated(Op)) deallocate(Op)
       allocate(Op(multip))
       do i=1,multip
         call Allocate_Operator(d,Op(i))
       end do
    End Subroutine Allocate_Operators

    Subroutine Allocate_Group(d,multip,Gr)
       integer,        intent(in)     :: d,multip
       type(Spg_Type), intent(in out) :: Gr
       integer :: i
       if(allocated(Gr%Op)) deallocate(Gr%Op)
       allocate(Gr%Op(multip))
       Gr%d=d
       Gr%multip=multip
       do i=1,multip
         call Allocate_Operator(d,Gr%Op(i))
       end do
       if(allocated(Gr%Symb_Op)) deallocate(Gr%Symb_Op)
       allocate(Gr%Symb_Op(multip))
    End Subroutine Allocate_Group

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

    !!---- Subroutine Get_Multiplication_Table(Op,Table)
    !!----   type(Symm_Oper_Type), dimension(:), intent(in) :: Op
    !!----   integer, dimension(:,:),allocatable,intent(out):: Table
    !!----
    !!----   This subroutine construct the Cayley table of a Group
    !!----   defined by the rational operators Op. It is assumed that
    !!----   the first element is the identity.
    !!----
    Subroutine Get_Multiplication_Table(Op,Table)
      type(Symm_Oper_Type), dimension(:), intent(in) :: Op
      integer, dimension(:,:),allocatable,intent(out):: Table
      integer:: i,j,m,Multip
      type(Symm_Oper_Type) :: Opm
      multip=size(Op)
      allocate(Table(multip,multip))
      Table=0
      Table(1,:) = [(i,i=1,multip)] !It is supposed that the first operator is the identity
      Table(:,1) = [(i,i=1,multip)]
      do i=2,multip
        do j=2,multip
          Opm=Op(i)*Op(j)
          do m=1,multip
            if(Opm == Op(m)) then
               Table(i,j)=m
               exit
            end if
          end do
        end do
      end do
    End Subroutine Get_Multiplication_Table

    Pure Function Equal_sets(set1,set2) result(eq)    !To be tested
      integer, dimension(:),intent(in) :: set1,set2
      logical                          :: eq
      integer :: ord,i

      eq=.false.
      ord=size(set1)
      if(size(set2) /= ord) return
      eq=.true.
      do i=1,ord
        if( any(set2 == set1(i)) ) cycle
        eq=.false.
      end do
    End Function Equal_sets

    !!---- Pure Function included_in_set(set1,set2) result(eq)
    !!----   integer, dimension(:),intent(in) :: set1,set2
    !!----   logical                          :: eq
    !!----
    !!---- The value of the function is true if all the elements
    !!---- of set1 are present in set2
    !!----
    Pure Function included_in_set(set1,set2) result(eq)  !To be tested
      integer, dimension(:),intent(in) :: set1,set2
      logical                          :: eq
      integer :: ord,i

      eq=.false.
      ord=size(set1)
      if(size(set2) < ord) return
      eq=.true.
      do i=1,ord
        if( any(set2 == set1(i)) ) cycle
        eq=.false.
      end do
    End Function included_in_set

    Subroutine Get_group_from_Table(n,G,Table,ord)   !To be tested
      integer,                 intent(in)     :: n   !Number of initial elements in G
      integer, dimension(:),   intent(in out) :: G   !Vector containing the different elements of the groups
      integer, dimension(:,:), intent(in)     :: Table
      integer,                 intent(out)    :: ord !order or the final group
      integer :: i,j,k,m,nt,mult,ij
      logical, dimension(size(Table,1),size(Table,1)) :: done
      integer, dimension(:),allocatable :: ind
      mult=size(Table,1)
      done=.false.
      !done(1,:) = .true.
      !done(:,1) = .true.
      nt=n
      do_ext:do
        m=nt
        do i=1,m
          do_j:do j=1,m
            if(done(i,j)) cycle do_j
            ij=table(i,j)
            do k=1,nt
              if(ij == G(k)) then
                done(i,j)=.true.
                cycle do_j
              end if
            end do
            done(i,j)=.true.
            nt=nt+1
            G(nt)=ij
            if(nt > mult) then
              nt=nt-1
              exit do_ext
            end if
          end do do_j
        end do  !i
        if ( m == nt) exit do_ext
      end do do_ext
      ord=nt
      !Order the group by ascending order
      !allocate(ind(nt))
      !call sort(G(1:nt),nt,ind)
      !G(1:nt)=G(ind(1:nt))
    End Subroutine Get_group_from_Table

    Subroutine Get_subgroups_from_Table(Table,G,ord,ns)  !This is not currently working or it is too low
      integer, dimension(:,:),intent(in):: Table
      integer, dimension(:,:),allocatable, intent(out) :: G   !Vector containing the different elements of the groups
      integer, dimension(:),  allocatable, intent(out) :: ord !order or each group
      integer,                             intent(out) :: ns  !Number of subgroups
      integer:: i,j,k,n,m,Multip,L,mi,mj,mk,mij,mji,mik,mjk,mki,mkj,maxsub
      logical:: new
      multip=size(Table,1)
      !Estimated maximum number of subgroups
      maxsub=4096
      allocate(G(multip,maxsub),ord(maxsub))
      Ord=1
      G=0
      G(1,:) = 1
      L=0
      !One generator
      if(multip > 1) then
         do i=2,multip
           L=L+1
           G(:,L)=0
           G(1,L)=1
           G(2,L)=i
           !write(*,*) L,G(1:2,L)
           call Get_group_from_Table(2,G(:,L),table,ord(L))
           !Verify if this group already exists
           do n=1,L-1
             if(Equal_sets(G(1:Ord(n),n),G(1:Ord(L),L))) then
                L=L-1
                exit
             end if
           end do
             write(*,"(2(a,i4),a,192i4)") "Group:",L," Order:",Ord(L), " Elements:",G(1:Ord(L),L)
         end do
      else
        ns=1
        return
      end if
      ns=L
      !two generatora
      if(multip > 3) then
        do_ext:do i=2,multip-1
          do j=i+1,multip
             L=L+1
             if(L > maxsub) then
                L= maxsub
                exit do_ext
             end if
             G(:,L)=0
             G(1,L)=1
             G(2,L)=i
             G(3,L)=j
             call Get_group_from_Table(3,G(:,L),table,ord(L))
             !Verify if this group already exists
             do n=1,L-1
               if(Equal_sets(G(1:Ord(n),n),G(1:Ord(L),L))) then
                  L=L-1
                  exit
               end if
             end do
             write(*,"(2(a,i4),a,192i4)") "Group:",L," Order:",Ord(L), " Elements:",G(1:Ord(L),L)
          end do
        end do do_ext
      else
        ns=L
        return
      end if
      ns=L
      !if(multip > 3) then
      !  do_ext1:do i=2,multip-2  !three generatora
      !    do j=i+1,multip-1
      !      do k=j+1,multip
      !         L=L+1
      !         if(L > maxsub) then
      !            L= maxsub
      !            exit do_ext1
      !         end if
      !         G(:,L)=0
      !         G(1,L)=1
      !         G(2,L)=i
      !         G(3,L)=j
      !         G(4,L)=k
      !         !ord(L)=4   !{ 1, i, j, k }
      !         call Get_group_from_Table(4,G(:,L),table,ord(L))
      !         do n=1,L-1
      !           if(Equal_sets(G(1:Ord(n),n),G(1:Ord(L),L))) then
      !              L=L-1
      !              exit
      !           end if
      !         end do
      !       write(*,"(2(a,i4),a,192i4)") "Group:",L," Order:",Ord(L), " Elements:",G(1:Ord(L),L)
      !      end do
      !    end do
      !  end do do_ext1
      !end if
      !ns=L
    End Subroutine Get_subgroups_from_Table

    Subroutine Get_Group_From_Generators(ngen,Op,multip)
      integer,                               intent(in)     :: ngen
      type(Symm_Oper_Type), dimension(:),    intent(in out) :: Op
      integer,                               intent(out)    :: multip
      !--- Local variables ---!
      integer :: i,j,k,n,nt,max_op
      type(Symm_Oper_Type) :: Opt
      logical, dimension(size(Op),size(Op)) :: done
      integer, dimension(size(Op),size(Op)) :: tb

      max_op=size(Op)
      n=size(Op(1)%Mat,dim=1)
      call Allocate_Operator(n,Opt)
      done=.false.
      done(1,:) = .true.
      done(:,1) = .true.
      nt=ngen
      Err_group=.false.
      Err_group_mess=" "
      !Ensure that determinant of generators are calculated
      do i=1,ngen
        !Op(i)%dt=rational_determinant(Op(i)%Mat(1:3,1:3))
        Op(i)%dt=rdet(Op(i)%Mat(1:3,1:3))
      end do

      do_ext:do
        n=nt
        do i=1,n
          do_j:do j=1,n
            if(done(i,j)) cycle
            Opt=Op(i)*Op(j)
            do k=1,nt
              if(Opt == Op(k)) then
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
            Op(nt)=Opt
          end do do_j
        end do
        if ( n == nt) exit do_ext
      end do do_ext

      if(any(done(1:nt,1:nt) .eqv. .false. ) ) then
      	Err_group=.true.
      	Err_group_mess="Table of SSG operators not exhausted! Increase the expected order of the group!"
      end if
      if(nt == max_op) then
        Err_group=.true.
      	write(Err_group_mess,"(a,i5,a)") "Max_Order (",max_op,") reached! The provided generators may not form a group!"
      end if
      multip=nt
    End Subroutine Get_Group_From_Generators

    Subroutine Get_Symb_Op_from_Mat(Mat,Symb,x1x2x3_type,invt)
       !---- Arguments ----!
       type(rational),dimension(:,:), intent( in) :: Mat
       character (len=*),             intent(out) :: symb
       character(len=*), optional,    intent( in) :: x1x2x3_type
       integer, optional,             intent( in) :: invt

       !---- Local Variables ----!
       character(len=3),dimension(10)    :: x_typ
       character(len=15)                 :: car
       character(len=40)                 :: translation
       character(len=15),dimension(10)   :: sym
       integer                           :: i,j,Dd,d,k
       logical                           :: abc_type

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
         if(present(invt)) then
           write(unit=car,fmt="(i2)") invt !print_rational(Mat(Dd,Dd))
           car=adjustl(car)
           symb=trim(symb)//","//trim(car)
         end if
       end if
    End Subroutine Get_Symb_Op_from_Mat

    !!---- Subroutine Get_Mat_From_SSymSymb(Symb,Mat,invt)
    !!----   character(len=*),                intent(in)  :: Symb
    !!----   type(rational),dimension(:,:),   intent(out) :: Mat
    !!----   integer, optional                            :: invt
    !!----
    !!----  This subroutine provides the rational matrix Mat in standard
    !!----  form (all translation components positive) corresponding to the
    !!----  operator symbol Symb in Jone's faithful notation.
    !!----  The symbol is not modified but the matrix contain reduced translation.
    !!----  Some checking on the correctness of the symbol is performed
    !!----
    !!----  Created: June 2017 (JRC)
    !!----
    Subroutine Get_Mat_From_Symb_Op(Symb,Mat,invt)
      character(len=*),                intent(in)  :: Symb
      type(rational),dimension(:,:),   intent(out) :: Mat
      integer, optional                            :: invt
      !---- local variables ----!
      type(rational) :: det
      integer :: i,j,k,Dd, d, np,ns, n,m,inv,num,den,ind,ier
      character(len=3),dimension(10)                        :: x_typ
      character(len=len(Symb)), dimension(size(Mat,dim=1))  :: split
      character(len=len(Symb))                              :: string,pSymb,translation,transf
      character(len=10),        dimension(size(Mat,dim=1))  :: subst
      integer,                  dimension(size(Mat,dim=1)-1):: pos,pn
      logical                                               :: abc_transf
      !for checking
      character(len=40),dimension(size(Mat,dim=1),size(Mat,dim=1)) :: matrix
      character(len=5) :: forma
      err_group=.false.
      err_group_mess=" "
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
          Err_group=.true.
          Err_group_mess="Error in the origin coordinates"
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

      if(index(pSymb,"=") /= 0) then
          err_group=.true.
          err_group_mess="Error in the symbol of the operator: symbol '=' is forbidden!"
          return
      end if

      if(index(pSymb,";") /= 0) then
          err_group=.true.
          err_group_mess="Error in the symbol of the operator: symbol ';' is forbidden!"
          return
      end if

      if(index(pSymb,".") /= 0) then
          err_group=.true.
          err_group_mess="Error in the symbol of the operator: symbol '.' is forbidden!"
          return
      end if

      pos=0
      call Get_Separator_Pos(pSymb,",",pos,np)
      if(np /= d) then
      	if(np == d-1) then
      		do i=1,d
      			j=index(pSymb,trim(x_typ(i)))
      			if(j == 0) then !error in the symbol
              Err_group=.true.
              Err_group_mess="Error in the symbol of the operator: Missing ( "//trim(x_typ(i))//" )"
              return
      			end if
      		end do
      		pSymb=trim(pSymb)//",1"
      		np=np+1
      		pos(np)=len_trim(pSymb)-1
        else
          err_group=.true.
          write(err_group_mess,"(a,2i3,a)") "Error in the dimension of the symbol operator: "//trim(pSymb), np,d," for n_commas and dimension"
          return
        end if
      else !Check the presence of all symbols
      	do i=1,d
      		j=index(pSymb(1:pos(np)),trim(x_typ(i)))
      		if(j == 0) then !error in the symbol
            err_group=.true.
            err_group_mess="Error in the symbol of the operator: Missing ( "//trim(x_typ(i))//" )"
            return
      		end if
      	end do
      end if

      read(unit=pSymb(pos(np)+1:),fmt=*,iostat=ier) inv
      if(ier == 0) then
        Mat(Dd,Dd)=1//1
        if (present(invt)) invt = inv
      else
        Mat(Dd,Dd)=1//1
        if (present(invt)) invt = 1
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

      !Put the generator in standard form with positive translations
      call reduced_translation(Mat)

      !Final check that the determinant of the rotational matrix is integer
      !det=rational_determinant(Mat)
      det=rdet(Mat)
      if(det%denominator /= 1) then
         err_group=.true.
         err_group_mess="The determinant of the matrix is not integer! -> "//print_rational(det)
      end if
      if(det%numerator == 0) then
         err_group=.true.
         err_group_mess="The matrix of the operator is singular! -> det="//print_rational(det)
         matrix=print_rational(Mat)
         Write(*,*) "The matrix of the operator is singular!"
         forma="( a8)"
         write(forma(2:2),"(i1)") Dd
         do j=1,Dd
            write(unit=*,fmt=forma) (trim( Matrix(j,k))//" ",k=1,Dd)
         end do
      end if
    End Subroutine Get_Mat_From_Symb_Op

    subroutine sort_op(n, Op,cod)
        integer, intent(in) :: n
        type(Symm_Oper_Type) ,dimension(n), intent(in out):: Op
        character(len=*), intent(in) :: cod
        !
        type(Symm_Oper_Type) :: Ops
        integer :: i, j,opso
        integer, dimension(n) :: option
        if(cod == "tim") then
          option=Op(:)%time_inv
        else
          option=Op(:)%dt
        end if
        do i = 2, n
            Ops = Op(i)
            opso= Ops%dt
            if(cod == "tim") opso= Ops%time_inv
            j = i - 1
            do while (j >= 1)
                if (option(j) >= opso) exit
                Op(j + 1) = Op(j)
                option(j + 1) =  option(j)
                j = j - 1
            end do
            Op(j + 1) = Ops
            option(j+1) = opso
        end do
    end subroutine sort_op



    Subroutine Reorder_Operators(multip,Op, centred, centre_coord, Numops, num_lat, num_alat, Lat_tr, aLat_tr,mag_type)
      use CFML_Rational_Arithmetic, equal_matrix => equal_rational_matrix
      integer,                            intent(in)     :: multip
      type(Symm_Oper_Type), dimension(:), intent(in out) :: Op
      integer,                            intent(out)    :: num_lat,num_alat,Numops, centred, mag_type
      type(rational),dimension(:,:),      intent(out)    :: Lat_tr
      type(rational),dimension(:,:),      intent(out)    :: aLat_tr
      type(rational),dimension(:),        intent(out)    :: centre_coord
      !--- Local variables ---!
      integer :: i,j,L,n,m,Ng,invt,i_centre,d
      real(kind=cp), dimension(multip) :: tr   !Sum of absolute values of Translations components associated to the array of operators
      logical,       dimension(multip) :: nul  !Logical to control the exclusion of an operator from the independet list
      type(rational), dimension(:,:),allocatable:: identity,invers,imat !(d,d)
      type(Symm_Oper_Type), dimension(multip)  :: Opr,Op_Lat,Op_aLat
      type(Symm_Oper_Type)                     :: Op_aux,Op_aux1,Op_aux2, &
                                                  Op_centre,Op_identp
      real(kind=cp), parameter :: loc_eps=0.001
      real(kind=cp)            :: tmin
      type(rational)           :: ZERO, ONE, ONE_HALF

      ZERO=0_ik/1_ik;  ONE=1_ik/1_ik; ONE_HALF=1_ik/2_ik

      !Initializing
      n=size(Op(1)%Mat,1) !dimension of the full square matrices
      d=n-1               !here d is the dimension of the square matrix containing rotational operator
      allocate(identity(d,d),invers(d,d),imat(d,d))
      call Allocate_Operator(n,Op_identp)   ! {1|0}'
      identity=ZERO; nul=.false.; mag_type=1
      do i=1,d
        Op_identp%Mat(i,i)=ONE
        identity(i,i)=ONE
      end do
      Op_identp%time_inv=-1
      invers=-identity !Inversion
      centred=1 !Default value for non-centrosymmetric groups

      !Insertion sort putting the negative determinants at the bottom
      call sort_Op(multip,Op(1:multip),"det")
      do i=1,Multip
        tr(i)=sum(abs(Op(i)%Mat(1:d,n)))
        call Allocate_Operator(n,Opr(i))
      end do


      !Testing
      !write(*,*) " "
      !write(*,*) "List of operators after re-ordering by determinant: "
      !m=0
      !do i=1,Multip
      !  if(nul(i)) cycle
      !  m=m+1
      !  write(*,"(2i5,a30,f12.5,i6,tr2,L)") m,i,trim(Symbol_Operator(Op(i))),tr(i),Op(i)%dt,nul(i)
      !end do
      !end Testing


      !Check if the group is paramagnetic
      j=0
      do i=2,Multip
        if(Op(i)== Op_identp) then
          j=i
          exit
        end if
      end do
      if(j /= 0) Then
        do i=2,Multip   !Nullify all primed operators
          if(Op(i)%time_inv < 0) nul(i) = .true.
        end do
        mag_type=2
      end if
      !----End intial re-ordering

      !Look for centre of symmetry, and centring translations
      num_lat=0; num_alat=0; tmin=1.0e8; i_centre=0
      do j=2,Multip
         if(nul(j)) cycle
         invt= Op(j)%time_inv
         imat=Op(j)%Mat(1:d,1:d)
         if(equal_matrix(identity,imat) .and. invt == 1) then
            num_lat=num_lat+1
            Lat_tr(:,num_lat)=Op(j)%Mat(1:d,n)
            Op_Lat(num_lat)=Op(j)
            nul(j)=.true.   !Nullify centring translations
            cycle
         end if
         if(equal_matrix(imat,invers) .and. invt == 1 ) then
             nul(j) = .true.
             if(tr(j) < tmin) Then
               tmin= tr(j)
               i_centre=j
             end if
         end if
      end do

      if(i_centre /= 0) then
         Op_centre=Op(i_centre)
         centre_coord=ONE_HALF*Op(i_centre)%Mat(1:d,n)
         if(tr(i_centre) < loc_eps) then
           centred = 2
         else
           centred = 0
         end if
      end if

      do j=2,Multip-1   !Nullify operators deduced by lattice translations and centre of symmetry
         if(nul(j)) cycle

           if(mag_type ==2) then
             Op_aux=Op(j)*Op_identp
             do i=j+1,Multip
                if(nul(i)) cycle
                if(Op_aux == Op(i)) then
                   nul(i)=.true.
                end if
             end do
           end if

           if(num_lat > 0 .and. i_centre /= 0) then
              Op_aux=Op(j)*Op_centre
              do L=1,num_lat
                 Op_aux1=Op(j)*Op_Lat(L)
                 Op_aux2=Op_aux1*Op_centre
                 do i=j+1,Multip
                    if(nul(i)) cycle
                    if(Op_aux == Op(i) .or. Op_aux1 == Op(i) .or. Op_aux2 == Op(i)) then
                       nul(i)=.true.
                    end if
                 end do
              end do
              cycle
           end if

           if(num_lat > 0) then
              do L=1,num_lat
                 Op_aux=Op(j)*Op_Lat(L)
                 do i=j+1,Multip
                    if(nul(i)) cycle
                    if(Op_aux == Op(i)) then
                       nul(i)=.true.
                    end if
                 end do
              end do
           end if

           if(i_centre /= 0) then
              Op_aux=Op(j)*Op_centre
              do i=j+1,Multip
                 if(nul(i)) cycle
                 if(Op_aux == Op(i) ) then
                    nul(i)=.true.
                 end if
              end do
           end if
      end do

      !Determine the lattice anti-translations
      do_ext: do j=2,Multip
        if(nul(j)) cycle
        invt= Op(j)%time_inv
        imat=Op(j)%Mat(1:d,1:d)
        if(equal_matrix(identity,imat) .and. invt == -1) then
          num_alat=num_alat+1
          aLat_tr(:,num_alat)=Op(j)%Mat(1:d,n)
          Op_aLat(num_alat)=Op(j)
        end if
      end do  do_ext !j=2,Multip


      if(num_alat > 0) then
        if(mag_type /= 2) then
          mag_type=4
        end if
      else
        if(mag_type /= 2) then
          if(any(Op(:)%time_inv < 0)) mag_type=3
        end if
      end if

      ! => Determine the reduced set of symmetry operators"
      j=0
      do i=1,Multip
        if(nul(i)) cycle
        j=j+1
        Opr(j) = Op(i)
      end do
      Numops=j

      !Promote the reduced set of symmetry operators to the top of the list
      Op(1:j)=Opr(1:j)

      !Reorder the reduced set putting primed elements at the bottom
      call sort_op(Numops,Op(1:Numops),"tim")

      !Testing
      !write(*,*) " "
      !write(*,*) "List of reduced set of operators after re-ordering by time inversion: "
      !m=0
      !do i=1,Numops
      !  if(nul(i)) cycle
      !  m=m+1
      !  write(*,"(2i5,a30,f12.5,i6,tr2,L)") m,i,trim(Symbol_Operator(Op(i))),tr(i),Op(i)%time_inv,nul(i)
      !end do
      !end Testing

      if(i_centre /= 0) then
        m=j*2*(num_lat+1)
      else
        m=j*(num_lat+1)
      end if
      if(mag_type == 2) m=m*2
      if( m /= Multip) then !Check that it is OK
        write(unit=Err_group_mess,fmt="(2(a,i4))") " Warning! Multip=",Multip, " Calculated Multip: ",m
        Err_group=.true.
        return
      end if

      !Re-Construct, in an ordered way, all the symmetry operators
      !starting with the reduced set
      m=Numops
      ng=m
      if(i_centre /= 0) then   !First apply the centre of symmetry
        do i=1,Numops
          m=m+1
          Op(m) = Op(i) * Op_centre
        end do
      end if

      ng=m  ! Number or symmetry operators including centre of symmetry

      if(Num_Lat > 0) then  !Fourth apply the lattice centring translations
        do L=1,Num_Lat
           do i=1,ng
             m=m+1
             Op(m)=Op_Lat(L)*Op(i)
           end do
        end do
      end if

      if(mag_type == 2) then
         ng=m
         do i=1,ng
           m=m+1
           Op(m)=Op_identp*Op(i)
         end do
      end if

      !Normally here the number of operators should be equal to multiplicity
      !Test that everything is OK
      ng=m
      if(ng /= Multip) then
        Err_group=.true.
        write(unit=Err_group_mess,fmt="(2(a,i3))") " => Problem! the multiplicity ",Multip," has not been recovered, value of ng=",ng
      end if
    End Subroutine Reorder_Operators



    !!----  Subroutine Group_Constructor_ext(n,Opx,Grp)
    !!----     integer,                                    intent(in)      :: n
    !!----     type(Symm_Oper_Type), dimension(:),intent(in)      :: Opx
    !!----     type(rational_spg_Type),                    intent(in out)  :: Grp
    !!----
    !!----  This subroutine completes Grp win "n" additional generators
    !!----
    !!----
    !Subroutine Group_Constructor_ext(n,Opx,Grp)
    !   integer,                                    intent(in)      :: n
    !   type(Symm_Oper_Type), dimension(:),intent(in)      :: Opx
    !   type(Spg_Type),                    intent(in out)  :: Grp
    !   !--- Local variables ---!
    !   type(Symm_Oper_Type), dimension(:),allocatable :: Op
    !   integer :: d,i,j,ngen,invt,multip,centred,Numops,num_lat,num_alat,mag_type,n
    !   type(rational), dimension(:)  ,allocatable  :: centre_coord
    !   type(rational), dimension(:,:),allocatable  :: Lat_tr, aLat_tr
    !   type(rational), dimension(:,:),allocatable  :: Mat
    !   character(len=80) :: Symb_Op
    !
    !   ngen=Grp%Numops-1
    !   ! Get dimension of the generator matrices
    !   d=Grp%d
    !   allocate(Op(maxnum_op))
    !   do i=1,maxnum_op
    !     call Allocate_Operator(d,Op(i))
    !   end do
    !   allocate(Mat(d,d))
    !
    !   !Construct the list of the generators on top of Op. The identity is always the first operator
    !   do i=1,n
    !     Op(i+1)%Mat=Opx(i)%Mat
    !     Op(i+1)%time_inv=Opx(i)%time_inv
    !   end do
    !   ngen=ngen+1
    !
    !   include "CFML_group_constructor_template_inc.f90"
    !End Subroutine Group_Constructor_ext


    Subroutine Group_Constructor_gen(gen,Grp)
       character(len=*),dimension(:),intent(in)     :: gen
       type(Spg_Type),               intent(in out) :: Grp
       !--- Local variables ---!
       type(Symm_Oper_Type), dimension(:),  allocatable :: Op
       type(rational),       dimension(:),  allocatable :: centre_coord
       type(rational),       dimension(:,:),allocatable :: Lat_tr, aLat_tr
       type(rational),       dimension(:,:),allocatable :: Mat
       character(len=80) :: Symb_Op
       integer :: d,i,j,n,ngen,invt,multip,centred,Numops,num_lat,num_alat,mag_type

       ngen=size(gen)
       ! Get dimension of the generator matrices
       d=get_dimension(gen(1))
       include "CFML_group_constructor_template_inc.f90"
       !do n=1,Multip
       !   write(*,"(2(a,i3))") "  Operator #",n," Time inversion: ",Grp%Op(n)%time_inv
       !   do i=1,Grp%d
       !      write(unit=*,fmt="(a,10a)") "      [ ",(trim(print_rational(Grp%Op(n)%Mat(j,i)))//" ",j=1,Grp%d),"]"
       !   end do
       !end do
    End Subroutine Group_Constructor_gen

    Subroutine Group_Constructor_string(generatorList,Grp)
       character(len=*), intent(in)     :: generatorList
       type(Spg_Type),   intent(in out) :: Grp
       !--- Local variables ---!
       character(len=40),    dimension(:),  allocatable :: gen
       type(Symm_Oper_Type), dimension(:),  allocatable :: Op
       type(rational),       dimension(:),  allocatable :: centre_coord
       type(rational),       dimension(:,:),allocatable :: Lat_tr, aLat_tr
       type(rational),       dimension(:,:),allocatable :: Mat
       integer :: d,i,j,ngen,n,invt,multip,centred,Numops,num_lat,num_alat,mag_type
       character(len=80) :: Symb_Op

       allocate(gen(maxnum_op))
       call Get_Operators_From_String(generatorList,d,ngen,gen)
       include "CFML_group_constructor_template_inc.f90"
       !do n=1,Multip
       !   write(*,"(2(a,i3))") "  Operator #",n," Time inversion: ",Grp%Op(n)%time_inv
       !   do i=1,Grp%d
       !      write(unit=*,fmt="(a,10a)") "      [ ",(trim(print_rational(Grp%Op(n)%Mat(j,i)))//" ",j=1,Grp%d),"]"
       !   end do
       !end do
    End Subroutine Group_Constructor_string

    Subroutine Get_Operators_From_String(generatorList,d,ngen,gen)
       character(len=*),              intent(in)  :: generatorList
       integer,                       intent(out) :: d,ngen
       character(len=*),dimension(:), allocatable, intent(out) :: gen
       !--- Local variables ---!
       character(len=:), allocatable   :: symbol
       integer,dimension(10) :: Pos
       integer :: i,j,k,np
       logical :: timerev_provided

       !Determine the dimension from the generatorList (first operator)

       i=index(generatorList,";")
       if(i == 0) then !there is only one generator and ";" is not given
         Symbol=generatorList
       else
         Symbol=generatorList(1:i-1)
       end if

       d=Get_dimension(Symbol)

       call Get_Separator_Pos(generatorList,";",pos,np)

       !Verify if there is a final ";" not followed by a generator
       if(np == 0) then !there is only one generator and ";" is not given
         np=1
         Pos(np)=len_trim(generatorList)+2 !add artificially a position for ";"
       end if
       if(len_trim(generatorList(Pos(np)+1:)) == 0) then
         ngen=np  !final ;
       else
         ngen=np+1
       end if
      allocate(gen(nGen))
       j = 1
       do i = 1,np
         k = pos(i)
         gen(i) = generatorList(j:k-1)
         j = k + 1
       end do
       if(ngen > np) then
         gen(ngen) =generatorList(j:)
         i=index(gen(ngen),";")
         if(i /= 0) gen(ngen)(i:i) = " "
       end if

       timerev_provided=.false.
       do i=1,ngen
         call Get_Separator_Pos(gen(i),",",pos,np)
         if(np < d-1) cycle
         timerev_provided=.true.
       end do

       !Add time inversion in those operators that have not been Read
       if(timerev_provided) then
         do i=1,ngen
           call Get_Separator_Pos(gen(i),",",pos,np)
           if(np < d-1) gen(i)=trim(gen(i))//",1"
         end do
       end if

    End Subroutine Get_Operators_From_String

    !!---- Subroutine Get_SubGroups(SpG,SubG,nsg,point)
    !!---- !   !---- Arguments ----!
    !!----     Class (Group_Type) ,              intent( in) :: SpG
    !!----     Class (Group_Type) ,dimension(:), intent(out) :: SubG
    !!----     integer,                          intent(out) :: nsg
    !!----     logical, dimension(:,:), optional,intent(out) :: point
    !!----
    !!----   This subroutine provides all subgroups of a generic space group
    !!----   It is supposed that the operators have been already ordered, so
    !!----   that if the group is centrosymmetric the position of the corresponding
    !!----   operator is Numops+1, moreover if it is centred the symmetry operators
    !!----   are orderd as:
    !!----  [1...Numops] {-1|t}*[1...Numops] {1|t1}*[1...Numops] {-1|t}*[1 ... Numops]}
    !!----  {1|t2}*[1...Numops] {-1|t}*[1 ... Numops]} ....
    !!----  where t1,t2, ... are the centring translations
    !!----  This ordering facilitates the calculation of the major part of subgroups
    !!----  (excluding those for which less centring translations are considered)
    !!----
    Subroutine Get_SubGroups(SpG,SubG,nsg,point)
    !!   !---- Arguments ----!
        type(Spg_Type),                      intent( in) :: SpG
        type(Spg_Type),dimension(:),         intent(out) :: SubG
        integer,                             intent(out) :: nsg
        logical, dimension(:,:), optional,   intent(out) :: point
       !--- Local variables ---!
       integer  :: i,L,j,k,m,d, nc, mp,maxg,ng,kb, nla, i1,i2,nop,n,ns_1,ns_2,ns_3,n_nc_group
       logical  :: newg, cen_added
       character (len=40), dimension(:),allocatable :: gen
       character (len=40), dimension(30)            :: gen_lat
       character (len=40)                           :: gen_cent
       type(Symm_Oper_Type)                         :: Op_cent
       type(Symm_Oper_Type), dimension(30)          :: Op_lat


       maxg=size(SubG)
       allocate(gen(SpG%multip))
       d=SpG%d
       !---- Construct first the generators of centring translations ----!
       ng=0; nc=0
       nop=SpG%numops !number of symmetry operators excluding lattice centrings & centre of symmetry
       if (SpG%centred /= 1) then
          nop=nop*2 !!number of symmetry operators excluding lattice centrings
          nc=SpG%Numops+1  !Position of the centre of symmetry if it exist
          gen_cent=SpG%Symb_Op(nc)
          call Allocate_Operator(SpG%d,Op_cent)
          Op_cent=SpG%Op(nc)
       end if
       if(SpG%num_lat > 0) then
         do i=1,SpG%num_lat
            ng=ng+1
            gen_lat(ng)= SpG%Symb_Op(1+nop*i)
            Op_lat(ng)= SpG%Op(1+nop*i)
         end do
       end if
       !First work with the Numops operators to determine the subgroups, the other subgroups
       !will be obtained adding progressively the rest of generators (centre of symmetry and
       !lattice centrings.
       L=0
       !---- Determine first the groups with only one rotational generator
       ng=1
       do i=2,SpG%numops
          gen(1) = SpG%Symb_Op(i)
          L=L+1
          if (L > maxg) then
             nsg=maxg
             return
          end if
          newg=.true.
          !write(*,*) (trim(gen(j))//" ; ",j=1,ng)
          call Group_Constructor(gen(1:ng),SubG(L))
          do k=1,L-1
            if (SubG(L) == SubG(k)) then
               newg=.false.
               exit
            end if
          end do
          if (.not. newg) L=L-1
       end do
       ns_1=L
       !---- Determine now the groups with two rotational generators
       if(SpG%numops > 2) then
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
              do k=1,L-1
                if (SubG(L) == SubG(k)) then
                   newg=.false.
                   exit
                end if
              end do
              if (.not. newg) L=L-1
            end do
         end do
         ns_2=L-ns_1
       end if
       nsg=L
       n_nc_group=L
       !write(*,*) " Number of subgroups of the first Numops elements: ",n_nc_group
       !---- Determine now the new groups adding a centre of symmetry if it exists
       if (SpG%centred /= 1) then !This doubles the number of groups
         do i=1,n_nc_group
           L=L+1
           call Allocate_Group(SpG%d,2*SubG(i)%multip,SubG(L))
           if(allocated(SubG(L)%centre_coord)) deallocate(SubG(L)%centre_coord)
           allocate(SubG(L)%centre_coord(d-1))
           if(SubG(i)%num_alat /= 0) then
               if(allocated(SubG(L)%aLat_tr)) deallocate(SubG(L)%aLat_tr)
               allocate(SubG(L)%aLat_tr(d-1,SubG(i)%num_alat))
               SubG(L)%aLat_tr=SubG(i)%aLat_tr
               SubG(L)%num_alat=SubG(i)%num_alat
           end if
           k=SubG(i)%numops
           do j=1,SubG(i)%numops
             SubG(L)%Op(j)=SubG(i)%Op(j)
             SubG(L)%Symb_Op(j)=SubG(i)%Symb_Op(j)
             k=k+1
             SubG(L)%Op(k)=SubG(i)%Op(j)*Op_cent
             SubG(L)%Symb_Op(k)=Symbol_Operator(SubG(L)%Op(k))
           end do
           SubG(L)%Numops=SubG(i)%Numops
           SubG(L)%mag_type=SubG(i)%mag_type
           SubG(L)%centred=SpG%centred
           SubG(L)%centre_coord=SpG%centre_coord
         end do
       end if
       nsg=L
       n_nc_group=L
       !if (SpG%centred /= 1) write(*,*) " Number of subgroups of adding a centre of symmetry: ",n_nc_group

       !Determine now the rest of groups adding the lattice translations if they exist in the
       !original space group
       if(SpG%num_lat > 0) then
         Select Case (SpG%num_lat)
           case(1)
             do i=1,n_nc_group
               L=L+1
               call Allocate_Group(SpG%d,2*SubG(i)%multip,SubG(L))
               if(allocated(SubG(L)%centre_coord)) deallocate(SubG(L)%centre_coord)
               allocate(SubG(L)%centre_coord(d-1))
               if(allocated(SubG(L)%Lat_tr)) deallocate(SubG(L)%Lat_tr)
               allocate(SubG(L)%Lat_tr(d-1,1))
               SubG(L)%Lat_tr(:,1)=SpG%Lat_tr(:,1)
               if(SubG(i)%num_alat /= 0) then
                   if(allocated(SubG(L)%aLat_tr)) deallocate(SubG(L)%aLat_tr)
                   allocate(SubG(L)%aLat_tr(d-1,SubG(i)%num_alat))
                   SubG(L)%aLat_tr=SubG(i)%aLat_tr
                   SubG(L)%num_alat=SubG(i)%num_alat
               end if
               SubG(L)%num_lat=1
               if(SubG(i)%centred /= 1) then
                  k=SubG(i)%numops*2
               else
                  k=SubG(i)%numops
               end if
               kb=k
               do j=1,kb
                 SubG(L)%Op(j)=SubG(i)%Op(j)
                 SubG(L)%Symb_Op(j)=SubG(i)%Symb_Op(j)
                 k=k+1
                 SubG(L)%Op(k)=SubG(i)%Op(j)*Op_lat(1)
                 SubG(L)%Symb_Op(k)=Symbol_Operator(SubG(L)%Op(k))
               end do
               SubG(L)%Multip=2*SubG(i)%Multip
               SubG(L)%Numops=SubG(i)%Numops
               SubG(L)%mag_type=SubG(i)%mag_type
               SubG(L)%centred=SubG(i)%centred
               if(SubG(L)%centred /= 1) SubG(L)%centre_coord=SpG%centre_coord
             end do
           case(2)
             do i=1,n_nc_group
               L=L+1
               call Allocate_Group(SpG%d,3*SubG(i)%multip,SubG(L))
               if(allocated(SubG(L)%centre_coord)) deallocate(SubG(L)%centre_coord)
               allocate(SubG(L)%centre_coord(d-1))
               if(allocated(SubG(L)%Lat_tr)) deallocate(SubG(L)%Lat_tr)
               allocate(SubG(L)%Lat_tr(d-1,2))
               SubG(L)%Lat_tr(:,:)=SpG%Lat_tr(:,:)
               SubG(L)%num_lat=2
               if(SubG(i)%num_alat /= 0) then
                   if(allocated(SubG(L)%aLat_tr)) deallocate(SubG(L)%aLat_tr)
                   allocate(SubG(L)%aLat_tr(d-1,SubG(i)%num_alat))
                   SubG(L)%aLat_tr=SubG(i)%aLat_tr
                   SubG(L)%num_alat=SubG(i)%num_alat
               end if
               if(SubG(i)%centred /= 1) then
                  k=SubG(i)%numops*2
               else
                  k=SubG(i)%numops
               end if
               kb=k
               do j=1,kb
                 SubG(L)%Op(j)=SubG(i)%Op(j)
                 SubG(L)%Symb_Op(j)=SubG(i)%Symb_Op(j)
                 k=k+1
                 SubG(L)%Op(k)=SubG(i)%Op(j)*Op_lat(1)
                 SubG(L)%Symb_Op(k)=Symbol_Operator(SubG(L)%Op(k))
               end do
               do j=1,kb
                 k=k+1
                 SubG(L)%Op(k)=SubG(i)%Op(j)*Op_lat(2)
                 SubG(L)%Symb_Op(k)=Symbol_Operator(SubG(L)%Op(k))
               end do
               SubG(L)%Multip=3*SubG(i)%Multip
               SubG(L)%Numops=SubG(i)%Numops
               SubG(L)%mag_type=SubG(i)%mag_type
               SubG(L)%centred=SubG(i)%centred
               if(SubG(L)%centred /= 1) SubG(L)%centre_coord=SpG%centre_coord
             end do
           case(3)
             do i=1,n_nc_group
               L=L+1
               call Allocate_Group(SpG%d,4*SubG(i)%multip,SubG(L))
               if(allocated(SubG(L)%centre_coord)) deallocate(SubG(L)%centre_coord)
               allocate(SubG(L)%centre_coord(d-1))
               if(allocated(SubG(L)%Lat_tr)) deallocate(SubG(L)%Lat_tr)
               allocate(SubG(L)%Lat_tr(d-1,3))
               SubG(L)%Lat_tr(:,:)=SpG%Lat_tr(:,:)
               SubG(L)%num_lat=3
               if(SubG(i)%num_alat /= 0) then
                   if(allocated(SubG(L)%aLat_tr)) deallocate(SubG(L)%aLat_tr)
                   allocate(SubG(L)%aLat_tr(d-1,SubG(i)%num_alat))
                   SubG(L)%aLat_tr=SubG(i)%aLat_tr
                   SubG(L)%num_alat=SubG(i)%num_alat
               end if
               if(SubG(i)%centred /= 1) then
                  k=SubG(i)%numops*2
               else
                  k=SubG(i)%numops
               end if
               kb=k
               do j=1,kb
                 SubG(L)%Op(j)=SubG(i)%Op(j)
                 SubG(L)%Symb_Op(j)=SubG(i)%Symb_Op(j)
                 k=k+1
                 SubG(L)%Op(k)=SubG(i)%Op(j)*Op_lat(1)
                 SubG(L)%Symb_Op(k)=Symbol_Operator(SubG(L)%Op(k))
               end do
               do j=1,kb
                 k=k+1
                 SubG(L)%Op(k)=SubG(i)%Op(j)*Op_lat(2)
                 SubG(L)%Symb_Op(k)=Symbol_Operator(SubG(L)%Op(k))
               end do
               do j=1,kb
                 k=k+1
                 SubG(L)%Op(k)=SubG(i)%Op(j)*Op_lat(3)
                 SubG(L)%Symb_Op(k)=Symbol_Operator(SubG(L)%Op(k))
               end do
               SubG(L)%Multip=4*SubG(i)%Multip
               SubG(L)%Numops=SubG(i)%Numops
               SubG(L)%mag_type=SubG(i)%mag_type
               SubG(L)%centred=SubG(i)%centred
               if(SubG(L)%centred /= 1) SubG(L)%centre_coord=SpG%centre_coord
             end do
         End Select

       end if

       nsg=L
       if(present(point)) then
         point=.false.
         do j=1,nsg
           L=1
           do i=1,SpG%multip
              do k=L,SubG(j)%multip
               if(SubG(j)%Symb_Op(k) == SpG%Symb_Op(i)) then
                  point(i,j) = .true.
                  L=k+1
                  exit
               end if
              end do
           end do
         end do
       end if
       !include "CFML_subgroups_template_inc.f90"

    End Subroutine Get_SubGroups


    subroutine print_Group(Grp,lun)
      type(Spg_Type),    intent(in)   :: Grp
      integer, optional, intent(in)   :: lun
      integer :: iout,i,j
      iout=6 !To be replaced by Fortran environment value
      if(present(lun)) iout=lun
      write(unit=iout,fmt="(a)")        "    General Space Group"
      write(unit=iout,fmt="(a)")        "    -------------------"
      write(unit=iout,fmt="(a,i4)")     "     Op-Dimension: ",Grp%d
      write(unit=iout,fmt="(a,i4)")     "  Space-Dimension: ",Grp%d-1
      write(unit=iout,fmt="(a,i4)")     "     Multiplicity: ",Grp%multip
      write(unit=iout,fmt="(a,i4)")     "          MagType: ",Grp%mag_type
      write(unit=iout,fmt="(a,i4)")     "           NumOps: ",Grp%numops
      write(unit=iout,fmt="(a,i4)")     "          Centred: ",Grp%centred
      write(unit=iout,fmt="(a,i4)")     "          Num_Lat: ",Grp%num_lat
      write(unit=iout,fmt="(a,i4)")     "         Num_aLat: ",Grp%num_alat
      if(Grp%centred == 1) then
         write(unit=iout,fmt="(a)")     "     Centre_coord: none!"
      else
         !write(unit=iout,fmt="(a,10f8.3)") "     Centre_coord: ",Grp%centre_coord
         write(unit=iout,fmt="(a,10a)") "     Centre_coord: [ ",(trim(print_rational(Grp%centre_coord(i)))//" ",i=1,Grp%d-1),"]"
      end if
      if(Grp%num_lat > 0) then
        write(unit=iout,fmt="(/a)")      "  Centring translations:"
        do i=1,Grp%num_lat
           !write(unit=iout,fmt="(i3,tr4,10f8.3)") i,Grp%Lat_tr(:,i)
         write(unit=iout,fmt="(a,10a)") "      [ ",(trim(print_rational(Grp%Lat_tr(j,i)))//" ",j=1,Grp%d-1),"]"
        end do
      end if
      if(Grp%num_alat > 0) then
        write(unit=iout,fmt="(/a)")      "  Anti-translations:"
        do i=1,Grp%num_alat
          ! write(*,"(i3,tr4,10f8.3)") i,Grp%aLat_tr(:,i)
         write(unit=iout,fmt="(a,10a)") "      [ ",(trim(print_rational(Grp%aLat_tr(j,i)))//" ",j=1,Grp%d-1),"]"
        end do
      end if
      write(unit=iout,fmt="(/a)")      "  Complete list of symmetry operators:"
      do i=1,Grp%Multip
        write(unit=iout,fmt="(i5,a)") i,"  ->  "//trim(Grp%Symb_Op(i))
      end do
    end subroutine print_Group

  End Module CFML_Rational_Groups