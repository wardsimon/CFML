  Module CFML_Groups
    Use CFML_GlobalDeps,   only : cp,sp,dp,pi
    use CFML_Math_General, only : equal_matrix,equal_vector,Zbelong,modulo_lat, sort, Set_Epsg
    use CFML_Math_3D,      only : determ_A
    use CFML_Rational_Arithmetic
    use CFML_String_Utilities, only : pack_string, Get_Separator_Pos
    Use CFML_Crystallographic_Symmetry, only: Set_SpaceGroup, Space_Group_Type,Get_Transl_Symbol, &
                                              Get_SymSymb, NS_Space_Group_Type,Sym_Oper_Type,&
                                              applyso, Lattice_trans,Get_Trasfm_Symbol,SpGr_Equal

    implicit none
    Private
    Public  :: Allocate_Operator, Allocate_Operators, Init_Group, Get_Operators_From_String
    Public  :: Group_Constructor, Get_Group_From_Generators, is_Lattice_vec, Get_dimension, &
               Get_Mat_From_Symb_Op, Get_Symb_Op_from_Mat, Reorder_Operators,  &
               Get_SubGroups

    integer, private :: maxnum_op=2048
    logical,           public :: Err_group
    character(len=256),public :: Err_group_mess

   !-------Type declarations

    Type, public :: Symm_Oper_Type
       integer       :: time_inv=1
       integer       :: dt=1  !determinant of the submatrix (1:3,1:3), it should be 1 or -1
       real(kind=cp), allocatable, dimension(:,:) :: Mat
    End Type Symm_Oper_Type

    Type, public :: rational_Symm_Oper_Type
       integer       :: time_inv=1
       integer       :: dt=1  !determinant of the submatrix (1:3,1:3), it should be 1 or -1
       type(rational), allocatable, dimension(:,:) :: Mat
    End Type rational_Symm_Oper_Type

    Type, public :: Group_type
      integer :: Multip
      integer :: d  !Dimension of operator matrices (for 2D, d=3, for 3D, d=4, etc.)
      type(Symm_Oper_Type),dimension(:),allocatable:: Op
      Character(len=80),   dimension(:),allocatable:: Symb_Op !It doesn't work ":" in gfortran
    End Type Group_type

    Type, extends(Group_type),public :: rational_Group_type
      type(rational_Symm_Oper_Type),dimension(:),allocatable:: rOp
    End Type rational_Group_type

    Type, extends(Group_type), public :: raw_spg_type
      integer :: Numops
      integer :: centred !=0 Centric(-1 no at origin), =1 Acentric, =2 Centric(-1 at origin)
      integer :: mag_type
      integer :: num_lat
      integer :: num_alat
      real(kind=cp),dimension(:),   allocatable :: centre_coord
      real(kind=cp),dimension(:,:), allocatable :: Lat_tr
      real(kind=cp),dimension(:,:), allocatable :: aLat_tr
    End Type raw_spg_type

    Type, extends(rational_Group_type), public :: rational_spg_type
      integer :: Numops
      integer :: centred !=0 Centric(-1 no at origin), =1 Acentric, =2 Centric(-1 at origin)
      integer :: mag_type
      integer :: num_lat
      integer :: num_alat
      type(rational),dimension(:),   allocatable :: centre_coord
      type(rational),dimension(:,:), allocatable :: Lat_tr
      type(rational),dimension(:,:), allocatable :: aLat_tr
    End Type rational_spg_type

   !------ Interfaces declarations (overloaded procedures)
    public :: operator (*)
    interface operator (*)
      module procedure multiply_Symm_Oper
      module procedure multiply_Symm_Oper_rational
    end interface
    private :: multiply_Symm_Oper

    public :: operator (==)
    interface operator (==)
      module procedure equal_Symm_Oper
      module procedure equal_rational_Symm_Oper
      module procedure equal_Group
    end interface
    private :: equal_Symm_Oper,equal_Group

    interface Get_Group_From_Generators
      module procedure Get_Group_From_Generators_real
      module procedure Get_Group_From_Generators_rational
    end interface Get_Group_From_Generators
    private :: Get_Group_From_Generators_real, Get_Group_From_Generators_rational

    interface Group_Constructor
      module procedure Group_Constructor_gen
      module procedure Group_Constructor_string
      module procedure Group_Constructor_rational
    end interface Group_Constructor
    Private :: Group_Constructor_gen,Group_Constructor_string, &
               Group_Constructor_rational

    interface Allocate_Operator
      module procedure Allocate_Operator_real
      module procedure Allocate_Operator_rational
    end interface Allocate_Operator
    private :: Allocate_Operator_real, Allocate_Operator_rational

    interface Allocate_Operators
      module procedure Allocate_Operators_real
      module procedure Allocate_Operators_rational
    end interface Allocate_Operators
    private :: Allocate_Operators_real, Allocate_Operators_rational

    interface Symbol_Operator
      module procedure Symbol_Operator_rational
      module procedure Symbol_Operator_real
    end interface Symbol_Operator
    private :: Symbol_Operator_real, Symbol_Operator_rational

    interface sort_op
      module procedure sort_op_real
      module procedure sort_op_rational
    end interface sort_op

    interface Reorder_Operators
      module procedure Reorder_Operators_real
      module procedure Reorder_Operators_rational
    end interface Reorder_Operators


    !Private symbols for parsing operators in both senses: Symbol <-> Operator
    character(len=*),dimension(10),parameter :: xyz=(/"x","y","z","t","u","v","w","p","q","r"/)
    character(len=*),dimension(10),parameter :: x1x2x3=(/"x1 ","x2 ","x3 ","x4 ","x5 ","x6 ","x7 ","x8 ","x9 ","x10"/)
    character(len=*),dimension(10),parameter :: abc=(/"a","b","c","d","e","f","g","h","i","j"/)

    Contains

    Pure function multiply_Symm_Oper(Op1,Op2) result (Op3)
      type(Symm_Oper_Type), intent(in) :: Op1,Op2
      type(Symm_Oper_Type)             :: Op3
      integer :: n,d,i
      n=size(Op1%Mat,dim=1)
      allocate(Op3%Mat(n,n))
      d=n-1
      Op3%Mat=matmul(Op1%Mat,Op2%Mat) !automatic allocation in f2003
      do i=1,d     !Put translations in the interval [0.0,1.0)
         if(Op3%Mat(i,n) < 0.0_cp) then
         	  Op3%Mat(i,n) = Op3%Mat(i,n) + 1.0_cp
         else if(Op3%Mat(i,n) >= 1.0_cp) then
         	  Op3%Mat(i,n) = Op3%Mat(i,n) - 1.0_cp
         end if
      end do
      Op3%time_inv=Op1%time_inv*Op2%time_inv
      Op3%dt=Op1%dt*Op2%dt
    End Function multiply_Symm_Oper

    Pure function multiply_Symm_Oper_rational(Op1,Op2) result (Op3)
      type(rational_Symm_Oper_Type), intent(in) :: Op1,Op2
      type(rational_Symm_Oper_Type)             :: Op3
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
    End Function multiply_Symm_Oper_rational

    Pure function equal_Symm_Oper(Op1,Op2) result (info)
      type(Symm_Oper_Type), intent(in) :: Op1,Op2
      logical                          :: info
      info=.false.
      if(Op1%time_inv == Op2%time_inv) then
        if(equal_matrix(Op1%Mat,Op2%Mat)) info=.true.
      end if
    end function equal_Symm_Oper

    Pure function equal_rational_Symm_Oper(Op1,Op2) result (info)
      type(rational_Symm_Oper_Type), intent(in) :: Op1,Op2
      logical                                   :: info
      info=.false.
      if(Op1%time_inv == Op2%time_inv) then
        if(equal_rational_matrix(Op1%Mat,Op2%Mat)) info=.true.
      end if
    end function equal_rational_Symm_Oper

    Pure Function equal_Group(Gr1,Gr2) result (info)
      class(Group_Type), intent(in) :: Gr1,Gr2
      logical                       :: info
      integer :: i,j
      info=.false.
      if(Gr1%multip /= Gr2%multip) return
      do i=2,Gr1%multip
        do j=2,Gr2%multip
           if(equal_Symm_Oper(Gr1%Op(i),Gr2%Op(j))) then
             info=.true.
             exit
           end if
        end do
        if(.not. info) return
      end do
    End Function equal_Group

    Pure Function is_Lattice_vec(V,Ltr,nlat) Result(Lattice_Transl)
       !---- Argument ----!
       real(kind=cp), dimension(:),   intent( in) :: v
       real(kind=cp), dimension(:,:), intent( in) :: Ltr
       integer,                       intent( in) :: nlat
       logical                                    :: Lattice_Transl

       !---- Local variables ----!
       real(kind=cp)   , dimension(size(v)) :: vec
       integer                              :: i

       Lattice_Transl=.false.

       if (Zbelong(v)) then       ! if v is an integral vector =>  v is a lattice vector
          Lattice_Transl=.true.
       else                       ! if not look for lattice type
          do i=1,nlat
            vec=Ltr(:,i)-v
            if (Zbelong(vec)) then
              Lattice_Transl=.true.
              return
            end if
          end do
       end if
    End Function is_Lattice_vec

    Function Symbol_Operator_real(Op) result(symb)
      type(Symm_Oper_Type) :: Op
      Character(len=80)    :: symb
      type(rational),dimension(:,:),allocatable :: Mat
      integer :: n
      n=size(Op%Mat,1)
      allocate (Mat(n,n))
      Mat=Op%Mat
      call Get_Symb_Op_from_Mat(Mat,Symb,"xyz",Op%time_inv)
    End Function Symbol_Operator_real

    Function Symbol_Operator_rational(Op) result(symb)
      type(rational_Symm_Oper_Type) :: Op
      Character(len=80)    :: symb
      call Get_Symb_Op_from_Mat(Op%Mat,Symb,"xyz",Op%time_inv)
    End Function Symbol_Operator_rational

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
         d=np+2    !time_inv not-provided
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

    Subroutine Allocate_Operator_real(d,Op)
       integer,              intent(in)     :: d
       type(Symm_Oper_Type), intent(in out) :: Op
       integer :: i
       if(allocated(Op%Mat)) deallocate(Op%Mat)
       allocate(Op%Mat(d,d))
       !Inititalize to identity matrix
       Op%Mat=0.0
       do i=1,d
         Op%Mat(i,i)=1.0
       end do
       Op%time_inv=1
       Op%dt=1
    End Subroutine Allocate_Operator_real

    Subroutine Allocate_Operator_rational(d,Op)
       integer,                       intent(in)     :: d
       type(rational_Symm_Oper_Type), intent(in out) :: Op
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
    End Subroutine Allocate_Operator_rational

    !!---- Subroutine Allocate_Operators(d,multip,Op)
    !!----    integer,              intent(in)     :: d,multip
    !!----    type(Symm_Oper_Type), dimension(:), allocatable, intent(in out) :: Op
    !!----
    !!----  Multip is the expected maximum number of operators
    !!----
    Subroutine Allocate_Operators_real(d,multip,Op)
       integer,              intent(in)     :: d,multip
       type(Symm_Oper_Type), dimension(:), allocatable, intent(in out) :: Op
       integer :: i
       if(allocated(Op)) deallocate(Op)
       allocate(Op(multip))
       do i=1,multip
         call Allocate_Operator(d,Op(i))
       end do
    End Subroutine Allocate_Operators_real

    Subroutine Allocate_Operators_rational(d,multip,Op)
       integer,                                                  intent(in)     :: d,multip
       type(rational_Symm_Oper_Type), dimension(:), allocatable, intent(in out) :: Op
       integer :: i
       if(allocated(Op)) deallocate(Op)
       allocate(Op(multip))
       do i=1,multip
         call Allocate_Operator(d,Op(i))
       end do
    End Subroutine Allocate_Operators_rational

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

    Subroutine Get_Group_From_Generators_real(ngen,Op,multip,table)
      integer,                                        intent(in)     :: ngen
      type(Symm_Oper_Type), dimension(:),             intent(in out) :: Op
      integer,                                        intent(out)    :: multip
      integer, dimension(:,:), allocatable, optional, intent(out)    :: table
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
      tb(1,:) = [(i,i=1,max_op)] !It is supposed that the first generator is the identity
      tb(:,1) = [(i,i=1,max_op)]
      nt=ngen
      Err_group=.false.
      Err_group_mess=" "
      !Ensure that determinant of generators are calculated
      do i=1,ngen
        Op(i)%dt=nint(determ_A(Op(i)%Mat(1:3,1:3)))
      end do

      do_ext:do
        n=nt
        do i=1,n
          do_j:do j=1,n
            if(done(i,j)) cycle
            Opt=Op(i)*Op(j)
            do k=1,nt
              if(Opt == Op(k)) then
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
      	Err_group=.true.
      	Err_group_mess="Table of SSG operators not exhausted! Increase the expected order of the group!"
      end if
      if(nt == max_op) then
        Err_group=.true.
      	write(Err_group_mess,"(a,i5,a)") "Max_Order (",max_op,") reached! The provided generators may not form a group!"
      end if

      multip=nt
      if(present(table)) then
        allocate(Table(multip,multip))
        Table=tb
      end if
    End Subroutine Get_Group_From_Generators_real

    Subroutine Get_Group_From_Generators_rational(ngen,Op,multip,table)
      integer,                                        intent(in)     :: ngen
      type(rational_Symm_Oper_Type), dimension(:),    intent(in out) :: Op
      integer,                                        intent(out)    :: multip
      integer, dimension(:,:), allocatable, optional, intent(out)    :: table
      !--- Local variables ---!
      integer :: i,j,k,n,nt,max_op
      type(rational_Symm_Oper_Type) :: Opt
      logical, dimension(size(Op),size(Op)) :: done
      integer, dimension(size(Op),size(Op)) :: tb

      max_op=size(Op)
      n=size(Op(1)%Mat,dim=1)
      call Allocate_Operator(n,Opt)
      done=.false.
      done(1,:) = .true.
      done(:,1) = .true.
      tb(1,:) = [(i,i=1,max_op)] !It is supposed that the first generator is the identity
      tb(:,1) = [(i,i=1,max_op)]
      nt=ngen
      Err_group=.false.
      Err_group_mess=" "
      !Ensure that determinant of generators are calculated
      do i=1,ngen
        Op(i)%dt=rational_determinant(Op(i)%Mat(1:3,1:3))
      end do

      do_ext:do
        n=nt
        do i=1,n
          do_j:do j=1,n
            if(done(i,j)) cycle
            Opt=Op(i)*Op(j)
            do k=1,nt
              if(Opt == Op(k)) then
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
      	Err_group=.true.
      	Err_group_mess="Table of SSG operators not exhausted! Increase the expected order of the group!"
      end if
      if(nt == max_op) then
        Err_group=.true.
      	write(Err_group_mess,"(a,i5,a)") "Max_Order (",max_op,") reached! The provided generators may not form a group!"
      end if

      multip=nt
      if(present(table)) then
        allocate(Table(multip,multip))
        Table=tb
      end if
    End Subroutine Get_Group_From_Generators_rational

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
    !!----  operator symbol Symb in Jone's faithfull notation.
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
      det=rational_determinant(Mat)
      if(det%denominator /= 1) then
         err_group=.true.
         err_group_mess="The determinant of the matrix is not integer! -> "//print_rational(det)
      end if
      if(det%numerator == 0) then
         err_group=.true.
         err_group_mess="The matrix of the operator is singular! -> det="//print_rational(det)
      end if
    End Subroutine Get_Mat_From_Symb_Op

    subroutine sort_op_real(n, Op,cod)
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
    end subroutine sort_op_real

    subroutine sort_op_rational(n, Op,cod)
        integer, intent(in) :: n
        type(rational_Symm_Oper_Type) ,dimension(n), intent(in out):: Op
        character(len=*), intent(in) :: cod
        !
        type(rational_Symm_Oper_Type) :: Ops
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
    end subroutine sort_op_rational

    Subroutine Reorder_Operators_real(multip,Op, centred, centre_coord, Numops, num_lat, num_alat, Lat_tr, aLat_tr,mag_type)
      integer,                            intent(in)     :: multip
      type(Symm_Oper_Type), dimension(:), intent(in out) :: Op
      integer,                            intent(out)    :: num_lat,num_alat,Numops, centred, mag_type
      real(kind=cp),dimension(:,:),       intent(out)    :: Lat_tr
      real(kind=cp),dimension(:,:),       intent(out)    :: aLat_tr
      real(kind=cp),dimension(:),         intent(out)    :: centre_coord
      !--- Local variables ---!
      integer :: i,j,L,n,m,Ng,invt,i_centre,d
      real(kind=cp), dimension(multip) :: tr   !Sum of absolute values of Translations components associated to the array of operators
      logical,       dimension(multip) :: nul  !Logical to control the exclusion of an operator from the independet list
      real(kind=cp), dimension(:,:),allocatable:: identity,invers,imat !(d,d)
      type(Symm_Oper_Type), dimension(multip)  :: Opr,Op_Lat,Op_aLat
      type(Symm_Oper_Type)                     :: Op_aux,Op_aux1,Op_aux2, &
                                                  Op_centre,Op_identp
      real(kind=cp), parameter :: loc_eps=0.001
      real(kind=cp)            :: tmin

      !Initializing
      n=size(Op(1)%Mat,1) !dimension of the full square matrices
      d=n-1               !here d is the dimension of the square matrix containing rotational operator
      allocate(identity(d,d),invers(d,d),imat(d,d))
      call Allocate_Operator(n,Op_identp)   ! {1|0}'
      identity=0.0; nul=.false.; mag_type=1
      do i=1,d
        Op_identp%Mat(i,i)=1.0_cp
        identity(i,i)=1.0_cp
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

      !write(*,*) " "
      !write(*,*) "List of operators generated by Get_Group_From_Generators: "
      !do i=1,Multip
      !  write(*,"(i5,a50,f12.5,2i6)") i,trim(Symbol_Operator(Op(i))),tr(i),Op(i)%time_inv,Op(i)%dt
      !end do


      !Check if the group is paramagnetic
      j=0
      do i=2,Multip
        if(Op(i)== Op_identp) then
          j=i
          nul(i)=.true.
          exit
        end if
      end do
      if(j /= 0) mag_type=2
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
         centre_coord=0.5*Op(i_centre)%Mat(1:d,n)
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

      !write(*,*) " "
      !write(*,*) "List of operators after removing centre of symmetry and lattice translations: "
      !m=0
      !do i=1,Multip
      !  if(nul(i)) cycle
      !  m=m+1
      !  write(*,"(2i5,a30,f12.5,i6,tr2,L)") m,i,trim(Symbol_Operator(Op(i))),tr(i),Op(i)%dt,nul(i)
      !end do

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

      !write(*,*) " Reduced set of Symmetry operators ({-1|t} and centrings removed)"
      !do i=1,Numops
      !  write(*,"(i5,a50,2i6)") i,trim(Symbol_Operator(Op(i))),Op(i)%time_inv,Op(i)%dt
      !end do
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
      !Final reordering for paramagnetic groups: all primed operators are put at the bottom
      if(mag_type == 2) then
        m=0
        do i=1,multip
          if(Op(i)%time_inv == -1) Cycle
          m=m+1
          Opr(m)=Op(i)
        end do
        do i=1,multip
          if(Op(i)%time_inv == 1) Cycle
          m=m+1
          Opr(m)=Op(i)
        end do
        Op(1:multip) = Opr(1:multip)
        !write(*,*)  " Ordered operators for a paramagnetic group:"
        !do i=1,multip
        !  write(*,"(i5,a50,2i6)") i,trim(Symbol_Operator(Op(i))),Op(i)%time_inv,Op(i)%dt
        !end do
      end if
    End Subroutine Reorder_Operators_real

    Subroutine Reorder_Operators_rational(multip,Op, centred, centre_coord, Numops, num_lat, num_alat, Lat_tr, aLat_tr,mag_type)
      integer,                            intent(in)     :: multip
      type(rational_Symm_Oper_Type), dimension(:), intent(in out) :: Op
      integer,                             intent(out)    :: num_lat,num_alat,Numops, centred, mag_type
      type(rational),dimension(:,:),       intent(out)    :: Lat_tr
      type(rational),dimension(:,:),       intent(out)    :: aLat_tr
      type(rational),dimension(:),         intent(out)    :: centre_coord
      !--- Local variables ---!
      integer :: i,j,L,n,m,Ng,invt,i_centre,d
      real(kind=cp), dimension(multip) :: tr   !Sum of absolute values of Translations components associated to the array of operators
      logical,       dimension(multip) :: nul  !Logical to control the exclusion of an operator from the independet list
      type(rational), dimension(:,:),allocatable:: identity,invers,imat !(d,d)
      type(rational_Symm_Oper_Type), dimension(multip)  :: Opr,Op_Lat,Op_aLat
      type(rational_Symm_Oper_Type)                     :: Op_aux,Op_aux1,Op_aux2, &
                                                           Op_centre,Op_identp
      real(kind=cp), parameter :: loc_eps=0.001
      real(kind=cp)            :: tmin

      !Initializing
      n=size(Op(1)%Mat,1) !dimension of the full square matrices
      d=n-1               !here d is the dimension of the square matrix containing rotational operator
      allocate(identity(d,d),invers(d,d),imat(d,d))
      call Allocate_Operator(n,Op_identp)   ! {1|0}'
      identity=0_ik/1_ik; nul=.false.; mag_type=1
      do i=1,d
        Op_identp%Mat(i,i)=1_ik//1_ik
        identity(i,i)=1_ik//1_ik
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
      !Check if the group is paramagnetic
      j=0
      do i=2,Multip
        if(Op(i)== Op_identp) then
          j=i
          nul(i)=.true.
          exit
        end if
      end do
      if(j /= 0) mag_type=2
      !----End intial re-ordering

      !Look for centre of symmetry, and centring translations
      num_lat=0; num_alat=0; tmin=1.0e8; i_centre=0
      do j=2,Multip
         if(nul(j)) cycle
         invt= Op(j)%time_inv
         imat=Op(j)%Mat(1:d,1:d)
         if(equal_rational_matrix(identity,imat) .and. invt == 1) then
            num_lat=num_lat+1
            Lat_tr(:,num_lat)=Op(j)%Mat(1:d,n)
            Op_Lat(num_lat)=Op(j)
            nul(j)=.true.   !Nullify centring translations
            cycle
         end if
         if(equal_rational_matrix(imat,invers) .and. invt == 1 ) then
             nul(j) = .true.
             if(tr(j) < tmin) Then
               tmin= tr(j)
               i_centre=j
             end if
         end if
      end do

      if(i_centre /= 0) then
         Op_centre=Op(i_centre)
         centre_coord=Op(i_centre)%Mat(1:d,n)/2_ik
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
        if(equal_rational_matrix(identity,imat) .and. invt == -1) then
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
      !Final reordering for paramagnetic groups: all primed operators are put at the bottom
      if(mag_type == 2) then
        m=0
        do i=1,multip
          if(Op(i)%time_inv == -1) Cycle
          m=m+1
          Opr(m)=Op(i)
        end do
        do i=1,multip
          if(Op(i)%time_inv == 1) Cycle
          m=m+1
          Opr(m)=Op(i)
        end do
        Op(1:multip) = Opr(1:multip)
      end if
    End Subroutine Reorder_Operators_rational

    Subroutine Group_Constructor_string(generatorList,Grp)
       character(len=*),   intent(in)     :: generatorList
       Class (Group_Type), intent(in out) :: Grp
       !--- Local variables ---!
       character(len=40),dimension(:),allocatable :: gen
       type(Symm_Oper_Type), dimension(maxnum_op) :: Op
       integer :: d,i,j,ngen,invt,multip,centred,Numops,num_lat,num_alat,mag_type
       real(kind=cp),  dimension(:)  ,allocatable  :: centre_coord
       real(kind=cp),  dimension(:,:),allocatable  :: Lat_tr, aLat_tr
       character(len=80) :: Symb_Op
       type(rational), dimension(:,:),allocatable :: Mat
       !
       allocate(gen(maxnum_op))
       call Get_Operators_From_String(generatorList,d,ngen,gen)
       allocate(Mat(d,d))
       do i=1,maxnum_op
         call Allocate_Operator(d,Op(i))
       end do
       !Construct the list of the generators on top of Op
       !
       do i=1,ngen
         call Get_Mat_From_Symb_Op(gen(i),Mat,invt)
         if(Err_group) return
         Op(i+1)%Mat=Mat
         Op(i+1)%time_inv=invt
       end do
       ngen=ngen+1

       !Construct the raw Group
       call Get_Group_From_Generators(ngen,Op,multip)
       if(Err_group) return

       !Allocate provisionally to 1/4 Multip the lattice translations and anti-Translations
       allocate(Lat_tr(d-1,multip/4), aLat_tr(d-1,multip/4))
       allocate(centre_coord(d-1))

       call Reorder_Operators(multip, Op, centred, centre_coord, Numops, num_lat, num_alat, Lat_tr, aLat_tr,mag_type)
       if(Err_group) return

       Grp%multip=multip
       Grp%d=d
       Grp%Op=Op(1:multip)
       !allocate(character(len=40) :: Grp%Symb_Op(multip)) !it doesn't work in gfortran
       if(allocated(Grp%Symb_Op)) Deallocate(Grp%Symb_Op)
       allocate(Grp%Symb_Op(multip))
       do i=1,multip
          Grp%Symb_Op(i)=trim(Symbol_Operator(Op(i)))
       end do

       Select Type(Grp)

         type is (raw_spg_type)
           if(num_lat > 0) then
             if(allocated(Grp%Lat_tr)) Deallocate(Grp%Lat_tr)
             allocate(Grp%Lat_tr(1:d-1,1:Num_Lat))
           end if
           if(num_alat > 0) then
             if(allocated(Grp%aLat_tr)) Deallocate(Grp%aLat_tr)
             allocate(Grp%aLat_tr(1:d-1,1:Num_aLat))
           end if
           if(allocated(Grp%centre_coord)) Deallocate(Grp%centre_coord)
           allocate(Grp%centre_coord(1:d-1))
           Grp%Numops      = Numops
           Grp%centred     = centred
           Grp%mag_type    = mag_type
           Grp%num_lat     = num_lat
           Grp%num_alat    = num_alat
           Grp%centre_coord= centre_coord
           if(num_lat  > 0)  Grp%Lat_tr = Lat_tr(1:d-1,1:Num_Lat)
           if(num_alat > 0)  Grp%aLat_tr=aLat_tr(1:d-1,1:Num_aLat)

       End Select

    End Subroutine Group_Constructor_string

    Subroutine Group_Constructor_gen(gen,Grp)
       character(len=*),dimension(:),intent(in) :: gen
       Type(raw_spg_Type),       intent(in out) :: Grp
       !--- Local variables ---!

       type(Symm_Oper_Type), dimension(maxnum_op)   :: Op
       integer :: d,i,j,ngen,invt,multip,centred,Numops,num_lat,num_alat,mag_type
       real(kind=cp),  dimension(:)  ,allocatable  :: centre_coord
       real(kind=cp),  dimension(:,:),allocatable  :: Lat_tr, aLat_tr
       character(len=80) :: Symb_Op
       type(rational), dimension(:,:),allocatable :: Mat
       !
       ngen=size(gen)
       ! Get dimension of the generator matrices
       d=get_dimension(gen(1))

       allocate(Mat(d,d))
       do i=1,maxnum_op
         call Allocate_Operator(d,Op(i))
       end do
       !Construct the list of the generators on top of Op
       !
       do i=1,ngen
         call Get_Mat_From_Symb_Op(gen(i),Mat,invt)
         if(Err_group) return
         Op(i+1)%Mat=Mat
         Op(i+1)%time_inv=invt
       end do
       ngen=ngen+1

       !Construct the raw Group
       call Get_Group_From_Generators(ngen,Op,multip)
       if(Err_group) return

       !Allocate provisionally to 1/4 Multip the lattice translations and anti-Translations
       allocate(Lat_tr(d-1,multip/4), aLat_tr(d-1,multip/4))
       allocate(centre_coord(d-1))

       call Reorder_Operators(multip, Op, centred, centre_coord, Numops, num_lat, num_alat, Lat_tr, aLat_tr,mag_type)
       if(Err_group) return

       Grp%multip=multip
       Grp%d=d
       Grp%Op=Op(1:multip)
       !allocate(character(len=40) :: Grp%Symb_Op(multip)) !it doesn't work in gfortran
       if(allocated(Grp%Symb_Op)) Deallocate(Grp%Symb_Op)
       allocate(Grp%Symb_Op(multip))
       do i=1,multip
          Grp%Symb_Op(i)=trim(Symbol_Operator(Op(i)))
       end do

       if(num_lat > 0) then
         if(allocated(Grp%Lat_tr)) Deallocate(Grp%Lat_tr)
         allocate(Grp%Lat_tr(1:d-1,1:Num_Lat))
       end if
       if(num_alat > 0) then
         if(allocated(Grp%aLat_tr)) Deallocate(Grp%aLat_tr)
         allocate(Grp%aLat_tr(1:d-1,1:Num_aLat))
       end if
       if(allocated(Grp%centre_coord)) Deallocate(Grp%centre_coord)
       allocate(Grp%centre_coord(1:d-1))
       Grp%Numops      = Numops
       Grp%centred     = centred
       Grp%mag_type    = mag_type
       Grp%num_lat     = num_lat
       Grp%num_alat    = num_alat
       Grp%centre_coord= centre_coord
       if(num_lat  > 0)  Grp%Lat_tr = Lat_tr(1:d-1,1:Num_Lat)
       if(num_alat > 0)  Grp%aLat_tr=aLat_tr(1:d-1,1:Num_aLat)

    End Subroutine Group_Constructor_gen

    Subroutine Group_Constructor_rational(generatorList,Grp,mode)
       character(len=*),        intent(in)     :: generatorList
       type(rational_spg_Type), intent(in out) :: Grp
       character(len=*),        intent(in)     :: mode
       !--- Local variables ---!
       character(len=80),dimension(:),allocatable :: gen
       type(rational_Symm_Oper_Type), dimension(:),allocatable :: Op
       integer :: d,i,j,ngen,invt,multip,centred,Numops,num_lat,num_alat,mag_type
       type(rational), dimension(:)  ,allocatable  :: centre_coord
       type(rational), dimension(:,:),allocatable  :: Lat_tr, aLat_tr
       type(rational), dimension(:,:),allocatable  :: Mat
       character(len=80) :: Symb_Op
       !
       allocate(gen(maxnum_op))
       call Get_Operators_From_String(generatorList,d,ngen,gen)
       if(mode == "direct") then
         call Allocate_Operators(d,maxnum_op,Op)
       else
         allocate(Op(maxnum_op))
         do i=1,maxnum_op
           call Allocate_Operator(d,Op(i))
         end do
       end if
       allocate(Mat(d,d))
       !Construct the list of the generators on top of Op
       !
       do i=1,ngen
         call Get_Mat_From_Symb_Op(gen(i),Mat,invt)
         if(Err_group) return
         Op(i+1)%Mat=Mat
         Op(i+1)%time_inv=invt
       end do
       ngen=ngen+1

       !Construct the raw Group
       call Get_Group_From_Generators(ngen,Op,multip)
       if(Err_group) return

       !Allocate provisionally to 1/4 Multip the lattice translations and anti-Translations
       allocate(Lat_tr(d-1,multip/4), aLat_tr(d-1,multip/4))
       allocate(centre_coord(d-1))

       call Reorder_Operators(multip, Op, centred, centre_coord, Numops, num_lat, num_alat, Lat_tr, aLat_tr,mag_type)
       if(Err_group) return

       Grp%multip=multip
       Grp%d=d
       Grp%rOp=Op(1:multip)
       !allocate(character(len=40) :: Grp%Symb_Op(multip)) !it doesn't work in gfortran
       if(allocated(Grp%Symb_Op)) Deallocate(Grp%Symb_Op)
       allocate(Grp%Symb_Op(multip))
       do i=1,multip
          Grp%Symb_Op(i)=trim(Symbol_Operator(Op(i)))
       end do

       if(num_lat > 0) then
         if(allocated(Grp%Lat_tr)) Deallocate(Grp%Lat_tr)
         allocate(Grp%Lat_tr(1:d-1,1:Num_Lat))
       end if
       if(num_alat > 0) then
         if(allocated(Grp%aLat_tr)) Deallocate(Grp%aLat_tr)
         allocate(Grp%aLat_tr(1:d-1,1:Num_aLat))
       end if
       if(allocated(Grp%centre_coord)) Deallocate(Grp%centre_coord)
       allocate(Grp%centre_coord(1:d-1))
       Grp%Numops      = Numops
       Grp%centred     = centred
       Grp%mag_type    = mag_type
       Grp%num_lat     = num_lat
       Grp%num_alat    = num_alat
       Grp%centre_coord= centre_coord
       if(num_lat  > 0)  Grp%Lat_tr = Lat_tr(1:d-1,1:Num_Lat)
       if(num_alat > 0)  Grp%aLat_tr=aLat_tr(1:d-1,1:Num_aLat)


    End Subroutine Group_Constructor_rational

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
       Symbol=generatorList(1:i-1)

       d=Get_dimension(Symbol)

       call Get_Separator_Pos(generatorList,";",pos,np)

       !Verify if there is a final ";" not followed by a generator
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
    !!----   operator is Numops+1
    !!----
    Subroutine Get_SubGroups(SpG,SubG,nsg,point)
    !   !---- Arguments ----!
        Class(Group_Type) ,              intent( in) :: SpG
        Class(Group_Type) ,dimension(:), intent(out) :: SubG
        integer,                          intent(out) :: nsg
        logical, dimension(:,:), optional,intent(out) :: point
       !--- Local variables ---!
       integer                            :: i,L,j,k, nc, maxg,ng , nla, i1,i2,nop
       character (len=40), dimension(192) :: gen
       logical                            :: newg, cen_added

       maxg=size(SubG)
    !   !---- Construct first the generators of centring translations ----!
    !   ng=0
    !   nop=SpG%numops !number of symmetry operators excluding lattice centrings
    !   if (SpG%centred /= 1) nop=nop*2
    !   do i=2,SpG%numlat
    !      ng=ng+1
    !      gen(ng)= SpG%SymopSymb(1+nop*(i-1))
    !   end do
    !
    !   nla=ng
    !   nc=SpG%Numops+1  !Position of the centre of symmetry if it exist
    !   L=0
    !   !---- Determine first the triclinic subgroups
    !   cen_added=.false.
    !   do
    !       L=L+1
    !       newg=.true.
    !       call set_spacegroup(" ",SubG(L),gen,ng,"gen")
    !       do j=1,L-1
    !          if (SpGr_Equal(SubG(L), SubG(j))) then
    !             newg=.false.
    !             exit
    !          end if
    !       end do
    !       if (newg) then
    !          call get_HallSymb_from_gener(SubG(L))
    !       else
    !          L=L-1
    !       end if
    !       if (SpG%centred /= 1 .and. newg .and. .not. cen_added) then !add the centre of symmetry if needed
    !          ng=ng+1
    !          gen(ng)=SpG%SymopSymb(nc)
    !          cen_added=.true.
    !       else
    !          exit
    !       end if
    !   end do
    !
    !   !---- Determine first the groups with only one rotational generator
    !   do i=2,nop
    !      ng=nla+1
    !      gen(ng) = SpG%SymopSymb(i)
    !      cen_added=.false.
    !      do
    !         L=L+1
    !         if (L > maxg) then
    !            nsg=maxg
    !            return
    !         end if
    !         newg=.true.
    !         call set_spacegroup(" ",SubG(L),gen,ng,"gen")
    !         if(SubG(L)%multip == 0) then
    !           L=L-1
    !           newg=.false.
    !         else
    !           do j=1,L-1
    !              if (SpGr_Equal(SubG(L), SubG(j))) then
    !                 newg=.false.
    !                 exit
    !              end if
    !           end do
    !           if (newg) then
    !              call get_HallSymb_from_gener(SubG(L))
    !           else
    !              L=L-1
    !           end if
    !         end if
    !         if (SpG%centred /= 1 .and. newg .and. .not. cen_added) then !add the centre of symmetry if needed
    !            ng=ng+1
    !            gen(ng)=SpG%SymopSymb(nc)
    !            cen_added=.true.
    !         else
    !            exit
    !         end if
    !      end do
    !   end do
    !
    !   !---- Determine now the groups with two rotational generator ----!
    !
    !   do i1=2,nop-1
    !      gen(nla+1) = SpG%SymopSymb(i1)
    !      do i2 = i1+1,nop
    !         gen(nla+2) = SpG%SymopSymb(i2)
    !         ng=nla+2
    !         cen_added=.false.
    !         do
    !            L=L+1
    !            if (L > maxg) then
    !               nsg=maxg
    !               return
    !            end if
    !            newg=.true.
    !            call set_spacegroup(" ",SubG(L),gen,ng,"gen")
    !            if(mod(nop,SubG(L)%Numops) /= 0 .or. SubG(L)%multip == 0) then
    !              L=L-1
    !              newg=.false.
    !            else
    !              do j=1,L-1
    !                 if (SpGr_Equal(SubG(L), SubG(j))) then
    !                    newg=.false.
    !                    exit
    !                 end if
    !              end do
    !              if (newg) then
    !                 call get_HallSymb_from_gener(SubG(L))
    !              else
    !                 L=L-1
    !              end if
    !            end if
    !            if (SpG%centred /= 1 .and. newg .and. .not. cen_added) then !add the centre of symmetry if needed
    !               ng=ng+1
    !               gen(ng)=SpG%SymopSymb(nc)
    !               cen_added=.true.
    !            else
    !               exit
    !            end if
    !         end do
    !      end do
    !   end do
    !   nsg=L
    !   if(present(point)) then
    !     point=.false.
    !     do j=1,nsg
    !       L=1
    !       do i=1,SpG%multip
    !          do k=L,SubG(j)%multip
    !           if(SubG(j)%SymopSymb(k) == SpG%SymopSymb(i)) then
    !              point(i,j) = .true.
    !              L=k+1
    !              exit
    !           end if
    !          end do
    !       end do
    !     end do
    !   end if
    !
    !   return
    End Subroutine Get_SubGroups

  End Module CFML_Groups

   Program test_groups
     use CFML_Groups
     use CFML_Rational_Arithmetic
     character(len=256)   :: generatorList
     !type(raw_spg_type)        :: Grp
     type(rational_spg_type)   :: Grp
     integer :: i,j
     call Init_Group(2048) !Maximum admissible multiplicity
     do
         write(*,'(/,a)',advance='no') "Introduce generators: "
         read(*,'(a)') generatorList
         if (len_trim(generatorList) == 0) exit
         call Group_Constructor(generatorList,Grp,"dummy")
         if (Err_group) then
            write(*,'(/,4x,a)') Err_group_Mess
         else
            write(*,"(a)")        "    General Space Group"
            write(*,"(a)")        "    -------------------"
            write(*,"(a,i3)")     "     Op-Dimension: ",Grp%d
            write(*,"(a,i3)")     "  Space-Dimension: ",Grp%d-1
            write(*,"(a,i3)")     "     Multiplicity: ",Grp%multip
            write(*,"(a,i3)")     "          MagType: ",Grp%mag_type
            write(*,"(a,i3)")     "           NumOps: ",Grp%numops
            write(*,"(a,i3)")     "          Centred: ",Grp%centred
            write(*,"(a,i3)")     "          Num_Lat: ",Grp%num_lat
            write(*,"(a,i3)")     "         Num_aLat: ",Grp%num_alat
            if(Grp%centred == 1) then
               write(*,"(a)")     "     Centre_coord: none!"
            else
               !write(*,"(a,10f8.3)") "     Centre_coord: ",Grp%centre_coord
               write(*,"(a,10a)") "     Centre_coord: [",(trim(print_rational(Grp%centre_coord(i)))//" ",i=1,Grp%d-1),"]"
            end if
            if(Grp%num_lat > 0) then
              write(*,"(/a)")      "  Centring translations:"
              do i=1,Grp%num_lat
                 !write(*,"(i3,tr4,10f8.3)") i,Grp%Lat_tr(:,i)
               write(*,"(a,10a)") "      [",(trim(print_rational(Grp%Lat_tr(j,i)))//" ",j=1,Grp%d-1),"]"
             end do
            end if
            if(Grp%num_alat > 0) then
              write(*,"(/a)")      "  Anti-translations:"
              do i=1,Grp%num_alat
                ! write(*,"(i3,tr4,10f8.3)") i,Grp%aLat_tr(:,i)
               write(*,"(a,10a)") "      [",(trim(print_rational(Grp%aLat_tr(j,i)))//" ",j=1,Grp%d-1),"]"
              end do
            end if
            write(*,"(/a)")      "  Complete list of symmetry operators:"
            do i=1,Grp%Multip
              write(*,"(i5,a)") i,"  ->  "//trim(Grp%Symb_Op(i))
            end do
         end if
     end do
   End Program test_groups