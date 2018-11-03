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
               print_group,Get_SubGroups,Allocate_Group

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
      Character(len=80),   dimension(:),allocatable:: Symb_Op !It doesn't work properly ":" in gfortran
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
      module procedure equal_Group_real
      module procedure equal_Group_rational
    end interface
    private :: equal_Symm_Oper,equal_group_real, equal_group_rational,equal_rational_Symm_Oper

    interface Get_Group_From_Generators
      module procedure Get_Group_From_Generators_real
      module procedure Get_Group_From_Generators_rational
    end interface Get_Group_From_Generators
    private :: Get_Group_From_Generators_real, Get_Group_From_Generators_rational

    interface Group_Constructor
      module procedure Group_Constructor_gen
      module procedure Group_Constructor_string
      module procedure Group_Constructor_rational_gen
      module procedure Group_Constructor_rational_string
    end interface Group_Constructor
    Private :: Group_Constructor_gen,Group_Constructor_string, &
               Group_Constructor_rational_gen,Group_Constructor_rational_string

    interface Allocate_Operator
      module procedure Allocate_Operator_real
      module procedure Allocate_Operator_rational
    end interface Allocate_Operator
    private :: Allocate_Operator_real, Allocate_Operator_rational

    interface Allocate_Group
      module procedure Allocate_Group_real
      module procedure Allocate_Group_rational
    end interface Allocate_Group
    private :: Allocate_Group_real, Allocate_Group_rational

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
    private :: Reorder_Operators_real, Reorder_Operators_rational

    interface print_group
      module procedure print_rational_group
      module procedure print_real_group
    end interface print_group
    private :: print_rational_group, print_real_group

    interface Get_SubGroups
      module procedure Get_SubGroups_real
      module procedure Get_SubGroups_rational
    end interface
    private :: Get_SubGroups_real, Get_SubGroups_rational


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

    Pure Function equal_Group_rational(Gr1,Gr2) result (info)
      type(rational_spg_type), intent(in) :: Gr1,Gr2
      logical                              :: info
      integer :: i,j
      info=.false.
      if(Gr1%multip /= Gr2%multip) return
      do i=2,Gr1%multip
        do j=2,Gr2%multip
           if(Gr1%rOp(i) == Gr2%rOp(j)) then
             info=.true.
             exit
           end if
        end do
        if(.not. info) return
      end do
    End Function equal_Group_rational

    Pure Function equal_Group_real(Gr1,Gr2) result (info)
      type(raw_spg_type), intent(in) :: Gr1,Gr2
      logical                        :: info
      integer :: i,j
      info=.false.
      if(Gr1%multip /= Gr2%multip) return
      do i=2,Gr1%multip
        do j=2,Gr2%multip
           if( Gr1%Op(i) == Gr2%Op(j)) then
             info=.true.
             exit
           end if
        end do
        if(.not. info) return
      end do
    End Function equal_Group_real

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

    Subroutine Allocate_Group_rational(d,multip,Gr)
       integer,                 intent(in)     :: d,multip
       type(rational_spg_Type), intent(in out) :: Gr
       integer :: i
       if(allocated(Gr%Op)) deallocate(Gr%Op)
       if(allocated(Gr%rOp)) deallocate(Gr%rOp)
       allocate(Gr%Op(multip),Gr%rOp(multip))
       do i=1,multip
         call Allocate_Operator(d,Gr%Op(i))
         call Allocate_Operator(d,Gr%rOp(i))
       end do
    End Subroutine Allocate_Group_rational

    Subroutine Allocate_Group_real(d,multip,Gr)
       integer,                 intent(in)     :: d,multip
       type(raw_spg_Type), intent(in out) :: Gr
       integer :: i
       if(allocated(Gr%Op)) deallocate(Gr%Op)
       allocate(Gr%Op(multip))
       do i=1,multip
         call Allocate_Operator(d,Gr%Op(i))
       end do
    End Subroutine Allocate_Group_real

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
       use CFML_Math_3D, only : determinant_typed => determ_A
      integer,                                        intent(in)     :: ngen
      type(Symm_Oper_Type), dimension(:),             intent(in out) :: Op
      integer,                                        intent(out)    :: multip
      integer, dimension(:,:), allocatable, optional, intent(out)    :: table
      !--- Local variables ---!
      integer :: i,j,k,n,nt,max_op
      type(Symm_Oper_Type) :: Opt
      logical, dimension(size(Op),size(Op)) :: done
      integer, dimension(size(Op),size(Op)) :: tb

      include "CFML_get_group_template_inc.f90"

    End Subroutine Get_Group_From_Generators_real

    Subroutine Get_Group_From_Generators_rational(ngen,Op,multip,table)
      use CFML_Rational_Arithmetic, determinant_typed => rational_determinant
      integer,                                        intent(in)     :: ngen
      type(rational_Symm_Oper_Type), dimension(:),    intent(in out) :: Op
      integer,                                        intent(out)    :: multip
      integer, dimension(:,:), allocatable, optional, intent(out)    :: table
      !--- Local variables ---!
      integer :: i,j,k,n,nt,max_op
      type(rational_Symm_Oper_Type) :: Opt
      logical, dimension(size(Op),size(Op)) :: done
      integer, dimension(size(Op),size(Op)) :: tb

      include "CFML_get_group_template_inc.f90"

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
         err_group_mess="The determinant of the matrix is not integer! -> "//print_rational(det)//" -> "//trim(Symb)
      end if
      if(det%numerator == 0) then
         err_group=.true.
         err_group_mess="The matrix of the operator is singular! -> det="//print_rational(det)//" -> "//trim(Symb)
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
      real(kind=cp)            :: ZERO,ONE,ONE_HALF
      ZERO=0.0_cp; ONE=1.0_cp; ONE_HALF=0.5_cp
      include "CFML_reorder_operators_template_inc.f90"

    End Subroutine Reorder_Operators_real

    Subroutine Reorder_Operators_rational(multip,Op, centred, centre_coord, Numops, num_lat, num_alat, Lat_tr, aLat_tr,mag_type)
      use CFML_Rational_Arithmetic, equal_matrix => equal_rational_matrix
      integer,                            intent(in)      :: multip
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
      type(rational)           :: ZERO, ONE, ONE_HALF

      ZERO=0_ik/1_ik;  ONE=1_ik/1_ik; ONE_HALF=1_ik/2_ik

      include "CFML_reorder_operators_template_inc.f90"
    End Subroutine Reorder_Operators_rational

    Subroutine Group_Constructor_gen(gen,Grp)
       character(len=*),dimension(:),intent(in) :: gen
       Type(raw_spg_Type),       intent(in out) :: Grp
       !--- Local variables ---!

       type(Symm_Oper_Type), dimension(:),allocatable  :: Op
       integer :: d,i,j,ngen,invt,multip,centred,Numops,num_lat,num_alat,mag_type
       real(kind=cp),  dimension(:)  ,allocatable :: centre_coord
       real(kind=cp),  dimension(:,:),allocatable :: Lat_tr, aLat_tr
       character(len=80) :: Symb_Op
       type(rational), dimension(:,:),allocatable :: Mat
       !
       ngen=size(gen)
       ! Get dimension of the generator matrices
       d=get_dimension(gen(1))
       allocate(Op(maxnum_op))
       do i=1,maxnum_op
         call Allocate_Operator(d,Op(i))
       end do
       allocate(Mat(d,d))

       !Construct the list of the generators on top of Op. The identity is always the first operator
       do i=1,ngen
         call Get_Mat_From_Symb_Op(gen(i),Mat,invt)
         if(Err_group) return
         Op(i+1)%Mat=Mat
         Op(i+1)%time_inv=invt
       end do
       ngen=ngen+1

       include "CFML_group_constructor_template_inc.f90"

    End Subroutine Group_Constructor_gen

    !!----  Subroutine Group_Constructor_rational_ext(n,Opx,Grp)
    !!----     integer,                                    intent(in)      :: n
    !!----     type(rational_Symm_Oper_Type), dimension(:),intent(in)      :: Opx
    !!----     type(rational_spg_Type),                    intent(in out)  :: Grp
    !!----
    !!----  This subroutine completes Grp win "n" additional generators
    !!----
    !!----
    !Subroutine Group_Constructor_rational_ext(n,Opx,Grp)
    !   integer,                                    intent(in)      :: n
    !   type(rational_Symm_Oper_Type), dimension(:),intent(in)      :: Opx
    !   type(rational_spg_Type),                    intent(in out)  :: Grp
    !   !--- Local variables ---!
    !   type(rational_Symm_Oper_Type), dimension(:),allocatable :: Op
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
    !   call Allocate_Operators(d,multip,Grp%rOp)
    !   Grp%rOp(1:multip)=Op(1:multip)
    !End Subroutine Group_Constructor_rational_ext


    Subroutine Group_Constructor_rational_gen(gen,Grp)
       character(len=*),dimension(:),intent(in) :: gen
       type(rational_spg_Type), intent(in out)  :: Grp
       !--- Local variables ---!
       type(rational_Symm_Oper_Type), dimension(:),allocatable :: Op
       integer :: d,i,j,ngen,invt,multip,centred,Numops,num_lat,num_alat,mag_type,n
       type(rational), dimension(:)  ,allocatable  :: centre_coord
       type(rational), dimension(:,:),allocatable  :: Lat_tr, aLat_tr
       type(rational), dimension(:,:),allocatable  :: Mat
       character(len=80) :: Symb_Op

       ngen=size(gen)
       ! Get dimension of the generator matrices
       d=get_dimension(gen(1))
       allocate(Op(maxnum_op))
       do i=1,maxnum_op
         call Allocate_Operator(d,Op(i))
       end do
       allocate(Mat(d,d))

       !Construct the list of the generators on top of Op. The identity is always the first operator
       do i=1,ngen
         call Get_Mat_From_Symb_Op(gen(i),Mat,invt)
         if(Err_group) return
         Op(i+1)%Mat=Mat
         Op(i+1)%time_inv=invt
       end do
       ngen=ngen+1

       include "CFML_group_constructor_template_inc.f90"
       call Allocate_Operators(d,multip,Grp%rOp)
       Grp%rOp(1:multip)=Op(1:multip)
        !do n=1,Multip
        !   write(*,"(2(a,i3))") "  Operator #",n," Time inversion: ",Grp%rOp(n)%time_inv
        !   do i=1,Grp%d
        !      write(unit=*,fmt="(a,10a)") "      [ ",(trim(print_rational(Grp%rOp(n)%Mat(j,i)))//" ",j=1,Grp%d),"]"
        !   end do
        !end do
    End Subroutine Group_Constructor_rational_gen


    Subroutine Group_Constructor_string(generatorList,Grp)
       character(len=*),   intent(in)     :: generatorList
       type(raw_spg_type), intent(in out) :: Grp
       !--- Local variables ---!
       character(len=40),dimension(:),allocatable :: gen
       type(Symm_Oper_Type), dimension(:), allocatable:: Op
       integer :: d,i,j,ngen,invt,multip,centred,Numops,num_lat,num_alat,mag_type
       real(kind=cp),  dimension(:)  ,allocatable  :: centre_coord
       real(kind=cp),  dimension(:,:),allocatable  :: Lat_tr, aLat_tr
       character(len=80) :: Symb_Op
       type(rational), dimension(:,:),allocatable :: Mat
       !
       allocate(gen(maxnum_op))
       call Get_Operators_From_String(generatorList,d,ngen,gen)

       include "CFML_group_constructor_template_inc.f90"

    End Subroutine Group_Constructor_string

    Subroutine Group_Constructor_rational_string(generatorList,Grp) !,mode)
       character(len=*),        intent(in)     :: generatorList
       type(rational_spg_Type), intent(in out) :: Grp
       !character(len=*),        intent(in)     :: mode
       !--- Local variables ---!
       character(len=40),dimension(:),allocatable :: gen
       type(rational_Symm_Oper_Type), dimension(:),allocatable :: Op
       integer :: d,i,j,ngen,invt,multip,centred,Numops,num_lat,num_alat,mag_type
       type(rational), dimension(:)  ,allocatable  :: centre_coord
       type(rational), dimension(:,:),allocatable  :: Lat_tr, aLat_tr
       type(rational), dimension(:,:),allocatable  :: Mat
       character(len=80) :: Symb_Op

       allocate(gen(maxnum_op))
       call Get_Operators_From_String(generatorList,d,ngen,gen)
       allocate(Op(maxnum_op))
       do i=1,maxnum_op
         call Allocate_Operator(d,Op(i))
       end do
       allocate(Mat(d,d))

       !Construct the list of the generators on top of Op. The identity is always the first operator
       do i=1,ngen
         call Get_Mat_From_Symb_Op(gen(i),Mat,invt)
         if(Err_group) return
         Op(i+1)%Mat=Mat
         Op(i+1)%time_inv=invt
       end do
       ngen=ngen+1
       !
       include "CFML_group_constructor_template_inc.f90"
       Grp%rOp=Op(1:multip)
    End Subroutine Group_Constructor_rational_string

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
    Subroutine Get_SubGroups_rational(SpG,SubG,nsg,point)
    !!   !---- Arguments ----!
        type(rational_spg_type),             intent( in) :: SpG
        type(rational_spg_type),dimension(:),intent(out) :: SubG
        integer,                             intent(out) :: nsg
        logical, dimension(:,:), optional,   intent(out) :: point
       !--- Local variables ---!
       integer  :: i,L,j,k,m, nc, maxg,ng , nla, i1,i2,nop,n,ns_1,ns_2,ns_3,n_nc_group
       logical  :: newg, cen_added
       character (len=40), dimension(:),allocatable :: gen
       character (len=40), dimension(30)            :: gen_lat
       character (len=40)                           :: gen_cent
       type(rational_Symm_Oper_Type)                :: Op_cent
       type(rational_Symm_Oper_Type), dimension(30) :: Op_lat


       maxg=size(SubG)
       allocate(gen(SpG%multip))
       !---- Construct first the generators of centring translations ----!
       ng=0; nc=0
       nop=SpG%numops !number of symmetry operators excluding lattice centrings & centre of symmetry
       if (SpG%centred /= 1) then
          nop=nop*2 !!number of symmetry operators excluding lattice centrings
          nc=SpG%Numops+1  !Position of the centre of symmetry if it exist
          gen_cent=SpG%Symb_Op(nc)
          Op_cent=SpG%rOp(nc)
       end if
       do i=2,SpG%num_lat
          ng=ng+1
          gen_lat(ng)= SpG%Symb_Op(1+nop*(i-1))
          Op_lat(ng)= SpG%rOp(1+nop*(i-1))
       end do
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
       !---- Determine now the groups with three rotational generators (probably this is not needed!)
       if(SpG%numops > 3) then
         ng=3
         do i=2,SpG%numops-1
            gen(1) = SpG%Symb_Op(i)
            do j=i+1,SpG%numops-1
              gen(2)=SpG%Symb_Op(j)
              do m=j+1,SpG%numops
                 gen(3)=SpG%Symb_Op(m)
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
         end do
         ns_3=L-(ns_2+ns_1)
       end if
       n_nc_group=L
       !---- Determine now the new groups adding a centre of symmetry if it exists
       if (SpG%centred /= 1) then !This doubles the number of groups
         do i=1,n_nc_group
           L=L+1

           call Allocate_Group(SpG%d,2*SubG(i)%multip,SubG(L))
           k=SubG(i)%numops
           do j=1,SubG(i)%numops
             SubG(L)%Op(j)=SubG(i)%Op(j)
             SubG(L)%Symb_Op(j)=SubG(i)%Symb_Op(j)
             k=k+1
             SubG(L)%rOp(k)=SubG(i)%rOp(j)*Op_cent
             !SubG(L)%Op(k)=(SubG(L)%rOp(k))
             SubG(L)%Symb_Op(k)=Symbol_Operator(SubG(L)%rOp(k))
           end do
           !call Extends_Group(SubG(L),gen_cent)
         end do
       end if


       write(*,*) "Number of subgroups with up to three rotational generators: ",L
       do i=1,L
          write(*,*) "Subgroup # ",i
          call print_group(SubG(i))
       end do

       !include "CFML_subgroups_template_inc.f90"
    End Subroutine Get_SubGroups_rational

    Subroutine Get_SubGroups_real(SpG,SubG,nsg,point)
    !   !---- Arguments ----!
        type(raw_spg_type) ,              intent( in) :: SpG
        type(raw_spg_type) ,dimension(:), intent(out) :: SubG
        integer,                          intent(out) :: nsg
        logical, dimension(:,:), optional,intent(out) :: point
       !--- Local variables ---!
       integer  :: i,L,j,k, nc, maxg,ng , nla, i1,i2,nop
       logical  :: newg, cen_added
       character (len=40), dimension(:),allocatable :: gen

       maxg=size(SubG)
       allocate(gen(SpG%multip))
       ng=0
       nop=SpG%numops !number of symmetry operators excluding lattice centrings
       if (SpG%centred /= 1) nop=nop*2
       do i=2,SpG%num_lat
          ng=ng+1
          gen(ng)= SpG%Symb_Op(1+nop*(i-1))
       end do
       include "CFML_subgroups_template_inc.f90"
    End Subroutine Get_SubGroups_real

    subroutine print_rational_Group(Grp,lun)
      type(rational_spg_type), intent(in)   :: Grp
      integer, optional,       intent(in)   :: lun
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
    end subroutine print_rational_Group

    subroutine print_real_Group(Grp,lun)
      type(raw_spg_type), intent(in)   :: Grp
      integer, optional,  intent(in)   :: lun
      integer :: iout,i,j
      character(len=12) :: forma="(a,  f8.4,a)"
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
      write(unit=forma(4:5),fmt="(i2)") Grp%d-1
      if(Grp%centred == 1) then
         write(unit=iout,fmt="(a)")     "     Centre_coord: none!"
      else
         !write(unit=iout,fmt="(a,10f8.3)") "     Centre_coord: ",Grp%centre_coord
         write(unit=iout,fmt=forma) "     Centre_coord: [",(Grp%centre_coord(i),i=1,Grp%d-1),"]"
      end if
      if(Grp%num_lat > 0) then
        write(unit=iout,fmt="(/a)")      "  Centring translations:"
        do i=1,Grp%num_lat
           !write(unit=iout,fmt="(i3,tr4,10f8.3)") i,Grp%Lat_tr(:,i)
         write(unit=iout,fmt=forma) "      [",(Grp%Lat_tr(j,i),j=1,Grp%d-1),"]"
       end do
      end if
      if(Grp%num_alat > 0) then
        write(unit=iout,fmt="(/a)")      "  Anti-translations:"
        do i=1,Grp%num_alat
          ! write(*,"(i3,tr4,10f8.3)") i,Grp%aLat_tr(:,i)
         write(unit=iout,fmt=forma) "      [",(Grp%aLat_tr(j,i),j=1,Grp%d-1),"]"
        end do
      end if
      write(unit=iout,fmt="(/a)")      "  Complete list of symmetry operators:"
      do i=1,Grp%Multip
        write(unit=iout,fmt="(i5,a)") i,"  ->  "//trim(Grp%Symb_Op(i))
      end do
    end subroutine print_real_Group

  End Module CFML_Groups