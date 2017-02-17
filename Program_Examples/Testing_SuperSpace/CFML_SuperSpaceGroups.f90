    Module CFML_SuperSpaceGroups
      use CFML_ssg_datafile
      use CFML_String_Utilities, only: pack_string, Get_Separator_Pos
      use CFML_Rational_Arithmetic

      Implicit None
      private
      public :: Allocate_SSG_SymmOps, Set_SSG_Reading_Database, Write_SSG, Gen_Group, &
                Get_SSymSymb_from_Mat,Get_Mat_From_SSymSymb

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
        !type(rational),      allocatable,dimension(:,:,:):: Om            !Operator matrices (3+d+1,3+d+1,Multip) common denominator at (4+d,4+d)
        type(SSym_Oper_Type),allocatable, dimension(:)   :: SymOp         ! Crystallographic symmetry operators
        character(len=80),   allocatable,dimension(:)    :: SymopSymb     ! Alphanumeric Symbols for SYMM
      End Type SuperSpaceGroup_Type

      logical, public :: Err_ssg
      character(len=80), public :: Err_ssg_mess

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
      d=n-1
      Op3%Mat=matmul(Op1%Mat,Op2%Mat)
      Op3%Mat(1:d,n)=mod(Op3%Mat(1:d,n),1)
       do i=1,d
        if(Op3%Mat(i,n) < 0//1) Op3%Mat(i,n) = Op3%Mat(i,n) + 1
      end do
    end function multiply_ssg_symop

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

      do_ext:do
        n=nt
        do i=1,n
          do_j:do j=1,n
            if(done(i,j)) cycle
            Opt=Op(i)*Op(j)
            do k=1,nt
              if(equal_rational_matrix(Opt%Mat,Op(k)%Mat)) then
                tb(i,j)=k
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
      multip=nt
      if(present(table)) then
        allocate(Table(multip,multip))
        Table=tb
      end if
    End Subroutine Gen_Group

    Subroutine Set_SSG_Reading_Database(num,ssg,ok,Mess)
      integer,                    intent(in)  :: num
      type(SuperSpaceGroup_Type), intent(out) :: ssg
      Logical,                    intent(out) :: ok
      character(len=*),           intent(out) :: Mess
      !
      integer :: i,j,nmod,Dp,D,iclass,m
      type(rational), dimension(:,:), allocatable :: Inv
      type(SSym_Oper_Type)                        :: transla
      character(len=15) :: forma
      logical :: inv_found

      if(.not. ssg_database_allocated)  then
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
      Allocate(ssg%SymopSymb(ssg%Multip))

      ssg%kv=0.0; ssg%Latt_trans=0//1; ssg%aLatt_trans=0//1; ssg%time_rev=1; ssg%Centre_coord=0//1
      nmod=iclass_nmod(iclass)
      !forma="(a,i3, i4)"
      !write(forma(7:7),"(i1)") Dp
      do j=1,ssg%Num_Lat
         !write(*,forma) " Vector #",j, iclass_centering(1:Dp,j,iclass)
         ssg%Latt_trans(1:D,j)= rational_simplify(iclass_centering(1:D,j,iclass)//iclass_centering(Dp,j,iclass))
      end do
      do i=1,ssg%NumOps
         ssg%SymOp(i)%Mat= rational_simplify(igroup_ops(1:Dp,1:Dp,i,num)//igroup_ops(Dp,Dp,i,num))
      end do

      !Look for a centre of symmetry
      inv_found=.false.
      allocate(Inv(D,D),transla%Mat(Dp,Dp))
      Inv=0//1; transla%Mat=0//1
      do i=1,D
        Inv(i,i) = -1//1
      end do
      do i=1,Dp
        transla%Mat(i,i) = 1//1
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
      !Extend the symmetry operator to the whole set of lattice centring and
      !set the symmetry operators symbols
      m=ssg%NumOps
      do i=2,ssg%Num_Lat
        transla%Mat(1:D,Dp)=ssg%Latt_trans(1:D,i)
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
      do i=1,ssg%Multip
        call Get_SSymSymb_from_Mat(ssg%SymOp(i)%Mat,ssg%SymOpSymb(i),"xyz")
      end do

    End Subroutine Set_SSG_Reading_Database


    Subroutine Get_Mat_From_SSymSymb(Symb,Mat)
      character(len=*),                intent(in)  :: Symb
      type(rational),dimension(:,:),   intent(out) :: Mat
      !---- local variables ----!
      integer :: i,j,k,Dp, d, np,ns, n,m,inv,num,den,ind
      character(len=*),dimension(10),parameter :: xyz=(/"x","y","z","t","u","v","w","p","q","r"/)
      character(len=*),dimension(10),parameter :: x1x2x3=(/"x1 ","x2 ","x3 ","x4 ","x5 ","x6 ","x7 ","x8 ","x9 ","x10"/)
      character(len=*),dimension(10),parameter :: abc=(/"a","b","c","d","e","f","g","h","i","j"/)
      character(len=3),dimension(10)           :: x_typ
      character(len=len(Symb)), dimension(size(Mat,dim=1))  :: split
      character(len=len(Symb))                              :: string,pSymb,translation,transf
      character(len=10),        dimension(size(Mat,dim=1))  :: subst
      integer,                  dimension(size(Mat,dim=1)-1):: pos,pn
      logical                                               :: abc_transf

      Dp=size(Mat,dim=1)
      d=Dp-1
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
            Mat(i,Dp)= num//den
        end do
        split=" "

      else

        pSymb=pack_string(Symb)

      end if

      pos=0
      call Get_Separator_Pos(pSymb,",",pos,np)
      if(np /= d) then
        err_ssg=.true.
        err_ssg_mess="Error in the symbol of the operator"
        return
      end if

      read(unit=pSymb(pos(np)+1:),fmt=*) inv
      Mat(Dp,Dp)=inv//1
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
            Mat(i,Dp)= num//den
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
       integer                          :: i,j,Dp,d,k
       logical                          :: abc_type
       Dp=size(Mat,dim=1)
       d=Dp-1
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
         End Select
       end if
       !---- Main ----!
       symb=" "; translation=" "
       do i=1,d
          sym(i)=" "
          do j=1,d
             if(Mat(i,j) == 1) then
                sym(i) = trim(sym(i))//"+"//trim(x_typ(j))
             else if(Mat(i,j) == -1) then
                sym(i) =  trim(sym(i))//"-"//trim(x_typ(j))
             else if(Mat(i,j) /= 0) then
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
               if(Mat(i,j) > 0) car="+"//trim(car)
               sym(i)=trim(sym(i))//pack_string(car)
             end if
          end do
          !Write here the translational part for each component
          if (Mat(i,Dp) /= 0) then
            car=adjustl(print_rational(Mat(i,Dp)))
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
          car=print_rational(Mat(Dp,Dp))
          symb=trim(symb)//","//trim(car)
       end if
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


  End Module CFML_SuperSpaceGroups
