   Module CFML_Groups
    Use CFML_GlobalDeps,   only : cp,sp,dp,pi
    use CFML_Math_General, only : equal_matrix,equal_vector,Zbelong,modulo_lat
    use CFML_Math_3D,      only : determ_A
    Use CFML_Crystallographic_Symmetry, only: Set_SpaceGroup, Space_Group_Type,Get_Transl_Symbol, &
                                              Get_SymSymb, NS_Space_Group_Type,Sym_Oper_Type,&
                                              applyso, Lattice_trans,Get_Trasfm_Symbol,SpGr_Equal

    implicit none

    Type, public :: Symm_Oper_Type
       integer       :: time_rev=1
       integer       :: dt=1  !determinant of the submatrix (1:3,1:3), it should be 1 or -1
       real(kind=cp), allocatable, dimension(:,:) :: Mat
    End Type Symm_Oper_Type

    Type, public :: Group
      integer :: order
      integer :: d  !Dimension of operator matrices (for 2D, d=3, for 3D, d=4, etc.)
      type(Symm_Oper_Type),dimension(:),allocatable:: Op
    End Type Group

    public :: operator (*)
    interface operator (*)
      module procedure multiply_Symm_Oper
    end interface
    private :: multiply_Symm_Oper

    public :: operator (==)
    interface operator (==)
      module procedure equal_Symm_Oper
    end interface
    private :: equal_Symm_Oper

    logical :: Err_group
    character(len=:), allocatable :: Err_group_mess

    Contains

    Pure function multiply_Symm_Oper(Op1,Op2) result (Op3)
      type(Symm_Oper_Type), intent(in) :: Op1,Op2
      type(Symm_Oper_Type)             :: Op3
      integer :: n,d,i
      n=size(Op1%Mat,dim=1)
      d=n-1
      Op3%Mat=matmul(Op1%Mat,Op2%Mat) !automatic allocation in f2003
      do i=1,d     !Put translations in the interval [0.0,1.0)
         if(Op3%Mat(i,n) < 0.0_cp) then
         	  Op3%Mat(i,n) = Op3%Mat(i,n) + 1.0_cp
         else if(Op3%Mat(i,n) >= 1.0_cp) then
         	  Op3%Mat(i,n) = Op3%Mat(i,n) - 1.0_cp
         end if
      end do
      Op3%time_rev=Op1%time_rev*Op2%time_rev
      Op3%dt=Op1%dt*Op2%dt
    end function multiply_Symm_Oper

    Pure function equal_Symm_Oper(Op1,Op2) result (info)
      type(Symm_Oper_Type), intent(in) :: Op1,Op2
      logical                          :: info
      info=.false.
      if(Op1%time_rev == Op2%time_rev) then
        if(equal_matrix(Op1%Mat,Op2%Mat)) info=.true.
      end if
    end function equal_Symm_Oper

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

    Subroutine Allocate_Operator(d,Op)
       integer,              intent(in)     :: d
       type(Symm_Oper_Type), intent(in out) :: Op
       if(allocated(Op%Mat)) deallocate(Op%Mat)
       allocate(Op%Mat(d,d))
       Op%Mat=0.0
    End Subroutine Allocate_Operator

    !!---- Subroutine Allocate_Operators(d,multip,Op)
    !!----    integer,              intent(in)     :: d,multip
    !!----    type(Symm_Oper_Type), dimension(:), allocatable, intent(in out) :: Op
    !!----
    !!----  Multip is the expected maximum number of operators
    !!----
    Subroutine Allocate_Operators(d,multip,Op)
       integer,              intent(in)     :: d,multip
       type(Symm_Oper_Type), dimension(:), allocatable, intent(in out) :: Op
       integer :: i
       if(allocated(Op)) deallocate(Op)
       allocate(Op(multip))
       do i=1,multip
         call Allocate_Operator(d,Op(i))
       end do
    End Subroutine Allocate_Operator

    Subroutine Get_Group_From_Generators(ngen,Op,multip,table)
      integer,                                        intent(in)     :: ngen
      type(Symm_Oper_Type), dimension(:),             intent(in out) :: Op
      integer,                                        intent(out)    :: multip
      integer, dimension(:,:), allocatable, optional, intent(out)    :: table
      !--- Local variables ---!
      integer :: i,j,k,n, nt,max_op,Dd,D
      type(Symm_Oper_Type) :: Opt
      logical, dimension(size(Op),size(Op)) :: done
      integer, dimension(size(Op),size(Op)) :: tb
      max_op=size(Op)
      done=.false.
      done(1,:) = .true.
      done(:,1) = .true.
      tb(1,:) = [(i,i=1,max_op)]
      tb(:,1) = [(i,i=1,max_op)]
      nt=ngen
      Err_group=.false.
      Err_group_mess=" "

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
    End Subroutine Get_Group_From_Generators

    Subroutine Reorder_Operators(Op, centre, centre_coord, Numops, num_lat, num_alat, Lat_tr, aLat_tr,mag_type)
      type(Symm_Oper_Type), dimension(:), intent(in out) :: Op
      integer,                            intent(out)    :: num_lat,num_alat,Numops, centre, mag_type
      real(kind=cp),dimension(:,:),       intent(out)    :: Lat_tr
      real(kind=cp),dimension(:,:),       intent(out)    :: aLat_tr
      real(kind=cp),dimension(:),         intent(out)    :: centre_coord
      !--- Local variables ---!
      integer :: i,j,k,kp,L,n,m,Ng,invt,i_centre,d,multip
      real(kind=cp), dimension(size(Op(1)%Mat-1,1))      :: v    !Translation associated to an operator
      real(kind=cp), dimension(size(Op))                 :: tr   !Sum of absolute values of Translations components associated to the array of operators
      logical,       dimension(size(Op))                 :: nul  !Logical to control
      real(kind=cp), dimension(size(Op(1)%Mat-1,1),size(Op(1)%Mat-1,1)):: identity,zerom,invers,imat
      type(Symm_Oper_Type), dimension(size(Op)):: Opr,Op_Lat,Op_aLat
      type(Symm_Oper_Type)                     :: Op_aux
      integer,              dimension(size(Op)):: ip,inp
      logical                  :: centrosymm
      real(kind=cp), parameter :: loc_eps=0.0005
      real(kind=cp)            :: tmin

      !Initializing
      multip=size(Op)
      n=size(Op(1)%Mat,1) !dimension of the full square matrices
      d=n-1               !dimension of the square matrix containing rotational operator
      identity=0.0; nul=.false.; mag_type=1; centrosymm=.false.; zerom=0.0
      do i=1,d
        identity(i,i)=1.0_cp
      end do
      invers=-identity !Inversion
      centre=1 !Default value for non-centrosymmetric groups
      do i=1,Multip
        tr(i)=sum(abs(Op(i)%Mat(1:d,n)))
      end do
      ip=0
      call sort(tr,Multip,ip)
      Op=Op(ip) !re-ordered by the smallest translations
      tr=tr(ip)
      !Reorder again the operators in case the identity is not given as the first operator
      j=0
      imat=Op(1)%Mat(1:d,1:d)
      if(.not. ( equal_matrix(imat,identity) .and.  sum(abs(Opr(1)%Mat(1:d,n)))  < loc_eps) ) then !Look for the identity
        do i=2,Multip
          imat=Op(i)%Mat(1:d,1:d)
          if(equal_matrix(imat,identity) .and. tr(i)  < loc_eps) then
            j=i  !Position of the identity
            exit
          end if
        end do
        if(j == 0) then
          Err_group=.true.
          Err_group_mess="The identity operator is not provided in the list of operators"
          return
        end if
        Op(j)=Op(1)
        Op(1)%Mat=0.0_cp
        Op(1)%Mat(1:d,1:d)=identity
      end if
      !Now look for centre of symmetry associated with time inversion and promote
      !it to the second position
      j=0
      do i=2,Multip
        imat=Op(i)%Mat(1:d,1:d)
        if(equal_matrix(imat,invers) .and. tr(i)  < loc_eps .and. Op(i)%time_rev == -1) then
          j=i
          exit
        end if
      end do
      if(j /= 0) then
        Op(j)=Op(2)
        Op(2)%Mat=0.0_cp
        Op(2)%Mat(1:d,1:d)=invers
      end if
      !----End intial re-ordering
      ip=0

      !Determine the lattice translations and anti-translations
      !Eliminate lattice translations and anti-translations
      num_lat=0; num_alat=0
      do j=2,Multip
        if(nul(j)) cycle
        invt= Op(j)%time_rev
        if(invt < 0) mag_type=3
        imat=Op(j)%Mat(1:d,1:d)
        if(equal_matrix(identity,imat)) then
           if(invt == 1) then
              num_lat=num_lat+1
              Lat_tr(:,num_lat)=Op(j)%Mat(1:d,n)
              Op_Lat(num_lat)=Op(j)
              nul(j)=.true.   !Nullify centring translations
           else
              num_alat=num_alat+1
              aLat_tr(:,num_alat)=Op(j)%Mat(1:d,n)
              Op_aLat(num_alat)=Op(j)
              nul(j)=.true.  !Nullify anti-centring translations
           end if
        end if
      end do  !j=2,Multip

      !Eliminate centre of symmetry {-1|t} and select that having
      !t=0 if it exist
      k=0; kp=0; ip=0; inp=0
      do j=2,Multip
          invt= Op(j)%time_rev
          imat=Op(j)%Mat(1:d,1:d)
          if(equal_matrix(imat,invers)) then
            if(invt == 1) then
              kp=kp+1  !c.o.symm without time reversal: -1
              ip(kp)=j
            else
              k=k+1    !c.o.symm with time reversal: -1'
              inp(k)=j
            end if
          end if
      end do
      i_centre=0
      if(kp > 0) then  !Localize Centres of symmetry, select that with minimal or no translations and nullify them
         centrosymm=.true.
         i_centre=ip(1)
         tmin=1.0e8
         do j=1,kp
           i=ip(j)
           nul(i)=.true.
           if(tr(i)  < tmin) then
             tmin=tr(i)
             i_centre=i
           end if
         end do
         if(tmin < loc_eps) then  !localization of the -x,-y,-z,+1 operator within the list
             centre=2             !Now this value indicates that the operator -x,-y,-z,+1 exists
         Else
             centre=0             !c.o.symm not at the origin
         end if

      end if
      !Nullify the operators of inversion centres associated with time inversion
      !and have a translation corresponding to a centring or anticentring vector
      do i=1,k
         j=inp(i)
         v=Op(j)%Mat(1:d,n)
         !if(tr(j) < loc_eps) cycle !Maintain the operator -x,-y,-z,-1'
         if(is_Lattice_vec(V,Lat_tr,num_lat)) then
            nul(j)=.true.
            cycle
         end if
         if(is_Lattice_vec(V,aLat_tr,num_alat)) then
            nul(j)=.true.
         end if
      end do
      !Nullify the operators that can be deduced from others by applying translations,
      !anti-translations and centre of symmetry

      ip=0
      do j=2,Multip-1
         if(nul(j)) cycle
         invt=Op(j)%time_rev
         do i=j+1,Multip
           if(nul(i)) cycle
           imat=Op(i)%Mat(1:d,1:d)-Op(j)%Mat(1:d,1:d)

           if(equal_matrix(imat,zerom) ) then  !Pure lattice translation or antitranslation
              v=Op(i)%Mat(1:d,n)-Op(j)%Mat(1:d,n)
              if(is_Lattice_vec(V,Lat_tr,num_lat)) then  !lattice translation
                 nul(i)=.true.
                 cycle
              end if
              if(is_Lattice_vec(V,aLat_tr,num_alat)) then  !lattice antitranslation
                 nul(i)=.true.
                 cycle
              end if
           end if

           if(centrosymm) then
              imat=Op(i)%Mat(1:d,1:d)+Op(j)%Mat(1:d,1:d)
              k=Op(i)%time_rev
              if(equal_matrix(imat,zerom) .and. k == invt) then
                 v=Op(i_centre)%Mat(1:d,n)-Op(i)%Mat(1:d,n)-Op(j)%Mat(1:d,n)
                 if(is_Lattice_vec(V,Lat_tr,num_lat)) then
                    nul(i)=.true.
                    cycle
                 end if
              end if

              if(equal_matrix(imat,zerom) .and. k /= invt) then
                 if(is_Lattice_vec(V,aLat_tr,num_alat)) then
                    nul(i)=.true.
                    cycle
                 end if
              end if

           end if
         end do
      end do

      ! => Determine the reduced set of symmetry operators"
      j=0
      do i=1,Multip
        if(nul(i)) cycle
        j=j+1
        Opr(j) = Op(i)
      end do
      m=j*max(centre,1)*(num_alat+num_lat+1)
      if( m /= Multip) then !Check that it is OK
        write(unit=Err_group_mess,fmt="(2(a,i4))") " Warning! Multip=",Multip, " Calculated Multip: ",m
        Err_group=.true.
        return
      end if

      !Promote the reduced set of symmetry operators to the top of the list
      Op(1:j)=Opr(1:j)
      Numops=j

      !Calculate the determinant of the rotational 3x3 submatrix
      do i=1,Numops
        Op(i)%dt=nint(determ_A(Op(i)%Mat(1:3,1:3)))
      end do

      !Re-Construct, in an ordered way, all the symmetry operators
      !starting with the reduced set
      m=Numops

      Centre_coord = 0.5_cp*Op(i_centre)%Mat(1:d,n)

      if(centrosymm) then   !First apply the centre of symmetry
        do i=1,Numops
          m=m+1
          Op(m) = Op(i) * Op(i_centre)
          if(Op(i)%dt == -1) then  !Swap the operators to put on top direct rotations
            Op_aux=Op(i)
            Op(i) = Op(m)
            Op(m) = Op_aux
          end if
        end do
      end if

      ng=m  ! Number or symmetry operators including centre of symmetry

      if(Num_aLat > 0) then   !Second apply the lattice centring anti-translations
        do L=1,Num_aLat
           do i=1,ng
             m=m+1
             Op(m)=Op_aLat(L)*Op(i)
           end do
        end do
      end if
      if(Num_Lat > 1) then  !Third apply the lattice centring translations
        do L=1,Num_Lat
           do i=1,ng
             m=m+1
             Op(m)=Op_Lat(L)*Op(i)
           end do
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



    Subroutine Get_T_SubGroups(SpG,SubG,nsg,point)
       !---- Arguments ----!
       type (Space_Group_Type) ,             intent( in) :: SpG
       type (Space_Group_Type) ,dimension(:),intent(out) :: SubG
       integer,                              intent(out) :: nsg
       logical, dimension(:,:), optional,    intent(out) :: point
       !--- Local variables ---!
       integer                            :: i,L,j,k, nc, maxg,ng , nla, i1,i2,nop
       character (len=40), dimension(192) :: gen
       logical                            :: newg, cen_added

       maxg=size(SubG)
       !---- Construct first the generators of centring translations ----!
       ng=0
       nop=SpG%numops !number of symmetry operators excluding lattice centrings
       if (SpG%centred /= 1) nop=nop*2
       do i=2,SpG%numlat
          ng=ng+1
          gen(ng)= SpG%SymopSymb(1+nop*(i-1))
       end do

       nla=ng
       nc=SpG%Numops+1  !Position of the centre of symmetry if it exist
       L=0
       !---- Determine first the triclinic subgroups
       cen_added=.false.
       do
           L=L+1
           newg=.true.
           call set_spacegroup(" ",SubG(L),gen,ng,"gen")
           do j=1,L-1
              if (SpGr_Equal(SubG(L), SubG(j))) then
                 newg=.false.
                 exit
              end if
           end do
           if (newg) then
              call get_HallSymb_from_gener(SubG(L))
           else
              L=L-1
           end if
           if (SpG%centred /= 1 .and. newg .and. .not. cen_added) then !add the centre of symmetry if needed
              ng=ng+1
              gen(ng)=SpG%SymopSymb(nc)
              cen_added=.true.
           else
              exit
           end if
       end do

       !---- Determine first the groups with only one rotational generator
       do i=2,nop
          ng=nla+1
          gen(ng) = SpG%SymopSymb(i)
          cen_added=.false.
          do
             L=L+1
             if (L > maxg) then
                nsg=maxg
                return
             end if
             newg=.true.
             call set_spacegroup(" ",SubG(L),gen,ng,"gen")
             if(SubG(L)%multip == 0) then
               L=L-1
               newg=.false.
             else
               do j=1,L-1
                  if (SpGr_Equal(SubG(L), SubG(j))) then
                     newg=.false.
                     exit
                  end if
               end do
               if (newg) then
                  call get_HallSymb_from_gener(SubG(L))
               else
                  L=L-1
               end if
             end if
             if (SpG%centred /= 1 .and. newg .and. .not. cen_added) then !add the centre of symmetry if needed
                ng=ng+1
                gen(ng)=SpG%SymopSymb(nc)
                cen_added=.true.
             else
                exit
             end if
          end do
       end do

       !---- Determine now the groups with two rotational generator ----!

       do i1=2,nop-1
          gen(nla+1) = SpG%SymopSymb(i1)
          do i2 = i1+1,nop
             gen(nla+2) = SpG%SymopSymb(i2)
             ng=nla+2
             cen_added=.false.
             do
                L=L+1
                if (L > maxg) then
                   nsg=maxg
                   return
                end if
                newg=.true.
                call set_spacegroup(" ",SubG(L),gen,ng,"gen")
                if(mod(nop,SubG(L)%Numops) /= 0 .or. SubG(L)%multip == 0) then
                  L=L-1
                  newg=.false.
                else
                  do j=1,L-1
                     if (SpGr_Equal(SubG(L), SubG(j))) then
                        newg=.false.
                        exit
                     end if
                  end do
                  if (newg) then
                     call get_HallSymb_from_gener(SubG(L))
                  else
                     L=L-1
                  end if
                end if
                if (SpG%centred /= 1 .and. newg .and. .not. cen_added) then !add the centre of symmetry if needed
                   ng=ng+1
                   gen(ng)=SpG%SymopSymb(nc)
                   cen_added=.true.
                else
                   exit
                end if
             end do
          end do
       end do
       nsg=L
       if(present(point)) then
         point=.false.
         do j=1,nsg
           L=1
           do i=1,SpG%multip
              do k=L,SubG(j)%multip
               if(SubG(j)%SymopSymb(k) == SpG%SymopSymb(i)) then
                  point(i,j) = .true.
                  L=k+1
                  exit
               end if
              end do
           end do
         end do
       end if

       return
    End Subroutine Get_T_SubGroups

   End Module CFML_Groups
