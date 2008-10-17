!!----
!!---- Copyleft(C) 1999-2005,              Version: 3.0
!!---- Juan Rodriguez-Carvajal & Javier Gonzalez-Platas
!!----
!!---- MODULE: MAGNETIC_SYMMETRY
!!----   INFO: Series of procedures handling operations with Magnetic Symmetry
!!----         and Magnetic Structures
!!----
!!---- HISTORY
!!----    Update: April - 2005
!!----
!!----            April - 2005. Created by JRC
!!----
!!---- DEPENDENCIES
!!--++    Use Math_gen,                  only: Sp, tpi, Modulo_Lat
!!--++    Use Symmetry_Tables,           only: ltr_a,ltr_b,ltr_c,ltr_i,ltr_r,ltr_f
!!--++    Use Crystallographic_Symmetry, only: Space_Group_Type, Read_Xsym, Get_SymSymb, &
!!--++                                         Sym_Oper_Type, Set_SpaceGroup,read_msymm, symmetry_symbol, &
!!--++                                         err_mess_symm,err_symm
!!--++    Use String_Utilities,          only: u_case, l_case, Frac_Trans_1Dig
!!--++    Use IO_Formats,                only: file_list_type
!!--++    Use Atom_Module,               only: Allocate_mAtom_list, mAtom_List_Type
!!----
!!---- VARIABLES
!!--..    Types
!!----    MSYM_OPER_TYPE
!!----    MAGNETIC_DOMAIN_TYPE
!!----    MAGNETIC_GROUP_TYPE
!!----    MAGSYMM_K_TYPE
!!--..
!!----    ERR_MAGSYM
!!----    ERR_MESS_MAGSYM
!!----
!!---- PROCEDURES
!!----    Functions:
!!----       APPLYMSO
!!----
!!----    Subroutines:
!!----       INIT_ERR_MAGSYM
!!--++       INIT_MAGSYMM_K_TYPE             [Private]
!!----       READN_SET_MAGNETIC_STRUCTURE
!!----       SET_SHUBNIKOV_GROUP
!!----       WRITE_MAGNETIC_STRUCTURE
!!----       WRITE_SHUBNIKOV_GROUP
!!----
!!
 Module Magnetic_Symmetry

    !---- Use Modules ----!
    Use Math_gen,                  only: Sp, tpi, Modulo_Lat
    Use Symmetry_Tables,           only: ltr_a,ltr_b,ltr_c,ltr_i,ltr_r,ltr_f
    Use Crystallographic_Symmetry, only: Space_Group_Type, Read_Xsym, Get_SymSymb, &
                                         Sym_Oper_Type, Set_SpaceGroup,read_msymm, symmetry_symbol, &
                                         err_mess_symm,err_symm, set_SpG_Mult_Table,ApplySO
    Use String_Utilities,          only: u_case, l_case, Frac_Trans_1Dig
    Use IO_Formats,                only: file_list_type
    Use Atom_Module,               only: Allocate_mAtom_list, mAtom_List_Type
    Use Scattering_Chemical_Tables,only: Set_Magnetic_Form, Remove_Magnetic_Form, num_mag_form, &
                                         Magnetic_Form
    Use Propagation_vectors,       only: K_Equiv_Minus_K

    !---- Variables ----!
    implicit none

    private

    !---- List of public functions ----!
    public :: ApplyMSO

    !---- List of public subroutines ----!
    public :: Readn_Set_Magnetic_Structure, Write_Magnetic_Structure, Set_Shubnikov_Group, &
              Write_Shubnikov_Group

    private:: Init_MagSymm_k_Type

    !---- Definitions ----!

    !!----
    !!---- TYPE :: MSYM_OPER_TYPE
    !!--..
    !!---- Type, public :: MSym_Oper_Type
    !!----    integer, dimension(3,3) :: Rot     !  Rotational Part of Symmetry Operator
    !!----    real(kind=sp)           :: Phas    !  Phase in fraction of 2pi
    !!---- End Type  MSym_Oper_Type
    !!----
    !!----  Definition of Magnetic symmetry operator type
    !!----
    !!---- Update: April - 2005
    !!
    Type, public :: MSym_Oper_Type
       integer, dimension(3,3) :: Rot
       real(kind=sp)           :: Phas
    End Type MSym_Oper_Type

    !!----
    !!---- TYPE :: MAGNETIC_DOMAIN_TYPE
    !!--..
    !!---- Type, public :: Magnetic_Domain_type
    !!----    integer                         :: nd=0          !Number of rotational domains (not counting chiral domains)
    !!----    logical                         :: Chir=.false.  !True if chirality domains exist
    !!----    integer,dimension(3,3,24)       :: DMat=0        !Domain matrices to be applied to Fourier Coefficients
    !!----    real(kind=sp), dimension (2,24) :: pop=0.0       !Populations of domains (sum=1,
    !!----                                                     !the second value is /=0 for chir=.true.)
    !!----    real(kind=sp), dimension (2,24) :: Lpop=0        !Number of the refined parameter
    !!----    real(kind=sp), dimension (2,24) :: Mpop=0.0      !Refinement codes for populations
    !!---- End type Magnetic_Domain_type
    !!----
    !!----
    !!--<<
    !!----  Magnetic S-domains corresponds to a different magnetic structure obtained from
    !!----  the domain 1 (actual model) by applying a rotational operator to the Fourier
    !!----  coefficients of magnetic moments. This rotational operator corresponds to a
    !!----  symmetry operator of the paramagnetic group that is lost in the ordered state.
    !!----  Chirality domains are simply obtained by changing the sign of the imaginary
    !!----  components of the Fourier coefficients. For each rotational domain two chiralities
    !!----  domains exist.
    !!-->>
    !!----
    !!---- Update: October - 2006
    !!
    Type, public :: Magnetic_Domain_type
       integer                                  :: nd=0          !Number of rotational domains (not counting chiral domains)
       logical                                  :: Chir=.false.  !True if chirality domains exist
       integer,dimension(3,3,24)                :: DMat=0        !Domain matrices to be applied to Fourier Coefficients
       real(kind=sp), dimension (2,24)          :: pop=0.0       !Populations of domains (sum=1,
                                                                 !the second value is /=0 for chir=.true.)
       real(kind=sp), dimension (2,24)          :: Lpop=0        !Number of the refined parameter
       real(kind=sp), dimension (2,24)          :: Mpop=0.0      !Refinement codes for populations
    End type Magnetic_Domain_type

    !!----
    !!---- TYPE :: MAGNETIC_GROUP_TYPE
    !!--..
    !!---- Type, Public :: Magnetic_Group_Type
    !!----    Character(len=30)           :: Shubnikov !Shubnikov symbol (Hermman-Mauguin + primes)
    !!----    type(Space_Group_Type)      :: SpG       !Crystallographic space group
    !!----    integer, dimension(192)     :: tinv      !When a component is +1 no time inversion is associated
    !!---- End Type Magnetic_Group_Type                !If tinv(i)=-1, the time inversion is associated to operator "i"
    !!----
    !!--<<
    !!----    A magnetic group type is adequate when k=(0,0,0). It contains as the second
    !!----    component the crystallographic space group. The first component is
    !!----    the Shubnikov Group symbol and the third component is an integer vector with
    !!----    values -1 or 1 when time inversion is associated (-1) with the corresponding
    !!----    crystallographic symmetry operator o not (1).
    !!-->>
    !!----
    !!---- Update: April - 2005
    !!
    Type, Public :: Magnetic_Group_Type
       Character(len=30)           :: Shubnikov
       type(Space_Group_Type)      :: SpG
       integer, dimension(192)     :: tinv
    End Type Magnetic_Group_Type

    !!----
    !!---- TYPE :: MAGSYMM_K_TYPE
    !!--..
    !!---- Type, Public :: MagSymm_k_Type
    !!----    character(len=31)                        :: MagModel   ! Name to characterize the magnetic symmetry
    !!----    character(len=1)                         :: Latt       ! Symbol of the crystallographic lattice
    !!----    integer                                  :: nirreps    ! Number of irreducible representations (max=4, if nirreps /= 0 => nmsym=0)
    !!----    integer                                  :: nmsym      ! Number of magnetic operators per crystallographic operator (max=8)
    !!----    integer                                  :: centred    ! =0 centric centre not at origin, =1 acentric, =2 centric (-1 at origin)
    !!----    integer                                  :: mcentred   ! =1 Anti/a-centric Magnetic symmetry, = 2 centric magnetic symmetry
    !!----    integer                                  :: nkv        ! Number of independent propagation vectors
    !!----    real(kind=sp),dimension(3,12)            :: kvec       ! Propagation vectors
    !!----    integer                                  :: NumLat     ! Number of centring lattice vectors
    !!----    real(kind=sp), dimension(3,4)            :: Ltr        ! Centring translations
    !!----    integer                                  :: Numops     ! Reduced number of crystallographic Symm. Op.
    !!----    integer                                  :: Multip     ! General multiplicity of the space group
    !!----    integer,             dimension(4)        :: nbas       ! Number of basis functions per irrep (if nbas < 0, the corresponding basis is complex).
    !!----    integer,             dimension(12,4)     :: icomp      ! Indicator (0 pure real/ 1 pure imaginary) for coefficients of basis fucntions
    !!----    character(len=40),   dimension(48)       :: SymopSymb  ! Alphanumeric Symbols for SYMM
    !!----    type( Sym_Oper_Type),dimension(48)       :: SymOp      ! Crystallographic symmetry operators
    !!----    character(len=40),   dimension(48,8)     :: MSymopSymb ! Alphanumeric Symbols for MSYMM
    !!----    type(MSym_Oper_Type),dimension(48,8)     :: MSymOp     ! Magnetic symmetry operators
    !!----    Complex(kind=sp),    dimension(3,12,48,4):: basf       ! Basis functions of the irreps of Gk
    !!---- End Type MagSymm_k_Type
    !!----
    !!----  Definition of the MagSymm_k_type derived type, encapsulating the information
    !!----  concerning the crystallographic symmetry, propagation vectors and magnetic matrices.
    !!----  Needed for calculating magnetic structure factors.
    !!----
    !!---- Update: April - 2005
    !!
    Type, Public :: MagSymm_k_Type
       character(len=31)                        :: MagModel
       character(len=1)                         :: Latt
       integer                                  :: nirreps
       integer                                  :: nmsym
       integer                                  :: centred
       integer                                  :: mcentred
       integer                                  :: nkv
       real(kind=sp),dimension(3,12)            :: kvec
       integer                                  :: NumLat
       real(kind=sp), dimension(3,4)            :: Ltr
       integer                                  :: Numops
       integer                                  :: Multip
       integer,             dimension(4)        :: nbas
       integer,             dimension(12,4)     :: icomp
       character(len=40),   dimension(48)       :: SymopSymb
       type( Sym_Oper_Type),dimension(48)       :: SymOp
       character(len=40),   dimension(48,8)     :: MSymopSymb
       type(MSym_Oper_Type),dimension(48,8)     :: MSymOp
       Complex(kind=sp),    dimension(3,12,48,4):: basf
    End Type MagSymm_k_Type

    !!----
    !!---- ERR_MAGSYM
    !!----    logical, public :: err_MagSym
    !!----
    !!----    Logical Variable indicating an error in Magnetic_Symmetry
    !!----
    !!---- Update: April - 2005
    !!
    logical, public :: err_MagSym

    !!----
    !!---- ERR_MESS_MAGSYM
    !!----    character(len=150), public :: err_mess_MagSym
    !!----
    !!----    String containing information about the last error
    !!----
    !!---- Update: April - 2005
    !!
    character(len=150), public :: err_mess_MagSym

 Contains

    !-------------------!
    !---- Functions ----!
    !-------------------!

    !!----
    !!---- Function ApplyMso(Op,V) Result(Skp)
    !!----    Type(MSym_Oper_Type),         intent(in) :: Op        !  Magnetic Symmetry Operator Type
    !!----    real(kind=sp), dimension(3) , intent(in) :: Sk        !  Complex vector
    !!----    real(kind=sp), dimension(3)              :: Skp       !  Transformed complex vector
    !!----
    !!----    Apply a magnetic symmetry operator to a complex vector:  Skp = ApplyMSO(Op,Sk)
    !!----
    !!---- Update: April - 2005
    !!
    Function ApplyMSO(Op,Sk) Result(Skp)
       !---- Arguments ----!
       Type(MSym_Oper_Type), intent(in) :: Op
       Complex, dimension(3),intent(in) :: Sk
       Complex, dimension(3)            :: Skp

       Skp = matmul(Op%Rot,Sk) * cmplx(cos(tpi*Op%Phas),sin(tpi*Op%Phas))

       return
    End Function ApplyMSO


    !---------------------!
    !---- Subroutines ----!
    !---------------------!

    !!----
    !!---- Subroutine Init_Err_MagSym()
    !!----
    !!----    Initialize the errors flags in this Module
    !!----
    !!---- Update: March - 2005
    !!
    Subroutine Init_Err_MagSym()

       err_magsym=.false.
       err_mess_magsym=" "

       return
    End Subroutine Init_Err_MagSym

    !!--++
    !!--++ subroutine Init_MagSymm_k_Type(MGp)
    !!--++   type(MagSymm_k_Type),  intent (in out) :: MGp
    !!--++
    !!--++   (Private)
    !!--++   Subroutine to initialize the MagSymm_k_Type variable MGp.
    !!--++   It is called inside Readn_set_Magnetic_Structure
    !!--++
    !!--++  Update: April 2005
    !!
    Subroutine Init_MagSymm_k_Type(MGp)
       !---- Arguments ----!
       type(MagSymm_k_Type),  intent (in out) :: MGp

       !---- Local variables ----!
       integer :: i,j


       MGp%MagModel="Unnamed Model"
       MGp%Latt="P"
       MGp%nmsym=0
       MGp%nirreps=0
       MGp%centred=1
       MGp%mcentred=1
       MGp%nkv=0
       MGp%kvec=0.0
       MGp%NumLat=1
       MGp%Ltr=0.0
       MGp%Numops=0
       MGp%Multip=0
       MGp%nbas=0
       MGp%icomp=0
       MGp%basf=cmplx(0.0,0.0)
       do i=1,48
          MGp%SymopSymb(i)=" "
          MGp%SymOp(i)%Rot(:,:)=0
          MGp%SymOp(i)%tr(:)=0.0
          do j=1,8
             MGp%MSymopSymb(i,j)=" "
             MGp%MSymOp(i,j)%Rot(:,:)=0
             MGp%MSymOp(i,j)%Phas=0.0
          end do
       end do

       return
    End Subroutine Init_MagSymm_k_Type

    !!----
    !!---- Subroutine Readn_Set_Magnetic_Structure_k(file_cfl,n_ini,n_end,MGp,Am,SGo,Mag_dom)
    !!----    type(file_list_type),                intent (in)     :: file_cfl
    !!----    integer,                             intent (in out) :: n_ini, n_end
    !!----    type(MagSymm_k_Type),                intent (out)    :: MGp
    !!----    type(mAtom_List_Type),               intent (out)    :: Am
    !!----    type(Magnetic_Group_Type), optional, intent (out)    :: SGo
    !!----    type(Magnetic_Domain_type),optional, intent (out)    :: Mag_dom
    !!----
    !!----    Subroutine for reading and construct the MagSymm_k_Type variable MGp.
    !!----    It is supposed that the CFL file is included in the file_list_type
    !!----    variable file_cfl. On output n_ini, n_end hold the lines with the
    !!----    starting and ending lines with information about a magnetic phase.
    !!----    Optionally the Magnetig space group (Shubnikov group) may be obtained
    !!----    separately for further use.
    !!----    Magnetic S-domains are also read in case of providing the optional variable Mag_dom.
    !!----
    !!---- Update: November 2006
    !!
    Subroutine Readn_Set_Magnetic_Structure(file_cfl,n_ini,n_end,MGp,Am,SGo,Mag_dom)
       !---- Arguments ----!
       type(file_list_type),                intent (in)     :: file_cfl
       integer,                             intent (in out) :: n_ini, n_end
       type(MagSymm_k_Type),                intent (out)    :: MGp
       type(mAtom_List_Type),               intent (out)    :: Am
       type(Magnetic_Group_Type), optional, intent (out)    :: SGo
       type(Magnetic_Domain_type),optional, intent (out)    :: Mag_dom

       !---- Local Variables ----!
       integer :: i,no_iline,no_eline, num_k, num_xsym, num_irrep, num_dom, &
                  num_msym, ier, j, m, n, num_matom, num_skp, ik,im, ip
       real                 :: ph
       real,dimension(3)    :: rsk,isk
       real,dimension(3,12) :: br,bi
       complex, dimension(3):: Sk
       real,dimension(12)   :: coef
       character(len=132)   :: lowline,line
       character(len=30)    :: magmod, shubk
       character(len=2)     :: lattice
       character(len=4)     :: symbcar
       character(len=30)    :: msyr
       logical              :: msym_begin, kvect_begin, skp_begin, shub_given, irreps_given, &
                               irreps_begin, bfcoef_begin, magdom_begin
       type(Magnetic_Group_Type)  :: SG

       call init_err_MagSym()

       if(n_ini == 0) n_ini=1
       if(n_end == 0) n_end= file_cfl%nlines

       no_iline=0
       no_eline=0

       do i=n_ini,n_end
          ! Read comment
          if (index(file_cfl%line(i)(1:1),"!")/=0 .or. index(file_cfl%line(i)(1:1),"#")/=0) cycle
          lowline=adjustl(l_case(file_cfl%line(i)))

          if (lowline(1:13) == "mag_structure") then
             no_iline=i
          end if
          if (lowline(1:7) =="end_mag" ) then
             no_eline=i
             exit
          end if
       end do
       n_ini=no_iline
       n_end=no_eline

       if (n_ini == 0 .or. n_end == 0) then
          err_magsym=.true.
          err_mess_magsym=" No magnetig phase found in file!"
          return
       end if

       num_matom=0
       do i=n_ini,n_end
          lowline=l_case(adjustl(file_cfl%line(i)))
          if (index(lowline(1:5),"matom") ==0 ) cycle
          num_matom=num_matom+1
       end do

       Call Allocate_mAtom_list(num_matom,Am)  !Am contains Am%natoms = num_matom
       num_matom=0

       num_k=0
       num_dom=0
       num_xsym=0
       kvect_begin=.true.
       magdom_begin=.true.
       call Init_MagSymm_k_Type(MGp)
       i=n_ini
       shub_given  =.false.
       irreps_given=.false.
       irreps_begin=.false.
       msym_begin  =.false.
       skp_begin   =.false.
       bfcoef_begin=.false.
       if (present(mag_dom)) then  !Initialise Mag_dom
      	  Mag_dom%nd=0
          Mag_dom%Chir=.false.
          Mag_dom%DMat=0
          Mag_dom%pop=0.0
          Mag_dom%Lpop=0
          Mag_dom%Mpop=0.0
       end if

       do
          i=i+1
          if(i >= n_end) exit

          ! Read comment
          if( len_trim(file_cfl%line(i)) == 0) cycle
          lowline=adjustl(l_case(file_cfl%line(i)))
          if (lowline(1:1) == "!" .or. lowline(1:1)=="#") cycle

          ! Detect keywords

          ! Read magnetic model
          ! write(unit=*,fmt="(i6,a)") i,"  -> "//trim(lowline)
          if (lowline(1:6) == "magmod") then
             read(unit=lowline(8:),fmt=*,iostat=ier) magmod
             if (ier /= 0) then
                err_magsym=.true.
                err_mess_magsym=" Error reading magnetic model name in magnetic phase"
                return
             end if
             MGp%MagModel= adjustl(magmod)
             cycle
          end if

          ! Read lattice
          if (lowline(1:7) == "lattice") then
             read(unit=lowline(9:),fmt=*,iostat=ier) lattice
             if (ier /= 0) then
                err_magsym=.true.
                err_mess_magsym=" Error reading lattice type in magnetic phase"
                return
             end if
             lattice=adjustl(lattice)
             if (lattice(1:1)=="-") then
                MGp%centred = 2
                MGp%latt=u_case(lattice(2:2))
             else
                MGp%centred = 1
                MGp%Latt= u_case(lattice(1:1))
             end if
             cycle
          end if

          ! Read magnetic centrig
          if (lowline(1:7) == "magcent") then
             MGp%mcentred = 2
             cycle
          end if

          ! Read propagation vectors
          if (lowline(1:5) == "kvect" .and. kvect_begin) then
             num_k=num_k+1
             read(unit=lowline(6:),fmt=*,iostat=ier) MGp%kvec(:,num_k)
             if (ier /= 0) then
                err_magsym=.true.
                err_mess_magsym=" Error reading propagation vectors"
                return
             end if
             do !repeat reading until continuous KVECT lines are exhausted
                i=i+1
                lowline=adjustl(l_case(file_cfl%line(i)))
                ! write(unit=*,fmt="(i6,a)") i,"  -> "//trim(lowline)
                if (lowline(1:1) == "!" .or. lowline(1:1) == "#") cycle
                if (lowline(1:5) == "kvect") then
                   num_k=num_k+1
                   read(unit=lowline(6:),fmt=*,iostat=ier) MGp%kvec(:,num_k)
                   if (ier /= 0) then
                      err_magsym=.true.
                      err_mess_magsym=" Error reading propagation vectors"
                      return
                   end if
                else
                   i=i-1
                   kvect_begin=.false.
                   exit
                end if
             end do
             cycle
          end if

          ! Read magnetic S-domains
          if (present(mag_dom)) then
             if (lowline(1:6) == "magdom" .and. magdom_begin) then
                num_dom=num_dom+1
                ip=index(lowline,":")
                msyr=lowline(8:ip-1)
                call read_msymm(msyr,Mag_Dom%Dmat(:,:,num_dom),ph)
                if (ph > 0.001) Mag_Dom%chir=.true.
                if (Mag_Dom%chir) then
                   read(unit=lowline(ip+1:),fmt=*, iostat=ier) Mag_Dom%Pop(1:2,num_dom)!, Mag_Dom%MPop(1:2,num_dom)
                else
                   read(unit=lowline(ip+1:),fmt=*, iostat=ier) Mag_Dom%Pop(1,num_dom)  !, Mag_Dom%MPop(1,num_dom)
                end if
                if (ier /= 0) then
                   err_magsym=.true.
                   err_mess_magsym=" Error reading magnetic S-domains"
                   return
                end if
                Mag_Dom%nd = num_dom

                do  !repeat reading until continuous MAGDOM lines are exhausted
                   i=i+1
                   lowline=adjustl(l_case(file_cfl%line(i)))
                   ! write(unit=*,fmt="(i6,a)") i,"  -> "//trim(lowline)
                   if (lowline(1:1) == "!" .or. lowline(1:1) == "#") cycle
                   if (lowline(1:6) == "magdom") then
                      num_dom=num_dom+1
                      ip=index(lowline,":")
                      msyr=lowline(8:ip-1)
                      call read_msymm(msyr,Mag_Dom%Dmat(:,:,num_dom),ph)
                      if (ph > 0.001) Mag_Dom%chir=.true.
                      if (Mag_Dom%chir) then
                         read(unit=lowline(ip+1:),fmt=*, iostat=ier) Mag_Dom%Pop(1:2,num_dom) !, Mag_Dom%MPop(1:2,num_dom)
                      else
                         read(unit=lowline(ip+1:),fmt=*, iostat=ier) Mag_Dom%Pop(1,num_dom) !, Mag_Dom%MPop(1,num_dom)
                      end if
                      if (ier /= 0) then
                         err_magsym=.true.
                         err_mess_magsym=" Error reading magnetic S-domains"
                         return
                      end if
                      Mag_Dom%nd = num_dom
                   else
                      i=i-1
                      magdom_begin=.false.
                      exit
                   end if
                end do
                cycle
             end if
          end if

          ! Read number of irreducible representations and number of basis functions for each
          if (lowline(1:6) == "irreps") then
             read(unit=lowline(7:),fmt=*,iostat=ier) MGp%nirreps
             if (ier /= 0) then
                err_magsym=.true.
                err_mess_magsym=" Error reading number of irreducible representations"
                return
             end if
             read(unit=lowline(7:),fmt=*,iostat=ier) n, (MGp%nbas(j),j=1,MGp%nirreps)
             if (ier /= 0) then
                err_magsym=.true.
                err_mess_magsym=" Error reading number of basis functions of irreducible representations"
                return
             end if
             irreps_given=.true.
             cycle
          end if

          ! Read the indicator real(0)/imaginary(1) of coefficients for basis functions of
          ! irreducible representations
          if (lowline(1:5) == "icomp" .and. irreps_given) then
             num_irrep=1
             n=MGp%nbas(num_irrep)
             read(unit=lowline(6:),fmt=*,iostat=ier) MGp%icomp(1:abs(n),num_irrep)
             if (ier /= 0) then
                err_magsym=.true.
                err_mess_magsym=" Error reading real/imaginary indicators of BF coeff. of irreducible representations"
                return
             end if
             do  !repeat reading until continuous icoebf lines are exhausted
                i=i+1
                lowline=adjustl(l_case(file_cfl%line(i)))
                ! write(unit=*,fmt="(i6,a)") i,"  -> "//trim(lowline)
                if (lowline(1:1) == "!" .or. lowline(1:1) == "#") cycle
                if (lowline(1:5) == "icomp") then
                   num_irrep=num_irrep+1
                   n=MGp%nbas(num_irrep)
                   read(unit=lowline(6:),fmt=*,iostat=ier) MGp%icomp(1:abs(n),num_irrep)
                   if (ier /= 0) then
                      err_magsym=.true.
                      err_mess_magsym=" Error reading real/imaginary indicators of BF coeff. of irreducible representations"
                      return
                   end if
                else
                   i=i-1
                   irreps_given=.false.
                   exit
                end if
             end do
             cycle
          end if

          ! Read Shubnikov group
          if (lowline(1:9) == "shubnikov") then
             shubk=adjustl(file_cfl%line(i)(10:))
             Call Set_Shubnikov_Group(shubk,SG,MGp)
             if (err_magsym) return
             shub_given=.true.
          end if

          ! Read SYMM operators
          if (lowline(1:4) == "symm" .and. .not. shub_given) then
             num_xsym=num_xsym+1
             num_msym=0
             num_irrep=0
             read(unit=lowline(5:),fmt="(a)") MGp%SymopSymb(num_xsym)
             msym_begin=.true.
             irreps_begin=.true.
          end if

          ! Read MSYM operators
          if (lowline(1:4) == "msym" .and. msym_begin .and. .not. shub_given) then
             num_msym=num_msym+1
             read(unit=lowline(5:),fmt="(a)") MGp%MSymopSymb(num_xsym,num_msym)
             do  !repeat reading until continuous MSYM lines are exhausted
                i=i+1
                lowline=adjustl(l_case(file_cfl%line(i)))
                ! write(unit=*,fmt="(i6,a)") i,"  -> "//trim(lowline)
                if (lowline(1:1) == "!" .or. lowline(1:1) == "#") cycle
                if (lowline(1:4) == "msym") then
                   num_msym=num_msym+1
                   read(unit=lowline(5:),fmt="(a)") MGp%MSymopSymb(num_xsym,num_msym)
                else
                   i=i-1
                   msym_begin=.false.
                   exit
                end if
             end do
             cycle
          end if

          ! Read basis functions of irreducible representations
          if (lowline(1:4) == "basr" .and. irreps_begin .and. .not. shub_given) then
             num_irrep=num_irrep+1
             n=MGp%nbas(num_irrep)
             br=0.0; bi=0.0
             read(unit=lowline(5:),fmt=*,iostat=ier) (br(:,j),j=1,abs(n))
             if (ier /= 0) then
                err_magsym=.true.
                write(unit=err_mess_magsym,fmt="(2(a,i3))")" Error reading basis fuctions (BASR) of irrep ",num_irrep,&
                                                           " for symmetry operator # ",num_xsym
                return
             end if
             if (n < 0) then  !Read the imaginary part of the basis functions
                i=i+1
                lowline=adjustl(l_case(file_cfl%line(i)))
                !write(unit=*,fmt="(i6,a)") i,"  -> "//trim(lowline)
                if (lowline(1:4) == "basi") then
                   read(unit=lowline(5:),fmt=*,iostat=ier) (bi(:,j),j=1,abs(n))
                   if (ier /= 0) then
                      err_magsym=.true.
                      write(unit=err_mess_magsym,fmt="(2(a,i3))")" Error reading basis fuctions (BASI) of irrep ",num_irrep,&
                                                                 " for symmetry operator # ",num_xsym
                      return
                   end if
                else
                   err_magsym=.true.
                   write(unit=err_mess_magsym,fmt="(2(a,i3))")" Lacking BASI keyword of irrep ",num_irrep,&
                                                               " for symmetry operator # ",num_xsym
                   return
                end if
             end if
             do j=1,abs(n)
                MGp%basf(:,j,num_xsym,num_irrep)=cmplx( br(:,j),bi(:,j) )
             end do

             do  !repeat reading until continuous BASR or BASI lines are exhausted
                i=i+1
                lowline=adjustl(l_case(file_cfl%line(i)))
                !write(unit=*,fmt="(i6,a)") i,"  -> "//trim(lowline)
                if (lowline(1:1) == "!" .or. lowline(1:1) == "#") cycle
                if (lowline(1:4) == "basr") then
                   num_irrep=num_irrep+1
                   n=MGp%nbas(num_irrep)
                   br=0.0; bi=0.0
                   read(unit=lowline(5:),fmt=*,iostat=ier) (br(:,j),j=1,abs(n))
                   if (ier /= 0) then
                      err_magsym=.true.
                      write(unit=err_mess_magsym,fmt="(2(a,i3))")" Error reading basis fuctions (BASR) of irrep ",num_irrep,&
                                                                 " for symmetry operator # ",num_xsym
                      return
                   end if
                   if (n < 0) then  !Read the imaginary part of the basis functions
                      i=i+1
                      lowline=adjustl(l_case(file_cfl%line(i)))
                      !write(unit=*,fmt="(i6,a)") i,"  -> "//trim(lowline)
                      if (lowline(1:4) == "basi") then
                         read(unit=lowline(5:),fmt=*,iostat=ier) (bi(:,j),j=1,abs(n))
                         if (ier /= 0) then
                            err_magsym=.true.
                            write(unit=err_mess_magsym,fmt="(2(a,i3))")" Error reading basis fuctions (BASI) of irrep ",num_irrep,&
                                                                       " for symmetry operator # ",num_xsym
                            return
                         end if
                      else
                         err_magsym=.true.
                         write(unit=err_mess_magsym,fmt="(2(a,i3))")" Lacking BASI keyword of irrep ",num_irrep,&
                                                                    " for symmetry operator # ",num_xsym
                         return
                      end if
                   end if
                   do j=1,abs(n)
                      MGp%basf(:,j,num_xsym,num_irrep)=cmplx( br(:,j),bi(:,j) )
                   end do
                else
                   i=i-1
                   irreps_begin=.false.
                   exit
                end if
             end do
             cycle
          end if

          ! Read magnetic atoms:  label, magnetic form factor label,x,y,z,Biso,occ
          if (lowline(1:5) == "matom") then
             num_matom=num_matom+1
             num_skp=0
             line=adjustl(file_cfl%line(i))
             read(unit=line(6:),fmt=*,iostat=ier) Am%atom(num_matom)%lab,      & !Label
                                                  Am%atom(num_matom)%SfacSymb, & !Formfactor label
                                                  Am%atom(num_matom)%x,        & !Fract. coord.
                                                  Am%atom(num_matom)%Biso,     & !Is. Temp. Fact.
                                                  Am%atom(num_matom)%occ         !occupation
             if (ier /= 0) then
                err_magsym=.true.
                write(unit=err_mess_magsym,fmt="(a,i4)")" Error reading magnetic atom #",num_matom
                return
             end if
             skp_begin=.true.
             bfcoef_begin=.true.
             cycle
          end if

          ! Read Fourier coefficients in cryst. axes and phase
          if (lowline(1:3) == "skp" .and. skp_begin) then
             num_skp=num_skp+1
             read(unit=lowline(4:),fmt=*,iostat=ier) ik,im,rsk,isk,ph
             if (ier /= 0) then
                err_magsym=.true.
                write(unit=err_mess_magsym,fmt="(a,i3)") " Error reading Fourier Coefficient #", num_skp
                return
             end if
             Am%atom(num_matom)%nvk= num_skp
             Am%atom(num_matom)%imat(ik)= im
             Am%atom(num_matom)%Skr(:,ik)= rsk(:)
             Am%atom(num_matom)%Ski(:,ik)= isk(:)
             Am%atom(num_matom)%mphas(ik)= ph

             do  !repeat reading until continuous SPK lines are exhausted
                i=i+1
                lowline=adjustl(l_case(file_cfl%line(i)))
                !write(unit=*,fmt="(i6,a)") i,"  -> "//trim(lowline)
                if (lowline(1:1) == "!" .or. lowline(1:1) == "#") cycle
                if (lowline(1:3) == "skp") then
                   num_skp=num_skp+1
                   if (num_skp > 12) then
                      err_magsym=.true.
                      err_mess_magsym= " Too many Fourier Coefficients, the maximum allowed is 12! "
                      return
                   end if
                   read(unit=lowline(4:),fmt=*,iostat=ier) ik,im,rsk,isk,ph
                   if (ier /= 0) then
                      err_magsym=.true.
                      write(unit=err_mess_magsym,fmt="(a,i3)") " Error reading Fourier Coefficient #", num_skp
                      return
                   end if
                   Am%atom(num_matom)%nvk= num_skp
                   Am%atom(num_matom)%imat(ik)= im
                   Am%atom(num_matom)%Skr(:,ik)= rsk(:)
                   Am%atom(num_matom)%Ski(:,ik)= isk(:)
                   Am%atom(num_matom)%mphas(ik)= ph
                else
                   i=i-1
                   skp_begin=.false.
                   Am%atom(num_matom)%nvk= num_skp
                   exit
                end if
             end do
          end if

          if (lowline(1:6) == "bfcoef" .and. bfcoef_begin) then
             num_skp=num_skp+1
             read(unit=lowline(7:),fmt=*,iostat=ier) ik,im
             n=abs(MGp%nbas(im))
             read(unit=lowline(7:),fmt=*,iostat=ier) ik,im,coef(1:n),ph
             if (ier /= 0) then
                err_magsym=.true.
                write(unit=err_mess_magsym,fmt="(a,i3)") " Error reading Coefficient of Basis Functions #", num_skp
                return
             end if
             Am%atom(num_matom)%nvk= num_skp
             Am%atom(num_matom)%imat(ik)= im
             Am%atom(num_matom)%cbas(1:n,ik)= coef(1:n)
             Am%atom(num_matom)%mphas(ik)= ph

             do  !repeat reading until continuous bfcoef lines are exhausted
                i=i+1
                lowline=adjustl(l_case(file_cfl%line(i)))
                if (lowline(1:1) == "!" .or. lowline(1:1) == "#") cycle
                write(unit=*,fmt="(i6,a)") i,"  -> "//trim(lowline)
                if (lowline(1:6) == "bfcoef" ) then
                   num_skp=num_skp+1
                   if (num_skp > 12) then
                      err_magsym=.true.
                      err_mess_magsym= " Too many sets of Coefficients, the maximum allowed is 12! "
                      return
                   end if
                   read(unit=lowline(7:),fmt=*,iostat=ier) ik,im
                   n=abs(MGp%nbas(im))
                   read(unit=lowline(7:),fmt=*,iostat=ier) ik,im,coef(1:n),ph
                   if (ier /= 0) then
                      err_magsym=.true.
                      write(unit=err_mess_magsym,fmt="(a,i3)") " Error reading Coefficient of Basis Functions #", num_skp
                      return
                   end if
                   Am%atom(num_matom)%nvk= num_skp
                   Am%atom(num_matom)%imat(ik)= im
                   Am%atom(num_matom)%cbas(1:n,ik)= coef(1:n)
                   Am%atom(num_matom)%mphas(ik)= ph
                else
                   i=i-1
                   bfcoef_begin=.false.
                   Am%atom(num_matom)%nvk= num_skp
                   exit
                end if
             end do
          end if
       end do

       !Arriving here we have exhausted reading magnetic phase

       !Get pointers to the magnetic form factors
       !Stored for each atom in the component ind(1)
       call Set_Magnetic_Form()

       !---- Find Species in Magnetic_Form ----!
       do i=1,Am%natoms
          symbcar=u_case(Am%atom(i)%SfacSymb)
          do j=1,num_mag_form
             if (symbcar /= Magnetic_Form(j)%Symb) cycle
             Am%atom(i)%ind(1)=j
             exit
          end do
       end do

       !Now construct the rest of magnetic symmetry type variable MGp
       MGp%nmsym =num_msym
       MGp%Numops=num_xsym
       MGp%nkv   =num_k

       !Construct the numerical symmetry operators
       do i=1,MGp%Numops
          Call Read_Xsym(MGp%SymopSymb(i),1,MGp%Symop(i)%Rot,MGp%Symop(i)%tr)
          do j=1,MGp%nmsym
             Call Read_Msymm(MGp%MSymopSymb(i,j),MGp%MSymop(i,j)%Rot,MGp%MSymop(i,j)%Phas)
          end do
       end do
       if (err_symm) then
          err_magsym=.true.
          write(unit=err_mess_magsym,fmt="(a)") " Error reading symmetry: "//trim(err_mess_symm)
          return
       end if

       !Complete the set of symmetry operators with the centre of symmetry
       m=MGp%Numops
       if (MGp%centred == 2) then
          do i=1,MGp%Numops
             m=m+1
             MGp%Symop(m)%Rot(:,:) = -MGp%Symop(i)%Rot(:,:)
             MGp%Symop(m)%tr(:)    =  modulo_lat(-MGp%Symop(m)%tr(:))
             call Get_SymSymb(MGp%Symop(m)%Rot(:,:), &
                              MGp%Symop(m)%tr(:), MGp%SymopSymb(m))
             if (Mgp%mcentred == 1) then  !Anticentre in the magnetic structure
                do j=1,MGp%nmsym
                   MGp%MSymop(m,j)%Rot(:,:) = -MGp%MSymop(i,j)%Rot(:,:)
                   MGp%MSymop(m,j)%Phas     = -MGp%MSymop(i,j)%Phas
                end do
             else if(Mgp%mcentred == 2) then
                do j=1,MGp%nmsym
                   MGp%MSymop(m,j)%Rot(:,:) =  MGp%MSymop(i,j)%Rot(:,:)
                   MGp%MSymop(m,j)%Phas     =  MGp%MSymop(i,j)%Phas
                end do
             end if
          end do
       end if

       !Get the centring lattice translations of the crystallographic structure
       !and calculate the general multiplicity of the group.
       Mgp%NumLat=1
       MGp%Ltr(:,:) = 0.0
       Select Case(MGp%Latt)
          case ("A")
             Mgp%NumLat=2
             MGp%Ltr(:,1:2)=Ltr_a(:,1:2)
          case ("B")
             Mgp%NumLat=2
             MGp%Ltr(:,1:2)=Ltr_b(:,1:2)
          case ("C")
             Mgp%NumLat=2
             MGp%Ltr(:,1:2)=Ltr_c(:,1:2)
          case ("I")
             Mgp%NumLat=2
             MGp%Ltr(:,1:2)=Ltr_i(:,1:2)
          case ("R")
             Mgp%NumLat=3
             MGp%Ltr(:,1:3)=Ltr_r(:,1:3)
          case ("F")
             Mgp%NumLat=4
             MGp%Ltr(:,1:4)=Ltr_f(:,1:4)
       End Select

       select case (MGp%centred)
          case (1)
             MGp%Multip =   MGp%Numops * Mgp%NumLat
          case (2)
             MGp%Multip = 2 * MGp%Numops * Mgp%NumLat
       end select

       if (present(SGo)) then
          if (shub_given) then
             SGo=SG
          else
             err_magsym=.true.
             err_mess_magsym=" Shubnikov Group has not been provided "
          end if
       end if

       return
    End Subroutine Readn_Set_Magnetic_Structure

    !!----
    !!---- Subroutine Set_Shubnikov_Group(shubk,SG,MGp)
    !!----    character (len=*),         intent (in)    :: Shubk
    !!----    type(Magnetic_Group_Type), intent (out)   :: SG
    !!----    type(MagSymm_k_Type),      intent (in out):: MGp
    !!----
    !!---- Update: April 2008
    !!
    Subroutine Set_Shubnikov_Group(shubk,SG,MGp)
       !---- Arguments ----!
       character (len=*),         intent (in)    :: Shubk
       type(Magnetic_Group_Type), intent (out)   :: SG
       type(MagSymm_k_Type),      intent (in out):: MGp

       !---- Local Variables ----!
       character (len=132) :: line
       character (len=20)  :: symb
       character (len=4)   :: gn
       character (len=4),dimension(10) :: gen
       logical,          dimension(10) :: found
       integer :: i,j, ng, k,m,nbl,n
       integer,              dimension(3)   :: bl
       integer,              dimension(10)  :: syp, numop
       integer, allocatable, dimension(:,:) :: tab
       character(len=*),parameter, dimension(26) :: oper = &
       (/"1 ","-1","m ","2 ","21","3 ","31","32","-3","4 ","41","42","43",&
         "-4","6 ","61","62","63","64","65","-6","a ","b ","c ","d ","n "/)
       character(len=40),allocatable, dimension(:) :: ope


       SG%Shubnikov=" "
       SG%Shubnikov=adjustl(Shubk)
       gen = " "

       ! Generate the space group
       j=0
       bl=len_trim(SG%Shubnikov)
       numop=0
       do i=1,len_trim(SG%Shubnikov)
          if (SG%Shubnikov(i:i) == " ") then
             j=j+1
             bl(j)=i
          end if
       end do

       SG%Shubnikov(bl(1):) = l_case( Shubk(bl(1):))   !Ensures lower case for symmetry elements

       nbl=j
       j=0
       ng=0
       symb=" "
       syp=0
       do i=1,len_trim(SG%Shubnikov)
          j=j+1
          symb(j:j)=SG%Shubnikov(i:i)
          if (symb(j:j) == "'") then
             ng=ng+1
             k=5
             gn=" "
             do m=j-1,1,-1
                if (symb(m:m) == " ") exit
                if (symb(m:m) == "/") exit
                k=k-1
                gn(k:k)= symb(m:m)
             end do
             gen(ng)=adjustl(gn)
             if (i > bl(1)) syp(ng) = 1
             if (i > bl(2)) syp(ng) = 2
             if (i > bl(3)) syp(ng) = 3
             symb(j:j)=" "
             j= j-1
          end if
       end do
       i=index(symb," ")
       if ( i > 2) then
          symb=symb(1:1)//symb(i:)
       end if
       !write(*,*) " Space group symbol: ", trim(symb)
       !write(*,*) "  Primed Generators: ", (gen(i),i=1,ng), " in positions: ",(syp(i),i=1,ng)

       call Set_SpaceGroup(symb, SG%SpG)

                  !Determine the vector tinv from the information given for the generators
       SG%tinv=1  !by default the magnetic group is identical to the crystallographic group

       m=SG%SpG%Multip
       if (allocated(ope)) deallocate(ope)
       allocate(ope(m))

       found=.false.
       do j=1,ng
          if (gen(j) == "-1") then
             SG%tinv(SG%SpG%Numops+1) = -1
               found(SG%SpG%Numops+1) =.true.
             numop(j)= SG%SpG%Numops+1
          end if
       end do

       !         "Triclinic   ","Monoclinic  ","Orthorhombic","Tetragonal  ",    &
       !  "Rhombohedral","Hexagonal   ","Cubic       " /)
       n=1
       gn=SG%SpG%CrystalSys(1:4)
       do i=2,m  !over all symmetry operators of Space Group
          if(n == 0) exit  !all operators have been found
          call Symmetry_Symbol(SG%SpG%SymopSymb(i),ope(i))
          n=0
          do j=1,ng
             if (found(j)) cycle
             n=n+1
             if (gen(j)(1:1) == "-") then           !Search for roto-inversion axes
                k=index(ope(i),gen(j)(1:2)//"+")
                if (k /= 0) then    !Operator found
                   if (gn == "Cubi") then
                      k=index(ope(i),"x,x,x")
                      if (k /= 0) then
                         found(j)=.true.
                         SG%tinv(i)=-1
                         numop(j)= i
                         exit
                      end if
                   else
                      found(j)=.true.
                      SG%tinv(i)=-1
                      numop(j)= i
                      exit
                   end if
                end if
             else

                Select Case (gn)
                   case("Mono")
                      k=index(ope(i),gen(j)(1:1))    !Valid for all operators
                      if (k /= 0) then    !Operator found
                         found(j)=.true.
                         SG%tinv(i)=-1
                         numop(j)= i
                         exit
                      end if

                   case("Orth")
                      Select Case (gen(j))
                         Case("2 ","21")           ! Look for 2-fold axes
                            k=index(ope(i),gen(j)(1:1))
                            if (k /= 0) then
                               Select Case (syp(j))
                                  Case(1)
                                     k=index(ope(i),"x")
                                     if (k /= 0) then    !Operator found
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                                  Case(2)
                                     k=index(ope(i),"y")
                                     if (k /= 0) then    !Operator found
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                                  Case(3)
                                     k=index(ope(i),"z")
                                     if (k /= 0) then    !Operator found
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                               End Select
                            end if

                         Case("m ","a ","b ","c ","d ","n ")
                            k=index(ope(i),gen(j)(1:1))
                            if (k /= 0) then
                               Select Case (syp(j))
                                  Case(1)
                                     k=index(ope(i),"x")
                                     if (k == 0) then    !Operator found
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                                  Case(2)
                                     k=index(ope(i),"y")
                                     if (k == 0) then    !Operator found
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                                  Case(3)
                                     k=index(ope(i),"z")
                                     if (k == 0) then    !Operator found
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                               End Select
                            end if
                      End Select

                   case("Tetr")
                      Select Case (gen(j))
                         Case("4 ","41","42","43")           ! Look for 4-fold axes
                            k=index(ope(i),gen(j)(1:1)//"+")
                            if (k /= 0) then
                               found(j)=.true.
                               SG%tinv(i)=-1
                               numop(j)= i
                               exit
                            end if

                         Case("2 ","21")           ! Look for 2-fold axes
                            k=index(ope(i),gen(j)(1:1))
                            if (k /= 0) then
                               Select Case (syp(j))
                                  Case(2)                ! along [100]
                                     k=index(ope(i),"x")
                                     if (k /= 0) then    !Operator found
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                                  Case(3)                ! along [1-10]
                                     k=index(ope(i),"-x")
                                     if (k /= 0) then    !Operator found
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                               End Select
                            end if

                         Case("m ","a ","b ","c ","d ","n ")
                            k=index(ope(i),gen(j)(1:1))
                            if (k /= 0) then
                               Select Case (syp(j))
                                  Case(1)               ! Plane perp. to z (x,y,..) plane
                                     k=index(ope(i),"z")
                                     if (k == 0) then    !z should not be in the symbol
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                                  Case(2)                ! perp. to [100]
                                     k=index(ope(i),"x")
                                     if (k == 0) then    !x should not be in the symbol
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                                  Case(3)                ! perp. to [1-100]
                                     k=index(ope(i),"-x")
                                     if (k == 0) then    !-x should not be in the symbol
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                               End Select
                            end if
                      End Select

                   case("Rhom")
                      Select Case (gen(j))
                         Case("3 ","31","32")           ! Look for 3-fold axes
                            k=index(ope(i),gen(j)(1:1)//"+")
                            if (k /= 0) then
                               found(j)=.true.
                               SG%tinv(i)=-1
                               numop(j)= i
                               exit
                            end if

                         Case("2 ","21")           ! Look for 2-fold axes
                            k=index(ope(i),gen(j)(1:1))
                            if (k /= 0) then
                               Select Case (syp(j))
                                  Case(2)                ! along [100]
                                     k=index(ope(i),"x,0")
                                     if (k /= 0) then    !Operator found
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                                  Case(3)                ! along [1-10]
                                     k=index(ope(i),"-x")
                                     if (k /= 0) then    !Operator found
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                               End Select
                            end if

                         Case("m ","c ")
                            k=index(ope(i),gen(j)(1:1))
                            if (k /= 0) then
                               Select Case (syp(j))
                                  Case(1)               ! Plane perp. to z (x,y,..) plane
                                     k=index(ope(i),"z")
                                     if (k == 0) then    !z should not be in the symbol
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                                  Case(2)                ! perp. to [100]
                                     k=index(ope(i),"x")
                                     if (k == 0) then    !x should not be in the symbol
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                                  Case(3)                ! perp. to [1-100]
                                     k=index(ope(i),"-x")
                                     if (k == 0) then    !-x should not be in the symbol
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                               End Select
                            end if
                      End Select

                   case("Hexa")
                      Select Case (gen(j))
                         Case("6 ","61","62","63","64","65")    ! Look for 6-fold axes
                            k=index(ope(i),gen(j)(1:1)//"+")    !only along z
                            if (k /= 0) then
                               found(j)=.true.
                               SG%tinv(i)=-1
                               numop(j)= i
                               exit
                            end if

                         Case("2 ","21")           ! Look for 2-fold axes
                            k=index(ope(i),gen(j)(1:1))
                            if (k /= 0) then
                               Select Case (syp(j))
                                  Case(2)                ! along [100]
                                     k=index(ope(i),"x,0")
                                     if (k /= 0) then    !Operator found
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                                  Case(3)                ! along [1-10]
                                     k=index(ope(i),"-x")
                                     if (k /= 0) then    !Operator found
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                               End Select
                            end if

                         Case("m ","c ")
                            k=index(ope(i),gen(j)(1:1))
                            if (k /= 0) then
                               Select Case (syp(j))
                                  Case(1)               ! Plane perp. to z (x,y,..) plane
                                     k=index(ope(i),"z")
                                     if (k == 0) then    !z should not be in the symbol
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                                  Case(2)                ! perp. to [100]
                                     k=index(ope(i),"x")
                                     if (k == 0) then    !x should not be in the symbol
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                                  Case(3)                ! perp. to [1-100]
                                     k=index(ope(i),"-x")
                                     if (k == 0) then    !-x should not be in the symbol
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                               End Select
                            end if
                      End Select

                   case("Cubi")
                      Select Case (gen(j))
                         Case("4 ","41","42","43")    ! Look for 4-fold axes
                            k=index(ope(i),gen(j)(1:1)//"+")    !only along z
                            if (k /= 0) then
                               k=index(ope(i),"z")
                               if (k /= 0) then    !Operator found
                                  found(j)=.true.
                                  SG%tinv(i)=-1
                                  numop(j)= i
                                  exit
                               end if
                            end if

                         Case("3 ")    ! Look for 3-fold axes
                            k=index(ope(i),gen(j)(1:1)//"+")    !only along [111]
                            if (k /= 0) then
                               k=index(ope(i),"x,x,x")
                               if (k /= 0 ) then    !Operator found
                                  found(j)=.true.
                                  SG%tinv(i)=-1
                                  numop(j)= i
                                  exit
                               end if
                            end if

                         Case("2 ","21")           ! Look for 2-fold axes
                            k=index(ope(i),gen(j)(1:1))
                            if (k /= 0) then
                               Select Case (syp(j))
                                  Case(1)                ! along [001]
                                     k=index(ope(i),"z")
                                     if (k /= 0) then    !Operator found
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                                  Case(3)                ! along [110]
                                     k=index(ope(i),"-x")
                                     if (k /= 0) then    !Operator found
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                               End Select
                            end if

                         Case("m ","a ","b ","c ","d ","n ")
                            k=index(ope(i),gen(j)(1:1))
                            if (k /= 0) then
                               Select Case (syp(j))
                                  Case(1)               ! Plane perp. to z (x,y,..) plane
                                     k=index(ope(i),"x,y,0")
                                     if (k /= 0) then
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                                  Case(3)                ! perp. to [1-10]
                                     k=index(ope(i),"x,x,z")
                                     if (k /= 0) then
                                        found(j)=.true.
                                        SG%tinv(i)=-1
                                        numop(j)= i
                                        exit
                                     end if
                               End Select
                            end if
                     End Select

                End Select
             end if

          end do  !j=1,ng over all primed symmetry
       end do    !i=2,m over all symmetry operators of Space Group

       !write(*,*) "  Primed Generators: ", (gen(i),i=1,ng), " Correspond to operators: ",(numop(i),i=1,ng)

       ! if(allocated(tab)) deallocate(tab)
       ! allocate(tab(m,m))
       ! call  Set_SpG_Mult_Table(SG%SpG,tab,.true.)

       !Construct MGp from the Shubnikov group
       return
    End Subroutine Set_Shubnikov_Group

    !!----
    !!---- Subroutine Write_Magnetic_Structure(Ipr,MGp,Am,Mag_Dom)
    !!----    Integer,                    intent(in)           :: Ipr
    !!----    type(MagSymm_k_Type),       intent(in)           :: MGp
    !!----    type(mAtom_List_Type),      intent(in)           :: Am
    !!----    type(Magnetic_Domain_Type), intent(in), optional :: Mag_Dom
    !!----
    !!----    Subroutine to write out the information about the magnetic symmetry
    !!----    and mangnetic structure in unit Ipr.
    !!----
    !!---- Update: November 2006
    !!
    Subroutine Write_Magnetic_Structure(Ipr,MGp,Am,Mag_Dom)
       !---- Arguments ----!
       Integer,                    intent(in)           :: Ipr
       type(MagSymm_k_Type),       intent(in)           :: MGp
       type(mAtom_List_Type),      intent(in)           :: Am
       type(Magnetic_Domain_Type), intent(in), optional :: Mag_Dom

       !---- Local Variables ----!
       character (len=100), dimension( 4):: texto
       character (len=40)                :: aux
       integer :: i,j,k,l, nlines,n,m
       real                  :: x
       complex               :: ci
       real, dimension(3)    :: xp,xo
       complex, dimension(3) :: Sk


       Write(unit=ipr,fmt="(/,a)")  "==================================="
       Write(unit=ipr,fmt="(  a)")  "== Magnetic Symmetry Information =="
       Write(unit=ipr,fmt="(a,/)")  "==================================="

       write(unit=ipr,fmt="(a)")    " => Magnetic  model name: "//trim(MGp%MagModel)
       write(unit=ipr,fmt="(a)")    " => Crystal lattice type: "//MGp%Latt
       if (MGp%nirreps == 0) then
          write(unit=ipr,fmt="(a,i2)") " => Number of Magnetic operators/Crystallographic operator: ",MGp%nmsym
       else
          write(unit=ipr,fmt="(a,i2)") " => Number of Irreducible Representations: ",MGp%nirreps
          do i=1,MGp%nirreps
             write(unit=ipr,fmt="(2(a,i3),a,12i2)") " => Number of basis functions of Irreducible Representation #",i," :", &
                                 MGp%nbas(i),"  Indicators for real(0)/imaginary(1): ", MGp%icomp(1:abs(MGp%nbas(i)),i)
          end do
       end if
       if (MGp%Centred == 2) then
          write(unit=ipr,fmt="(a)")    " => The crystallographic structure is centric (-1 at origin) "
       else
          write(unit=ipr,fmt="(a)")    " => The crystallographic structure is acentric  "
       End if
       if (MGp%MCentred == 2) then
          write(unit=ipr,fmt="(a)")    " => The magnetic structure is centric "
       else
          if (MGp%Centred == 2) then
             write(unit=ipr,fmt="(a)")    " => The magnetic structure is anti-centric  "
          else
             write(unit=ipr,fmt="(a)")    " => The magnetic structure is acentric  "
          end if
       End if
       write(unit=ipr,fmt="(a,i2)") " => Number of propagation vectors: ",MGp%nkv
       do i=1,MGp%nkv
          write(unit=ipr,fmt="(a,i2,a,3f8.4,a)") " => Propagation vectors #",i," = (",MGp%kvec(:,i)," )"
       end do
       if (MGp%Numlat > 1) then
          write(unit=ipr,fmt="(a,i3)")  " => Centring vectors:",MGp%Numlat-1
          nlines=1
          texto(:) (1:100) = " "
          do i=2,MGp%Numlat
             call Frac_Trans_1Dig(MGp%Ltr(:,i),aux)
             if (mod(i-1,2) == 0) then
                write(unit=texto(nlines)(51:100),fmt="(a,i2,a,a)") " => Latt(",i-1,"): ",trim(aux)
                nlines=nlines+1
             else
                write(unit=texto(nlines)( 1:50),fmt="(a,i2,a,a)") " => Latt(",i-1,"): ",trim(aux)
             end if
          end do
          do i=1,nlines
             write(unit=ipr,fmt="(a)") texto(i)
          end do
       end if
       write(unit=ipr,fmt="(/,a,/)")        " => List of all Symmetry Operators and Symmetry Symbols"

       do i=1,MGp%Numops
          texto(1)=" "
          call Symmetry_Symbol(MGp%SymopSymb(i),texto(1))
          write(unit=ipr,fmt="(a,i3,2a,t50,2a)") " => SYMM(",i,"): ",trim(MGp%SymopSymb(i)), &
                                                          "Symbol: ",trim(texto(1))
          if (MGp%nirreps == 0) then
             do j=1,MGp%NMSym
                write(unit=ipr,fmt="(a,2(i2,a))")      "    MSYMM(",i,",",j,"): "//trim(MGp%MSymopSymb(i,j))
             end do
          else
             do j=1,MGp%nirreps
                write(unit=ipr,fmt="(a,2(i2,a),12(3f9.4,tr2))") "    BASR(",i,",",j,"): ",real(MGp%Basf(:,1:abs(MGp%nbas(j)),i,j))
                if (MGp%nbas(j) < 0) &
                write(unit=ipr,fmt="(a,2(i2,a),12(3f9.4,tr2))") "    BASI(",i,",",j,"): ",AIMAG(MGp%Basf(:,1:abs(MGp%nbas(j)),i,j))
             end do
          end if
       end do

       Write(unit=ipr,fmt="(/,a)")  "===================================="
       Write(unit=ipr,fmt="(  a)")  "== Magnetic Structure Information =="
       Write(unit=ipr,fmt="(a,/)")  "===================================="

       Write(unit=ipr,fmt="(a)")    " "
       Write(unit=ipr,fmt="(  a)")  "== Magnetic Asymmetric Unit Data =="
       Write(unit=ipr,fmt="(a,/)")  " "

       if (MGp%nirreps == 0) then
          Write(unit=ipr,fmt="(a)")  &
          "  The Fourier coefficients are of the form: Sk(j) = 1/2 { Rk(j) + i Ik(j) } exp {-2pi i Mphask(j)}"
          Write(unit=ipr,fmt="(a)")  &
          "  They are written for each atom j as Sk( j)= 1/2 {(Rx Ry Rz) + i ( Ix Iy Iz)} exp {-2pi i Mphask} -> MagMatrix # imat"
          Write(unit=ipr,fmt="(a)")  "  In case of k=2H (H reciprocal lattice vector) Sk(j)= (Rx Ry Rz)"

          do i=1,Am%Natoms
             Write(unit=ipr,fmt="(a,a,5f10.5)")  &
               "   Atom "//Am%Atom(i)%Lab, Am%Atom(i)%SfacSymb, Am%Atom(i)%x,Am%Atom(i)%Biso,Am%Atom(i)%occ
             do j=1,Am%Atom(i)%nvk
                if (K_Equiv_Minus_K(MGp%kvec(:,j),MGp%latt)) then
                   Write(unit=ipr,fmt="(a,i2,a,3f10.5,a,i4)")  &
                   "     Sk(",j,") =  (", Am%Atom(i)%Skr(:,j),")  -> MagMatrix #", Am%Atom(i)%imat(j)
                else
                   Write(unit=ipr,fmt="(a,i2,a,2(3f10.5,a),f9.5,a,i4)")  &
                   "     Sk(",j,") = 1/2 {(", Am%Atom(i)%Skr(:,j),") + i (",Am%Atom(i)%Ski(:,j),")}  exp { -2pi i ",&
                   Am%Atom(i)%MPhas(j),"}  -> MagMatrix #", Am%Atom(i)%imat(j)
                end if
             end do
          end do

       else

          Write(unit=ipr,fmt="(a)")  &
          "  The Fourier coefficients are of the form: Sk(j) = 1/2 Sum(i){Ci* Basf(i,imat)} exp {-2pi i Mphask(j)}"
          Write(unit=ipr,fmt="(a)")  &
          "  Where Ci are coefficients given below, Basf are the basis functions given above -> Irrep# imat"

          do i=1,Am%Natoms
             Write(unit=ipr,fmt="(a,a,5f10.5)")  &
               "   Atom "//Am%Atom(i)%Lab, Am%Atom(i)%SfacSymb, Am%Atom(i)%x,Am%Atom(i)%Biso,Am%Atom(i)%occ
             do j=1,Am%Atom(i)%nvk
                m=Am%Atom(i)%imat(j)
                n=abs(MGp%nbas(m))
                !1234567890123456789012345678
                aux="(a,i2,a,  f10.5,a,f9.5,a,i4)"
                write(unit=aux(9:10),fmt="(i2)") n
                Write(unit=ipr,fmt=aux)  &
                   "  Coef_BasF(",j,") = 1/2 {(", Am%Atom(i)%cbas(1:n,j),")}  exp { -2pi i ",&
                Am%Atom(i)%MPhas(j),"}  -> Irrep #", m
             end do
          end do
       end if

       ! Complete list of all atoms per primitive cell
       Write(unit=ipr,fmt="(/,a)")  " "
       Write(unit=ipr,fmt="(  a)")  "== List of all atoms and Fourier coefficients in the primitive cell =="
       Write(unit=ipr,fmt="(a,/)")  " "

       ! Construct the Fourier coefficients in case of Irreps
       if (MGp%nirreps /= 0 ) then
          do i=1,Am%natoms
             xo=Am%Atom(i)%x
             do k=1,MGp%NumOps
                xp=ApplySO(MGp%SymOp(k),xo)
                Write(unit=ipr,fmt="(a,i2,a,3f8.5)") " =>  Atom "//Am%Atom(i)%lab//"(",k,") :",xp
                do j=1,Am%Atom(i)%nvk
                   m=Am%Atom(i)%imat(j)
                   n=abs(MGp%nbas(m))
                   Sk(:) = cmplx(0.0,0.0)
                   do l=1,n !cannot be greater than 12 at present
                      x=real(MGp%icomp(l,m))
                      ci=cmplx(1.0-x,x)
                      Sk(:)=Sk(:)+ Am%atom(i)%cbas(l,m)*ci* MGp%basf(:,l,k,m)
                   end do
                   x=-tpi*Am%atom(i)%mphas(j)
                   Sk=Sk*cmplx(cos(x),sin(x))
                   Write(unit=ipr,fmt="(a,i2,a,2(3f10.5,a),f9.5,a)")  &
                    "     Sk(",j,") = 1/2 {(", real(Sk),") + i (",aimag(Sk),")}"
                end do
             end do  !Ops
             Write(unit=ipr,fmt="(a)") "  "
          end do  !atoms

       else
          do i=1,Am%natoms
             xo=Am%Atom(i)%x
             do k=1,MGp%NumOps
                xp=ApplySO(MGp%SymOp(k),xo)
                Write(unit=ipr,fmt="(a,i2,a,3f8.5)") " =>  Atom "//Am%Atom(i)%lab//"(",k,") :",xp
                do j=1,Am%Atom(i)%nvk
                   m=Am%Atom(i)%imat(j)
                   n=abs(MGp%nbas(m))
                   x=-tpi*Am%atom(i)%mphas(j)
                   Sk=cmplx(Am%Atom(i)%Skr(:,j),Am%Atom(i)%Ski(:,j))
                   Sk= ApplyMSO(MGp%MSymOp(k,j),Sk)*cmplx(cos(x),sin(x))
                   Write(unit=ipr,fmt="(a,i2,a,2(3f10.5,a),f9.5,a)")  &
                    "     Sk(",j,") = 1/2 {(", real(Sk),") + i (",aimag(Sk),")}"
                end do
             end do  !Ops
             Write(unit=ipr,fmt="(a)") "  "
          end do  !atoms
       end if

       ! Writing information about domains (like in FullProf)
       if (present(Mag_Dom)) then
          write(unit=ipr,fmt="(a)") " => Magnetic S-Domains are present"
          if (Mag_Dom%chir) write(unit=ipr,fmt="(a)")"    Chirality domains are also present                     Chir-1      Chir-2"
          do i=1,Mag_Dom%nd
             if (Mag_Dom%chir) then
                write(unit=ipr,fmt="(a,i2,1(a,2f12.4))")"      Matrix of Magnetic Domain #:",i, &
                   " -> Populations: ",Mag_Dom%Pop(1:2,i) !,'  Codes:',MagDom(iom)%MPop(1:2,i)
             else
                write(unit=ipr,fmt="(a,i2,1(a,f12.4))")"      Matrix of Magnetic Domain #:",i,  &
                   " -> Population: ",Mag_Dom%Pop(1,i) !,'  Code:',MagDom(iom)%MPop(1,i)
             end if
             do j=1,3
                write(unit=ipr,fmt="(a,3i4)")  "                    ",Mag_Dom%Dmat(j,:,i)
            end do
          end do
       end if

       return
    End Subroutine Write_Magnetic_Structure

    !!----
    !!---- Subroutine Write_Shubnikov_Group(SG,Iunit)
    !!----    type (Magnetic_Group_Type),intent(in) :: SG
    !!----    integer,   optional,       intent(in) :: iunit
    !!----
    !!----    Subroutine to write out the information about the Shubnikov_Group
    !!----
    !!---- Update: April 2008
    !!
    Subroutine Write_Shubnikov_Group(SG,Iunit)
       !---- Arguments ----!
       type (Magnetic_Group_Type),intent(in) :: SG
       integer,   optional,       intent(in) :: iunit

       !---- Local variables ----!
       character (len=100), dimension(24):: texto
       character (len=40)                :: aux
       integer                           :: lun
       integer                           :: i, nlines
       logical                           :: print_latt

       !---- Initializing variables ----!
       lun=6
       if (present(iunit)) lun=iunit
       print_latt=.true.

       !---- Printing ----!
       write(unit=lun,fmt="(/,/,a)")          "        Information on Space Group: "
       write(unit=lun,fmt="(a,/ )")           "        --------------------------- "
       write(unit=lun,fmt="(a,a )")          " =>       Shubnikov Symbol: ", SG%Shubnikov
       write(unit=lun,fmt="(a,i3)")          " =>  Number of Space group: ", SG%SpG%NumSpg
       write(unit=lun,fmt="(a,a)")           " => Hermann-Mauguin Symbol: ", SG%SpG%SPG_Symb
       write(unit=lun,fmt="(a,a)")           " =>            Hall Symbol: ", SG%SpG%Hall
       write(unit=lun,fmt="(a,a)")           " =>   Table Setting Choice: ", SG%SpG%info
       write(unit=lun,fmt="(a,a)")           " =>           Setting Type: ", SG%SpG%SG_setting
       write(unit=lun,fmt="(a,a)")           " =>         Crystal System: ", SG%SpG%CrystalSys
       write(unit=lun,fmt="(a,a)")           " =>             Laue Class: ", SG%SpG%Laue
       write(unit=lun,fmt="(a,a)")           " =>            Point Group: ", SG%SpG%Pg
       write(unit=lun,fmt="(a,a)")           " =>        Bravais Lattice: ", SG%SpG%SPG_Lat
       write(unit=lun,fmt="(a,a)")           " =>         Lattice Symbol: ", SG%SpG%SPG_Latsy
       write(unit=lun,fmt="(a,i3)")          " => Reduced Number of S.O.: ", SG%SpG%NumOps
       write(unit=lun,fmt="(a,i3)")          " =>   General multiplicity: ", SG%SpG%Multip
       write(unit=lun,fmt="(a,a)")           " =>         Centrosymmetry: ", SG%SpG%Centre
       write(unit=lun,fmt="(a,i3)")          " => Generators (exc. -1&L): ", SG%SpG%num_gen
       write(unit=lun,fmt="(a,f6.3,a,f6.3)") " =>        Asymmetric unit: ", SG%SpG%R_Asym_Unit(1,1), &
                                                                " <= x <= ", SG%SpG%R_Asym_Unit(1,2)
       write(unit=lun,fmt="(a,f6.3,a,f6.3)") "                            ", SG%SpG%R_Asym_Unit(2,1), &
                                                                " <= y <= ", SG%SpG%R_Asym_Unit(2,2)
       write(unit=lun,fmt="(a,f6.3,a,f6.3)") "                            ", SG%SpG%R_Asym_Unit(3,1), &
                                                                " <= z <= ", SG%SpG%R_Asym_Unit(3,2)

       if (SG%SpG%centred == 0) then
          call Frac_Trans_1Dig(SG%SpG%Centre_coord,texto(1))
          write(unit=lun,fmt="(a,a)")          " =>              Centre at: ", trim(texto(1))
       end if
       if (SG%SpG%SPG_Lat == "Z" .or. print_latt) then
          texto(:) (1:100) = " "
          if (SG%SpG%SPG_Lat == "Z") then
             write(unit=lun,fmt="(a,i3)")          " => Non-conventional Centring vectors:",SG%SpG%Numlat
          else
             write(unit=lun,fmt="(a,i3)")          " => Centring vectors:",SG%SpG%Numlat-1
          end if
          nlines=1
          do i=2,SG%SpG%Numlat
             call Frac_Trans_1Dig(SG%SpG%Latt_trans(:,i),aux)
             if (mod(i-1,2) == 0) then
                write(unit=texto(nlines)(51:100),fmt="(a,i2,a,a)") &
                                           " => Latt(",i-1,"): ",trim(aux)
                nlines=nlines+1
             else
                write(unit=texto(nlines)( 1:50),fmt="(a,i2,a,a)")  &
                                           " => Latt(",i-1,"): ",trim(aux)
             end if
          end do
          do i=1,nlines
             write(unit=lun,fmt="(a)") texto(i)
          end do
       end if

       !---- Symmetry Operators ----!
       write(unit=lun,fmt="(/,a,/)")        " => List of all Symmetry Operators and Symmetry Symbols"

       do i=1,SG%SpG%Multip
          texto(1)=" "
          call Symmetry_Symbol(SG%SpG%SymopSymb(i),texto(1))
          write(unit=lun,fmt="(a,i3,2a,t50,2a)") " => SYMM(",i,"): ",trim(SG%SpG%SymopSymb(i)), &
                                                    "Symbol: ",trim(texto(1))
       end do

       return
    End Subroutine Write_Shubnikov_Group

 End Module Magnetic_Symmetry

