module CFML_SpG_Standard_Representation

    !---- Used External Modules ----!
    use CFML_GlobalDeps,                only: cp
    use CFML_String_Utilities,          only: Get_Symb_From_Mat,Get_Separator_Pos,pack_string
    use CFML_Math_General,              only: Co_Linear,Rank,Determinant,Invert_Matrix,SmithNormalForm,&
                                              Equal_Matrix,Equal_Vector,RowEchelonForm, Trace,Zbelong
    use CFML_crystallographic_symmetry, only: Space_Group_Type,NS_Space_Group_Type,NS_Sym_Oper_Type,&
                                              Set_SpaceGroup,Setting_Change,Get_Transl_Symbol,&
                                              Sym_Oper_Type,Err_Symm,Err_Symm_Mess,Init_Err_Symm,&
                                              Get_SymSymb,Get_Setting_Info
     use CFML_Symmetry_Tables,           only: Spgr_Info,Set_Spgr_Info

    implicit none
    Private
    public:: Get_Standard_Representation

    !---- Global Variables ----!

    !---- Definitions ----!

    !!----
    !!---- ERR_STDREP_MESS
    !!----    character(len=256), public :: ERR_StdRep_Mess
    !!----
    !!----    String containing information about the last error
    !!----
    !!---- Update: May - 2018
    !!
    character(len=256), public :: Err_StdRep_Mess

    !!----
    !!---- ERR_STDREP
    !!----    logical, public :: ERR_StdRep
    !!----
    !!----    Logical Variable to indicate an error on this module.
    !!----
    !!---- Update: May - 2018
    !!
    logical, public :: Err_StdRep

    interface Get_Rotation_Id
        module procedure Get_Rotation_Id_Standard
        module procedure Get_Rotation_Id_Non_Standard
    end interface

contains

    subroutine Get_A_Matrices(LaueClass,A,n)

        !---- Arguments ----!
        character(len=5),                intent(in)  :: LaueClass
        real(kind=cp), dimension(3,3,6), intent(out) :: A
        integer,                         intent(out) :: n

        !---- Local variables ----!
        integer                         :: i
        logical                         :: singular
        integer,       dimension(3,3)   :: R1,R2,R3
        real(kind=cp), dimension(3,3)   :: B
        real(kind=cp), dimension(3,3,6) :: Ainv

        select case (trim(LaueClass))

            case ("2/m")
                n  = 6
                R1(1,1:3) = (/ 1, 0, 0 /)
                R1(2,1:3) = (/ 0, 1, 0 /)
                R1(3,1:3) = (/ 0, 0, 1 /)
                R2(1,1:3) = (/ 0, 0, 1 /)
                R2(2,1:3) = (/ 0,-1, 0 /)
                R2(3,1:3) = (/ 1, 0, 0 /)
                R3(1,1:3) = (/-1, 0, 1 /)
                R3(2,1:3) = (/ 0, 1, 0 /)
                R3(3,1:3) = (/-1, 0, 0 /)

                Ainv(:,:,1) = R1(:,:)
                Ainv(:,:,2) = R3(:,:)
                Ainv(:,:,3) = matmul(R3,R3)
                Ainv(:,:,4) = R2
                Ainv(:,:,5) = matmul(R2,R3)
                Ainv(:,:,6) = matmul(Ainv(:,:,5),R3)
            case ("mmm")
                n  = 6
                R1(1,1:3) = (/ 1, 0, 0 /)
                R1(2,1:3) = (/ 0, 1, 0 /)
                R1(3,1:3) = (/ 0, 0, 1 /)
                R2(1,1:3) = (/ 0, 1, 0 /)
                R2(2,1:3) = (/ 1, 0, 0 /)
                R2(3,1:3) = (/ 0, 0,-1 /)
                R3(1,1:3) = (/ 0, 0, 1 /)
                R3(2,1:3) = (/ 1, 0, 0 /)
                R3(3,1:3) = (/ 0, 1, 0 /)

                Ainv(:,:,1) = R1(:,:)
                Ainv(:,:,2) = R3(:,:)
                Ainv(:,:,3) = matmul(R3,R3)
                Ainv(:,:,4) = R2
                Ainv(:,:,5) = matmul(R2,R3)
                Ainv(:,:,6) = matmul(Ainv(:,:,5),R3)
            case ("m-3")
                n = 2
                R1(1,1:3) = (/ 1, 0, 0 /)
                R1(2,1:3) = (/ 0, 1, 0 /)
                R1(3,1:3) = (/ 0, 0, 1 /)
                R2(1,1:3) = (/ 0,-1, 0 /)
                R2(2,1:3) = (/ 1, 0, 0 /)
                R2(3,1:3) = (/ 0, 0, 1 /)
                Ainv(:,:,1) = R1(:,:)
                Ainv(:,:,2) = R2(:,:)
            case default
                n = 1
                R1(1,1:3) = (/ 1, 0, 0 /)
                R1(2,1:3) = (/ 0, 1, 0 /)
                R1(3,1:3) = (/ 0, 0, 1 /)
                Ainv(:,:,1) = R1(:,:)
        end select

        do i = 1 , n
            call Invert_Matrix(Ainv(:,:,i),A(:,:,i),singular)
        end do

    end subroutine Get_A_Matrices

    subroutine Get_Generators(spaceGroupNumber,symOp,nSymOp,G,nGen)

        ! Stores in Gt and Gx the generators for the spacegroups
        ! SGtarget and SGaux. ng is the number of generators (three
        ! as maximum)

        integer,                                intent(in)  :: spaceGroupNumber
        type(Sym_Oper_Type), dimension(nSymOp), intent(in)  :: symOp     ! symmetry operations
        integer,                                intent(in)  :: nSymOp    ! number of symmetry operations
        real(kind=cp),       dimension(4,4,3),  intent(out) :: G         ! generators
        integer,                                intent(out) :: nGen      ! number of generators

        integer                              :: inversion,ng_,n,i,j,ngaux
        real(kind=cp)                        :: d
        logical                              :: matched,positive
        integer,       dimension(3)          :: axis
        integer,       dimension(3,3)        :: identity
        real(kind=cp), dimension(4,4)        :: W
        integer, dimension(:,:), allocatable :: idd

        ERR_StdRep = .false.
        if (spaceGroupNumber < 1 .or. spaceGroupNumber > 230) then
            ERR_StdRep = .true.
            ERR_StdRep_Mess = "Wrong space group number"
            return
        end if

        ! Initialization
        nGen          = 0
        inversion     = 0
        identity(1,:) = (/ 1,0,0 /)
        identity(2,:) = (/ 0,1,0 /)
        identity(3,:) = (/ 0,0,1 /)
        allocate(idd(nSymOp,2))

        ! Look for an inversion center
        do i = 1 , nSymOp
            if (Equal_Matrix(-symOp(i)%Rot,identity,3)) then
                inversion = i
                exit
            end if
        end do

        if (spaceGroupNumber < 3) then ! Triclinic

            if (inversion == 0) then
                nGen         = 1
                ! Search for the onefold axis
                call Get_Rotation_Id(symOp(:),nSymOp,1,n,idd)
                call Get_Rotation_Id(symOp(:),nSymOp,2,n,idd)
                G(1:3,1:3,1) = symOp(idd(1,1))%Rot(1:3,1:3)
                G(1:3,4,1)   = symOp(idd(1,1))%Tr(1:3)
                G(4,1:4,1)   = (/ 0.,0.,0.,1. /)
            end if

        else if (spaceGroupNumber < 16)  then ! Monoclinic

            nGen = 1
            ! Search for a twofold axis
            call Get_Rotation_Id(symOp(:),nSymOp,2,n,idd)
            G(1:3,1:3,1) = symOp(idd(1,1))%Rot(1:3,1:3)
            G(1:3,4,1)   = symOp(idd(1,1))%Tr(1:3)
            G(4,1:4,1)   = (/ 0.,0.,0.,1. /)
            ! Choose proper rotations if the spacegroup is centrosymmetric
            if (inversion > 0) G(1:3,1:3,1) = idd(1,2) * G(1:3,1:3,1)

        else if (spaceGroupNumber < 75)  then ! Orthorhombic

            nGen = 2
            ! Search for the two fold axes along [001] and [010]
            call Get_Rotation_Id(symOP(:),nSymOp,2,n,idd)
            ngaux = 0
            do i = 1 , n
                call Get_Rotation_Axis(symOp(idd(i,1))%Rot,axis)
                if ((axis(1) == 0 .and. axis(2) == 0 .and. axis(3) == 1) .or. &
                    (axis(1) == 0 .and. axis(2) == 1 .and. axis(3) == 0)) then
                    ngaux = ngaux + 1
                    G(1:3,1:3,ngaux) = symOp(idd(i,1))%Rot(1:3,1:3)
                    G(1:3,4,ngaux)   = symOp(idd(i,1))%Tr(1:3)
                    G(4,1:4,ngaux)   = (/ 0.,0.,0.,1. /)
                    if (ngaux == 2) exit
                end if
            end do

        else if (spaceGroupNumber < 143) then ! Tetragonal

            ! Search for the fourfold axis along [001]
            call Get_Rotation_Id(symOp(:),nSymOp,4,n,idd)
            do i = 1 , n
                call Get_Rotation_Axis(symOp(idd(i,1))%Rot,axis)
                if (axis(1) == 0 .and. axis(2) == 0 .and. axis(3) == 1) then
                    call Positive_Sense_of_Rotation(symOp(idd(i,1))%Rot,axis,positive)
                    if (positive) then
                        G(1:3,1:3,1) = symOp(idd(i,1))%Rot(1:3,1:3)
                        G(1:3,4,1)   = symOp(idd(i,1))%Tr(1:3)
                        G(4,1:4,1)   = (/ 0.,0.,0.,1. /)
                        nGen         = 1
                        exit
                    end if
                end if
            end do
            ! Look for a possible twofold axis along [100]
            call Get_Rotation_Id(symOp(:),nSymOp,2,n,idd)
            do i = 1 , n
                call Get_Rotation_Axis(symOp(idd(i,1))%Rot,axis)
                if (axis(1) == 1 .and. axis(2) == 0 .and. axis(3) == 0) then
                    G(1:3,1:3,2) = symOp(idd(i,1))%Rot(1:3,1:3)
                    G(1:3,4,2)   = symOp(idd(i,1))%Tr(1:3)
                    G(4,1:4,2)   = (/ 0.,0.,0.,1. /)
                    nGen         = 2
                    exit
                end if
            end do

        else if (spaceGroupNumber < 168) then ! Trigonal

            ! Search for the threefold axis along [001]
            call Get_Rotation_Id(symOP(:),nSymOp,3,n,idd)
            do i = 1 , n
                call Get_Rotation_Axis(symOp(idd(i,1))%Rot,axis)
                if (axis(1) == 0 .and. axis(2) == 0 .and. axis(3) == 1) then
                    call Positive_Sense_of_Rotation(symOp(idd(i,1))%Rot,axis,positive)
                    if (positive) then
                        G(1:3,1:3,1) = symOp(idd(i,1))%Rot(1:3,1:3)
                        G(1:3,4,1)   = symOp(idd(i,1))%Tr(1:3)
                        G(4,1:4,1)   = (/ 0.,0.,0.,1. /)
                        nGen         = 1
                        exit
                    end if
                end if
            end do
            ! Search for a possible twofold axis along [110] or [-110]
            call Get_Rotation_Id(symOp(:),nSymOp,2,n,idd)
            do i = 1 , n
                call Get_Rotation_Axis(symOp(idd(i,1))%Rot,axis)
                if ((axis(1) == 1 .and. axis(2) == 1 .and. axis(3) == 0) .or. &
                    (axis(1) ==-1 .and. axis(2) == 1 .and. axis(3) == 0)) then
                    G(1:3,1:3,2) = symOp(idd(i,1))%Rot(1:3,1:3)
                    G(1:3,4,2)   = symOp(idd(i,1))%Tr(1:3)
                    G(4,1:4,2)   = (/ 0.,0.,0.,1. /)
                    nGen         = 2
                    exit
                end if
            end do

        else if (spaceGroupNumber < 195) then ! Hexagonal

            ! Search for the sixfold axis along [001] in SGtarget
            call Get_Rotation_Id(symOp(:),nSymOp,6,n,idd)
            do i = 1 , n
                call Get_Rotation_Axis(symOp(idd(i,1))%Rot,axis)
                if (axis(1) == 0 .and. axis(2) == 0 .and. axis(3) == 1) then
                    call Positive_Sense_of_Rotation(symOp(idd(i,1))%Rot,axis,positive)
                    if (positive) then
                        G(1:3,1:3,1) = symOp(idd(i,1))%Rot(1:3,1:3)
                        G(1:3,4,1)   = symOp(idd(i,1))%Tr(1:3)
                        G(4,1:4,1)   = (/ 0.,0.,0.,1. /)
                        nGen         = 1
                        exit
                    end if
                end if
            end do
            ! Look for a possible twofold axis along [-110] in SGtarget
            call Get_Rotation_Id(symOp(:),nSymOp,2,n,idd)
            do i = 1 , n
                call Get_Rotation_Axis(symOp(idd(i,1))%Rot,axis)
                if (axis(1) ==-1 .and. axis(2) == 1 .and. axis(3) == 0) then
                    G(1:3,1:3,2) = symOp(idd(i,1))%Rot(1:3,1:3)
                    G(1:3,4,2)   = symOp(idd(i,1))%Tr(1:3)
                    G(4,1:4,2)   = (/ 0.,0.,0.,1. /)
                    nGen         = 2
                    exit
                end if
            end do

        else

            ! Search for the fourfold axis along [001] in SGtarget
            call Get_Rotation_Id(symOp(:),nSymOp,4,n,idd)
            do i = 1 , n
                call Get_Rotation_Axis(symOp(idd(i,1))%Rot,axis)
                if (axis(1) == 0 .and. axis(2) == 0 .and. axis(3) == 1) then
                    call Positive_Sense_of_Rotation(symOp(idd(i,1))%Rot,axis,positive)
                    if (positive) then
                        G(1:3,1:3,1) = symOp(idd(i,1))%Rot(1:3,1:3)
                        G(1:3,4,1)   = symOp(idd(i,1))%Tr(1:3)
                        G(4,1:4,1)   = (/ 0.,0.,0.,1. /)
                        nGen         = 1
                        exit
                    end if
                end if
            end do
            if (nGen == 0) then
                ! Search for the twofold axis along [001] in SGtarget
                call Get_Rotation_Id(symOp(:),nSymOp,2,n,idd)
                do i = 1 , n
                    call Get_Rotation_Axis(symOp(idd(i,1))%Rot,axis)
                    if (axis(1) == 0 .and. axis(2) == 0 .and. axis(3) == 1) then
                        G(1:3,1:3,1) = symOp(idd(i,1))%Rot(1:3,1:3)
                        G(1:3,4,1)   = symOp(idd(i,1))%Tr(1:3)
                        G(4,1:4,1)   = (/ 0.,0.,0.,1. /)
                        nGen         = 1
                        exit
                    end if
                end do
            end if
            ! Search for a threefold axis along {111} in SGtarget
            call Get_Rotation_Id(symOp(:),nSymOp,3,n,idd)
            do i = 1 , n
                call Get_Rotation_Axis(symOp(idd(i,1))%Rot,axis)
                if (axis(1) == 1 .and. axis(2) == 1 .and. axis(3) == 1) then
                    call Positive_Sense_of_Rotation(symOp(idd(i,1))%Rot,axis,positive)
                    if (positive) then
                        G(1:3,1:3,2) = symOp(idd(i,1))%Rot(1:3,1:3)
                        G(1:3,4,2)   = symOp(idd(i,1))%Tr(1:3)
                        G(4,1:4,2)   = (/ 0.,0.,0.,1. /)
                        nGen         = 2
                        exit
                    end if
                end if
            end do

        end if

        if (inversion > 0) then
            ! Choose proper rotations if the spacegroup is centrosymmetric
            do i = 1 , nGen
                call Determinant(real(G(1:3,1:3,i)),3,d)
                G(1:3,1:3,i) = nint(d) * G(1:3,1:3,i)
            end do
            nGen = nGen + 1
            G(1,1:3,nGen)   = (/ -1.,0.,0. /)
            G(2,1:3,nGen)   = (/ 0.,-1.,0. /)
            G(3,1:3,nGen)   = (/ 0.,0.,-1. /)
            G(1:3,4,nGen)   = symOp(inversion)%Tr(:)
            G(4,1:4,nGen)   = (/ 0.,0.,0.,1. /)
        end if

    end subroutine Get_Generators

    subroutine Get_HM_Standard(sgNumber,symbolHM)

        ! Returns the Herman-Maugin symbol for the standard.
        implicit none

        integer,           intent(in)  :: sgNumber
        character(len=12), intent(out) :: symbolHM

        integer i,n
        integer, dimension(1) :: posSep

        if (sgNumber < 0 .or. sgNumber > 230) then
            ERR_StdRep = .true.
            ERR_StdRep_Mess = "Error in Get_HM_Standard subroutine: sgNumber out of range"
            return
        end if

        call Set_Spgr_Info()
        i = 1
        do
            if (spgr_info(i)%n == sgNumber) then
                call Get_Separator_Pos(spgr_info(i)%HM,':',posSep,n)
                if (n == 0) then
                    symbolHM = spgr_info(i)%HM
                else ! Origin choice 2
                    symbolHM = spgr_info(i+1)%HM
                end if
                exit
            end if
            i = i + 1
        end do

    end subroutine Get_HM_Standard

    subroutine Get_Lattice_Type(L, Latc, lattyp)

       !---- Arguments ----!
       integer,                        intent( in) :: L
       real(kind=cp), dimension(:,:),  intent( in) :: Latc
       character(len=*),               intent(out) :: lattyp

       !---- Local variables ----!
       logical :: latt_p, latt_a, latt_b, latt_c, latt_i, latt_r, latt_s, latt_h, latt_f, latt_z
       integer, dimension(10):: latt_given
       integer, dimension(3) :: tt
       integer               :: i, j
       integer, dimension(3,10), parameter :: lattice=reshape((/0,6,6, 6,0,6, &
                                                     6,6,0, 6,6,6, 8,4,4, 4,8,8, &
                                                     4,8,4, 8,4,8, 8,4,0, 4,8,0/),(/3,10/))

       if (l > 3) then  !non conventional centring
          lattyp="Z"
          return
       else if(l == 0) then !primitive lattice
          lattyp="P"
          return
       end if

       latt_p=.true.
       latt_a=.false.
       latt_b=.false.
       latt_c=.false.
       latt_i=.false.
       latt_r=.false.
       latt_s=.false. ! Rombohedral,reverse setting
       latt_h=.false. ! Hexagonally centred
       latt_f=.false.
       latt_z=.false.

       do i=1,L
          tt(1:3)=nint(12.0 * Latc(1:3,i))   ! Translations x 12

          !---- Compare the translation part of the operator with tabulated array ----!
          latt_given(:) = 0
          do j=1,10
             if (equal_vector(tt,lattice(:,j),3)) then
                latt_given(j) = 1
                select case (j)
                   case (1)
                      latt_a=.true.
                   case (2)
                      latt_b=.true.
                   case (3)
                      latt_c=.true.
                   case (4)
                      latt_i=.true.
                   case (5,6)
                      latt_r=.true.
                   case (7,8)
                      latt_s=.true.
                   case (9,10)
                      latt_h=.true.
                end select
                exit
             end if
          end do
          if (sum(latt_given) == 0) then
             latt_z = .true.
             exit
          end if
       end do

       !---- Lattice Type ----!
       if (latt_z) then
           lattyp="Z"
           return
       end if
       if ( (latt_a .and. latt_b .and. latt_c) .or. (latt_a .and. latt_b) .or. &
            (latt_a .and. latt_c) .or. (latt_b .and. latt_c) ) then
            latt_f=.true.
            latt_a=.false.
            latt_b=.false.
            latt_c=.false.
            latt_p=.false.
            latt_i=.false.
       end if
       if (latt_p) lattyp="P"
       if (latt_a) lattyp="A"
       if (latt_b) lattyp="B"
       if (latt_c) lattyp="C"
       if (latt_i) lattyp="I"
       if (latt_r) lattyp="R"
       if (latt_s) lattyp="S"
       if (latt_h) lattyp="H"
       if (latt_f) lattyp="F"

       return
    end subroutine Get_Lattice_Type

    subroutine Get_Lattice_Type_from_M(M,lattyp,info)

        ! M matrix transforms a primitive setting to a centred
        ! setting (standard basis). The columns of the inverse of
        ! the M matrix contain the primitive vectors of the lattice
        ! expressed in the standard basis

        !---- Arguments ----!
        real(kind=cp), dimension(3,3), intent(in)  :: M
        character,                     intent(out) :: lattyp
        logical, optional,             intent(in)  :: info

        !---- Local variables ----!
        integer   :: i,j,nLatticePoints,nCentringVectors
        real      :: d
        logical   :: singular,newt,output=.false.
        real(kind=cp), dimension(3)    :: t
        logical,       dimension(3)    :: isOrigin
        real(kind=cp), dimension(3,3)  :: Minv,Maux
        real(kind=cp), dimension(3,11) :: latc

        output = .false.
        if (present(info)) output = .true.
        call Invert_Matrix(M,Minv,singular)
        call Determinant(Minv,3,d)

        nLatticePoints = Nint(abs(1/d))

        if ( nLatticePoints > 1 ) then
            if (output) write(*,'(8x,a,i2,1x,a)') " => Standard lattice is non primitive.",&
                                        nLatticePoints,"lattice points in the unit cell"
            nCentringVectors = 0
            Maux             = Minv ! Maux stores lattice points inside the unit cell
            isOrigin         = .True.
            do i = 1 , 3
                do j = 1 , 3
                    if ( abs(Maux(i,j) - Nint(Maux(i,j))) < 0.001 ) then ! integer component
                        Maux(i,j) = 0.0
                    else                                                 ! fractional component
                        Maux(i,j) = Maux(i,j) - Int(Maux(i,j))
                        if (Maux(i,j) < 0.0) Maux(i,j) = Maux(i,j) + 1.0
                        isOrigin(j) = .False.
                    end if
                end do
            end do
            do i = 1 , 3
                if ( .not. isOrigin(i) ) then
                    t = Maux(:,i) ! t is a centring vector
                    do j = 1 , 3
                       if (t(j) < 0.0) t(j) = t(j) + 1.0
                    end do
                    newt = .true.
                    do j = 1 , nCentringVectors
                       if (Co_Linear(t,latc(:,j),3)) newt = .false.
                    end do
                    if (newt) then
                       nCentringVectors = nCentringVectors + 1
                       latc(:,nCentringVectors) = t
                    end if
                end if
            end do
            call Get_Lattice_Type(nCentringVectors,latc,lattyp)
            if (output) write(*,'(12x,2a)') "Lattice type: ",lattyp
        else
            if (output) write(*,'(8x,a)') " => Standard lattice is primitive "
            lattyp = "P"
        end if

    end subroutine Get_Lattice_Type_from_M

    subroutine Get_Laue_Class_from_Point_Group(pointGroup,laueClass)
        character(len=5), intent(in)  :: pointGroup
        character(len=5), intent(out) :: laueClass

        laueClass = " "

        select case (trim(pointGroup))

            case ("1","-1")
                laueClass = "-1"
            case ("2","m","2/m")
                laueClass = "2/m"
            case ("222","mm2","mmm")
                laueClass = "mmm"
            case ("4","-4","4/m")
                laueClass = "4/m"
            case ("422","4mm","-4m2","4/mmm")
                laueClass = "4/mmm"
            case ("3","-3")
                laueClass = "-3"
            case ("32","3m","-3m")
                laueClass = "-3m"
            case ("6","-6","6/m")
                laueClass = "6/m"
            case ("622","6mm","-6m2","6/mmm")
                laueClass = "6/mmm"
            case ("23","m-3")
                laueClass = "m-3"
            case ("432","-43m","m-3m")
                laueClass = "m-3m"
        end select

    end subroutine Get_Laue_Class_from_Point_Group

    subroutine Get_Mc_Matrix(LaueClass,Mp,Mc,info)

        ! Compute a correction matrix if necessary
        ! The correction needed in some cases for cubic
        ! groups with primitive lattices has been introduced
        ! in the subroutine Get_A_Matrices, for the case m-3

        !---- Arguments ----!
        character(len=5),              intent(in)  :: LaueClass
        real(kind=cp), dimension(3,3), intent(in)  :: Mp
        real(kind=cp), dimension(3,3), intent(out) :: Mc
        logical, optional,             intent(in)  :: info ! if present, write output

        !---- Local variables ----!
        integer           :: i
        character         :: lattyp
        character(len=60) :: symb
        logical           :: singular,output
        real(kind=cp), dimension(3,3) :: Mcinv

        output = .false.
        if (present(info)) output = .true.
        if (output) then
            write(*,'(8x,a)') " => Constructing (Mc,0) matrix...."
            write(*,'(8x,a)') "    This matrix corrects (M',0) if (M',0) does not &
                                   bring the system to a standard setting"
        end if

        select case (trim(LaueClass))

            case ("2/m")  ! Put the two fold axis along b

                Mcinv(1,:) = (/ 0.,1.,0. /)
                Mcinv(2,:) = (/ 0.,0.,1. /)
                Mcinv(3,:) = (/ 1.,0.,0. /)
                call Invert_Matrix(Mcinv,Mc,singular)

            case ("4/m","4/mmm")

                call Get_Lattice_Type_from_M(Mp,lattyp)
                if (lattyp == "C") then  ! C -> P
                    Mcinv(1,:) = (/ 1., 1., 0. /)
                    Mcinv(2,:) = (/ 1.,-1., 0. /)
                    Mcinv(3,:) = (/ 0., 0.,-1. /)
                    call Invert_Matrix(Mcinv,Mc,singular)
                else if (lattyp == "F") then ! F -> I
                    Mcinv(1,:) = (/ 1., 1., 0. /)
                    Mcinv(2,:) = (/-1., 1., 0. /)
                    Mcinv(3,:) = (/ 0., 0., 1. /)
                    call Invert_Matrix(Mcinv,Mc,singular)
                else
                    Mc(1,:) = (/ 1.,0.,0. /)
                    Mc(2,:) = (/ 0.,1.,0. /)
                    Mc(3,:) = (/ 0.,0.,1. /)
                end if

            case ("-3","-3 R")

                call Get_Lattice_Type_from_M(Mp,lattyp)
                if (lattyp == "S") then ! reverse -> obverse setting
                    Mc(1,:) = (/-1., 0., 0. /)
                    Mc(2,:) = (/ 0.,-1., 0. /)
                    Mc(3,:) = (/ 0., 0., 1. /)
                else
                    Mc(1,:) = (/ 1.,0.,0. /)
                    Mc(2,:) = (/ 0.,1.,0. /)
                    Mc(3,:) = (/ 0.,0.,1. /)
                end if

            case ("-3m","-3m R","-3m1","-31m")
                call Get_Lattice_Type_from_M(Mp,lattyp)
                if (lattyp == "S") then ! reverse -> obverse setting
                    Mc(1,:) = (/-1., 0., 0. /)
                    Mc(2,:) = (/ 0.,-1., 0. /)
                    Mc(3,:) = (/ 0., 0., 1. /)
                else if (lattyp == "H") then ! H -> P
                    Mcinv(1,:) = (/ 1., 1., 0. /)
                    Mcinv(2,:) = (/-1., 2., 0. /)
                    Mcinv(3,:) = (/ 0., 0., 1. /)
                    call Invert_Matrix(Mcinv,Mc,singular)
                else
                    Mc(1,:) = (/ 1.,0.,0. /)
                    Mc(2,:) = (/ 0.,1.,0. /)
                    Mc(3,:) = (/ 0.,0.,1. /)
                end if

            case ("6/mmm")

                call Get_Lattice_Type_from_M(Mp,lattyp)
                if (lattyp == "H") then ! H -> P
                    Mcinv(1,:) = (/ 1., 1., 0. /)
                    Mcinv(2,:) = (/-1., 2., 0. /)
                    Mcinv(3,:) = (/ 0., 0., 1. /)
                    call Invert_Matrix(Mcinv,Mc,singular)
                else
                    Mc(1,:) = (/ 1.,0.,0. /)
                    Mc(2,:) = (/ 0.,1.,0. /)
                    Mc(3,:) = (/ 0.,0.,1. /)
                end if

            case default

                Mc(1,:) = (/ 1.,0.,0. /)
                Mc(2,:) = (/ 0.,1.,0. /)
                Mc(3,:) = (/ 0.,0.,1. /)

        end select

        if (output) then
            write(*,'(12x,a)',advance='no') "Quasi-standard setting --> Standard setting transformation: "
            call Get_Symb_From_Mat(transpose(Mc(1:3,1:3)),symb,(/"a","b","c"/))
            write(*,'(a)') trim(symb)
        end if

    end subroutine Get_Mc_Matrix

    subroutine Get_Mprime_Matrix(SG,P,Mp,info)

        ! Given a space group SG and a matrix P -see Get_P_Matrix
        ! subroutine-, computes the matrix Mp to generate a
        ! representation where basis vectors are along the symmetry
        ! axes of the Laue class of the SG. In the notation used in
        ! this subroutine:
        !
        !           [ bx1 by1 bz1 ]
        !      Mp = [ bx2 by2 bz2 ]
        !           [ bx3 by3 bz3 ]

        !---- Arguments ----!
        type(NS_Space_Group_Type),     intent(in)  :: SG    ! space group
        real(kind=cp), dimension(3,3), intent(in)  :: P     ! P-matrix
        real(kind=cp), dimension(3,3), intent(out) :: Mp    ! Mprime matrix
        logical, optional,             intent(in)  :: info  ! if present, write output

        !---- Local variables ----!
        integer                           :: i,j,k,n,n_,tr,nCubicAxes
        real                              :: d
        logical                           :: colin,singular,integral,standard,output
        integer,           dimension(2)   :: order
        integer,           dimension(3)   :: bx,by,bz
        real(kind=cp),     dimension(3)   :: v
        character(len=20), dimension(6)   :: axisName
        character(len=60)                 :: symb
        integer,           dimension(3,3) :: W,U
        integer,           dimension(3,4) :: paxis,cubicAxes
        real(kind=cp),     dimension(3,3) :: Pinv,Waux,Uaux,PM,PMinv
        integer,           dimension(:,:), allocatable :: idd
        real(kind=cp),     dimension(:,:), allocatable :: T

        output = .false.
        if (present(info)) output = .true.
        axisName(2) = "twofold"
        axisName(3) = "threefold"
        axisName(4) = "fourfold"
        axisName(6) = "sixfold"
        call Invert_Matrix(P,Pinv,singular)

        if (output) then
            write(*,'(8x,a)') " => Constructing (M',0) matrix...."
            write(*,'(12x,a)') "This matrix transforms the primitive basis in a basis with &
            vectors along the symmetry axes of the Laue class"
            write(*,'(12x,2a)') "Laue Class: ", SG%Laue
        end if

        select case (trim(SG%Laue))
            case ("-1")
                Mp = reshape ( (/1.0,0.0,0.0, &
                                 0.0,1.0,0.0, &
                                 0.0,0.0,1.0/),(/3,3/) )
            case ("2/m","4/m","-3","-3 R","6/m")
                select case (trim(SG%Laue))
                    case ("2/m")
                        order(1) = 2
                    case ("4/m")
                        order(1) = 4
                    case ("-3","-3 R","6/m")
                        order(1) = 3
                end select
                if (output) write(*,'(12x,3a)') "Searching for the ",trim(axisName(order(1)))," axis..."
                allocate(idd(SG%Multip,2))
                call Get_Rotation_Id(SG%SymOP(:),SG%Multip,order(1),n,idd)
                if (Err_Symm) then
                    Err_StdRep = .true.
                    Err_StdRep_Mess = Err_Symm_Mess
                    return
                end if
                ! Put the symmetry operation in the primitive setting
                Waux = MatMul(SG%SymOP(idd(1,1))%Rot,P)
                Waux = MatMul(Pinv,Waux)
                W = NInt(Waux) ! In the primitive setting W is always integral
                               ! This has been checked in Get_P_Matrix, so
                               ! doing W = NInt(Waux) is safe
                call Get_Rotation_Axis(W,bz)
                if (Err_Symm) then
                    Err_StdRep = .true.
                    Err_StdRep_Mess = Err_Symm_Mess
                    return
                end if
                if (output) then
                    write(*,'(16x,a,i2,",",i2,",",i2,"]")') "Axis direction in the primitive setting: [", bz
                    write(*,'(12x,a)') "Searching for lattice vectors perpendicular to the rotation axis. &
                    Building the complete basis..."
                end if
                call Get_Vectors_Perpendicular_To_Rotation_Axis(W,idd(1,2),order(1),paxis)
                if (Err_Symm) then
                    Err_StdRep = .true.
                    Err_StdRep_Mess = Err_Symm_Mess
                    return
                end if
                call Get_Pseudo_Standard_Base(W,idd(1,2),order(1),paxis,bz,bx,by)
                Mp(:,1) = real(bx)
                Mp(:,2) = real(by)
                Mp(:,3) = real(bz)
                deallocate(idd)
            case ("4/mmm")
                if (output) write(*,'(12x,a)') "Searching for the fourfold axis..."
                allocate(idd(SG%Multip,2))
                call Get_Rotation_Id(SG%SymOP(:),SG%Multip,4,n,idd)
                if (Err_Symm) then
                    Err_StdRep = .true.
                    Err_StdRep_Mess = Err_Symm_Mess
                    return
                end if
                ! Put the symmetry operation in the primitive setting
                Waux = MatMul(SG%SymOP(idd(1,1))%Rot,P)
                Waux = MatMul(Pinv,Waux)
                W = NInt(Waux) ! In the primitive setting W is always integral
                               ! This has been checked in Get_P_Matrix, so
                               ! doing W = NInt(Waux) is safe
                call Get_Rotation_Axis(W,bz)
                if (Err_Symm) then
                    Err_StdRep = .true.
                    Err_StdRep_Mess = Err_Symm_Mess
                    return
                end if
                if (output) then
                    write(*,'(16x,a,i2,",",i2,",",i2,"]")') "Axis direction in the primitive setting: [", bz
                    write(*,'(12x,a)') "Searching for the twofold axis..."
                end if
                call Get_Rotation_Id(SG%SymOP(:),SG%Multip,2,n,idd)
                if (Err_Symm) then
                    Err_StdRep = .true.
                    Err_StdRep_Mess = Err_Symm_Mess
                    return
                end if
                colin = .true.
                do i = 1 , n
                    ! Put the symmetry operation in the primitive setting
                    Uaux = MatMul(SG%SymOP(idd(i,1))%Rot,P)
                    Uaux = MatMul(Pinv,Uaux)
                    U = NInt(Uaux) ! In the primitive setting W is always integral
                                   ! This has been checked in Get_P_Matrix, so
                                   ! doing W = NInt(Waux) is safe
                    call Get_Rotation_Axis(U,bx)
                    if (Err_Symm) then
                        Err_StdRep = .true.
                        Err_StdRep_Mess = Err_Symm_Mess
                        return
                    end if
                    if (.not. Co_Linear(bz,bx,3)) then
                        colin = .false.
                        exit
                    end if
                end do
                if (colin) then
                    ERR_StdRep = .true.
                    ERR_StdRep_Mess = "Error in Get_Mprime_Matrix. Unable to find a second rotation &
                    axis linearly independent from the first"
                    return
                end if
                if (output) write(*,'(16x,a,i2,",",i2,",",i2,"]")') "Axis direction in the primitive setting: [", bx
                by = matmul(W,bx)
                Mp(:,1) = real(bx)
                Mp(:,2) = real(by)
                Mp(:,3) = real(bz)
                deallocate(idd)
            case ("-3m","-3m R","-3m1","-31m","6/mmm")
                if (output) write(*,'(12x,a)') "Searching for the threefold axis..."
                allocate(idd(SG%Multip,2))
                call Get_Rotation_Id(SG%SymOP(:),SG%Multip,3,n,idd)
                if (Err_Symm) then
                    Err_StdRep = .true.
                    Err_StdRep_Mess = Err_Symm_Mess
                    return
                end if
                ! Put the symmetry operation in the primitive setting
                Waux = MatMul(SG%SymOP(idd(1,1))%Rot,P)
                Waux = MatMul(Pinv,Waux)
                W = NInt(Waux) ! In the primitive setting W is always integral
                               ! This has been checked in Get_P_Matrix, so
                               ! doing W = NInt(Waux) is safe
                if (idd(1,2) == -1)  W = -W
                call Get_Rotation_Axis(W,bz)
                if (Err_Symm) then
                    Err_StdRep = .true.
                    Err_StdRep_Mess = Err_Symm_Mess
                    return
                end if
                if (output) then
                    write(*,'(16x,a,i2,",",i2,",",i2,"]")') "Axis direction in the primitive setting: [", bz
                    write(*,'(12x,a)') "Searching for the twofold axis..."
                end if
                call Get_Rotation_Id(SG%SymOP(:),SG%Multip,2,n,idd)
                if (Err_Symm) then
                    Err_StdRep = .true.
                    Err_StdRep_Mess = Err_Symm_Mess
                    return
                end if
                colin = .true.
                do i = 1 , n
                    ! Put the symmetry operation in the primitive setting
                    Uaux = MatMul(SG%SymOP(idd(i,1))%Rot,P)
                    Uaux = MatMul(Pinv,Uaux)
                    U = NInt(Uaux) ! In the primitive setting W is always integral
                                   ! This has been checked in Get_P_Matrix, so
                                   ! doing W = NInt(Waux) is safe
                    call Get_Rotation_Axis(U,bx)
                    if (Err_Symm) then
                        Err_StdRep = .true.
                        Err_StdRep_Mess = Err_Symm_Mess
                        return
                    end if
                    if (.not. Co_Linear(bz,bx,3)) then
                        colin = .false.
                        exit
                    end if
                end do
                if (colin) then
                    ERR_StdRep = .true.
                    ERR_StdRep_Mess = "Error in Get_Mprime_Matrix. Unable to find a second rotation &
                    axis linearly independent from the first"
                    return
                end if
                if (output) write(*,'(16x,a,i2,",",i2,",",i2,"]")') "Axis direction in the primitive setting: [", bx
                by = matmul(W,bx)
                Mp(:,1) = real(bx)
                Mp(:,2) = real(by)
                Mp(:,3) = real(bz)
                deallocate(idd)
            case ("mmm")
                if (output) write(*,'(12x,a)') "Searching for the three twofold axis..."
                allocate(idd(SG%Multip,2))
                call Get_Rotation_Id(SG%SymOP(:),SG%Multip,2,n,idd)
                if (Err_Symm) then
                    Err_StdRep = .true.
                    Err_StdRep_Mess = Err_Symm_Mess
                    return
                end if
                n_ = 0
                do i = 1 , n
                    ! Put the symmetry operation in the primitive setting
                    Waux = MatMul(SG%SymOP(idd(i,1))%Rot,P)
                    Waux = MatMul(Pinv,Waux)
                    W = NInt(Waux) ! In the primitive setting W is always integral
                                   ! This has been checked in Get_P_Matrix, so
                                   ! doing W = NInt(Waux) is safe
                    call Get_Rotation_Axis(W,bz)
                    if (Err_Symm) then
                        Err_StdRep = .true.
                        Err_StdRep_Mess = Err_Symm_Mess
                        return
                    end if
                    colin = .false.
                    do j = 1 , n_
                        if (Co_Linear(bz,Nint(Mp(:,j)),3)) colin = .true.
                    end do
                    if (.not. colin) then
                        n_ = n_ + 1
                        Mp(:,n_) = real(bz)
                        if (output) write(*,'(16x,a,i2,",",i2,",",i2,"]")') "Axis direction in the primitive setting: [", bz
                    end if
                    if (n_ == 3) exit
                end do
                deallocate(idd)
            case ("m3","m-3","m3m","m-3m")
                if (output) write(*,'(12x,a)') "Searching for the four threefold axes..."
                allocate(idd(SG%Multip,2))
                call Get_Rotation_Id(SG%SymOP(:),SG%Multip,3,n,idd)
                if (Err_Symm) then
                    Err_StdRep = .true.
                    Err_StdRep_Mess = Err_Symm_Mess
                    return
                end if
                nCubicAxes = 0
                standard   = .true. ! Initially we assume P is a standard basis
                do i = 1 , n
                    Waux = MatMul(SG%SymOP(idd(i,1))%Rot,P)
                    Waux = MatMul(Pinv,Waux)
                    W = NInt(Waux) ! In the primitive setting W is always integral
                                   ! This has been checked in Get_P_Matrix, so
                                   ! doing W = NInt(Waux) is safe
                    call Get_Rotation_Axis(W,bz)
                    if (Err_Symm) then
                        Err_StdRep = .true.
                        Err_StdRep_Mess = Err_Symm_Mess
                        return
                    end if
                    colin = .false.
                    do j = 1 , nCubicAxes
                        if (Co_Linear(bz,cubicAxes(:,j),3)) colin = .true.
                    end do
                    if (.not. colin) then
                        nCubicAxes = nCubicAxes + 1
                        cubicAxes(:,nCubicAxes) = bz(:)
                        if (output) write(*,'(16x,a,i2,",",i2,",",i2,"]")') "Axis direction in the primitive setting: [", bz
                        if (abs(bz(1)) /= 1 .or. &
                            abs(bz(2)) /= 1 .or. &
                            abs(bz(3)) /= 1) standard = .false.
                    end if
                end do
                if (standard) then
                    Mp(:,1) = (/ 1.,0.,0. /)
                    Mp(:,2) = (/ 0.,1.,0. /)
                    Mp(:,3) = (/ 0.,0.,1. /)
                else
                    ! Find a combination of cubicAxes for which the threefold axes are along {111}
                    do i = 0 , 1
                        if (standard) exit
                        bx(:) = cubicAxes(:,2) - 2 * i * cubicAxes(:,2)
                        do j = 0 , 1
                            if (standard) exit
                            by(:) = cubicAxes(:,3) - 2 * j * cubicAxes(:,3)
                            do k = 0 , 1
                                if (standard) exit
                                bz(:) = cubicAxes(:,4) - 2 * k * cubicAxes(:,4)
                                Mp(:,1) = 0.5 * (cubicAxes(:,1) + bx)
                                Mp(:,2) = 0.5 * (cubicAxes(:,1) + by)
                                Mp(:,3) = 0.5 * (cubicAxes(:,1) + bz)
                                ! For I lattices, Mp is not integral (there is a lattice
                                ! point at the middle of each diagonal. We must multiply
                                ! by two in order to get the diagonal of the cube
                                integral= Zbelong(Mp)
                                if (.not. integral) Mp = 2 * Mp
                                PM = matmul(P,Mp)
                                call Invert_Matrix(PM,PMinv,singular)
                                if (.not. singular) then
                                    standard = .true.
                                    do n_ = 1 , n
                                        Waux = MatMul(SG%SymOP(idd(n_,1))%Rot,PM)
                                        Waux = MatMul(PMinv,Waux)
                                        if (Zbelong(Waux)) then
                                            W = NInt(Waux)
                                            call Get_Rotation_Axis(W,bz)
                                            if (abs(bz(1)) /= 1 .or. &
                                                abs(bz(2)) /= 1 .or. &
                                                abs(bz(3)) /= 1) then
                                                standard = .false.
                                                exit
                                            end if
                                        end if
                                    end do
                                end if
                            end do
                        end do
                    end do
                end if
                if (.not. standard) then
                    Err_StdRep = .true.
                    Err_StdRep_Mess = "Error in Get_Mprime_Matrix: &
                    Unable to build a standard setting for the cubic crystal system..."
                    return
                end if
        end select

        call Set_Right_Handedness(Mp)

        if (output) then
            write(*,'(12x,a)',advance='no') "Primitive setting --> Quasi-Standard setting transformation: "
            call Get_Symb_From_Mat(transpose(Mp(1:3,1:3)),symb,(/"a","b","c"/))
            write(*,'(a)') trim(symb)
        end if

    end subroutine Get_Mprime_Matrix

   subroutine Get_Origin_Shift(G,G_,ng,P,oShift,getShift)

        implicit none

        real(kind=cp), dimension(4,4,ng), intent(in)  :: G,G_
        integer,                          intent(in)  :: ng
        real(kind=cp), dimension(3,3),    intent(in)  :: P
        real(kind=cp), dimension(3),      intent(out) :: oShift
        logical,                          intent(out) :: getShift

        integer :: i,j,k,l,r,nr
        logical :: singular
        integer,       dimension(3,3)    :: identity
        real(kind=cp), dimension(4,4)    :: Pt,Ptinv
        real(kind=cp), dimension(4,4,ng) :: Gx,Gt
        integer,       dimension(:,:), allocatable :: U,D,T,V
        real(kind=cp), dimension(:,:), allocatable :: b,x

        getShift      = .true.
        Pt(1:3,1:3)   = P
        Pt(1:3,4)     = (/ 0.,0.,0. /)
        Pt(4,1:4)     = (/ 0.,0.,0.,1. /)
        identity(1,:) = (/ 1,0,0 /)
        identity(2,:) = (/ 0,1,0 /)
        identity(3,:) = (/ 0,0,1 /)
        Gx            = G(:,:,1:ng)
        Gt            = G_(:,:,1:ng)

        call Invert_Matrix(Pt,Ptinv,singular)

        ! Transform generators to a primitive setting
        do i = 1 , ng
            Gt(:,:,i) = matmul(Gt(:,:,i),Pt)
            Gt(:,:,i) = matmul(Ptinv,Gt(:,:,i))
            Gx(:,:,i) = matmul(Gx(:,:,i),Pt)
            Gx(:,:,i) = matmul(Ptinv,Gx(:,:,i))
        end do
        ! Build the equation to compute the origin shift
        !      U x cp = b (mod Z)
        if (.not. allocated(U)) then
            nr = 3 * ng
            allocate(U(nr,3))
            allocate(b(nr,1))
            allocate(x(3,1))
            allocate(D(nr,3))
            allocate(T(nr,nr))
            allocate(V(3,3))
        end if
        r = 0
        do j = 1 , ng
            do k = 1 , 3
                r = r + 1
                b(r,1) = Gt(k,4,j) - Gx(k,4,j)
                do l = 1 , 3
                    U(r,l) = Nint(Gt(k,l,j)) - identity(k,l)
                end do
            end do
        end do
        call SmithNormalForm(U,D,T,V)
        b = matmul(T,b)
        do j = 1 , 3
            if (D(j,j) == 0) then
                if (abs(b(j,1) - Nint(b(j,1))) > 0.0001) then
                    getShift = .false.
                    return
                else
                    x(j,1) = 0.
                end if
            else
                x(j,1) = b(j,1) / D(j,j)
            end if
        end do
        do j = 4 , nr
            if (abs(b(j,1) - Nint(b(j,1))) > 0.0001) then
                getShift = .false.
                return
            end if
        end do
        if (getShift) then
            x = matmul(V,x)
            x = matmul(Pt(1:3,1:3),x)
            oshift = x(:,1)
        end if

    end subroutine Get_Origin_Shift

    subroutine Get_P_Matrix(SG,P,info)

        ! If the spacegroup is not given in a primitive basis,
        ! computes a P matrix to change it to a primitive basis

        !---- Arguments ----!
        type(NS_Space_Group_Type),     intent(in)  :: SG   ! space group
        real(kind=cp), dimension(3,3), intent(out) :: P    ! P-matrix
        logical, optional,             intent(in)  :: info ! if present, write output

        !---- Local variables ----!
        integer                       :: n,i,j,k,l,m
        integer                       :: nCentringVectors
        real(kind=cp)                 :: d
        logical                       :: linearDependence,sorted,&
                                         singular,primitive,integral,&
                                         output
        character(len=60)             :: symb,tr_symb
        real(kind=cp), dimension(4)   :: v
        real(kind=cp), dimension(3,3) :: A,Pinv
        real(kind=cp), dimension(:,:), allocatable :: centringVectors

        output = .false.
        if (present(info)) output = .true.
        if (output) then
            write(*,'(8x,a)') " => Constructing (P,0) matrix...."
            write(*,'(12x,a)') "This matrix transforms the original basis in a primitive basis"
        end if
        P    = 0.
        Pinv = 0.
        primitive = .false.
        if (SG%NumLat == 1) then
            ! Basis is already primitive
            if (output) write(*,'(12x,a)') "The basis is already primitive"
            primitive = .true.
            do i = 1 , 3
                P(i,i)    = 1.
                Pinv(i,i) = 1.
            end do
            return
        else
            ! Build an expanded list of centring vectors
            Allocate(centringVectors(4,8*SG%NumLat))
            nCentringVectors = 0
            do n = 2 , SG%NumLat
                do i = 0 , 1
                    if (i == 1 .and. SG%Latt_Trans(1,n) == 0) cycle
                    v(1) = SG%Latt_Trans(1,n) - i
                    do j = 0 , 1
                        if (j == 1 .and. SG%Latt_Trans(2,n) == 0) cycle
                        v(2) = SG%Latt_Trans(2,n) - j
                        do k = 0 , 1
                            if (k == 1 .and. SG%Latt_Trans(3,n) == 0) cycle
                            v(3) = SG%Latt_Trans(3,n) - k
                            v(4) = dot_product(v(1:3),v(1:3))
                            ! Check if there is linear dependence with previous vectors
                            linearDependence = .False.
                            do m = 1 , nCentringVectors
                                if (Co_Linear(v(1:3),centringVectors(1:3,m),3)) Then
                                    linearDependence = .True.
                                    if (v(4) < centringVectors(4,m)) centringVectors(:,m) = v
                                    exit
                                end if
                            end do
                            if (.not. linearDependence) then
                                nCentringVectors = nCentringVectors + 1
                                centringVectors(:,nCentringVectors) = v
                            end if
                        end do
                    end do
                end do
            end do
            ! Sort vectors from shorter to longer
            do
                sorted = .True.
                do i = 1 , nCentringVectors - 1
                    if (centringVectors(4,i) > centringVectors(4,i+1)) Then
                        v = centringVectors(:,i)
                        centringVectors(:,i)   = centringVectors(:,i+1)
                        centringVectors(:,i+1) = v
                        sorted = .False.
                    end if
                end do
                if (sorted) exit
            end do
            ! Append the unit translations
            centringVectors(:,nCentringVectors+1) = (/ 1.,0.,0.,1. /)
            centringVectors(:,nCentringVectors+2) = (/ 0.,1.,0.,1. /)
            centringVectors(:,nCentringVectors+3) = (/ 0.,0.,1.,1. /)
            nCentringVectors = nCentringVectors + 3
            ! Combine three vectors until a primitive setting is found
            primitive = .false.
            do i = 1 , nCentringVectors
                if (primitive) exit
                do j = i + 1 , nCentringVectors
                    if (primitive) exit
                    do k = j + 1 , nCentringVectors
                        if (primitive) exit
                        P(:,1) = centringVectors(1:3,i)
                        P(:,2) = centringVectors(1:3,j)
                        P(:,3) = centringVectors(1:3,k)
                        call Invert_Matrix(P,Pinv,singular)
                        if (.not. singular) then
                            call Determinant(P,3,d)
                            if (abs(abs(1/d) - SG%NumLat) < 0.001) then
                                primitive = .true.
                                if (d < 0.) then
                                    P(:,1) = -P(:,1)
                                    call Invert_Matrix(P,Pinv,singular)
                                end if
                                do n = 1 , SG%Multip
                                    A = MatMul(SG%SymOP(n)%Rot,P)
                                    A = MatMul(Pinv,A)
                                    integral= Zbelong(A)
                                    if (.not. integral) then
                                        primitive = .false.
                                        exit
                                    end if
                                end do
                            end if
                        end if
                    end do
                end do
            end do
        end if

        if (.not. primitive) then
            ERR_StdRep = .true.
            ERR_StdRep_Mess = "Error in Get_P_Matrix subroutine. A primitive basis cannot be found"
            return
        else if (output) then
            write(*,'(12x,a)',advance='no') "Original setting --> Primitive setting transformation: "
            call Get_Symb_From_Mat(transpose(P(1:3,1:3)),symb,(/"a","b","c"/))
            write(*,'(a)') trim(symb)
        end if

    end subroutine Get_P_Matrix

    subroutine Get_Point_Group(SG)

        !---- Arguments ----!
        type(Space_Group_Type), intent(inout) :: SG ! Space Group

        !---- Local variables ----!
        integer       :: i,j,n,t
        integer       :: nSelected
        logical       :: selected
        real(kind=cp) :: d
        integer,             dimension(6)                :: nRot
        type(Sym_Oper_Type), dimension(:),   allocatable :: repSymOp ! representative operations
        integer,             dimension(:,:), allocatable :: idd

        nRot(:) = 0  ! number of selected rotations of order 1,2,....6
        SG%PG   = "" ! point group
        SG%NumOps = SG%Multip / SG%NumLat
        ERR_Symm = .false.
        if (SG%Centred /= 1) SG%NumOps = SG%NumOps / 2
        allocate(repSymOp(SG%NumOps))

        ! Count all the rotations in the set of representative matrices
        nSelected = 0
        do i = 1, SG%Multip
            call Determinant(real(SG%SymOP(i)%Rot),3,d)
            t = Trace(SG%SymOp(i)%Rot)
            selected = .false.
            do j = 1 , nSelected
                if (Equal_Matrix(SG%SymOp(i)%Rot,repSymOp(j)%Rot,3) .or. &
                    Equal_Matrix(SG%SymOp(i)%Rot,-repSymOp(j)%Rot,3)) then
                    selected = .true.
                    exit
                end if
            end do
            if (.not. selected) then
                nSelected           = nSelected + 1
                if (nSelected > SG%NumOps) then
                    ERR_Symm = .True.
                    ERR_Symm_Mess = "Error in Get_Point_Group. &
                    nSelected > SG%NumOps"
                    return
                end if
                repSymOp(nSelected) = SG%SymOp(i)
                select case (abs(t))
                    case (0)
                        nRot(3) = nRot(3) + 1
                    case(1)
                        if (Nint(d * t) ==  1) then
                            nRot(4) = nRot(4) + 1
                        else
                            nRot(2) = nRot(2) + 1
                        end if
                    case(2)
                        nRot(6) = nRot(6) + 1
                    case(3)
                        nRot(1) = nRot(1) + 1
                end select
            end if
        end do
        if (nSelected < SG%NumOps) then
            ERR_Symm = .true.
            ERR_Symm_Mess = "Error in Get_Point_Group. &
            nSelected < SG%NumOps"
            return
        end if

        ! Get the point group
        allocate(idd(SG%NumOps,2))
        if (nRot(3) == 8) then ! Cubic
            if (SG%NumOps == 12) then
                if (SG%Centred == 1) then
                    SG%PG = "23"
                else
                    SG%PG = "m-3"
                end if
            else if (SG%NumOps == 24) then
                if (SG%Centred == 1) then
                    call Get_Rotation_Id(repSymOp,SG%NumOps,4,n,idd)
                    if (n == 6 .and. idd(1,2) == 1) then
                        SG%PG = "432"
                    else if (n == 6 .and. idd(1,2) == -1) then
                        SG%PG = "-43m"
                    end if
                else
                    SG%PG = "m-3m"
                end if
            end if
        else if (nRot(6) == 2) then ! Hexagonal
            if (SG%NumOps == 6) then
                if (SG%Centred == 1) then
                    call Get_Rotation_Id(repSymOp,SG%NumOps,6,n,idd)
                    if (idd(1,2) == 1) then
                        SG%PG = "6"
                    else
                        SG%PG = "-6"
                    end if
                else
                    SG%PG = "6/m"
                end if
            else if (SG%NumOps == 12) then
                if (SG%Centred == 1) then
                    call Get_Rotation_Id(repSymOp,SG%NumOps,6,n,idd)
                    if (idd(1,2) == 1) then
                        call Get_Rotation_Id(repSymOp,SG%NumOps,2,n,idd)
                        if (n == 7 .and. idd(1,2) == 1) then
                            SG%PG = "622"
                        else
                            SG%PG = "6mm"
                        end if
                    else if (idd(1,2) == -1) then
                        SG%PG = "-6m2"
                    end if
                else
                    SG%PG = "6/mmm"
                end if
            end if
        else if (nRot(3) == 2) then ! Trigonal
            if (SG%NumOps == 3) then
                if (SG%Centred == 1) then
                    SG%PG = "3"
                else
                    SG%PG = "-3"
                end if
            else if (SG%NumOps == 6) then
                if (SG%Centred == 1) then
                    call Get_Rotation_Id(repSymOp,SG%NumOps,2,n,idd)
                    if (n == 3 .and. idd(1,2) == 1) then
                        SG%PG = "32"
                    else if (n == 3 .and. idd(1,2) == -1) then
                        SG%PG = "3m"
                    end if
                else
                    SG%PG = "-3m"
                end if
            end if
        else if (nRot(4) == 2) then ! Tetragonal
            if (SG%NumOps == 4) then
                if (SG%Centred == 1) then
                    call Get_Rotation_Id(repSymOp,SG%NumOps,4,n,idd)
                    if (n == 2 .and. idd(1,2) == 1) then
                        SG%PG = "4"
                    else if (n == 2 .and. idd(1,2) == -1) then
                        SG%PG = "-4"
                    end if
                else
                    SG%PG = "4/m"
                end if
            else if (SG%NumOps == 8) then
                if (SG%Centred == 1) then
                    call Get_Rotation_Id(repSymOp,SG%NumOps,4,n,idd)
                    if (n == 2 .and. idd(1,2) == 1) then
                        call Get_Rotation_Id(repSymOp,SG%NumOps,2,n,idd)
                        if (n == 5 .and. idd(1,2) == 1) then
                            SG%PG = "422"
                        else
                            SG%PG = "4mm"
                        end if
                    else if (n == 2 .and. idd(1,2) == -1) then
                        SG%PG = "-4m2"
                    end if
                else
                    SG%PG = "4/mmm"
                end if
            end if
        else if (nRot(2) == 3) then ! Orthorhombic
            if (SG%Centred == 1) then
                call Get_Rotation_Id(repSymOp,SG%NumOps,2,n,idd)
                if (n == 3 .and. idd(1,2) == 1) then
                    SG%PG = "222"
                else
                    SG%PG = "mm2"
                end if
            else
                SG%PG = "mmm"
            end if
        else if (nRot(2) == 1) then ! Monoclinic
            if (SG%Centred == 1) then
                call Get_Rotation_Id(repSymOp,SG%NumOps,2,n,idd)
                if (idd(1,2) == 1) then
                    SG%PG = "2"
                else
                    SG%PG = "m"
                end if
            else
                SG%PG = "2/m"
            end if
        else
            if (SG%Centred == 1) then ! Triclinic
                SG%PG = "1"
            else
                SG%PG = "-1"
            end if
        end if

    end subroutine Get_Point_Group

    subroutine Get_Pseudo_Standard_Base(W,d,order,paxis,bz,bx,by)

        ! Given the rotation matrix of the principal axis W, the shortest
        ! lattice vector along the axis, bz, and the four shortest lattice
        ! vectors perpendicular to that axis, paxis, computes the basis
        ! bx,by,bz which gives the smallest unit cell

        !---- Arguments ----!
        integer, dimension(3,3), intent(in)  :: W
        integer,                 intent(in)  :: d
        integer,                 intent(in)  :: order
        integer, dimension(3,4), intent(in)  :: paxis
        integer, dimension(3),   intent(in)  :: bz
        integer, dimension(3),   intent(out) :: bx
        integer, dimension(3),   intent(out) :: by

        !---- Local variables ----!
        integer                           :: i,j,n,imin
        integer, dimension(3,3)           :: A,B
        real                              :: mind
        real, dimension(:,:), allocatable :: trials,byaux

        B(:,3) = bz(:)

        select case (order)

            case (2)

                Allocate(trials(6,3))
                n = 0
                do i = 1 , 4
                    do j = i + 1 , 4
                        n = n + 1
                        trials(n,1) = i
                        trials(n,2) = j
                        B(:,1)      = paxis(:,i)
                        B(:,2)      = paxis(:,j)
                        call Determinant(real(B),3,trials(n,3))
                    end do
                end do

                ! Select the combination with the smallest determinant

                imin = 1
                mind = abs(trials(1,3))
                do i = 2 , 6
                    if (abs(trials(i,3)) < mind) then
                        imin = i
                        mind = abs(trials(i,3))
                    end if
                end do

                bx(:) = paxis(:,int(trials(imin,1)))
                by(:) = paxis(:,int(trials(imin,2)))

            case (3,4,6)

                A = (d/abs(d)) * W

                allocate(trials(4,1))
                allocate(byaux(3,4))
                do i = 1 , 4
                    byaux(:,i) = matmul(A,paxis(:,i))
                    B(:,1) = paxis(:,i)
                    B(:,2) = byaux(:,i)
                    call Determinant(real(B),3,trials(i,1))
                end do

                imin = 1
                mind = abs(trials(1,1))
                do i = 2 , 4
                    if (abs(trials(i,1)) < mind) then
                        imin = i
                        mind = abs(trials(i,1))
                    end if
                end do

                bx(:) = paxis(:,imin)
                by(:) = byaux(:,imin)

        end select

    end subroutine Get_Pseudo_Standard_Base

    subroutine Get_Rotation_Axis(W,axis)

        ! Computes the shortest lattice vector along
        ! the rotation axis for the symmetry operation W

        !---- Arguments ----!
        integer, dimension(3,3), intent(in)  :: W     !rotation matrix
        integer, dimension(3),   intent(out) :: axis  !shortest vector along the rotation axis

        !---- Local variables ----!
        integer                 :: i,r,nrow=3,ncolumn=3,nzeros_aux,d
        real                    :: dr
        logical                 :: integral,ordered
        real,    dimension(3)   :: raxis,axisaux
        logical, dimension(3)   :: fix
        integer, dimension(3,3) :: A,U
        integer, dimension(:), allocatable :: nzeros,row

        call Determinant(real(W),3,dr)
        d = Nint(dr)

        A = (d / abs(d)) * W
        do i = 1 , 3
            A(i,i) = A(i,i) - 1
        end do

        !call Get_Row_Echelon_Form(A,U)
        U = A
        call RowEchelonForm(U)

        ! Check that the rank of the U matrix is two

        call Rank(real(U),0.001,r)

        if ( r /= 2 ) then
            ERR_Symm = .true.
            ERR_Symm_Mess = "Error in Get_Rotation_Axis subroutine. &
            Rank of matrix U is not two"
            return
        end if

        ! If there is a row with all zero entries,
        ! we put it in the last row

        allocate(nzeros(nrow),row(ncolumn))

        nzeros = 0
        do r = 1 , nrow
            do i = 1 , ncolumn
                if ( U(r,i) == 0 ) then
                    nzeros(r) = nzeros(r) + 1
                else
                    exit
                end if
            end do
        end do

        do
            r = 1
            ordered = .true.
            do r = 1 , nrow - 1
                if (nzeros(r+1) < nzeros(r)) then
                    nzeros_aux  = nzeros(r)
                    row(:)      = U(r,:)
                    U(r,:)      = U(r+1,:)
                    U(r+1,:)    = row(:)
                    nzeros(r)   = nzeros(r+1)
                    nzeros(r+1) = nzeros_aux
                    ordered      = .false.
                end if
            end do
            if ( ordered ) exit
        end do

        ! Compute the axis direction

        if ( U(3,3) /= 0 ) then
            raxis(3) = 0.
        else if ( U(2,3) /= 0 .and. U(2,2) == 0 ) then
            raxis(3) = 0.
        else
            raxis(3) = 1.      ! Free choice
        end if

        if ( U(2,2) /= 0 ) then
            raxis(2) = - U(2,3) * raxis(3) / U(2,2)
            if (U(1,1) == 0) then
                raxis(1) = 1. ! Free choice
            else
                raxis(1) = - ( U(1,2) * raxis(2) + U(1,3) * raxis(3) ) / U(1,1)
            end if
        else ! raxis(3) must be zero because row(2) cannot be zero
            if ( U(1,2) == 0 ) then
                raxis(2) = 1. ! Free choice
                raxis(1) = - ( U(1,2) * raxis(2) + U(1,3) * raxis(3) ) / U(1,1)
            else
                raxis(1) = 1. ! Free choice
                raxis(2) = - ( U(1,1) * raxis(1) + U(1,3) * raxis(3) ) / U(1,2)
            end if
        end if

        ! Find the vector with the smallest set of
        ! integers which is along the axis

        do i = 1 , 10
            raxis = i * raxis
            integral=Zbelong(raxis)
            if (integral) exit
        end do

        if (integral) then
            do i = 2, 10
                axisaux = raxis / i
                if (Zbelong(axisaux)) raxis = axisaux
            end do
            axis = Nint(raxis)
        else
            ERR_Symm = .true.
            ERR_Symm_Mess = "Error in Get_Rotation_Axis subroutine. Unable to find a set of integers"
            return
        end if

        ! Choose the eigenvector raxis with axis(3) > 0. If axis(3) = 0, choose axis with
        ! axis(2) > 0. If axis(2) = 0, choose raxis with axis(1) > 0

        do i = 3, 1, -1
            if (axis(i) /= 0) then
                if (axis(i) < 0) axis = -axis
                exit
            end if
        end do

    end subroutine Get_Rotation_Axis

    subroutine Get_Rotation_Id_Non_Standard(SymOP,nSymOP,order,n,idd)

        ! Identifies the symmetry operations in the array SymOP
        ! which correspond to rotations of order 'order'. It returns
        ! the number of selected operations, and an integer array
        ! with the associated id and determinant values. The id
        ! is the index of the symmetry operation in the array
        ! SymOP.

        !---- Arguments ----!
        type(NS_Sym_Oper_Type), dimension(nSymOP), intent(in)  :: SymOP
        integer,                                   intent(in)  :: nSymOP
        integer,                                   intent(in)  :: order
        integer,                                   intent(out) :: n
        integer, dimension(nSymOP,2),              intent(out) :: idd

        !---- Local variables ----!
        integer                           :: i,tr
        real(kind=cp)                     :: d
        integer,           dimension(6)   :: traces
        character(len=20), dimension(6)   :: axisName

        axisName(2) = "twofold"
        axisName(3) = "threefold"
        axisName(4) = "fourfold"
        axisName(6) = "sixfold"

        traces(2) = -1
        traces(3) =  0
        traces(4) =  1
        traces(6) =  2

        n = 0
        do i = 1 , nSymOP
            call Determinant(SymOP(i)%Rot,3,d)
            tr = Nint(Trace(SymOP(i)%Rot))
            if ( (nint(d) == 1 .and. tr == traces(order)) .or. (nint(d) == -1 .and. tr == -traces(order))) then
                n = n + 1
                idd(n,1) = i
                idd(n,2) = nint(d)
            end if
        end do
        if (n == 0) then
            Err_StdRep = .true.
            Err_StdRep_Mess = "Error in Get_Rotation_Id. Unable to find the "&
                //trim(axisName(order))//" axis"
        end if

    end subroutine Get_Rotation_Id_Non_Standard

    subroutine Get_Rotation_Id_Standard(SymOP,nSymOP,order,n,idd)

        ! Identifies the symmetry operations in the array SymOP
        ! which correspond to rotations of order 'order'. It returns
        ! the number of selected operations, and an integer array
        ! with the associated id and determinant values. The id
        ! is the index of the symmetry operation in the array
        ! SymOP.

        !---- Arguments ----!
        type(Sym_Oper_Type), dimension(nSymOP), intent(in)  :: SymOP
        integer,                                intent(in)  :: nSymOP
        integer,                                intent(in)  :: order
        integer,                                intent(out) :: n
        integer, dimension(nSymOP,2),           intent(out) :: idd

        !---- Local variables ----!
        integer                           :: i,tr
        real(kind=cp)                     :: d
        integer,           dimension(6)   :: traces
        character(len=20), dimension(6)   :: axisName

        axisName(1) = "onefold"
        axisName(2) = "twofold"
        axisName(3) = "threefold"
        axisName(4) = "fourfold"
        axisName(6) = "sixfold"

        traces(1) =  3
        traces(2) = -1
        traces(3) =  0
        traces(4) =  1
        traces(6) =  2

        n = 0
        do i = 1 , nSymOP
            call Determinant(real(SymOP(i)%Rot),3,d)
            tr = Trace(SymOP(i)%Rot)
            if ( (nint(d) == 1 .and. tr == traces(order)) .or. (nint(d) == -1 .and. tr == -traces(order))) then
                n = n + 1
                idd(n,1) = i
                idd(n,2) = nint(d)
            end if
        end do
        if (n == 0) then
            ERR_Symm = .true.
            ERR_Symm_Mess = "Error in Get_Rotation_Id. Unable to find the "&
                //trim(axisName(order))//" axis"
        end if

    end subroutine Get_Rotation_Id_Standard

    subroutine Get_S_Matrix(W,order,S)

        ! Given a rotation matrix W, computes the matrix S defined as:
        !           S = W + W^2 + W^3 + ... + W^n
        ! where n is the order of the rotation axis. This matrix is
        ! used to find vectors perpendicular to the rotation axis. For
        ! a vector x perpendicular to the rotation axis, Sx = 0

        !---- Arguments ----!
        integer, dimension(3,3), intent(in)  :: W
        integer,                 intent(in)  :: order
        integer, dimension(3,3), intent(out) :: S

        !---- Local variables ----!
        integer i
        integer, dimension(3,3) :: Waux

        S    = W
        Waux = W

        do i = 2 , order
            Waux = MatMul(Waux,W)
            S    = S + Waux
        end do

    end subroutine Get_S_Matrix

    subroutine Get_Standard_Representation(SG,SGstd,C,symb_set,outp)

        ! Given an arbitrary representation of a space group,
        ! it finds its standard representation. It follows
        ! Acta Cryst. A55 383-395 (1999). P,M,A and C matrices
        ! in the paper correspond to Pinv,Minv,Ainv and Cinv in
        ! this code.

        !---- Arguments ----!
        type(NS_Space_Group_Type),  intent(in out) :: SG     ! space group in the original setting (the name is changed to standard+trans)
        type(Space_Group_Type),        intent(out) :: SGstd  ! space group in the standard setting
        real(kind=cp), dimension(4,4), intent(out) :: C      ! transformation matrix
        character(len=*),              intent(out) :: symb_set !a,b,c;1/2,0,0
        logical, optional,             intent(in)  :: outp

        !---- Local variables ----!
        integer                         :: n,sgNumber,i,j
        character                       :: lattyp
        character(len=90)               :: setting
        real(kind=cp), dimension(3,3)   :: P,Mp,Mc,M
        real(kind=cp), dimension(3,3,6) :: A
        integer, dimension(3) :: axis

        call Init_Err_StdRep()
        call Init_Err_Symm()
        if(present(outp)) then
          call Get_P_Matrix(SG,P,info=.true.)
          if (Err_StdRep) return
          call Get_Mprime_Matrix(SG,P,Mp,info=.true.)
          if (Err_StdRep) return
          call Get_Mc_Matrix(SG%Laue,Mp,Mc,info=.true.)
        else
          call Get_P_Matrix(SG,P)
          if (Err_StdRep) return
          call Get_Mprime_Matrix(SG,P,Mp)
          if (Err_StdRep) return
          call Get_Mc_Matrix(SG%Laue,Mp,Mc)
        end if
        M = matmul(Mp,Mc)
        call Get_Lattice_Type_from_M(M,lattyp)
        call Get_A_Matrices(SG%Laue,A,n)
        call Space_Group_Matching(SG,P,M,A(:,:,1:n),n,C,SGstd)

        if (Err_StdRep) then
          return
        else
          SG%NumSpg=SGstd%NumSpg
          SG%SPG_Symb=SGstd%SPG_Symb
          call Get_Setting_Info(C(1:3,1:3),C(1:3,4),setting)
          SG%SG_setting=trim(SGstd%SPG_Symb)//"  <-- "//trim(setting)
          i=index(setting,"a'=")
          if(i /= 0) setting(i:i+2)=" "
          i=index(setting,"b'=")
          if(i /= 0) setting(i:i+2)=" "
          i=index(setting,"c'=")
          if(i /= 0) setting(i:i+2)=" "
          i=index(setting,"-> Origin: (")
          if(i /= 0) setting(i:i+11)=";           "
          setting=pack_string(setting)
          i=len_trim(setting)
          symb_set=setting(1:i-1)
        end if

    end subroutine Get_Standard_Representation

    subroutine Get_Vectors_Perpendicular_To_Rotation_Axis(W,d,order,paxis)

        ! Given a rotation matrix W, computes the shortest
        ! four lattice vectors perpendicular to the rotation
        ! axis

        !---- Arguments ----!
        integer, dimension(3,3), intent(in)  :: W
        integer,                 intent(in)  :: d
        integer,                 intent(in)  :: order
        integer, dimension(3,4), intent(out) :: paxis

        !---- Local variables ----!
        integer                 :: i,j,r,nzeros,row
        integer, dimension(3,3) :: A,U,S
        real, dimension(3)      :: v
        real, dimension(3,4)    :: paxisaux
        logical                 :: info,zerorow

        A = (d/abs(d)) * W
        call Get_S_Matrix(A,order,S)
        call RowEchelonForm(S)
        U = S

        ! Check that the rank of the U matrix is one

        call Rank(real(U),0.001,r)

        if ( r /= 1 ) then
            ERR_Symm = .true.
            ERR_Symm_Mess = "Error in Get_Vectors_Perpendicular_To_Rotation_Axis subroutine. &
            Rank of matrix U is not one"
            return
        end if

        ! Find a row different from zeros

        zerorow = .true.
        do i = 1 , 3
            do j = 1 , 3
                if (U(i,j) /= 0) then
                    zerorow = .false.
                    exit
                end if
            end do
            if (.not. zerorow) exit
        end do
        row = i

        ! Count the number of zeros of the first row of the echelon matrix

        nzeros = 0
        do i = 1 , 3
            if (U(row,i) == 0) nzeros = nzeros + 1
        end do

        ! Build the four shortest vectors perpendicular to the rotation axis

        select case (nzeros)
            case (0)
                    paxisaux(1,1) =  1
                    paxisaux(2,1) =  0
                    paxisaux(1,2) =  0
                    paxisaux(2,2) =  1
                    paxisaux(1,3) =  1
                    paxisaux(2,3) =  1
                    paxisaux(1,4) =  1
                    paxisaux(2,4) = -1
                    do i = 1 , 4
                        paxisaux(3,i) = -(paxisaux(1,i) * U(row,1) + paxisaux(2,i) * U(row,2)) / U(row,3)
                    end do
            case (1)
                do i = 1 , 3
                    if (U(row,i) == 0) exit
                end do
                select case (i)
                    case (1)
                        paxisaux(1,1) = 1
                        paxisaux(2,1) = 0
                        paxisaux(3,1) = 0
                        paxisaux(1,2) = 0
                        paxisaux(2,2) = 1
                        paxisaux(3,2) = - U(row,2) / U(row,3)
                        paxisaux(1,3) = 1
                        paxisaux(2,3) = 1
                        paxisaux(3,3) = - U(row,2) / U(row,3)
                        paxisaux(1,4) = -1
                        paxisaux(2,4) = 1
                        paxisaux(3,4) = - U(row,2) / U(row,3)
                    case (2)
                        paxisaux(1,1) = 0
                        paxisaux(2,1) = 1
                        paxisaux(3,1) = 0
                        paxisaux(1,2) = 1
                        paxisaux(2,2) = 0
                        paxisaux(3,2) = - U(row,1) / U(row,3)
                        paxisaux(1,3) = 1
                        paxisaux(2,3) = 1
                        paxisaux(3,3) = - U(row,1) / U(row,3)
                        paxisaux(1,4) = 1
                        paxisaux(2,4) = -1
                        paxisaux(3,4) = - U(row,1) / U(row,3)
                    case (3)
                        paxisaux(1,1) = 0
                        paxisaux(2,1) = 0
                        paxisaux(3,1) = 1
                        paxisaux(1,2) = 1
                        paxisaux(2,2) = -U(row,1) / U(row,2)
                        paxisaux(3,2) = 0
                        paxisaux(1,3) = 1
                        paxisaux(2,3) = -U(row,1) / U(row,2)
                        paxisaux(3,3) = 1
                        paxisaux(1,4) = 1
                        paxisaux(2,4) = -U(row,1) / U(row,2)
                        paxisaux(3,4) = -1
                end select
            case (2)
                do i = 1 , 3
                    if (U(row,i) /= 0) exit
                end do
                select case (i)
                    case (1)
                        paxisaux(1,:) =  0
                        paxisaux(2,1) =  1
                        paxisaux(3,1) =  0
                        paxisaux(2,2) =  0
                        paxisaux(3,2) =  1
                        paxisaux(2,3) =  1
                        paxisaux(3,3) =  1
                        paxisaux(2,4) =  1
                        paxisaux(3,4) = -1
                    case (2)
                        paxisaux(2,:) =  0
                        paxisaux(1,1) =  1
                        paxisaux(3,1) =  0
                        paxisaux(1,2) =  0
                        paxisaux(3,2) =  1
                        paxisaux(1,3) =  1
                        paxisaux(3,3) =  1
                        paxisaux(1,4) =  1
                        paxisaux(3,4) = -1
                    case (3)
                        paxisaux(3,:) =  0
                        paxisaux(1,1) =  1
                        paxisaux(2,1) =  0
                        paxisaux(1,2) =  0
                        paxisaux(2,2) =  1
                        paxisaux(1,3) =  1
                        paxisaux(2,3) =  1
                        paxisaux(1,4) =  1
                        paxisaux(2,4) = -1
                end select
        end select

        ! For each solution, determine the smallest set of
        ! integers

        do i = 1 , 4
            v(:) = paxisaux(:,i)
            do j = 1 , 10
                paxisaux (:,i)= j * paxisaux(:,i)
                info= Zbelong(paxisaux(:,i))
                if (info) exit
            end do
            if (info) then
                paxis(:,i) = Nint(paxisaux(:,i))
            else
                ERR_Symm = .true.
                ERR_Symm_Mess = "Error in Get_Vectors_Perpendicular_To_Rotation_Axis &
                subroutine. Unable to find a set of integers"
                return
            end if
        end do

    end subroutine Get_Vectors_Perpendicular_To_Rotation_Axis

    !!----
    !!---- Subroutine Init_Err_StdRep()
    !!----
    !!----    Initialize the errors flags in this Module
    !!----
    !!---- Update: May - 2018
    !!
    subroutine Init_Err_StdRep()

       Err_StdRep=.false.
       ERR_StdRep_Mess=" "

       return

    end subroutine Init_Err_StdRep


    subroutine Init_Matrices(P,M,A)

        !---- Arguments ----!
        real(kind=cp), dimension(3,3),   intent(out) :: P,M
        real(kind=cp), dimension(3,3,6), intent(out) :: A

        !---- Local variables ----'
        integer i

        P = 0.0
        M = 0.0
        A = 0.0

        do i = 1 , 3
            P(i,i)   = 1.0
            M(i,i)   = 1.0
            A(i,i,1) = 1.0
        end do

    end subroutine Init_Matrices

    subroutine Positive_Sense_of_Rotation(W,axis,positive)

        !---- Arguments ----!
        integer, dimension(3,3), intent(in)  :: W
        integer, dimension(3),   intent(in)  :: axis
        logical,                 intent(out) :: positive

        !---- Local variables ----!
        real    :: d
        integer, dimension(3,3) :: U

        call Determinant(real(W),3,d)
        if (d < 0.) then
            U = -W
        else
            U = W
        end if

        if (axis(2) == 0 .and. axis(3) == 0 .and. (axis(1) * U(3,2)) > 0) then
            positive = .true.
        else if ( (U(2,1) * axis(3) - U(3,1) * axis(2)) > 0 ) then
            positive = .true.
        else
            positive = .false.
        end if

    end subroutine Positive_Sense_of_Rotation

    subroutine Set_Right_Handedness(A)

        !---- Arguments ----!
        real(kind=cp), dimension(3,3), intent(inout) :: A

        !---- Local variables ----!
        real                        :: d
        real(kind=cp), dimension(3) :: row

        call Determinant(A,3,d)

        if (d < 0.) then
            row(:) = A(:,1)
            A(:,1) = A(:,2)
            A(:,2) = row(:)
        end if

    end subroutine Set_Right_Handedness


    subroutine Space_Group_Matching(SG,P,M,A,n,C,SGstd,info)

        !---- Arguments ----!
        type(NS_Space_Group_Type),       intent(in)  :: SG       ! space group in the original setting
        integer,                         intent(in)  :: n        ! number of A matrices (six as maximum)
        real(kind=cp), dimension(3,3),   intent(in)  :: P        ! P matrix   -see Get_P_Matrix-
        real(kind=cp), dimension(3,3),   intent(in)  :: M        ! M matrix   -see Get_M_Matrix-
        real(kind=cp), dimension(3,3,n), intent(in)  :: A        ! A matrices -see Get_A_Matrices-
        real(kind=cp), dimension(4,4),   intent(out) :: C        ! transformation matrix
        type(Space_Group_Type),          intent(out) :: SGstd    ! space group in the standard setting
        logical, optional,               intent(in)  :: info

        !---- Local variables ----!
        integer                   :: i,j,k,l,r,s,ng,ng_,nr,firstSpaceGroup,lastSpaceGroup,sgNumber
        real                      :: dx
        character(len=12)         :: sgString
        character(len=60)         :: symb,tr_symb
        character(len=256)        :: generators
        logical                   :: output,singular,matching,getShift
        type(NS_Space_Group_Type) :: SGtargetaux
        type(Space_Group_Type)    :: SGtarget
        character(len=1),          dimension(n) :: lattyp
        type(NS_Space_Group_Type), dimension(n) :: SGaux

        real,          dimension(3)         :: nullVector,oShift
        integer,       dimension(3,3)       :: identity
        real(kind=cp), dimension(3,3)       :: MA,Pt
        real(kind=cp), dimension(4,4)       :: Cinv,W,Wp
        real(kind=cp), dimension(4,4,3)     :: Gt,Gx
        real(kind=cp), dimension(4,4,n)     :: C_

        identity = 0
        identity(1,1) = 1
        identity(2,2) = 1
        identity(3,3) = 1
        nullVector    = 0.

        select case (trim(SG%Laue))

            case ("-1")

                firstSpaceGroup = 1
                lastSpaceGroup  = 2
                if(present(info)) write(*,'(8x,a)') " => Matching representation against standard triclinic space groups...."

            case ("2/m")

                firstSpaceGroup = 3
                lastSpaceGroup  = 15
                if(present(info)) write(*,'(8x,a)') " => Matching representation against standard monoclinic space groups...."

            case ("mmm")
                firstSpaceGroup = 16
                lastSpaceGroup  = 74
                if(present(info)) write(*,'(8x,a)') " => Matching representation against standard orthorhombic space groups...."

            case ("4/m","4/mmm")
                firstSpaceGroup = 75
                lastSpaceGroup  = 142
                if(present(info)) write(*,'(8x,a)') " => Matching representation against standard tetragonal space groups...."

            case ("-3","-3 R","-3m","-3m R","-3m1","-31m")
                firstSpaceGroup = 143
                lastSpaceGroup  = 167
                if(present(info)) write(*,'(8x,a)') " => Matching representation against standard trigonal space groups...."

            case ("6/m","6/mmm")
                firstSpaceGroup = 168
                lastSpaceGroup  = 194
                if(present(info)) write(*,'(8x,a)') " => Matching representation against standard hexagonal space groups...."

            case ("m3","m-3","m3m","m-3m")
                firstSpaceGroup = 195
                lastSpaceGroup  = 230
                if(present(info)) write(*,'(8x,a)') " => Matching representation against standard cubic space groups...."

            case default
                firstSpaceGroup = 0
                lastSpaceGroup  = -1

        end select

        do s = 1 , n  ! Loop over C matrices
            SGaux(s) = SG
            MA            = matmul(M,A(:,:,s))
            C_(1:3,1:3,s) = matmul(P,MA)
            C_(1:3,4,s)   = (/ 0.,0.,0. /)
            C_(4,1:4,s)   = (/ 0.,0.,0.,1. /)
            call Invert_Matrix(C_(:,:,s),Cinv,singular)
            ! Compute symmetry operations in the new basis
            do j = 1 , SG%Multip
                W(1:3,1:3) = SG%SymOP(j)%Rot(1:3,1:3)
                W(1:3,4)   = SG%SymOP(j)%Tr(1:3)
                W(4,1:4)   = (/ 0.,0.,0.,1. /)
                Wp         = matmul(W,C_(:,:,s))
                Wp         = matmul(Cinv,Wp)
                SGaux(s)%SymOP(j)%Rot = Wp(1:3,1:3)
                SGaux(s)%SymOP(j)%Tr  = Wp(1:3,4)
                !call Get_SymSymb(SGaux(s)%SymOp(j)%Rot,SG%SymOp(j)%tr,symb)
                !write(*,"(12x,i2,a)") j, " -> "//trim(symb)
            end do
            ! Get Lattice Type
            call Get_Lattice_Type_from_M(MA,lattyp(s))
        end do
        do sgNumber = firstSpaceGroup , lastSpaceGroup
            call Get_HM_Standard(sgNumber,sgString)
            call Set_Spacegroup(sgString,SGtarget)
            do s = 1 , n
                if (SGtarget%NumOPS == SGaux(s)%NumOPS .and. &
                    SGtarget%Spg_Lat == lattyp(s) .and. &
                    (SGtarget%Centred == 1 .and. SGaux(s)%Centred == 1 .or. &
                     SGtarget%Centred /= 1 .and. SGaux(s)%Centred /= 1)) then
                    ! Get generators from the standard space group
                    call Get_Generators(SGtarget%NumSpg,SGtarget%SymOp,SGtarget%Multip,Gt,ng)
                    if (ERR_StdRep) return
                    ! Try to get these generators from SGaux(s)
                    ng_ = 0
                    do i = 1 , ng
                        do j = 1 , SGaux(s)%Multip
                            if (Equal_Matrix(Gt(1:3,1:3,i),SGaux(s)%SymOp(j)%Rot,3)) then
                                ng_ = ng_ + 1
                                Gx(1:3,1:3,ng_) = Gt(1:3,1:3,i)
                                Gx(1:3,4,ng_)   = SGaux(s)%SymOp(j)%Tr
                                Gx(4,:,ng_)     = (/ 0.0,0.0,0.0,1.0 /)
                                exit
                            end if
                        end do
                    end do
                    if (ng /= ng_) cycle
                     if(present(info)) write(*,'(tr12,3a,i2)',advance = 'no') 'Trying to match space group ', &
                                                  SGtarget%Spg_Symb, 'setting ', s
                    ! Get a primitive setting
                    call Setting_Change(real(identity),nullVector,SGtarget,SGtargetaux,"IT")
                    call Get_P_Matrix(SGtargetaux,Pt(1:3,1:3))
                    ! Try to match the standard by an origin shift
                    call Get_Origin_Shift(Gx(:,:,1:ng),&
                                         Gt(:,:,1:ng),ng,Pt,oShift,getShift)
                    if (getShift) then
                        if(present(info)) write(*,'(tr8,a)') 'Done!'
                        C(1:3,1:3) = C_(1:3,1:3,s)
                        C(4,:)     = C_(4,:,s)
                        oShift = matmul(C(1:3,1:3),oShift)
                        C(1:3,4) = oShift(1:3)
                        SGstd = SGtarget
                        return
                    else
                        if(present(info)) write(*,'(tr8,a)') 'Failed'
                        cycle
                    end if
                end if
            end do
        end do

        Err_StdRep = .true.
        Err_StdRep_Mess = "Unable to indentify the space group"
        SGstd%NumSpg = 0

    end subroutine Space_Group_Matching

end module CFML_SpG_Standard_Representation