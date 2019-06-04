!!---- MODULE: CFML_Standard_Settings
!!----   INFO: This module contains everything needed for the identification of
!!----         Shubnikov groups, and for the calculation of the associated
!!----         transformation matrices to standard settings.
!!----
!!---- HISTORY
!!----    Update: 09/01/2019
!!----
!!---- DEPENDENCIES
!!----
!!----    Use CFML_Rational_Arithmetic
!!----    Use CFML_Rational_Groups
!!----    Use CFML_String_Utilities,          only: Get_Separator_Pos,Pack_String
!!----    Use CFML_Symmetry_Tables,           only: Spgr_Info,Set_Spgr_Info
!!----    Use CFML_Crystallographic_Symmetry, only: Space_Group_Type,Set_SpaceGroup
!!----    Use CFML_Magnetic_Groups
!!----
!!---- VARIABLES
!!----    err_std
!!----    err_std_mess
!!----
!!---- PROCEDURES
!!----    Functions:
!!----       POSITIVE_SENSE_OF_ROTATION
!!----
!!----    Subroutines:
!!----       GET_A_MATRICES_CRYS
!!----       GET_A_MATRICES_SHUB
!!----       GET_GENERATORS
!!----       GET_LATTICE_TYPE
!!----       GET_LATTICE_TYPE_FROM_M
!!----       GET_MAGNETIC_LATTICE_TYPE
!!----       GET_Mc_MATRIX
!!----       GET_Mp_MATRIX
!!----       GET_P_MATRIX
!!----       GET_PSEUDO_STANDARD_BASE
!!----       GET_ROTATION_AXIS
!!----       GET_ROTATION_ORDER
!!----       GET_ROTATIONS
!!----       GET_S_MATRIX
!!----       GET_VECTORS_PERPENDICULAR_TO_ROTATION_AXIS
!!----       IDENTIFY_CRYSTALLOGRAPHIC_POINT_GROUP
!!----       IDENTIFY_CRYSTALLOGRAPHIC_SPACE_GROUP
!!----       IDENTIFY_GROUP
!!----       IDENTIFY_LAUE_CLASS
!!----       IDENTIFY_SHUBNIKOV_GROUP
!!----       MATCH_SHUBNIKOV_GROUP
!!----       MATCH_CRYSTALLOGRAPHIC_SPACE_GROUP
!!----       SET_RIGHT_HANDEDNESS
!!----       SMALLEST_INTEGRAL_VECTOR
!!----
!!
module CFML_Standard_Settings

    !---- Used External Modules ----!
    use CFML_Rational_Arithmetic
    use CFML_Rational_Groups
    use CFML_String_Utilities,          only: Get_Separator_Pos,Pack_String
    use CFML_Symmetry_Tables,           only: Spgr_Info,Set_Spgr_Info
    use CFML_Crystallographic_Symmetry, only: Space_Group_Type,Set_SpaceGroup
    use CFML_Magnetic_Groups

    implicit none

    private

    !---- List of public variables and types ----!
    public :: err_std, err_std_mess

    !---- List of public overloaded operators ----!

    !---- List of public functions ----!

    !---- List of public subroutines ----!
    public  :: Identify_Group, Identify_Laue_Class, Identify_Crystallographic_Point_Group, Get_Lattice_Type

    !!----
    !!---- err_std
    !!----    logical, public :: err_std
    !!----
    !!----    Logical Variable to indicate an error on this module.
    !!----
    !!---- Update: January - 2019
    !!
    logical :: err_std

    !!----
    !!---- err_std_mess
    !!----    character(len=:), public :: err_std_mess
    !!----
    !!----    String containing information about the last error
    !!----
    !!---- Update: January - 2019
    !!
    character(len=:), allocatable :: err_std_mess

contains

    !---- Functions ----!

    !!---- function Positive_Sense_of_Rotation(W,axis) result(positive)
    !!----      type(rational), dimension(3,3), intent(in)  :: W
    !!----      type(rational), dimension(3),   intent(in)  :: axis
    !!----      logical,                        intent(out) :: positive
    !!----
    !!---- Evaluates if axis corresponds with the positive sense of rotation
    !!---- of the rotation matrix W.
    !!----
    !!---- Updated: September - 2018
    !!----

    logical function Positive_Sense_of_Rotation(W,axis) result(positive)

        !---- Arguments ----!
        type(rational), dimension(3,3), intent(in)  :: W
        type(rational), dimension(3),   intent(in)  :: axis

        !---- Local variables ----!
        type(rational)                 :: det
        type(rational), dimension(3,3) :: U

        det = rational_determinant(W)
        if (det < (0//1)) then
            U = -W
        else
            U = W
        end if

        if (axis(2) == (0//1) .and. axis(3) == (0//1) .and. (axis(1) * U(3,2)) > (0//1)) then
            positive = .true.
        else if ( (U(2,1) * axis(3) - U(3,1) * axis(2)) > (0//1) ) then
            positive = .true.
        else
            positive = .false.
        end if

    end function Positive_Sense_of_Rotation

    !---- Subroutines ----!

    !!---- Subroutine Get_A_Matrices_Crys(LaueClass,A,n)
    !!----      character(len=5),                 intent(in)  :: LaueClass
    !!----      type(rational), dimension(3,3,6), intent(out) :: A
    !!----      integer,                          intent(out) :: n
    !!----
    !!---- Build A matrices -see Acta Cryst. A55 383-395 (1999) -.
    !!----
    !!---- Updated: September - 2018
    !!

    subroutine Get_A_Matrices_Crys(LaueClass,A,n)

        !---- Arguments ----!
        character(len=*),                 intent(in)  :: LaueClass
        type(rational), dimension(3,3,6), intent(out) :: A
        integer,                          intent(out) :: n

        !---- Local variables ----!
        integer                          :: i
        type(rational), dimension(3,3)   :: R1,R2,R3
        type(rational), dimension(3,3,6) :: Ainv

        select case (trim(LaueClass))

            case ("2/m")
                n  = 6
                R1(1,1:3) = [ 1//1, 0//1, 0//1 ]
                R1(2,1:3) = [ 0//1, 1//1, 0//1 ]
                R1(3,1:3) = [ 0//1, 0//1, 1//1 ]
                R2(1,1:3) = [ 0//1, 0//1, 1//1 ]
                R2(2,1:3) = [ 0//1,-1//1, 0//1 ]
                R2(3,1:3) = [ 1//1, 0//1, 0//1 ]
                R3(1,1:3) = [-1//1, 0//1, 1//1 ]
                R3(2,1:3) = [ 0//1, 1//1, 0//1 ]
                R3(3,1:3) = [-1//1, 0//1, 0//1 ]

                Ainv(:,:,1) = R1(:,:)
                Ainv(:,:,2) = R3(:,:)
                Ainv(:,:,3) = matmul(R3,R3)
                Ainv(:,:,4) = R2
                Ainv(:,:,5) = matmul(R2,R3)
                Ainv(:,:,6) = matmul(Ainv(:,:,5),R3)
            case ("mmm")
                n  = 6
                R1(1,1:3) = [ 1//1, 0//1, 0//1 ]
                R1(2,1:3) = [ 0//1, 1//1, 0//1 ]
                R1(3,1:3) = [ 0//1, 0//1, 1//1 ]
                R2(1,1:3) = [ 0//1, 1//1, 0//1 ]
                R2(2,1:3) = [ 1//1, 0//1, 0//1 ]
                R2(3,1:3) = [ 0//1, 0//1,-1//1 ]
                R3(1,1:3) = [ 0//1, 0//1, 1//1 ]
                R3(2,1:3) = [ 1//1, 0//1, 0//1 ]
                R3(3,1:3) = [ 0//1, 1//1, 0//1 ]

                Ainv(:,:,1) = R1(:,:)
                Ainv(:,:,2) = R3(:,:)
                Ainv(:,:,3) = matmul(R3,R3)
                Ainv(:,:,4) = R2
                Ainv(:,:,5) = matmul(R2,R3)
                Ainv(:,:,6) = matmul(Ainv(:,:,5),R3)
            case ("m-3")
                n = 2
                R1(1,1:3) = [ 1//1, 0//1, 0//1 ]
                R1(2,1:3) = [ 0//1, 1//1, 0//1 ]
                R1(3,1:3) = [ 0//1, 0//1, 1//1 ]
                R2(1,1:3) = [ 0//1,-1//1, 0//1 ]
                R2(2,1:3) = [ 1//1, 0//1, 0//1 ]
                R2(3,1:3) = [ 0//1, 0//1, 1//1 ]
                Ainv(:,:,1) = R1(:,:)
                Ainv(:,:,2) = R2(:,:)
            case default
                n = 1
                R1(1,1:3) = [ 1//1, 0//1, 0//1 ]
                R1(2,1:3) = [ 0//1, 1//1, 0//1 ]
                R1(3,1:3) = [ 0//1, 0//1, 1//1 ]
                Ainv(:,:,1) = R1(:,:)
        end select

        do i = 1 , n
            call Rational_Inv_Matrix(Ainv(:,:,i),A(:,:,i))
        end do

    end subroutine Get_A_Matrices_Crys

    !!---- Subroutine Get_A_Matrices_Shub(laueClass,A,n)
    !!----      character(len=*),                 intent(in)  :: laueClass
    !!----      type(rational), dimension(3,3,6), intent(out) :: A
    !!----      integer,                          intent(out) :: n
    !!----
    !!---- Similar to Get_A_Matrices_Crys, but adapted to Shubnikov groups
    !!----
    !!---- Updated: September - 2018
    !!

    subroutine Get_A_Matrices_Shub(laueClass,A,n)

        ! Matrices for different settings

        !---- Arguments ----!
        character(len=*),                 intent(in)  :: laueClass
        type(rational), dimension(3,3,6), intent(out) :: A
        integer,                          intent(out) :: n

        n = 1
        ! setting a,b,c
        A(1,1:3,1) = [ 1, 0, 0 ]
        A(2,1:3,1) = [ 0, 1, 0 ]
        A(3,1:3,1) = [ 0, 0, 1 ]

        select case (trim(laueClass))
            case ("2/m") ! Monoclinic
                n = 6
                ! setting c,b,-a-c
                A(1,1:3,2) = [ 0, 0,-1 ]
                A(2,1:3,2) = [ 0, 1, 0 ]
                A(3,1:3,2) = [ 1, 0,-1 ]
                ! setting -a-c,b,a
                A(1,1:3,3) = [-1, 0, 1 ]
                A(2,1:3,3) = [ 0, 1, 0 ]
                A(3,1:3,3) = [-1, 0, 0 ]
                ! setting c,-b,a
                A(1,1:3,4) = [ 0, 0, 1 ]
                A(2,1:3,4) = [ 0,-1, 0 ]
                A(3,1:3,4) = [ 1, 0, 0 ]
                ! setting -a-c,-b,c
                A(1,1:3,5) = [-1, 0, 0 ]
                A(2,1:3,5) = [ 0,-1, 0 ]
                A(3,1:3,5) = [-1, 0, 1 ]
                ! setting a,-b,-a-c
                A(1,1:3,6) = [ 1, 0,-1 ]
                A(2,1:3,6) = [ 0,-1, 0 ]
                A(3,1:3,6) = [ 0, 0,-1 ]
            case ("mmm") ! Orthorhombic
                n = 6
                ! setting b,c,a
                A(1,1:3,2) = [ 0, 0, 1 ]
                A(2,1:3,2) = [ 1, 0, 0 ]
                A(3,1:3,2) = [ 0, 1, 0 ]
                ! setting c,a,b
                A(1,1:3,3) = [ 0, 1, 0 ]
                A(2,1:3,3) = [ 0, 0, 1 ]
                A(3,1:3,3) = [ 1, 0, 0 ]
                ! setting a,c,-b
                A(1,1:3,4) = [ 1, 0, 0 ]
                A(2,1:3,4) = [ 0, 0,-1 ]
                A(3,1:3,4) = [ 0, 1, 0 ]
                ! setting c,b,-a
                A(1,1:3,5) = [ 0, 0,-1 ]
                A(2,1:3,5) = [ 0, 1, 0 ]
                A(3,1:3,5) = [ 1, 0, 0 ]
                ! setting b,a,-c
                A(1,1:3,6) = [ 0, 1, 0 ]
                A(2,1:3,6) = [ 1, 0, 0 ]
                A(3,1:3,6) = [ 0, 0,-1 ]
            case ("-3","-3 R","-3m","-3m R","-3m1","-31m") ! Trigonal
                n = 2
                ! reverse -> obverse setting
                A(1,1:3,2) = [-1, 0, 0 ]
                A(2,1:3,2) = [ 0,-1, 0 ]
                A(3,1:3,2) = [ 0, 0, 1 ]
        end select

    end subroutine Get_A_Matrices_Shub

    !!---- Subroutine Get_Generators(spaceGroupNumber,symOp,nSymOp,G,nGen)
    !!----      character(len=*),                        intent(in)  :: laueClass
    !!----      type(Symm_Oper_Type), dimension(nSymOp), intent(in)  :: symOp
    !!----      integer,                                 intent(in)  :: nSymOp
    !!----      type(Symm_Oper_Type), dimension(3),      intent(out) :: G
    !!----      integer,                                 intent(out) :: nGen
    !!----
    !!---- Returns the generators for the space group = spaceGroupNumber
    !!----
    !!---- Updated: September - 2018
    !!

    subroutine Get_Generators(laueClass,symOp,nSymOp,G,nGen)

        character(len=*),                        intent(in)  :: laueClass ! Laue class
        type(Symm_Oper_Type), dimension(nSymOp), intent(in)  :: symOp     ! symmetry operations
        integer,                                 intent(in)  :: nSymOp    ! number of symmetry operations
        type(Symm_Oper_Type), dimension(3),      intent(out) :: G         ! generators
        integer,                                 intent(out) :: nGen      ! number of generators

        integer                              :: i,n,ngaux,inversion
        type(rational), dimension(3)         :: axis
        type(rational), dimension(3,3)       :: identity
        integer, dimension(:,:), allocatable :: idd

        err_std = .false.

        ! Initialization
        nGen          = 0
        inversion     = 0
        call Rational_Identity_Matrix(3,identity)
        allocate(idd(nSymOp,2))
        do i = 1 , 3
            allocate(G(i)%Mat(4,4))
        end do

        ! Look for an inversion center
        do i = 1 , nSymOp
            if (Equal_Rational_Matrix(-symOp(i)%Mat(1:3,1:3),identity)) then
                inversion = i
                exit
            end if
        end do

        select case(trim(laueClass))

            case ("-1") ! Triclinic

                if (inversion == 0) then
                    nGen = 1
                    ! Search for the onefold axis
                    call Get_Rotations(symOp(:),nSymOp,1,n,idd)
                    G(1) = symOp(idd(1,1))
                end if

            case ("2/m") ! Monoclinic

                nGen = 1
                ! Search for a twofold axis
                call Get_Rotations(symOp(:),nSymOp,2,n,idd)
                G(1) = symOp(idd(1,1))

            case ("mmm")! Orthorhombic

                nGen = 2
                ! Search for the two fold axes along [001] and [010]
                call Get_Rotations(symOP(:),nSymOp,2,n,idd)
                ngaux = 0
                do i = 1 , n
                    call Get_Rotation_Axis(symOp(idd(i,1))%Mat(1:3,1:3),axis)
                    if ((axis(1) == (0//1) .and. axis(2) == (0//1) .and. axis(3) == (1//1)) .or. &
                        (axis(1) == (0//1) .and. axis(2) == (1//1) .and. axis(3) == (0//1))) then
                        ngaux = ngaux + 1
                        G(ngaux) = symOp(idd(i,1))
                        if (ngaux == 2) exit
                    end if
                end do

            case ("4/m","4/mmm") ! Tetragonal

                ! Search for the fourfold axis along [001]
                call Get_Rotations(symOp(:),nSymOp,4,n,idd)
                do i = 1 , n
                    call Get_Rotation_Axis(symOp(idd(i,1))%Mat(1:3,1:3),axis)
                    if (axis(1) == (0//1) .and. axis(2) == (0//1) .and. axis(3) == (1//1)) then
                        if (Positive_Sense_of_Rotation(symOp(idd(i,1))%Mat(1:3,1:3),axis)) then
                            G(1) = symOp(idd(i,1))
                            nGen = 1
                            exit
                        end if
                    end if
                end do
                ! Look for a possible twofold axis along [100]
                call Get_Rotations(symOp(:),nSymOp,2,n,idd)
                do i = 1 , n
                    call Get_Rotation_Axis(symOp(idd(i,1))%Mat(1:3,1:3),axis)
                    if (axis(1) == (1//1) .and. axis(2) == (0//1) .and. axis(3) == (0//1)) then
                        G(2) = symOp(idd(i,1))
                        nGen = 2
                        exit
                    end if
                end do

            case ("-3","-3 R","-3m","-3m R","-3m1","-31m") ! Trigonal

                ! Search for the threefold axis along [001]
                call Get_Rotations(symOP(:),nSymOp,3,n,idd)
                do i = 1 , n
                    call Get_Rotation_Axis(symOp(idd(i,1))%Mat(1:3,1:3),axis)
                    if (axis(1) == (0//1) .and. axis(2) == (0//1) .and. axis(3) == (1//1)) then
                        if (Positive_Sense_of_Rotation(symOp(idd(i,1))%Mat(1:3,1:3),axis)) then
                            G(1) = symOp(idd(i,1))
                            nGen     = 1
                            exit
                        end if
                    end if
                end do
                ! Search for a possible twofold axis along [110] or [-110]
                call Get_Rotations(symOp(:),nSymOp,2,n,idd)
                do i = 1 , n
                    call Get_Rotation_Axis(symOp(idd(i,1))%Mat(1:3,1:3),axis)
                    if ((axis(1) == (1//1) .and. axis(2) == (1//1) .and. axis(3) == (0//1)) .or. &
                        (axis(1) == (-1//1) .and. axis(2) == (1//1) .and. axis(3) == (0//1))) then
                        G(2) = symOp(idd(i,1))
                        nGen     = 2
                        exit
                    end if
                end do

            case ("6/m","6/mmm") ! Hexagonal

                ! Search for the sixfold axis along [001] in SGtarget
                call Get_Rotations(symOp(:),nSymOp,6,n,idd)
                do i = 1 , n
                    call Get_Rotation_Axis(symOp(idd(i,1))%Mat(1:3,1:3),axis)
                    if (axis(1) == (0//1) .and. axis(2) == (0//1) .and. axis(3) == (1//1)) then
                        if (Positive_Sense_of_Rotation(symOp(idd(i,1))%Mat(1:3,1:3),axis)) then
                            G(1) = symOp(idd(i,1))
                            nGen     = 1
                            exit
                        end if
                    end if
                end do
                ! Look for a possible twofold axis along [-110] in SGtarget
                call Get_Rotations(symOp(:),nSymOp,2,n,idd)
                do i = 1 , n
                    call Get_Rotation_Axis(symOp(idd(i,1))%Mat(1:3,1:3),axis)
                    if (axis(1) == (-1//1) .and. axis(2) == (1//1) .and. axis(3) == (0//1)) then
                        G(2) = symOp(idd(i,1))
                        nGen         = 2
                        exit
                    end if
                end do

            case ("m3","m-3","m3m","m-3m") ! Cubic

                ! Search for the fourfold axis along [001] in SGtarget
                call Get_Rotations(symOp(:),nSymOp,4,n,idd)
                do i = 1 , n
                    call Get_Rotation_Axis(symOp(idd(i,1))%Mat(1:3,1:3),axis)
                    if (axis(1) == (0//1) .and. axis(2) == (0//1) .and. axis(3) == (1//1)) then
                        if (Positive_Sense_of_Rotation(symOp(idd(i,1))%Mat,axis)) then
                            G(1) = symOp(idd(i,1))
                            nGen = 1
                            exit
                        end if
                    end if
                end do
                if (nGen == 0) then
                    ! Search for the twofold axis along [001] in SGtarget
                    call Get_Rotations(symOp(:),nSymOp,2,n,idd)
                    do i = 1 , n
                        call Get_Rotation_Axis(symOp(idd(i,1))%Mat(1:3,1:3),axis)
                        if (axis(1) == (0//1) .and. axis(2) == (0//1) .and. axis(3) == (1//1)) then
                            G(1) = symOp(idd(i,1))
                            nGen = 1
                            exit
                        end if
                    end do
                end if
                ! Search for a threefold axis along {111} in SGtarget
                call Get_Rotations(symOp(:),nSymOp,3,n,idd)
                do i = 1 , n
                    call Get_Rotation_Axis(symOp(idd(i,1))%Mat(1:3,1:3),axis)
                    if (axis(1) == (1//1) .and. axis(2) == (1//1) .and. axis(3) == (1//1)) then
                        if (Positive_Sense_of_Rotation(symOp(idd(i,1))%Mat(1:3,1:3),axis)) then
                            G(2) = symOp(idd(i,1))
                            nGen = 2
                            exit
                        end if
                    end if
                end do

            case default
                err_std = .true.
                err_std_mess = "Error in Match_Shubnikov_Group. Unknown Laue class"
                return

        end select

        if (inversion > 0) then
            ! Choose proper rotations if the spacegroup is centrosymmetric
            !do i = 1 , nGen
            !    det = rational_determinant(G(i)%Mat(1:3,1:3))
            !    G(i)%Mat(1:3,1:3) = det * G(i)%Mat(1:3,1:3)
            !    !G(i)%time_inv     = symOp(inversion)%time_inv * G(i)%time_inv
            !    write(*,*) "Hola ",symOp(inversion)%time_inv,G(i)%time_inv,G(i)%time_inv*symOp(inversion)%time_inv
            !end do
            nGen = nGen + 1
            G(nGen) = symOp(inversion)
        end if

    end subroutine Get_Generators

    !!---- Subroutine Get_HM_Standard(numSpg,symbolHM)
    !!----      integer,           intent(in)  :: numSpg
    !!----      character(len=12), intent(out) :: symbolHM
    !!----
    !!---- Returns the Herman-Maugin symbol for the standard
    !!----
    !!---- Updated: September - 2018
    !!----

    subroutine Get_HM_Standard(numSpg,symbolHM)

        !---- Arguments ----!
        integer,           intent(in)  :: numSpg
        character(len=12), intent(out) :: symbolHM

        !---- Local variables ----!
        integer i,n
        integer, dimension(1) :: posSep

        if (numSpg < 0 .or. numSpg > 230) then
            write(*,'(a)') "Error in Get_HM_Standard subroutine: numSpg out of range"
            return
        end if

        call Set_Spgr_Info()
        i = 1
        do
            if (spgr_info(i)%n == numSpg) then
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

    !!---- Subroutine Get_Lattice_Type(L,Latc,lattyp)
    !!----      integer,                        intent( in) :: L
    !!----      type(rational), dimension(:,:), intent( in) :: Latc
    !!----      character(len=*),               intent(out) :: lattyp
    !!----
    !!---- Returns the lattice type symbol from lattice centring vectors
    !!----
    !!---- Updated: September - 2018
    !!

    subroutine Get_Lattice_Type(L,Latc,lattyp)

        !---- Arguments ----!
        integer,                        intent( in) :: L
        type(rational), dimension(:,:), intent( in) :: Latc
        character(len=*),               intent(out) :: lattyp

        !---- Local variables ----!
        integer :: i,j,nlat
        logical :: latt_p,latt_a,latt_b,latt_c,latt_i,latt_r,latt_s,latt_h,latt_f,latt_z
        integer, dimension(10):: latt_given

        type(rational), dimension(3,10) :: lattice

        lattice(:,1)  = [ 0//1,1//2,1//2 ]
        lattice(:,2)  = [ 1//2,0//1,1//2 ]
        lattice(:,3)  = [ 1//2,1//2,0//1 ]
        lattice(:,4)  = [ 1//2,1//2,1//2 ]
        lattice(:,5)  = [ 2//3,1//3,1//3 ]
        lattice(:,6)  = [ 1//3,2//3,2//3 ]
        lattice(:,7)  = [ 1//3,2//3,1//3 ]
        lattice(:,8)  = [ 2//3,1//3,2//3 ]
        lattice(:,9)  = [ 2//3,1//3,0//1 ]
        lattice(:,10) = [ 1//3,2//3,0//1 ]

        ! Check for primitive and non conventional centrings
        nlat   = 0
        latt_p = .true.
        do i = 1 , L
            if (.not. IsInteger(latc(1:3,i))) then
                latt_p = .false.
                nlat   = nlat + 1
            end if
        end do
        if (latt_p) then    ! primitive
            lattyp="P"
            return
        end if
        if (nlat > 3) then  !non conventional centring
            lattyp="Z"
            err_std = .true.
            err_std_mess = "Error in Get_Lattice_Type. Number of centring vectors > 3."
            return
        end if

        latt_a=.false.
        latt_b=.false.
        latt_c=.false.
        latt_i=.false.
        latt_r=.false.
        latt_s=.false. ! Rombohedral,reverse setting
        latt_h=.false. ! Hexagonally centred
        latt_f=.false.
        latt_z=.false.

        do i = 1 , L
            latt_given(:) = 0
            do j = 1 , 10
                if (equal_rational_vector(latc(1:3,i),lattice(1:3,j))) then
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
            lattyp  = "Z"
            err_std = .true.
            err_std_mess = "Error in Get_Lattice_Type. Unable to identify lattice type."
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

        if (latt_a) lattyp="A"
        if (latt_b) lattyp="B"
        if (latt_c) lattyp="C"
        if (latt_i) lattyp="I"
        if (latt_r) lattyp="R"
        if (latt_s) lattyp="S"
        if (latt_h) lattyp="H"
        if (latt_f) lattyp="F"

    end subroutine Get_Lattice_Type

    !!---- Subroutine Get_Lattice_Type_from_M(M,lattyp,info)
    !!----      type(rational), dimension(3,3), intent(in)  :: M
    !!----      character,                      intent(out) :: lattyp
    !!----      logical, optional,              intent(in)  :: output
    !!----
    !!---- The M matrix transforms a primitive basis to a standard
    !!---- basis. The columns of the inverse of the M matrix contain
    !!---- the primitive vectors of the lattice expressed in the
    !!---- standard basis.
    !!----
    !!---- Updated: September - 2018
    !!----

    subroutine Get_Lattice_Type_from_M(M,lattyp,output)

        !---- Arguments ----!
        type(rational), dimension(3,3), intent(in)  :: M
        character,                      intent(out) :: lattyp
        logical, optional,              intent(in)  :: output

        !---- Local variables ----!
        integer   :: i,j,n,nLatticePoints,nCentringVectors
        logical   :: newt
        type(rational) :: det
        logical,        dimension(3)   :: isOrigin
        type(rational), dimension(3)   :: t
        type(rational), dimension(3,3) :: Minv,Maux
        type(rational), dimension(3,3) :: latc
        logical :: pout
        pout = .false.
        if(present(output)) pout=output
        err_std = .false.
        call Rational_Inv_Matrix(M,Minv)
        det = rational_determinant(Minv)
        det = (1//1) / det

        if (mod(det%Numerator,det%Denominator) /= 0) then
            err_std = .true.
            err_std_mess = "Error in Get_Lattice_Type_from_M. &
            Inverse of the determinant is not integer"
            return
        end if

        nLatticePoints = det%Numerator / det%Denominator

        if ( nLatticePoints > 1 ) then
            if (pout) write(*,'(8x,a,i2,1x,a)') " => Standard lattice is non primitive.",&
                                        nLatticePoints,"lattice points in the unit cell"
            nCentringVectors = 0
            Maux             = Minv ! Maux stores lattice points inside the unit cell
            isOrigin         = .True.
            do i = 1 , 3
                do j = 1 , 3
                    if (mod(Maux(i,j)%Numerator,Maux(i,j)%Denominator) == 0) then ! integer component
                        Maux(i,j) = 0 // 1
                    else ! fractional component
                        n = Maux(i,j)%Numerator / Maux(i,j)%Denominator
                        Maux(i,j) = Maux(i,j) - (n//1)
                        if (Maux(i,j) < (0//1)) Maux(i,j) = Maux(i,j) + (1//1)
                        isOrigin(j) = .False.
                    end if
                end do
            end do
            do i = 1 , 3
                if ( .not. isOrigin(i) ) then
                    t = Maux(:,i) ! t is a centring vector
                    newt = .true.
                    do j = 1 , nCentringVectors
                       if (Rational_Colinear(t,latc(:,j),3)) newt = .false.
                    end do
                    if (newt) then
                       nCentringVectors = nCentringVectors + 1
                       latc(:,nCentringVectors) = t
                    end if
                end if
            end do
            if (nCentringVectors > 0) then
                call Get_Lattice_Type(nCentringVectors,latc,lattyp)
            else
                lattyp = "P"
            end if
            if (err_std) return
            if (pout) write(*,'(12x,2a)') "Lattice type: ",lattyp
        else
            if (pout) write(*,'(8x,a)') " => Standard lattice is primitive "
            lattyp = "P"
        end if

    end subroutine Get_Lattice_Type_from_M

    !!---- Subroutine Get_Magnetic_Lattice_Type(nLat,latV,naLat,alatV,lattyp)
    !!----      integer,                            intent(in)  :: nLat
    !!----      type(rational), dimension(3,nLat),  intent(in)  :: latV
    !!----      integer,                            intent(in)  :: naLat
    !!----      type(rational), dimension(3,naLat), intent(in)  :: alatV
    !!----      character(len=1), dimension(2),     intent(out) :: lattyp
    !!----
    !!---- Returns the magnetic type symbol from lattice and antilattice translations.
    !!----
    !!---- Updated: September - 2018
    !!

    subroutine Get_Magnetic_Lattice_Type(G)!nLat,latV,naLat,alatV,lattyp)

        !---- Arguments ----!
        !integer,                              intent(in)  :: nLat
        !type(rational),   dimension(3,nLat),  intent(in)  :: latV
        !integer,                              intent(in)  :: naLat
        !type(rational),   dimension(3,naLat), intent(in)  :: alatV
        !character(len=1), dimension(2),       intent(out) :: lattyp
        class(spg_type), intent(in out) :: G

        !---- Local variables ----!
        integer :: i,j
        logical :: latt_p,latt__a,latt__b,latt__c,latt_a,latt_b,latt_c,latt_i,latt_r
        integer,          dimension(9)   :: latt_given
        type(rational),   dimension(3,9) :: lattice

        err_std = .false.
        G%shu_lat(:)=" "
        if (G%num_lat > 0) then
            call Get_Lattice_Type(G%num_Lat,G%Lat_tr,G%shu_lat(1))
        else
            G%shu_lat(1) = "P"
        end if
        if (err_std)    return
        if (G%num_aLat == 0) return

        lattice(:,1)  = [ 1//2,0//1,0//1 ]
        lattice(:,2)  = [ 0//1,1//2,0//1 ]
        lattice(:,3)  = [ 0//1,0//1,1//2 ]
        lattice(:,4)  = [ 0//1,1//2,1//2 ]
        lattice(:,5)  = [ 1//2,0//1,1//2 ]
        lattice(:,6)  = [ 1//2,1//2,0//1 ]
        lattice(:,7)  = [ 1//2,1//2,1//2 ]
        lattice(:,8)  = [ 2//3,1//3,5//6 ]
        lattice(:,9)  = [ 1//3,2//3,1//6 ]

        latt_p  = .false.
        latt__a = .false.
        latt__b = .false.
        latt__c = .false.
        latt_a  = .false.
        latt_b  = .false.
        latt_c  = .false.
        latt_i  = .false.
        latt_r  = .false.

        do i = 1 , G%num_aLat
            latt_given(:) = 0
            do j = 1 , 9
                if (equal_rational_vector(G%aLat_tr(1:3,i),lattice(1:3,j))) then
                    latt_given(j) = 1
                    select case (j)
                    case (1)
                        latt__a = .true.
                    case (2)
                        latt__b = .true.
                    case (3)
                        latt__c = .true.
                    case (4)
                        latt_a  = .true.
                    case (5)
                        latt_b  = .true.
                    case (6)
                        latt_c  = .true.
                    case (7)
                        latt_i  = .true.
                    case (8,9)
                        latt_r  = .true.
                    end select
                    exit
                end if
            end do
            if (sum(latt_given) == 0) then
                G%shu_lat(2) = "Z"
                exit
            end if
        end do

        if (G%shu_lat(2) == " ") then
            G%shu_lat(2) = "Z"
            select case (G%shu_lat(1))
                case ("P")
                    if (sum(latt_given) > 1) then
                        err_std      = .true.
                        err_std_mess = "Error in Get_Magnetic_Lattice_Type. Primitive lattice with more than one anti-translation."
                        return
                    end if
                    if (latt__a) then
                        G%shu_lat(2) = "a"
                        return
                    else if (latt__b) then
                        G%shu_lat(2) = "b"
                        return
                    else if (latt__c) then
                        G%shu_lat(2) = "c"
                        return
                    else if (latt_a) then
                        G%shu_lat(2) = "A"
                        return
                    else if (latt_b) then
                        G%shu_lat(2) = "B"
                    else if (latt_c) then
                        G%shu_lat(2) = "C"
                        return
                    else if (latt_i) then
                        G%shu_lat(2) = "I"
                    end if
                case ("A")
                    if (latt__a .or. latt_i) then
                        G%shu_lat(2) = "a"
                    else if (latt__b .or. latt__c) then
                        G%shu_lat(2) = "b"
                    else if (latt_b .or. latt_c) then
                        G%shu_lat(2) = "B"
                    end if
                case ("C")
                    if (latt__a .or. latt__b) then
                        G%shu_lat(2) = "a"
                    else if (latt__c .or. latt_i) then
                        G%shu_lat(2) = "c"
                    else if (latt_a .or. latt_b) then
                        G%shu_lat(2) = "A"
                    end if
                case ("I")
                    if (latt__a .or. latt_a) then
                        G%shu_lat(2) = "a"
                    else if (latt__b .or. latt_b) then
                        G%shu_lat(2) = "b"
                    else if (latt__c .or. latt_c) then
                        G%shu_lat(2) = "c"
                    end if
                case ("F")
                    if (latt__a .or. latt__b .or. latt__c .or. latt_i) G%shu_lat(2) = "S"
                case ("R")
                    if (latt__c .or. latt_r) G%shu_lat(2) = "I"
            end select
        end if

        if (G%shu_lat(2) == "Z") then
            err_std = .true.
            err_std_mess = "Error in Get_Magnetic_Lattice_Type. Unable to identify lattice type."
        end if

    end subroutine Get_Magnetic_Lattice_Type

    !!---- Subroutine Get_Mc_Matrix
    !!----      character(len=5),               intent(in)  :: LaueClass
    !!----      type(rational), dimension(3,3), intent(in)  :: Mp
    !!----      type(rational), dimension(3,3), intent(out) :: Mc
    !!----      logical, optional,              intent(in)  :: output
    !!----
    !!---- Computes a correction matrix if necessary. A correction needed
    !!---- in some cases for cubic groups with primitive lattices has been
    !!---- introduced in the subroutine Get_A_Matrices_Crys, for the case m-3.
    !!----
    !!---- Updated: September - 2018
    !!

    subroutine Get_Mc_Matrix(LaueClass,Mp,Mc,output)

        !---- Arguments ----!
        character(len=*),               intent(in)  :: LaueClass
        type(rational), dimension(3,3), intent(in)  :: Mp
        type(rational), dimension(3,3), intent(out) :: Mc
        logical, optional,              intent(in)  :: output

        !---- Local variables ----!
        character                      :: lattyp
        character(len=60)              :: symb
        type(rational), dimension(3,3) :: Mcinv
        type(rational), dimension(4,4) :: McAux
        logical :: pout

        pout=.false.
        if(present(output)) pout=output

        if (pout) then
            write(*,'(8x,a)') " => Constructing (Mc,0) matrix...."
            write(*,'(8x,a)') "    This matrix corrects (M',0) if (M',0) does not &
                                   bring the system to a standard setting"
        end if

        select case (trim(LaueClass))

            case ("2/m")  ! Put the two fold axis along b

                Mcinv(1,:) = [ 0//1, 1//1, 0//1 ]
                Mcinv(2,:) = [ 0//1, 0//1, 1//1 ]
                Mcinv(3,:) = [ 1//1, 0//1, 0//1 ]
                call Rational_Inv_Matrix(Mcinv,Mc)

            case ("4/m","4/mmm")

                call Get_Lattice_Type_from_M(Mp,lattyp)
                if (err_std) return
                if (lattyp == "C") then  ! C -> P
                    Mcinv(1,:) = [ 1//1, 1//1, 0//1 ]
                    Mcinv(2,:) = [ 1//1,-1//1, 0//1 ]
                    Mcinv(3,:) = [ 0//1, 0//1,-1//1 ]
                    call Rational_Inv_Matrix(Mcinv,Mc)
                else if (lattyp == "F") then ! F -> I
                    Mcinv(1,:) = [ 1//1, 1//1, 0//1 ]
                    Mcinv(2,:) = [-1//1, 1//1, 0//1 ]
                    Mcinv(3,:) = [ 0//1, 0//1, 1//1 ]
                    call Rational_Inv_Matrix(Mcinv,Mc)
                else
                    Mc(1,:) = [ 1//1, 0//1, 0//1 ]
                    Mc(2,:) = [ 0//1, 1//1, 0//1 ]
                    Mc(3,:) = [ 0//1, 0//1, 1//1 ]
                end if

            case ("-3","-3 R")

                call Get_Lattice_Type_from_M(Mp,lattyp)
                if (err_std) return
                if (lattyp == "S") then ! reverse -> obverse setting
                    Mc(1,:) = [-1//1, 0//1, 0//1 ]
                    Mc(2,:) = [ 0//1,-1//1, 0//1 ]
                    Mc(3,:) = [ 0//1, 0//1, 1//1 ]
                else
                    Mc(1,:) = [ 1//1, 0//1, 0//1 ]
                    Mc(2,:) = [ 0//1, 1//1, 0//1 ]
                    Mc(3,:) = [ 0//1, 0//1, 1//1 ]
                end if

            case ("-3m","-3m R","-3m1","-31m")
                call Get_Lattice_Type_from_M(Mp,lattyp)
                if (err_std) return
                if (lattyp == "S") then ! reverse -> obverse setting
                    Mc(1,:) = [-1//1, 0//1, 0//1 ]
                    Mc(2,:) = [ 0//1,-1//1, 0//1 ]
                    Mc(3,:) = [ 0//1, 0//1, 1//1 ]
                else if (lattyp == "H") then ! H -> P
                    Mcinv(1,:) = [ 1//1, 1//1, 0//1 ]
                    Mcinv(2,:) = [-1//1, 2//1, 0//1 ]
                    Mcinv(3,:) = [ 0//1, 0//1, 1//1 ]
                    call Rational_Inv_Matrix(Mcinv,Mc)
                else
                    Mc(1,:) = [ 1//1, 0//1, 0//1 ]
                    Mc(2,:) = [ 0//1, 1//1, 0//1 ]
                    Mc(3,:) = [ 0//1, 0//1, 1//1 ]
                end if

            case ("6/mmm")

                call Get_Lattice_Type_from_M(Mp,lattyp)
                if (err_std) return
                if (lattyp == "H") then ! H -> P
                    Mcinv(1,:) = [ 1//1, 1//1, 0//1 ]
                    Mcinv(2,:) = [-1//1, 2//1, 0//1 ]
                    Mcinv(3,:) = [ 0//1, 0//1, 1//1 ]
                    call Rational_Inv_Matrix(Mcinv,Mc)
                else
                    Mc(1,:) = [ 1//1, 0//1, 0//1 ]
                    Mc(2,:) = [ 0//1, 1//1, 0//1 ]
                    Mc(3,:) = [ 0//1, 0//1, 1//1 ]
                end if

            case default

                Mc(1,:) = [ 1//1, 0//1, 0//1 ]
                Mc(2,:) = [ 0//1, 1//1, 0//1 ]
                Mc(3,:) = [ 0//1, 0//1, 1//1 ]

        end select

        if (pout) then
            write(*,'(12x,a)',advance='no') "Quasi-standard setting --> Standard setting transformation: "
            call set_identity_matrix(4)
            McAux = identity_matrix
            McAux(1:3,1:3) = Mc
            call Get_Symb_Op_from_Mat(transpose(McAux),symb,"abc")
            write(*,'(a)') trim(symb)
        end if

    end subroutine Get_Mc_Matrix

    !!---- Subroutine Get_Mp_Matrix(G,P,Mp,output)
    !!----     class(spg_type),                 intent(in)  :: G
    !!----     type(rational), dimension(3,3), intent(in)  :: P
    !!----     type(rational), dimension(3,3), intent(out) :: Mp
    !!----     logical, optional,              intent(in)  :: output
    !!----
    !!---- It returns a 3x3 matrix Mp which transforms the
    !!---- original setting in a setting where basis vectors
    !!---- are along the symmetry axes of the Laue class.
    !!
    !!----           [ bx1 by1 bz1 ]
    !!----      Mp = [ bx2 by2 bz2 ]
    !!----           [ bx3 by3 bz3 ]
    !!----
    !!---- Updated January - 2019
    !!

    subroutine Get_Mp_Matrix(G,P,Mp,output)

        !---- Arguments ----!
        class(spg_type),                 intent(in)  :: G
        type(rational), dimension(3,3), intent(in)  :: P
        type(rational), dimension(3,3), intent(out) :: Mp
        logical, optional,              intent(in)  :: output

        !---- Local variables ----!
        integer            :: i,j,k,n,n_,nCubicAxes
        character(len=256) :: symb
        logical            :: colinear,standard
        integer,           dimension(2)   :: order
        character(len=20), dimension(6)   :: axisName
        type(rational),    dimension(3)   :: bx,by,bz
        type(rational),    dimension(3,3) :: Pinv,W,U,PM,PMinv
        type(rational),    dimension(3,4) :: vPerp,cubicAxes
        type(rational),    dimension(4,4) :: MpAux
        integer,           dimension(:,:), allocatable :: idd
        logical :: pout

        pout        = .false.
        if (present(output)) pout=output
        axisName(2) = "twofold"
        axisName(3) = "threefold"
        axisName(4) = "fourfold"
        axisName(6) = "sixfold"

        call Rational_Inv_Matrix(P,Pinv)

        if (pout) then
            write(*,'(8x,a)') " => Constructing (M',0) matrix...."
            write(*,'(12x,a)') "This matrix transforms the primitive basis in a basis with &
            vectors along the symmetry axes of the Laue class"
            write(*,'(12x,2a)') "Laue Class: ", G%laue
        end if

        select case (trim(G%laue))
            case ("-1")
                Mp = reshape ( [ 1//1,0//1,0//1,&
                                  0//1,1//1,0//1,&
                                  0//1,0//1,1//1 ],(/3,3/) )
            case ("2/m","4/m","-3","-3 R","6/m")
                select case (trim(G%laue))
                    case ("2/m")
                        order(1) = 2
                    case ("4/m")
                        order(1) = 4
                    case ("-3","-3 R","6/m")
                        order(1) = 3
                end select
                if (pout) write(*,'(12x,3a)') "Searching for the ",trim(axisName(order(1)))," axis..."
                allocate(idd(G%Multip,2))
                call Get_Rotations(G%op(:),G%Multip,order(1),n,idd)
                if (err_std) return
                ! Put the symmetry operation in the primitive setting
                W = MatMul(G%op(idd(1,1))%Mat(1:3,1:3),P)
                W = MatMul(Pinv,W)
                call Get_Rotation_Axis(W,bz)
                if (err_std) return
                if (pout) then
                    write(*,'(16x,a,i2,",",i2,",",i2,"]")') "Axis direction in the primitive setting: [", bz(1:3)%Numerator
                    write(*,'(12x,a)') "Searching for lattice vectors perpendicular to the rotation axis. &
                    Building the complete basis..."
                end if
                call Get_Vectors_Perpendicular_To_Rotation_Axis(W,vPerp)
                if (err_std) return
                call Get_Pseudo_Standard_Base(W,vPerp,bz,bx,by)
                Mp(:,1) = bx
                Mp(:,2) = by
                Mp(:,3) = bz
                deallocate(idd)
            case ("4/mmm")
                if (pout) write(*,'(12x,a)') "Searching for the fourfold axis..."
                allocate(idd(G%Multip,2))
                call Get_Rotations(G%op(:),G%Multip,4,n,idd)
                if (err_std) return
                ! Put the symmetry operation in the primitive setting
                W = MatMul(G%op(idd(1,1))%Mat(1:3,1:3),P)
                W = MatMul(Pinv,W)
                call Get_Rotation_Axis(W,bz)
                if (err_std) return
                if (pout) then
                    write(*,'(16x,a,i2,",",i2,",",i2,"]")') "Axis direction in the primitive setting: [", bz(:)%Numerator
                    write(*,'(12x,a)') "Searching for the twofold axis..."
                end if
                call Get_Rotations(G%op(:),G%Multip,2,n,idd)
                if (err_std) return
                colinear = .true.
                do i = 1 , n
                    ! Put the symmetry operation in the primitive setting
                    U = MatMul(G%op(idd(i,1))%Mat(1:3,1:3),P)
                    U = MatMul(Pinv,U)
                    call Get_Rotation_Axis(U,bx)
                    if (err_std) return
                    if (.not. rational_colinear(bz,bx,3)) then
                        colinear = .false.
                        exit
                    end if
                end do
                if (colinear) then
                    err_std = .true.
                    err_std_mess = "Error in Get_Mprime_Matrix. Unable to find a second rotation &
                    axis linearly independent from the first"
                    return
                end if
                if (pout) write(*,'(16x,a,i2,",",i2,",",i2,"]")') "Axis direction in the primitive setting: [", bx(:)%Numerator
                by = matmul(W,bx)
                Mp(:,1) = bx
                Mp(:,2) = by
                Mp(:,3) = bz
                deallocate(idd)
            case ("-3m","-3m R","-3m1","-31m","6/mmm")
                if (pout) write(*,'(12x,a)') "Searching for the threefold axis..."
                allocate(idd(G%Multip,2))
                call Get_Rotations(G%op(:),G%Multip,3,n,idd)
                if (err_std) return
                ! Put the symmetry operation in the primitive setting
                W = MatMul(G%op(idd(1,1))%Mat(1:3,1:3),P)
                W = MatMul(Pinv,W)
                if (idd(1,2) == -1)  W = -W
                call Get_Rotation_Axis(W,bz)
                if (err_std) return
                if (pout) then
                    write(*,'(16x,a,i2,",",i2,",",i2,"]")') "Axis direction in the primitive setting: [", bz(:)%Numerator
                    write(*,'(12x,a)') "Searching for the twofold axis..."
                end if
                call Get_Rotations(G%op(:),G%Multip,2,n,idd)
                if (err_std) return
                colinear = .true.
                do i = 1 , n
                    ! Put the symmetry operation in the primitive setting
                    U = MatMul(G%op(idd(i,1))%Mat(1:3,1:3),P)
                    U = MatMul(Pinv,U)
                    call Get_Rotation_Axis(U,bx)
                    if (err_std) return
                    if (.not. rational_colinear(bz,bx,3)) then
                        colinear = .false.
                        exit
                    end if
                end do
                if (colinear) then
                    err_std = .true.
                    err_std_mess = "Error in Get_Mprime_Matrix. Unable to find a second rotation &
                    axis linearly independent from the first"
                    return
                end if
                if (pout) write(*,'(16x,a,i2,",",i2,",",i2,"]")') "Axis direction in the primitive setting: [", bx(:)%Numerator
                by = matmul(W,bx)
                Mp(:,1) = bx
                Mp(:,2) = by
                Mp(:,3) = bz
                deallocate(idd)
            case ("mmm")
                if (pout) write(*,'(12x,a)') "Searching for the three twofold axis..."
                allocate(idd(G%Multip,2))
                call Get_Rotations(G%op(:),G%Multip,2,n,idd)
                if (err_std) return
                n_ = 0
                do i = 1 , n
                    ! Put the symmetry operation in the primitive setting
                    W = MatMul(G%op(idd(i,1))%Mat(1:3,1:3),P)
                    W = MatMul(Pinv,W)
                    call Get_Rotation_Axis(W,bz)
                    if (err_std) return
                    colinear = .false.
                    do j = 1 , n_
                        if (rational_colinear(bz,Mp(:,j),3)) colinear = .true.
                    end do
                    if (.not. colinear) then
                        n_ = n_ + 1
                        Mp(:,n_) = bz
                        if (pout) write(*,'(16x,a,i2,",",i2,",",i2,"]")') "Axis direction in the primitive setting: [", bz(:)%Numerator
                    end if
                    if (n_ == 3) exit
                end do
                deallocate(idd)
            case ("m3","m-3","m3m","m-3m")
                if (pout) write(*,'(12x,a)') "Searching for the four threefold axes..."
                allocate(idd(G%Multip,2))
                call Get_Rotations(G%op(:),G%Multip,3,n,idd)
                if (err_std) return
                nCubicAxes = 0
                standard   = .true. ! Initially we assume P is a standard basis
                do i = 1 , n
                    W = MatMul(G%op(idd(i,1))%Mat(1:3,1:3),P)
                    W = MatMul(Pinv,W)
                    call Get_Rotation_Axis(W,bz)
                    if (err_std) return
                    colinear = .false.
                    do j = 1 , nCubicAxes
                        if (rational_colinear(bz,cubicAxes(:,j),3)) colinear = .true.
                    end do
                    if (.not. colinear) then
                        nCubicAxes = nCubicAxes + 1
                        cubicAxes(:,nCubicAxes) = bz(:)
                        if (pout) write(*,'(16x,a,i2,",",i2,",",i2,"]")') "Axis direction in the primitive setting: [", bz(:)%Numerator
                        if (abs(bz(1)) /= (1//1) .or. &
                            abs(bz(2)) /= (1//1) .or. &
                            abs(bz(3)) /= (1//1)) standard = .false.
                    end if
                end do
                if (standard) then
                    Mp(:,1) = [ 1//1,0//1,0//1 ]
                    Mp(:,2) = [ 0//1,1//1,0//1 ]
                    Mp(:,3) = [ 0//1,0//1,1//1 ]
                else
                    ! Find a combination of cubicAxes for which the threefold axes are along {111}
                    do i = 0 , 1
                        if (standard) exit
                        bx(:) = cubicAxes(:,2) - ((2 * i)//1) * cubicAxes(:,2)
                        do j = 0 , 1
                            if (standard) exit
                            by(:) = cubicAxes(:,3) - ((2 * j)//1) * cubicAxes(:,3)
                            do k = 0 , 1
                                if (standard) exit
                                bz(:) = cubicAxes(:,4) - ((2 * k)//1) * cubicAxes(:,4)
                                Mp(:,1) = (1//2) * (cubicAxes(:,1) + bx)
                                Mp(:,2) = (1//2) * (cubicAxes(:,1) + by)
                                Mp(:,3) = (1//2) * (cubicAxes(:,1) + bz)
                                ! For lattice I, Mp is not integral (there is a lattice
                                ! point at the middle of each diagonal. We must multiply
                                ! by two in order to get the diagonal of the cube
                                if (.not. IsInteger(Mp)) Mp = (2//1) * Mp
                                PM = matmul(P,Mp)
                                call Rational_Inv_Matrix(PM,PMinv)
                                if (.not. Err_Rational) then
                                    standard = .true.
                                    do n_ = 1 , n
                                        !write(*,*) G%OP(idd(n_,1))%Mat(1:3,1)
                                        !write(*,*) G%OP(idd(n_,1))%Mat(1:3,2)
                                        !write(*,*) G%OP(idd(n_,1))%Mat(1:3,3)
                                        W = MatMul(G%op(idd(n_,1))%Mat(1:3,1:3),PM)
                                        W = MatMul(PMinv,W)
                                        if (IsInteger(W)) then
                                            call Get_Rotation_Axis(W,bz)
                                            if (abs(bz(1)) /= (1//1) .or. &
                                                abs(bz(2)) /= (1//1) .or. &
                                                abs(bz(3)) /= (1//1)) then
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
                    err_std = .true.
                    err_std_mess = "Error in Get_Mp_Matrix: &
                    Unable to build a standard setting for the cubic crystal system..."
                    return
                end if
        end select

        call Set_Right_Handedness(Mp)

        if (pout) then
            write(*,'(tr12,a)',advance='no') "Primitive setting --> Quasi-Standard setting transformation: "
            call Set_Identity_Matrix(4)
            MpAux          = identity_matrix
            MpAux(1:3,1:3) = Mp
            call Get_Symb_Op_from_Mat(transpose(MpAux),symb,"abc")
            write(*,'(a)') trim(symb)
        end if

    end subroutine Get_Mp_Matrix

    !!---- Subroutine Get_Origin_Shift(G,G_,ng,P,origShift,shift)
    !!----      type(symm_oper_type), dimension(ng), intent(in)  :: G
    !!----      type(symm_oper_type), dimension(ng), intent(in)  :: G_
    !!----      integer,                             intent(in)  :: ng
    !!----      type(rational), dimension(3,3),      intent(in)  :: P
    !!----      type(rational), dimension(3),        intent(out) :: origShift
    !!----      logical,                             intent(out) :: shift
    !!----
    !!---- Tries to make G = G_ by an origin shift. If a solution is found,
    !!---- it returns shift = .true. P is used to express G and G_ in a
    !!---- primitive basis
    !!----
    !!---- Updated: September - 2018
    !!

    subroutine Get_Origin_Shift(G,G_,ng,P_,origShift,shift)

        type(symm_oper_type), dimension(ng), intent(in)  :: G
        type(symm_oper_type), dimension(ng), intent(in)  :: G_
        integer,                             intent(in)  :: ng
        type(rational), dimension(3,3),      intent(in)  :: P_
        type(rational), dimension(3),        intent(out) :: origShift
        logical,                             intent(out) :: shift

        integer :: i,j,k,l,r,nr
        type(rational), dimension(3,3)    :: identity
        type(rational), dimension(4,4)    :: P,Pinv
        type(rational), dimension(4,4,ng) :: Gx,Gt
        type(rational), dimension(:,:), allocatable :: U,D,T,V
        type(rational), dimension(:,:), allocatable :: b,x

        call Rational_Identity_Matrix(4,P)
        call Rational_Identity_Matrix(3,identity)
        shift      = .true.
        P(1:3,1:3) = P_
        do i = 1 , ng
            Gx(:,:,i)  = G(i)%Mat
            Gt(:,:,i)  = G_(i)%Mat
        end do
        call Rational_Inv_Matrix(P,Pinv)
        ! Transform generators to a primitive setting
        do i = 1 , ng
            Gt(:,:,i) = matmul(Gt(:,:,i),P)
            Gt(:,:,i) = matmul(Pinv,Gt(:,:,i))
            Gx(:,:,i) = matmul(Gx(:,:,i),P)
            Gx(:,:,i) = matmul(Pinv,Gx(:,:,i))
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
                    U(r,l) = Gt(k,l,j) - identity(k,l)
                end do
            end do
        end do
        call Rational_SmithNormalForm(U,nr,3,D,T,V)
        b = matmul(T,b)
        do j = 1 , 3
            if (D(j,j) == (0//1)) then
                if (mod(b(j,1)%Numerator,b(j,1)%Denominator) /= 0) then
                    shift = .false.
                    return
                else
                    x(j,1) = (0//1)
                end if
            else
                x(j,1) = b(j,1) / D(j,j)
            end if
        end do
        do j = 4 , nr
            if (mod(b(j,1)%Numerator,b(j,1)%Denominator) /= 0) then
                shift = .false.
                return
            end if
        end do
        if (shift) then
            x = matmul(V,x)
            x = matmul(P(1:3,1:3),x)
            origShift = x(:,1)
        end if

    end subroutine Get_Origin_Shift

    !!---- Subroutine Get_P_Matrix(G,P,nospin,output)
    !!----     class(spg_type),                 intent(in)  :: G
    !!----     type(rational), dimension(3,3), intent(out) :: P
    !!----     logical, optional,              intent(in)  :: nospin
    !!----     logical,        optional,       intent(in)  :: output
    !!----
    !!---- It returns a 3x3 matrix P which transforms the
    !!---- current setting to a primitive setting
    !!----
    !!---- Updated January - 2019
    !!

    subroutine Get_P_Matrix(G,P,nospin,output)

        !---- Arguments ----!
        class(spg_type),                 intent(in)  :: G
        type(rational), dimension(3,3), intent(out) :: P
        logical, optional,              intent(in)  :: nospin
        logical, optional,              intent(in)  :: output

        !---- Local variables ----!
        integer                                     :: i,j,k,n,m
        integer                                     :: nLatt,nAuxVec,nCentringVec
        logical                                     :: primitive,linearDependence,sorted
        character(len=256)                          :: symb
        type(rational)                              :: det,rNumLat
        type(rational), dimension(3)                :: nullVec
        type(rational), dimension(4)                :: v
        type(rational), dimension(3,3)              :: Pinv,A
        type(rational), dimension(4,4)              :: PAux
        type(rational), dimension(:,:), allocatable :: auxVec,centringVec
        logical :: pout

        err_std     = .false.
        primitive   = .false.
        nullVec(:)  = 0//1
        pout        = .false.
        if (present(output)) pout=output

        if (pout) then
            write(*,'(8x,a)') " => Constructing (P,0) matrix...."
            write(*,'(12x,a)') "This matrix transforms the original basis in a primitive basis"
        end if

        nLatt = 1
        do i = 1 , G%num_lat
            if (.not. equal_rational_vector(G%Lat_tr(:,i),nullVec))  nLatt = nLatt + 1
        end do
        if (present(nospin)) then
            do i = 1 , G%num_alat
                if (.not. equal_rational_vector(G%aLat_tr(:,i),nullVec)) nLatt = nLatt + 1
            end do
        end if

        if (nLatt == 1) then
            ! Basis is already primitive
            if (pout) write(*,'(12x,a)') "The basis is already primitive"
            call Rational_Identity_Matrix(3,P)
            primitive = .true.
            return
        else
            ! Build an expanded list of centring vectors
            nAuxVec = 0
            if (present(nospin)) then
                rNumLat = (G%num_lat+G%num_alat+1)//1
                allocate(auxVec(3,G%num_lat+G%num_alat))
            else
                rNumLat = (G%num_lat+1)//1
                allocate(auxVec(3,G%num_lat))
            end if
            do i = 1 , G%num_lat
                if (.not. equal_rational_vector(G%lat_tr(:,i),nullVec)) then
                    nAuxVec = nAuxVec + 1
                    auxVec(:,nAuxVec) = G%lat_tr(:,i)
                end if
            end do
            if (present(nospin)) then
                do i = 1 , G%num_alat
                    if (.not. equal_rational_vector(G%alat_tr(:,i),nullVec)) then
                        nAuxVec = nAuxVec + 1
                        auxVec(:,nAuxVec) = G%alat_tr(:,i)
                    end if
                end do
            end if
            allocate(centringVec(4,10*nAuxVec))
            nCentringVec = 0
            do n = 1 , nAuxVec
                do i = 0 , 1
                    if (i == 1 .and. auxVec(1,n)%numerator == 0) cycle
                    v(1)%numerator   = auxVec(1,n)%numerator - i*auxVec(1,n)%denominator
                    v(1)%denominator = auxVec(1,n)%denominator
                    do j = 0 , 1
                        if (j == 1 .and. auxVec(2,n)%numerator == 0) cycle
                        v(2)%numerator   = auxVec(2,n)%numerator - j*auxVec(2,n)%denominator
                        v(2)%denominator = auxVec(2,n)%denominator
                        do k = 0 , 1
                            if (k == 1 .and. auxVec(3,n)%numerator == 0) cycle
                            v(3)%numerator   = auxVec(3,n)%numerator - k*auxVec(3,n)%denominator
                            v(3)%denominator = auxVec(3,n)%denominator
                            v(4) = dot_product(v(1:3),v(1:3))
                            ! Check if there is linear dependence with previous vectors
                            linearDependence = .False.
                            do m = 1 , nCentringVec
                                if (rational_colinear(v(1:3),centringVec(1:3,m),3)) Then
                                    linearDependence = .True.
                                    if (v(4) < centringVec(4,m)) centringVec(1:4,m) = v(1:4)
                                    exit
                                end if
                            end do
                            if (.not. linearDependence) then
                                nCentringVec = nCentringVec + 1
                                centringVec(:,nCentringVec) = v
                            end if
                        end do
                    end do
                end do
            end do
        end if
        ! Sort vectors from shorter to longer
        do
            sorted = .True.
            do i = 1 , ncentringVec - 1
                if (centringVec(4,i) > centringVec(4,i+1)) Then
                    v = centringVec(:,i)
                    centringVec(:,i)   = centringVec(:,i+1)
                    centringVec(:,i+1) = v
                    sorted = .False.
                end if
            end do
            if (sorted) exit
        end do
        ! Append the unit translations
        centringVec(:,ncentringVec+1) = [ 1//1, 0//1, 0//1, 1//1 ]
        centringVec(:,ncentringVec+2) = [ 0//1, 1//1, 0//1, 1//1 ]
        centringVec(:,ncentringVec+3) = [ 0//1, 0//1, 1//1, 1//1 ]
        ncentringVec = ncentringVec + 3

        ! Combine three vectors until a primitive setting is found
        primitive = .false.
        do i = 1 , nCentringVec
            if (primitive) exit
            do j = i + 1 , nCentringVec
                if (primitive) exit
                do k = j + 1 , nCentringVec
                    if (primitive) exit
                    P(:,1) = centringVec(1:3,i)
                    P(:,2) = centringVec(1:3,j)
                    P(:,3) = centringVec(1:3,k)
                    call Rational_Inv_Matrix(P,Pinv)
                    if (.not. Err_Rational) then
                        det = rational_determinant(P)
                        if (abs((1//1)/det) == rNumLat) then
                            primitive = .true.
                            if (det < 0//1) then
                                P(:,1) = -P(:,1)
                                call Rational_Inv_Matrix(P,Pinv)
                            end if
                            do n = 1 , G%Multip
                                A = MatMul(G%Op(n)%Mat(1:3,1:3),P)
                                A = MatMul(Pinv,A)
                                if (.not. IsInteger(A)) then
                                    primitive = .false.
                                    exit
                                end if
                            end do
                        end if
                    end if
                end do
            end do
        end do

        if (.not. primitive) then
            err_std = .true.
            err_std_mess = "Error in Get_P_Matrix subroutine. A primitive basis cannot be found"
            return
        else if (pout) then
            write(*,'(12x,a)',advance='no') "Original setting --> Primitive setting transformation: "
            call Rational_Identity_Matrix(4,Paux)
            PAux(1:3,1:3) = P
            call Get_Symb_Op_from_Mat(transpose(PAux),symb,"abc")
            write(*,'(a)') trim(symb)
        end if

    end subroutine Get_P_Matrix

    !!---- Subroutine Get_Pseudo_Standard_Base(W,perpAxis,bz,bx,by)
    !!----      type(rational), dimension(3,3), intent(in)  :: W
    !!----      type(rational), dimension(3,4), intent(in)  :: perpAxis
    !!----      type(rational), dimension(3),   intent(in)  :: bz
    !!----      type(rational), dimension(3),   intent(out) :: bx
    !!----      type(rational), dimension(3),   intent(out) :: by
    !!----
    !!---- Given the rotation matrix of the principal axis W, the shortest
    !!---- lattice vector along the axis, bz, and the four shortest lattice
    !!---- vectors perpendicular to that axis, perpAxis, computes the basis
    !!---- bx,by,bz which gives the smallest unit cell
    !!
    !!---- Updated: September - 2018
    !!

    subroutine Get_Pseudo_Standard_Base(W,perpAxis,bz,bx,by)

        !---- Arguments ----!
        type(rational), dimension(3,3), intent(in)  :: W
        type(rational), dimension(3,4), intent(in)  :: perpAxis
        type(rational), dimension(3),   intent(in)  :: bz
        type(rational), dimension(3),   intent(out) :: bx
        type(rational), dimension(3),   intent(out) :: by

        !---- Local variables ----!
        integer                           :: i,j,n,imin,order
        type(rational)                    :: detW,minDet
        type(rational), dimension(3,3)    :: A,B
        type(rational), dimension(:),   allocatable :: det
        integer,        dimension(:,:), allocatable :: test_
        type(rational), dimension(:,:), allocatable :: byAux

        err_std  = .false.
        detW     = rational_determinant(W)

        if (detW%Numerator == detW%Denominator) then
            A = W
        else if (detW%Numerator == -detW%Denominator) then
            A = -W
        else
            err_std = .true.
            err_std_mess = "Error in Get_Vectors_Perpendicular_To_Rotation_Axis subroutine. &
            Determinant is not +-1"
            return
        end if

        B(:,3)   = bz(:)
        call Get_Rotation_Order(A,order)

        select case (order)

            case (2)

                Allocate(det(6),test_(6,2))
                n = 0
                do i = 1 , 4
                    do j = i + 1 , 4
                        n = n + 1
                        test_(n,1) = i
                        test_(n,2) = j
                        B(:,1)     = perpAxis(:,i)
                        B(:,2)     = perpAxis(:,j)
                        det(n)     = rational_determinant(B)
                    end do
                end do

                ! Select the combination with the smallest determinant

                imin   = 1
                minDet = abs(det(1))
                do i = 2 , 6
                    if (abs(det(i)) < minDet) then
                        imin   = i
                        minDet = abs(det(i))
                    end if
                end do

                bx(:) = perpAxis(:,test_(imin,1))
                by(:) = perpAxis(:,test_(imin,2))

            case (3,4,6)

                allocate(det(4),test_(4,1),byAux(3,4))
                do i = 1 , 4
                    byAux(:,i) = matmul(A,perpAxis(:,i))
                    B(:,1)     = perpAxis(:,i)
                    B(:,2)     = byAux(:,i)
                    det(i)     = rational_determinant(B)
                end do

                imin = 1
                minDet = abs(det(1))
                do i = 2 , 4
                    if (abs(det(i)) < minDet) then
                        imin = i
                        minDet = abs(det(i))
                    end if
                end do

                bx(:) = perpAxis(:,imin)
                by(:) = byAux(:,imin)

        end select

    end subroutine Get_Pseudo_Standard_Base

    !!---- Subroutine Get_Rotation_Axis(W,axis)
    !!----    type(rational), dimension(3,3), intent(in)  :: W     !rotation matrix
    !!----    type(rational), dimension(3),   intent(out) :: axis  !shortest vector along the rotation axisP
    !!----
    !!---- It computes the shortest lattice vector along the rotation
    !!---- axis of the symmetry operation W
    !!----
    !!---- Updated: January - 2019
    !!

    subroutine Get_Rotation_Axis(W,axis)

        !---- Arguments ----!
        type(rational), dimension(3,3), intent(in)  :: W
        type(rational), dimension(3),   intent(out) :: axis

        !---- Local variables ----!
        integer        :: i,j,rnk,nzeros_aux
        logical        :: ordered
        type(rational) :: det
        type(rational), dimension(3,3) :: A,U
        integer,        dimension(3) :: nzeros
        type(rational), dimension(3) :: row

        det = rational_determinant(W)

        if (det%Numerator == det%Denominator) then
            A = W
        else if (det%Numerator == -det%Denominator) then
            A = -W
        else
            err_std = .true.
            err_std_mess = "Error in Get_Rotation_Axis subroutine. &
            Determinant is not +-1"
            return
        end if

        do i = 1 , 3
            A(i,i) = A(i,i) - (1//1)
        end do

        U = A
        call Rational_RowEchelonForm(U)
        call Rational_Rank(U,rnk)
        if ( rnk /= 2 ) then
            err_std = .true.
            err_std_mess = "Error in Get_Rotation_Axis subroutine. &
            Rank of matrix U is not two"
            return
        end if

        ! If there is a row with all zero entries,
        ! we put it in the last row

        nzeros = 0
        do i = 1 , 3
            do j = 1 , 3
                if ( U(i,j) == (0//1) ) then
                    nzeros(i) = nzeros(i) + 1
                else
                    exit
                end if
            end do
        end do

        do
            i = 1
            ordered = .true.
            do i = 1 , 2
                if (nzeros(i+1) < nzeros(i)) then
                    nzeros_aux  = nzeros(i)
                    row(:)      = U(i,:)
                    U(i,:)      = U(i+1,:)
                    U(i+1,:)    = row(:)
                    nzeros(i)   = nzeros(i+1)
                    nzeros(i+1) = nzeros_aux
                    ordered      = .false.
                end if
            end do
            if ( ordered ) exit
        end do

        ! Compute the axis direction

        if ( U(3,3) /= (0//1) ) then
            axis(3) = (0//1)
        else if ( U(2,3) /= (0//1) .and. U(2,2) == (0//1) ) then
            axis(3) = (0//1)
        else
            axis(3) = (1//1)  ! Free choice
        end if

        if ( U(2,2) /= (0//1) ) then
            axis(2) = - U(2,3) * axis(3) / U(2,2)
            if (U(1,1) == (0//1)) then
                axis(1) = (1//1) ! Free choice
            else
                axis(1) = - ( U(1,2) * axis(2) + U(1,3) * axis(3) ) / U(1,1)
            end if
        else ! axis(3) must be zero because row(2) cannot be zero
            if ( U(1,2) == (0//1) ) then
                axis(2) = (1//1) ! Free choice
                axis(1) = - ( U(1,2) * axis(2) + U(1,3) * axis(3) ) / U(1,1)
            else
                axis(1) = (1//1) ! Free choice
                axis(2) = - ( U(1,1) * axis(1) + U(1,3) * axis(3) ) / U(1,2)
            end if
        end if

        call Smallest_Integral_Vector(axis)

        ! Choose the eigenvector axis with axis(3) > 0. If axis(3) = 0, choose axis with
        ! axis(2) > 0. If axis(2) = 0, choose axis with axis(1) > 0

        do i = 3, 1, -1
            if (axis(i) /= (0//1)) then
                if (axis(i) < (0//1)) axis = -axis
                exit
            end if
        end do

    end subroutine Get_Rotation_Axis

    !!---- Subroutine Get_Rotation_Order(W,order)
    !!----     type(rational), dimension(3,3), intent(in)  :: W
    !!----     integer,                        intent(out) :: order
    !!----
    !!---- Returns the rotation order of a symmetry operation
    !!----
    !!---- Updated: September - 2018
    !!----

    subroutine Get_Rotation_Order(W,order)

        !---- Arguments ----!
        type(rational), dimension(3,3), intent(in)  :: W
        integer,                        intent(out) :: order

        !---- Local variables ----!
        type(rational) :: tr

        err_std = .false.
        tr = rational_trace(W)
        if (mod(tr%Numerator,tr%Denominator) == 0) then
            select case (tr%Numerator / tr%Denominator)
                case (3)
                    order = 1
                case (-1)
                    order = 2
                case (0)
                    order = 3
                case (1)
                    order = 4
                case (2)
                    order = 6
                case default
                    err_std = .True.
                    err_std_mess = "Error in Get_Rotation_Order. Rotation matrix is not an allowed crystallographic rotation"
            end select
        else
            err_std = .True.
            err_std_mess = "Error in Get_Rotation_Order. Rotation matrix is not an allowed crystallographic rotation"
        end if

    end subroutine Get_Rotation_Order

    !!---- Subroutine Get_Rotations(symOP,nSymOP,order,n,idd)
    !!----    type(Symm_Oper_Type), dimension(nSymOP), intent(in)  :: symOP
    !!----    integer,                                 intent(in)  :: nSymOP
    !!----    integer,                                 intent(in)  :: n
    !!----    integer,                                 intent(out) :: nso
    !!----    integer, dimension(nSymOP,2),            intent(out) :: idd
    !!----
    !!---- It returns the number of symmetry operations in array symOp
    !!---- with rotational part of order n. The corresponding index
    !!---- and value of the determinant are returned in idd.
    !!----
    !!---- Updated: September - 2018
    !!

    subroutine Get_Rotations(symOP,nSymOP,n,nso,idd)

        !---- Arguments ----!
        type(Symm_Oper_Type), dimension(nSymOP), intent(in)  :: symOP
        integer,                                 intent(in)  :: nSymOP
        integer,                                 intent(in)  :: n
        integer,                                 intent(out) :: nso
        integer, dimension(nSymOP,2),            intent(out) :: idd

        !---- Local variables ----!
        integer                           :: i,d,t
        type(rational)                    :: det,tr
        integer,           dimension(6)   :: traces
        character(len=20), dimension(6)   :: axisName

        err_std = .false.
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

        nso = 0
        do i = 1 , nSymOP
            det = rational_determinant(symOp(i)%Mat(1:3,1:3))
            if (mod(det%Numerator,det%Denominator) /= 0) then
                err_std      = .true.
                err_std_mess = "Error in Get_Rotations. Determinant is not an integer."
                return
            end if
            d = det%Numerator / det%Denominator
            tr  = rational_trace(symOp(i)%Mat(1:3,1:3))
            if (mod(tr%Numerator,tr%Denominator) /= 0) then
                err_std      = .true.
                err_std_mess = "Error in Get_Rotations. Trace is not an integer."
                return
            end if
            t = tr%Numerator / tr%Denominator
            if ( d == 1 .and. t == traces(n) .or. d == -1 .and. t == -traces(n)) then
                nso = nso + 1
                idd(nso,1) = i
                idd(nso,2) = d
            end if
        end do

        !if (nso == 0) then
        !    err_std      = .true.
        !    err_std_mess = "Error in Get_Rotations. Unable to find the "&
        !        //trim(axisName(n))//" axis"
        !end if

    end subroutine Get_Rotations

    !!---- Subroutine Get_S_Matrix
    !!----      type(rational), dimension(3,3), intent(in)  :: W
    !!----      type(rational), dimension(3,3), intent(out) :: S
    !!----
    !!---- Given a rotation matrix W, computes the matrix S defined as:
    !!----           S = W + W^2 + W^3 + ... + W^n
    !!---- where n is the order of the rotation axis. This matrix is
    !!---- used to find vectors perpendicular to the rotation axis. For
    !!---- a vector x perpendicular to the rotation axis, since Sx = 0.
    !!----
    !!---- Updated: September - 2018
    !!----
    !!

    subroutine Get_S_Matrix(W,S)

        !---- Arguments ----!
        type(rational), dimension(3,3), intent(in)  :: W
        type(rational), dimension(3,3), intent(out) :: S


        !---- Local variables ----!
        integer                        :: i,order
        type(rational), dimension(3,3) :: Waux

        call Get_Rotation_Order(W,order)

        S    = W
        Waux = W

        do i = 2 , order
            Waux = MatMul(Waux,W)
            S    = S + Waux
        end do

    end subroutine Get_S_Matrix

    !!---- Subroutine Get_Vectors_Perpendicular_To_Rotation_Axis(W,vPerp)
    !!----    type(rational), dimension(3,3), intent(in)  :: W
    !!----    type(rational), dimension(3,4), intent(out) :: vPerp
    !!----
    !!---- Given a rotation matrix W, computes the shortest
    !!---- four lattice vectors perpendicular to the rotation
    !!---- axis.
    !!----
    !!---- Updated: September - 2018
    !!----
    !!

    subroutine Get_Vectors_Perpendicular_To_Rotation_Axis(W,vPerp)

        !---- Arguments ----!
        type(rational), dimension(3,3), intent(in)  :: W
        type(rational), dimension(3,4), intent(out) :: vPerp

        !---- Local variables ----!
        integer        :: i
        integer        :: rnk,row,nzeros
        type(rational) :: det
        type(rational), dimension(3)   :: nullVector
        type(rational), dimension(3,3) :: A,S,U

        nullVector(:) = 0 // 1
        det = rational_determinant(W)

        if (det%Numerator == det%Denominator) then
            A = W
        else if (det%Numerator == -det%Denominator) then
            A = -W
        else
            err_std = .true.
            err_std_mess = "Error in Get_Vectors_Perpendicular_To_Rotation_Axis subroutine. &
            Determinant is not +-1"
            return
        end if

        call Get_S_Matrix(A,S)
        call Rational_RowEchelonForm(S)
        U = S

        ! Check that the rank of the U matrix is one

        call Rational_Rank(U,rnk)

        if ( rnk /= 1 ) then
            err_std = .true.
            err_std_mess = "Error in Get_Vectors_Perpendicular_To_Rotation_Axis subroutine. &
            Rank of matrix U is not one"
            return
        end if

        ! Find a row different from zeros

        do i = 1 , 3
            if (.not. equal_rational_vector(U(i,:),nullVector)) exit
        end do
        row = i

        ! Count the number of zeros of the first row of the echelon matrix

        nzeros = 0
        do i = 1 , 3
            if (U(row,i) == (0//1)) nzeros = nzeros + 1
        end do

        ! Build the four shortest vectors perpendicular to the rotation axis

        select case (nzeros)
            case (0)
                    vPerp(1,1) =  1//1
                    vPerp(2,1) =  0//1
                    vPerp(1,2) =  0//1
                    vPerp(2,2) =  1//1
                    vPerp(1,3) =  1//1
                    vPerp(2,3) =  1//1
                    vPerp(1,4) =  1//1
                    vPerp(2,4) = -1//1
                    do i = 1 , 4
                        vPerp(3,i) = -(vPerp(1,i) * U(row,1) + vPerp(2,i) * U(row,2)) / U(row,3)
                    end do
            case (1)
                do i = 1 , 3
                    if (U(row,i) == (0//1)) exit
                end do
                select case (i)
                    case (1)
                        vPerp(1,1) = 1//1
                        vPerp(2,1) = 0//1
                        vPerp(3,1) = 0//1
                        vPerp(1,2) = 0//1
                        vPerp(2,2) = 1//1
                        vPerp(3,2) = - U(row,2) / U(row,3)
                        vPerp(1,3) = 1//1
                        vPerp(2,3) = 1//1
                        vPerp(3,3) = - U(row,2) / U(row,3)
                        vPerp(1,4) = -1//1
                        vPerp(2,4) = 1//1
                        vPerp(3,4) = - U(row,2) / U(row,3)
                    case (2)
                        vPerp(1,1) = 0//1
                        vPerp(2,1) = 1//1
                        vPerp(3,1) = 0//1
                        vPerp(1,2) = 1//1
                        vPerp(2,2) = 0//1
                        vPerp(3,2) = - U(row,1) / U(row,3)
                        vPerp(1,3) = 1//1
                        vPerp(2,3) = 1//1
                        vPerp(3,3) = - U(row,1) / U(row,3)
                        vPerp(1,4) = 1//1
                        vPerp(2,4) = -1//1
                        vPerp(3,4) = - U(row,1) / U(row,3)
                    case (3)
                        vPerp(1,1) = 0//1
                        vPerp(2,1) = 0//1
                        vPerp(3,1) = 1//1
                        vPerp(1,2) = 1//1
                        vPerp(2,2) = -U(row,1) / U(row,2)
                        vPerp(3,2) = 0//1
                        vPerp(1,3) = 1//1
                        vPerp(2,3) = -U(row,1) / U(row,2)
                        vPerp(3,3) = 1//1
                        vPerp(1,4) = 1//1
                        vPerp(2,4) = -U(row,1) / U(row,2)
                        vPerp(3,4) = -1//1
                end select
            case (2)
                do i = 1 , 3
                    if (U(row,i) /= (0//1)) exit
                end do
                select case (i)
                    case (1)
                        vPerp(1,:) =  0//1
                        vPerp(2,1) =  1//1
                        vPerp(3,1) =  0//1
                        vPerp(2,2) =  0//1
                        vPerp(3,2) =  1//1
                        vPerp(2,3) =  1//1
                        vPerp(3,3) =  1//1
                        vPerp(2,4) =  1//1
                        vPerp(3,4) = -1//1
                    case (2)
                        vPerp(2,:) =  0//1
                        vPerp(1,1) =  1//1
                        vPerp(3,1) =  0//1
                        vPerp(1,2) =  0//1
                        vPerp(3,2) =  1//1
                        vPerp(1,3) =  1//1
                        vPerp(3,3) =  1//1
                        vPerp(1,4) =  1//1
                        vPerp(3,4) = -1//1
                    case (3)
                        vPerp(3,:) =  0//1
                        vPerp(1,1) =  1//1
                        vPerp(2,1) =  0//1
                        vPerp(1,2) =  0//1
                        vPerp(2,2) =  1//1
                        vPerp(1,3) =  1//1
                        vPerp(2,3) =  1//1
                        vPerp(1,4) =  1//1
                        vPerp(2,4) = -1//1
                end select
        end select

        do i = 1 , 4
            call Smallest_Integral_Vector(vPerp(:,i))
        end do

    end subroutine Get_Vectors_Perpendicular_To_Rotation_Axis

    !!---- Subroutine Identify_Crystallographic_Point_Group(G)
    !!----     type(spg_type), intent(inout) :: G
    !!----
    !!----  Determines the crystallographic point group of the group G.
    !!----
    !!----  Updated: January - 2019
    !!

    subroutine Identify_Crystallographic_Point_Group(G)

        !---- Arguments ----!
        class(spg_type), intent(in out) :: G

        !---- Local variables ----!
        integer                                            :: i,j,n,d,t
        integer                                            :: numops,nRepSymOp
        logical                                            :: selected
        type(rational)                                     :: det,tr
        integer,               dimension(6)                :: nRot
        type(Symm_Oper_Type),  dimension(:),   allocatable :: repSymOp ! representative operations
        integer,               dimension(:,:), allocatable :: idd

        ! Initialization
        nRot(:)  = 0  ! number of selected rotations of order 1,2,...,6
        err_std  = .false.
        G%pg     = ""
        numops = G%multip / (G%num_lat+G%num_alat+1)
        if (G%centred /= 1 .or. G%anticentred /= 1) numops = numops / 2
        allocate(repSymOp(numops))
        ! Get the rotations of the representative matrices
        nRepSymOp = 0
        do i = 1 , G%Multip
            !if (Equal_Rational_Matrix(G%op(i)%Mat(1:3,1:3),identity)) cycle
            det = rational_determinant(G%op(i)%Mat(1:3,1:3))
            if (mod(det%numerator,det%denominator) /= 0) then
                err_std      = .true.
                err_std_mess = "Error in Identify_Crystallographic_Point_Group. Determinant is not an integer."
                return
            end if
            d = det%numerator / det%denominator
            selected = .false.
            do j = 1 , nRepSymOp
                if (Equal_Rational_Matrix(G%op(i)%Mat(1:3,1:3),repSymOp(j)%Mat(1:3,1:3)) .or. &
                    Equal_Rational_Matrix(G%op(i)%Mat(1:3,1:3),-repSymOp(j)%Mat(1:3,1:3))) then
                    selected = .true.
                    exit
                end if
            end do
            if (.not. selected) then
                nRepSymOp = nRepSymOp + 1
                if (nRepSymOp > numops) then
                    err_std      = .True.
                    err_std_mess = "Error in Get_Crystallographic_Point_Group. &
                    nRepSymOp > numops"
                    return
                end if
                repSymOp(nRepSymOp) = G%op(i)
                tr  = rational_trace(G%op(i)%Mat(1:3,1:3))
                if (mod(tr%numerator,tr%denominator) /= 0) then
                    err_std      = .true.
                    err_std_mess = "Error in Get_Point_Group. Trace is not an integer."
                    return
                end if
                t = tr%numerator / tr%denominator
                select case (abs(t))
                    case (0)
                        nRot(3) = nRot(3) + 1
                    case(1)
                        if (d * t ==  1) then
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
        ! Get the point group
        allocate(idd(nRepSymOp,2))
        if (nRot(3) == 8) then ! Cubic
            if (nRepSymOp == 12) then
                if (G%Centred == 1 .and. G%anticentred == 1) then
                    G%PG = "23"
                else
                    G%PG = "m-3"
                end if
            else if (nRepSymOp == 24) then
                if (G%Centred == 1 .and. G%anticentred == 1) then
                    call Get_Rotations(repSymOp,nRepSymOp,4,n,idd)
                    if (n == 6 .and. idd(1,2) == 1) then
                        G%PG = "432"
                    else if (n == 6 .and. idd(1,2) == -1) then
                        G%PG = "-43m"
                    end if
                else
                    G%PG = "m-3m"
                end if
            end if
        else if (nRot(6) == 2) then ! Hexagonal
            if (nRepSymOp == 6) then
                if (G%Centred == 1 .and. G%anticentred == 1) then
                    call Get_Rotations(repSymOp,nRepSymOp,6,n,idd)
                    if (idd(1,2) == 1) then
                        G%PG = "6"
                    else
                        G%PG = "-6"
                    end if
                else
                    G%PG = "6/m"
                end if
            else if (nRepSymOp == 12) then
                if (G%Centred == 1 .and. G%anticentred == 1) then
                    call Get_Rotations(repSymOp,nRepSymOp,6,n,idd)
                    if (idd(1,2) == 1) then
                        call Get_Rotations(repSymOp,nRepSymOp,2,n,idd)
                        if (n == 7 .and. (idd(1,2) == -1 .or. idd(2,2) == -1)) then
                            G%PG = "6mm"
                        else
                            G%PG = "622"
                        end if
                    else if (idd(1,2) == -1) then
                        G%PG = "-6m2"
                    end if
                else
                    G%PG = "6/mmm"
                end if
            end if
        else if (nRot(3) == 2) then ! Trigonal
            if (nRepSymOp == 3) then
                if (G%Centred == 1 .and. G%anticentred == 1) then
                    G%PG = "3"
                else
                    G%PG = "-3"
                end if
            else if (nRepSymOp == 6) then
                if (G%Centred == 1 .and. G%anticentred == 1) then
                    call Get_Rotations(repSymOp,nRepSymOp,2,n,idd)
                    if (n == 3 .and. idd(1,2) == 1) then
                        G%PG = "32"
                    else if (n == 3 .and. idd(1,2) == -1) then
                        G%PG = "3m"
                    end if
                else
                    G%PG = "-3m"
                end if
            end if
        else if (nRot(4) == 2) then ! Tetragonal
            if (nRepSymOp == 4) then
                if (G%Centred == 1 .and. G%anticentred == 1) then
                    call Get_Rotations(repSymOp,nRepSymOp,4,n,idd)
                    if (n == 2 .and. idd(1,2) == 1) then
                        G%PG = "4"
                    else if (n == 2 .and. idd(1,2) == -1) then
                        G%PG = "-4"
                    end if
                else
                    G%PG = "4/m"
                end if
            else if (nRepSymOp == 8) then
                if (G%Centred == 1 .and. G%anticentred == 1) then
                    call Get_Rotations(repSymOp,nRepSymOp,4,n,idd)
                    if (n == 2 .and. idd(1,2) == 1) then
                        call Get_Rotations(repSymOp,nRepSymOp,2,n,idd)
                        if (n == 5 .and. (idd(1,2) == -1 .or. idd(2,2) == -1)) then
                            G%PG = "4mm"
                        else
                            G%PG = "422"
                        end if
                    else if (n == 2 .and. idd(1,2) == -1) then
                        G%PG = "-4m2"
                    end if
                else
                    G%PG = "4/mmm"
                end if
            end if
        else if (nRot(2) == 3) then ! Orthorhombic
            if (G%Centred == 1 .and. G%anticentred == 1) then
                call Get_Rotations(repSymOp,nRepSymOp,2,n,idd)
                if (n == 3 .and. idd(1,2) == -1 .or. idd(2,2) == -1) then
                    G%PG = "mm2"
                else
                    G%PG = "222"
                end if
            else
                G%PG = "mmm"
            end if
        else if (nRot(2) == 1) then ! Monoclinic
            if (G%Centred == 1 .and. G%anticentred == 1) then
                call Get_Rotations(repSymOp,nRepSymOp,2,n,idd)
                if (idd(1,2) == 1) then
                    G%PG = "2"
                else
                    G%PG = "m"
                end if
            else
                G%PG = "2/m"
            end if
        else
            if (G%Centred == 1 .and. G%anticentred == 1) then ! Triclinic
                G%PG = "1"
            else
                G%PG = "-1"
            end if
        end if

        if (G%PG == "") then
            err_std = .true.
            err_std_mess = "Error in Identify_Crystallographic_Point_Group_Basic. Unable to identify point group"
        end if

    end subroutine Identify_Crystallographic_Point_Group

    !!---- Subroutine Identify_Crystallographic_Space_Group(G)
    !!----     class(spg_type), intent(inout)  :: G
    !!----
    !!---- For a given crystallographic group G in an arbitrary
    !!---- setting, identifies the space group and computes the
    !!---- transformation matrix to the standard setting.
    !!----
    !!---- It follows Acta Cryst. A55 383-395 (1999). P,M,A and
    !!---- C matrices in the paper correspond to Pinv, Minv,
    !!---- Ainv and Cinv in this subroutine.
    !!----
    !!---- Updated January - 2019
    !!

    subroutine Identify_Crystallographic_Space_Group(G,output)

        !---- Arguments ----!
        class(spg_type),    intent(in out) :: G
        logical, optional, intent(in)     :: output

        !---- Local variables ---!
        integer                          :: n
        character                        :: lattyp
        type(rational), dimension(3,3)   :: P,Mp,Mc,M
        type(rational), dimension(3,3,6) :: A
        logical :: pout

        pout=.false.
        if(present(output)) pout=output

        call Get_P_Matrix(G,P,nospin=.true.,output=pout)
        !call Get_P_Matrix(G,P,nospin=.true.)
        if (err_std) return

        call Get_Mp_Matrix(G,P,Mp,output=pout)
        !call Get_Mp_Matrix(G,P,Mp)
        if (err_std) return

        call Get_Mc_Matrix(G%laue,Mp,Mc,output=pout)
        !call Get_Mc_Matrix(G%laue,Mp,Mc)
        if (err_std) return

        M = matmul(Mp,Mc)
        call Get_Lattice_Type_from_M(M,lattyp)
        if (err_std) return
        call Get_A_Matrices_Crys(G%laue,A,n)
        call Match_Crystallographic_Space_Group(G,P,M,A(:,:,1:n),n,output=pout)
        !call Match_Crystallographic_Space_Group(G,P,M,A(:,:,1:n),n)

    end subroutine Identify_Crystallographic_Space_Group

    !!---- Subroutine Identify_Group(G)
    !!----     class(spg_type), intent(inout)  :: G
    !!----
    !!---- Initialize the identification of the group by calling
    !!---- the appropiate subroutine according to the nature  of
    !!---- the group -crystallographic, magnetic, superspace-.
    !!----
    !!---- Updated January - 2019
    !!

    subroutine Identify_Group(G,output)

        !---- Arguments ----!
        class(spg_type),    intent(in out) :: G
        logical, optional, intent(in)     :: output

        logical :: pout

        pout=.false.
        if(present(output)) pout=output

        err_std = .false.

        if (G%d < 4) then
            err_std = .true.
            err_std_mess = "Error in Identify_Group. Group dimension < 4."
            return

        else if (G%d == 4) then ! Shubnikov groups

            call Identify_Shubnikov_Group(G,pout)

        else ! Superspace groups
            err_std      = .true.
            err_std_mess = "Error in Identify_Group. Superspace groups not implemented yet"
            return
        end if

    end subroutine Identify_Group

    !!---- Subroutine Identify_Laue_Class(G)
    !!----     type(spg_type), intent(inout) :: G
    !!----
    !!---- Sets the Laue class from the crystallographic point group
    !!----
    !!---- Updated: January - 2019
    !!

    subroutine Identify_Laue_Class(G)

        !---- Arguments ----!
        class(spg_type), intent(in out) :: G

        select case (trim(G%pg))

            case ("1","-1")
                G%laue = "-1"
            case ("2","m","2/m")
                G%laue = "2/m"
            case ("222","mm2","mmm")
                G%laue = "mmm"
            case ("4","-4","4/m")
                G%laue = "4/m"
            case ("422","4mm","-4m2","4/mmm")
                G%laue = "4/mmm"
            case ("3","-3")
                G%laue = "-3"
            case ("32","3m","-3m")
                G%laue = "-3m"
            case ("6","-6","6/m")
                G%laue = "6/m"
            case ("622","6mm","-6m2","6/mmm")
                G%laue = "6/mmm"
            case ("23","m-3")
                G%laue = "m-3"
            case ("432","-43m","m-3m")
                G%laue = "m-3m"
            case default
                err_std = .true.
                err_std_mess = "Error in Identify_Laue_Class. Inconsistent crystallographic point group."
        end select

    end subroutine Identify_Laue_Class

    !!---- Subroutine Identify_Magnetic_Space_Group(G)
    !!----      type(spg_type), intent(inout) :: G
    !!
    !!---- Identifies the Shubnikov group of group G. Before
    !!---- using this subroutine, the subroutine Identify_
    !!---- _Crystallographic_Space_Group must be called
    !!----
    !!---- Updated January - 2019
    !!----

    subroutine Identify_Shubnikov_Group(G,output)

        !---- Arguments ----!
        class(spg_type),    intent(in out) :: G
        logical, optional, intent(in)     :: output

        !---- Local variables ---!
        character                        :: lattyp
        type(rational), dimension(3,3)   :: P,Mp,Mc,M
        logical :: pout
        !debug variables
        !character(len=20), dimension(3,3):: matrix
        !integer :: i,j,k

        pout=.false.
        if(present(output)) pout=output
        if(.not. database_allocated) then
          call Read_Magnetic_Data()
          if(err_magg) then
            write(*,"(a)") trim(err_magg_mess)
            return
          end if
        end if

        call Identify_Crystallographic_Point_Group(G)
        if (err_std) return

        call Identify_Laue_Class(G)
        if (err_std) return

        call Get_P_Matrix(G,P,output=pout)
        !call Get_P_Matrix(G,P)
        if (err_std) return

        call Get_Mp_Matrix(G,P,Mp,output=pout)
        !call Get_Mp_Matrix(G,P,Mp)
        if (err_std) return

        call Get_Mc_Matrix(G%laue,Mp,Mc,output=pout)
        !call Get_Mc_Matrix(G%laue,Mp,Mc)
        if (err_std) return

        M = matmul(Mp,Mc)
        call Get_Lattice_Type_from_M(M,lattyp)
        if (err_std) return

        call Match_Shubnikov_Group(G,P,M,output=pout)
        !call Match_Shubnikov_Group(G,P,M)

    end subroutine Identify_Shubnikov_Group

    !!---- Subroutine Match_Crystallographic_Space_Group(G,P,M,A,n,output)
    !!----      class(spg_type),                   intent(inout) :: G
    !!----      type(rational), dimension(3,3),   intent(in)    :: P
    !!----      type(rational), dimension(3,3),   intent(in)    :: M
    !!----      type(rational), dimension(3,3,n), intent(in)    :: A
    !!----      integer,                          intent(in)    :: n
    !!----      logical,        optional,         intent(in)    :: output
    !!----
    !!---- Tries to match the space group G against one of the 230
    !!---- standard crystallographic groups. It returns the space
    !!---- group number G%numspg, the space group symbol G%spg_symb
    !!---- and the transformation matrix to the standard G%mat2std
    !!----
    !!---- Updated: January - 2019
    !!

    subroutine Match_Crystallographic_Space_Group(G,P,M,A,n,output)

        !---- Arguments ----!
        class(spg_type),                  intent(inout) :: G        ! space group in the original setting
        type(rational), dimension(3,3),   intent(in)    :: P        ! P matrix   -see Get_P_Matrix-
        type(rational), dimension(3,3),   intent(in)    :: M        ! M matrix   -see Get_M_Matrix-
        type(rational), dimension(3,3,n), intent(in)    :: A        ! A matrices -see Get_A_Matrices_Crys-
        integer,                          intent(in)    :: n        ! number of A matrices (six as maximum)
        logical,        optional,         intent(in)    :: output

        !---- Local variables ----!
        integer                                         :: s,i,j,ng,ng_
        integer                                         :: firstSpaceGroup,lastSpaceGroup,numSpg
        character(len=12)                               :: sgString
        character(len=256)                              :: symb
        logical                                         :: shift
        type(spg_type)                                  :: G_target
        type(space_group_type)                          :: G_std
        type(rational),       dimension(3)              :: origShift
        type(spg_type),       dimension(n)              :: G_
        type(Symm_Oper_Type), dimension(3)              :: gen_std,gen_x
        type(rational),       dimension(3,3)            :: MA,P_target
        type(rational),       dimension(4,4)            :: C,Cinv,W
        type(rational),       dimension(4,4,n)          :: C_
        type(Symm_Oper_Type), dimension(:), allocatable :: op
        logical :: pout

        pout=.false.
        if(present(output)) pout=output

        select case (trim(G%laue))

            case ("-1")

                firstSpaceGroup = 1
                lastSpaceGroup  = 2
                if (pout) write(*,'(8x,a)') " => Matching representation against standard crystallographic triclinic space groups..."

            case ("2/m")

                firstSpaceGroup = 3
                lastSpaceGroup  = 15
                if (pout) write(*,'(8x,a)') " => Matching representation against standard crystallographic monoclinic space groups..."

            case ("mmm")
                firstSpaceGroup = 16
                lastSpaceGroup  = 74
                if (pout) write(*,'(8x,a)') " => Matching representation against standard crystallographic orthorhombic space groups..."

            case ("4/m","4/mmm")
                firstSpaceGroup = 75
                lastSpaceGroup  = 142
                if (pout) write(*,'(8x,a)') " => Matching representation against standard crystallographic tetragonal space groups..."

            case ("-3","-3 R","-3m","-3m R","-3m1","-31m")
                firstSpaceGroup = 143
                lastSpaceGroup  = 167
                if (pout) write(*,'(8x,a)') " => Matching representation against standard crystallographic trigonal space groups..."

            case ("6/m","6/mmm")
                firstSpaceGroup = 168
                lastSpaceGroup  = 194
                if (pout) write(*,'(8x,a)') " => Matching representation against standard crystallographic hexagonal space groups..."

            case ("m3","m-3","m3m","m-3m")
                firstSpaceGroup = 195
                lastSpaceGroup  = 230
                if (pout) write(*,'(8x,a)') " => Matching representation against standard crystallographic cubic space groups..."

            case default
                firstSpaceGroup = 0
                lastSpaceGroup  = -1

        end select

        do i = 1 , 3
            allocate(gen_std(i)%Mat(4,4),gen_x(i)%Mat(4,4))
        end do

        do s = 1 , n  ! Loop over C matrices
            G_(s)         = G
            G_(s)%numops  = G%multip / (G%num_lat + G%num_alat + 1)
            if (G_(s)%centred == 1) G_(s)%centred = G_(s)%anticentred
            if (G_(s)%centred /= 1) G_(s)%numops  = G_(s)%numops / 2
            MA            = matmul(M,A(:,:,s))
            C_(1:3,1:3,s) = matmul(P,MA)
            C_(1:3,4,s)   = [ 0//1,0//1,0//1 ]
            C_(4,1:4,s)   = [ 0//1,0//1,0//1,1//1 ]
            call Rational_Inv_Matrix(C_(:,:,s),Cinv)
            ! Compute symmetry operations in the new basis
            do j = 1 , G%Multip
                W(1:3,1:3) = G%op(j)%Mat(1:3,1:3)
                W(1:3,4)   = G%op(j)%Mat(1:3,G%d)
                W(4,1:4)   = [ 0//1,0//1,0//1,1//1 ]
                W          = matmul(W,C_(:,:,s))
                W          = matmul(Cinv,W)
                G_(s)%op(j)%Mat(1:3,1:3) = W(1:3,1:3)
                G_(s)%op(j)%Mat(1:3,G%d) = W(1:3,4)
            end do
            ! Compute vectors in the new basis
            do j = 1 , G%num_lat
                G_(s)%lat_tr(1:3,j) = matmul(Cinv(1:3,1:3),G%lat_tr(1:3,j))
                call PBC(G_(s)%lat_tr(1:3,j))
            end do
            do j = 1 , G%num_alat
                G_(s)%alat_tr(1:3,j) = matmul(Cinv(1:3,1:3),G%alat_tr(1:3,j))
                call PBC(G_(s)%alat_tr(1:3,j))
            end do
            ! Get Lattice Type
            call Get_Lattice_Type_from_M(MA,G_(s)%spg_lat)
        end do

        do numSpg = firstSpaceGroup , lastSpaceGroup
            call Get_HM_Standard(numSpg,sgString)
            call Set_Spacegroup(sgString,G_std)
            do s = 1 , n
                if (G_std%NumOps  == G_(s)%NumOps .and. &
                    G_std%Spg_Lat == G_(s)%Spg_Lat .and. &
                   (G_std%Centred == 1 .and. G_(s)%Centred == 1 .or. &
                    G_std%Centred /= 1 .and. G_(s)%Centred /= 1)) then
                    if (allocated(op)) deallocate(op)
                    allocate(op(G_std%Multip))
                    call set_identity_matrix(4)
                    do i = 1 , G_std%Multip
                        allocate(op(i)%Mat(4,4))
                        op(i)%Mat          = identity_matrix
                        op(i)%Mat(1:3,1:3) = G_std%symop(i)%Rot
                        op(i)%Mat(1:3,4)   = G_std%symop(i)%Tr
                    end do
                    ! Get generators from the standard space group
                    call Get_Generators(G_std%laue,op,G_std%Multip,gen_std,ng)
                    if (err_std) return
                    ! Try to get these generators from SGaux(s)
                    ng_ = 0
                    do i = 1 , ng
                        do j = 1 , G_(s)%Multip
                            if (Equal_Rational_Matrix(gen_std(i)%Mat(1:3,1:3),G_(s)%op(j)%Mat(1:3,1:3))) then
                                ng_ = ng_ + 1
                                gen_x(ng_) = G_(s)%op(j)
                                !gen_x(1:3,1:3,ng_) = gen_std(1:3,1:3,i)
                                !gen_x(1:3,4,ng_)   = G_(s)%op(j)%Mat(1:3,G%d)
                                !gen_x(4,:,ng_)     = [ 0//1, 0//1, 0//1, 1//1 ]
                                exit
                            end if
                        end do
                    end do
                    if (ng /= ng_) cycle
                    if (pout) write(*,'(12x,3a,i2)',advance = 'no') 'Trying to match space group ', G_std%Spg_Symb, 'setting ', s
                    ! Build a spg_type object from space_group_type object
                    G_target%num_alat = 0
                    G_target%num_lat  = G_std%NumLat - 1 ! In spg_type (000) is not included in lattice translations
                    G_target%Multip   = G_std%Multip
                    if (allocated(G_target%lat_tr)) deallocate(G_target%lat_tr)
                    if (G_target%num_lat > 0) then
                        allocate(G_target%lat_tr(3,G_target%num_lat))
                        G_target%lat_tr(:,:) = G_std%latt_trans(:,2:G_std%NumLat)
                    end if
                    if (allocated(G_target%op)) deallocate(G_target%op)
                    allocate(G_target%op(G_std%Multip))
                    do i = 1 , G_std%Multip
                        allocate (G_target%op(i)%Mat(3,3))
                        G_target%op(i)%Mat = G_std%symop(i)%Rot
                    end do
                    call Get_P_Matrix(G_target,P_target(1:3,1:3))
                    if (err_std) return
                    ! Try to match the standard by an origin shift
                    call Get_Origin_Shift(gen_x(1:ng),gen_std(1:ng),ng,P_target,origShift,shift)
                    if (shift) then
                        if (pout) write(*,'(8x,a)') 'Done!'
                        C(1:3,1:3) = C_(1:3,1:3,s)
                        C(4,:)     = C_(4,:,s)
                        origShift  = matmul(C(1:3,1:3),origShift)
                        C(4,1:3)   = origShift(1:3)
                        call Get_Symb_Op_from_Mat(transpose(C),symb,"abc")
                        if (pout) then
                            write(*,'(12x,a)',advance='no') "Original setting --> Standard crystallographic setting: "
                            write(*,'(a)') trim(symb)
                        end if
                        G%numspg   = G_std%numspg
                        G%spg_symb = G_std%spg_symb
                        G%mat2std  = trim(symb)
                        return
                    else
                        if (pout) write(*,'(tr8,a)') 'Failed'
                        cycle
                    end if
                end if
            end do
        end do

        err_std = .true.
        err_std_mess = "Unable to indentify the space group"

    end subroutine Match_Crystallographic_Space_Group

    !!---- Subroutine Match_Shubnikov_Group(G,P,M,output)
    !!----      class(spg_type),                   intent(inout) :: G
    !!----      type(rational), dimension(3,3),   intent(in)    :: P        ! P matrix   -see Get_P_Matrix-
    !!----      type(rational), dimension(3,3),   intent(in)    :: M        ! M matrix   -see Get_M_Matrix-
    !!----      logical,        optional,         intent(in)    :: output
    !!----
    !!---- Tries to match the space group G against Shubnikov groups. If the
    !!---- matching is successfull, it writes the group number G%numshu, the
    !!---- group symbol G%shu_symb and the transformation matrix to the
    !!---- standard G%mat2std_shu
    !!----
    !!---- Updated: February - 2019
    !!

    subroutine Match_Shubnikov_Group(G,P,M,output)

        !---- Arguments ----!
        class(spg_type),                  intent(in out):: G
        type(rational), dimension(3,3),   intent(in)    :: P        ! P matrix   -see Get_P_Matrix-
        type(rational), dimension(3,3),   intent(in)    :: M        ! M matrix   -see Get_M_Matrix-
        logical,        optional,         intent(in)    :: output

        !---- Local variables ----!
        integer                                                    :: i,j,k,n,ngen,ngen_
        integer                                                    :: iniG,endG
        integer                                                    :: nPointOper,nA
        logical                                                    :: hexagonal,shift
        character(len=256)                                         :: symb
        type(symm_oper_type)                                       :: Op
        character(len=1),          dimension(2)                    :: magLat
        type(symm_oper_type),      dimension(3)                    :: Gt
        type(rational),            dimension(3)                    :: origShift
        type(rational),            dimension(4,4)                  :: Caux
        type(symm_oper_type),      dimension(3,6)                  :: Gx
        type(rational),            dimension(3,3,6)                :: A
        logical,                   dimension(:),       allocatable :: doTest
        character(len=40),         dimension(:),       allocatable :: gener
        type(spg_type),            dimension(:),       allocatable :: G_aux
        integer,                   dimension(:,:),     allocatable :: idx,pointerToOper
        type(rational),            dimension(:,:,:),   allocatable :: C,Cinv,pointOper,Paux
        logical :: pout
        !character(len=40),         dimension(4,4)                  :: matrix

        pout=.false.
        if(present(output)) pout=output

        ! Initialize variables
        hexagonal = .false.

        ! Check dimension of the group is four
        if (G%d /= 4) then
            err_std = .true.
            err_std_mess = "Error in Match_Shubnikov_Group. Group dimension different from 4."
            return
        end if

        ! Load table for standard magnetic groups
        if(.not. database_allocated) then
          call Read_Magnetic_Data()
          if(err_magg) then
            write(*,"(a)") trim(err_magg_mess)
            return
          end if
        end if
        ! Set the range of possible magnetic groups from the laue class
        select case (trim(G%laue))
            case ("-1")
                iniG  =    1
                endG  =    7
                if (pout) write(*,'(8x,a)') " => Matching representation against standard triclinic magnetic groups..."
            case ("2/m")
                iniG =    8
                endG  =   98
                if (pout) write(*,'(8x,a)') " => Matching representation against standard monoclinic magnetic groups..."
            case ("mmm")
                iniG =   99
                endG  =  660
                if (pout) write(*,'(8x,a)') " => Matching representation against standard orthorhombic magnetic groups..."
            case ("4/m","4/mmm")
                iniG =  661
                endG  = 1230
                if (pout) write(*,'(8x,a)') " => Matching representation against standard tetragonal magnetic groups..."
            case ("-3","-3 R","-3m","-3m R","-3m1","-31m")
                iniG = 1231
                endG  = 1338
                hexagonal = .true.
                if (pout) write(*,'(8x,a)') " => Matching representation against standard trigonal magnetic groups..."
            case ("6/m","6/mmm")
                iniG = 1339
                endG  = 1502
                hexagonal = .true.
                if (pout) write(*,'(8x,a)') " => Matching representation against standard hexagonal magnetic groups..."
            case ("m3","m-3","m3m","m-3m")
                iniG = 1503
                endG  = 1651
                if (pout) write(*,'(8x,a)') " => Matching representation against standard cubic magnetic groups..."
            case default
                err_std = .true.
                err_std_mess = "Error in Match_Shubnikov_Group. Unknown Laue class"
                return
        end select

        ! Set pointers
        if (hexagonal) then
            ! Hexagonal systems
            allocate(pointOper(3,3,24))
            pointOper  = point_op_hex_matrix
            nPointOper = 24
        else
            ! Non hexagonal systems
            allocate(pointOper(3,3,48))
            pointOper  = point_op_matrix
            nPointOper = 48
        end if

        ! Get matrices needed for test different settings
        call Get_A_Matrices_Shub(G%laue,A,nA)

        ! Allocate memory
        allocate(G_aux(nA))
        allocate(doTest(nA))
        allocate(idx(3,nA))
        allocate(pointerToOper(2*G%multip,nA))
        allocate(C(4,4,nA),Cinv(4,4,nA),Paux(3,3,nA))

        ! Initialize arrays
        G_aux(:)           = G
        doTest(:)          = .true.
        idx(:,:)           = 0
        pointerToOper(:,:) = 0

        ! Get transformation matrices for every setting we will test
        Caux(1:3,1:3) = matmul(P,M)
        do n = 1 , nA
            C(1:3,1:3,n) = matmul(Caux(1:3,1:3),A(1:3,1:3,n))
            C(1:3,4,n)   = [ 0,0,0 ]
            C(4,1:4,n)   = [ 0,0,0,1 ]
            call Rational_Inv_Matrix(C(:,:,n),Cinv(:,:,n))
         end do

        ! Put G in the different settings
        allocate(gener(G%multip))
        !call Get_Operators_From_String(G%generators_list,n,ngen,gener)
        !allocate(Op(ngen,nA))
        !call Allocate_Operator(G%d,Op_)
        !do i = 1 , nGen
        !    call Get_Mat_From_Symb_Op(gener(i),Op_%Mat,Op_%time_inv)
        !    do n = 1 , nA
        !        Op(i,n) = Op_
        !        Op(i,n)%Mat = matmul(Op(i,n)%Mat,C(:,:,n))
        !        Op(i,n)%Mat = matmul(Cinv(:,:,n),Op(i,n)%Mat)
        !    end do
        !end do

        ! Build groups from generators
        call Allocate_Operator(G%d,Op)
        do n = 1 , nA
            !write(*,'("Setting ",i1)') n
            do i = 1 , G%multip
                Op%Mat      = matmul(G%Op(i)%Mat,C(:,:,n))
                Op%Mat      = matmul(Cinv(:,:,n),Op%Mat)
                Op%time_inv = G%Op(i)%time_inv
                call Get_Symb_Op_from_Mat(Op%Mat,gener(i),"xyz",Op%time_inv)
                !write(*,"(a)") trim(gener(i))
            end do
            call Group_Constructor(gener,G_aux(n),"xyz")
            G_aux(n)%laue = G%laue
            G_aux(n)%pg   = G%pg
        end do
        ! Put G in every setting we will test ---> [G_aux] = [Cinv][G_aux][C]
        !call Rational_Identity_Matrix(3,identity)
        !do n = 1 , nA
        !    G_aux(n)          = G
        !    G_aux(n)%num_lat  = 0
        !    G_aux(n)%num_alat = 0
        !    do i = 1 , G_aux(n)%multip
        !        G_aux(n)%Op(i)%Mat = matmul(G_aux(n)%Op(i)%Mat,C(:,:,n))
        !        G_aux(n)%Op(i)%Mat = matmul(Cinv(:,:,n),G_aux(n)%Op(i)%Mat)
        !        ! if the matrix is not integral, the setting is discarded
        !        if (.not. IsInteger(G_aux(n)%Op(i)%Mat(1:3,1:3))) then
        !            doTest(n) = .false.
        !            exit
        !        end if
        !        call reduced_translation(G_aux(n)%Op(i)%Mat)
        !        if (Equal_Rational_Matrix(G_aux(n)%Op(i)%Mat(1:3,1:3),identity) .and. &
        !            .not. IsInteger(G_aux(n)%Op(i)%Mat(1:3,4))) then
        !            if (G_aux(n)%Op(i)%time_inv == 1) then
        !                newt = .true.
        !                do j = 1 , G_aux(n)%num_lat
        !                    if (Equal_Rational_Vector(G_aux(n)%Op(i)%Mat(1:3,4),&
        !                                              G_aux(n)%Lat_tr(1:3,j))) then
        !                        newt = .false.
        !                        exit
        !                    end if
        !                end do
        !                if (newt) then
        !                    G_aux(n)%num_lat = G_aux(n)%num_lat + 1
        !                    G_aux(n)%Lat_tr(1:3,G_aux(n)%num_lat) = G_aux(n)%Op(i)%Mat(1:3,4)
        !                end if
        !            else
        !                newt = .true.
        !                do j = 1 , G_aux(n)%num_alat
        !                    if (Equal_Rational_Vector(G_aux(n)%Op(i)%Mat(1:3,4),&
        !                                              G_aux(n)%aLat_tr(1:3,j))) then
        !                        newt = .false.
        !                        exit
        !                    end if
        !                end do
        !                if (newt) then
        !                    G_aux(n)%num_alat = G_aux(n)%num_alat + 1
        !                    G_aux(n)%aLat_tr(1:3,G_aux(n)%num_alat) = G_aux(n)%Op(i)%Mat(1:3,4)
        !                end if
        !            end if
        !        end if
        !    end do
        !    ! If there are integer anti-translations, the setting is discarded
        !    do i = 1 , G_aux(n)%num_alat
        !        if (.not. isNullVector(G_aux(n)%aLat_Tr(1:3,i)) .and. isInteger(G_aux(n)%aLat_Tr(1:3,i))) then
        !            doTest(n) = .false.
        !            exit
        !        end if
        !    end do
        !    !if (doTest(n)) then
        !    !    write(*,'(a,i2)')  "Setting: ",n
        !    !    write(*,'(a)')     "    Translations: "
        !    !    do i = 1 , G_aux(n)%num_lat
        !    !        write(*,'(10x,6i2)') G_aux(n)%lat_tr(:,i)
        !    !    end do
        !    !    write(*,'(a)')     "    Anti-Translations: "
        !    !    do i = 1 , G_aux(n)%num_alat
        !    !        write(*,'(10x,6i2)') G_aux(n)%alat_tr(:,i)
        !    !    end do
        !    !end if
        !end do!;stop

        ! Map each symmetry operation in pointerToOper for every possible setting
        do n = 1 , nA
            !if (doTest(n)) then
                do i = 1 , G_aux(n)%multip
                    j = 1
                    do
                        if (Equal_Rational_Matrix(G_aux(n)%Op(i)%Mat(1:3,1:3),pointOper(:,:,j))) then
                            pointerToOper(i,n) = j
                            exit
                        end if
                        j = j + 1
                        if (j > nPointOper) then
                            err_std = .true.
                            err_std_mess = "Error in Magnetic_Space_Group_Matching. &
                            A symmetry operation of the group cannot be matched against tabulated symmetry operations"
                            return
                        end if
                    end do
                end do
            !end if
        end do

        ! Get the magnetic lattice type for every setting
        do n = 1 , nA
            if (doTest(n)) then
                call Get_Magnetic_Lattice_Type(G_aux(n))!%num_Lat,G_aux(n)%lat_tr,G_aux(n)%num_aLat,&
                                               !G_aux(n)%alat_tr,G_aux(n)%shu_lat)
                ! Correct symbol for triclinic systems with anti-translations
                if (trim(G_aux(n)%laue) == "-1" .and. G_aux(n)%shu_lat(2) /= " ") G_aux(n)%shu_lat(2) = "S"
                !write(*,*) n,G_aux(n)%shu_lat,G_aux(n)%num_alat,G_aux(n)%num_lat
            end if
        end do!;stop

        ! Get generators for every setting
        do n = 1 , nA
            if (doTest(n)) then
                call Get_Generators(G_aux(n)%laue,G_aux(n)%op,G_aux(n)%Multip,Gx(:,n),ngen)
                ! Map generators to G_aux(n)%Op
                do i = 1 , ngen
                    idx(i,n) = 0
                    do j = 1 , G_aux(n)%multip
                        if (G_aux(n)%Op(j) == Gx(i,n)) idx(i,n) = j
                    end do
                    if (idx(i,n) == 0) then
                        err_std = .true.
                        err_std_mess = "Error in Match_Shubnikov_Group. Generator cannot be mapped."
                        return
                    end if
                end do
                !write(*,'("Generators for setting ",i2)') n
                !do i = 1 , ngen
                !    write(*,'(5x,"Generator ",i1)') i
                !    do j = 1 , 4
                !        write(*,'(5x,4i3,"/",i1)') Gx(i,n)%Mat(j,:)%numerator,Gx(i,n)%Mat(j,4)%denominator
                !    end do
                !    write(*,'(5x,i3)') Gx(i,n)%time_inv
                !end do
                !write(*,'(5x,a)')     "Translations: "
                !do i = 1 , G_aux(n)%num_lat
                !    write(*,'(10x,6i2)') G_aux(n)%lat_tr(:,i)
                !end do
                !write(*,'(5x,a)')     "Anti-Translations: "
                !do i = 1 , G_aux(n)%num_alat
                !    write(*,'(10x,6i2)') G_aux(n)%alat_tr(:,i)
                !end do
            end if
        end do

        ! Compute the P matrix for every setting
        do n = 1 , nA
            if (doTest(n)) then
                call Get_P_Matrix(G_aux(n),Paux(:,:,n))
                if (err_std) doTest(n) = .false.
            end if
        end do

         ! Try to match G against one of the Shubnikov groups

        do i = iniG , endG
            if (magtype(i) == G%mag_type .and. ops_count(i) == (G%multip/(G%num_lat+1))) then
                ! Extract magnetic lattice type from symbol
                magLat(1) = spacegroup_label_bns(i)(1:1)
                magLat(2) = " "
                if (spacegroup_label_bns(i)(2:2) == "_") magLat(2) = spacegroup_label_bns(i)(3:3)
                do n = 1 , nA
                    if (doTest(n)) then
                        if (G_aux(n)%shu_lat(1) .ne. magLat(1) .or. &
                            G_aux(n)%shu_lat(2) .ne. magLat(2)) cycle
                        ! Try to find the same set of generators in the standard magnetic group
                        ngen_ = 0
                        do j = 1 , ngen
                            do k = 1 , ops_count(i)
                                if (ops_bns_point_op(k,i) == pointerToOper(idx(j,n),n) .and. &
                                    ops_bns_timeinv(k,i) == Gx(j,n)%time_inv) then
                                    Gt(j) = Gx(j,n)
                                    Gt(j)%Mat(1:3,4) = (ops_bns_trans(1:3,k,i) // ops_bns_trans_denom(k,i))
                                    ngen_ = ngen_ + 1
                                    exit
                                end if
                            end do
                        end do
                        if (ngen_ == ngen) then
                            if (pout) write(*,'(12x,3a,i1)',advance = 'no') "Trying to match magnetic space group ",&
                                trim(spacegroup_label_bns(i)),", setting ", n
                            call Get_Origin_Shift(Gx(1:nGen,n),Gt(1:nGen),nGen,Paux(1:3,1:3,n),origShift,shift)
                            if (shift) then
                                if (pout) write(*,'(8x,a)') "Done!"
                                origShift(1:3) = matmul(C(1:3,1:3,n),origShift)
                                C(4,1:3,n) = origShift(1:3)
                                call Get_Symb_Op_from_Mat(transpose(C(:,:,n)),symb,"abc")
                                if (pout) then
                                    write(*,'(12x,a)',advance='no') "Original setting --> Standard magnetic setting: "
                                    write(*,'(a)') trim(symb)
                                end if
                                G%shu_symb   = spacegroup_label_bns(i)
                                G%numshu      = i
                                G%mat2std_shu = trim(symb)
                                return
                            else
                                if (pout) write(*,'(8x,a)') "Failed"
                            end if
                        end if
                    end if
                end do
            end if
        end do

        err_std = .true.
        err_std_mess = "Unable to indentify the magnetic group"

    end subroutine Match_Shubnikov_Group

    !!---- Subroutine PBC(vector)
    !!----     type(rational), dimension(:), intent(inout) :: vector
    !!----
    !!---- Apply periodic boundary conditions to get the equivalent vector
    !!---- with components in the range [0-1)
    !!----
    !!---- Updated: October - 2018
    !!

    subroutine PBC(vector)

        !---- Arguments ----!
        type(rational), dimension(:), intent(inout) :: vector

        !---- Local variables ----!
        integer :: i

        do i = 1 , size(vector)
            if (vector(i) > (1//1)) then
                vector(i) = vector(i) - ((vector(i)%Numerator/vector(i)%Denominator)//1_ik)
            else if (vector(i) < (0//1)) then
                vector(i) = vector(i) - ((vector(i)%Numerator/vector(i)%Denominator)//1_ik) + (1//1)
            end if
        end do

    end subroutine PBC

    !!---- Subroutine Set_Right_Handedness
    !!----      type(rational), dimension(3,3), intent(inout) :: A
    !!
    !!---- Updated: September - 2018
    !!

    subroutine Set_Right_Handedness(A)

        !---- Arguments ----!
        type(rational), dimension(3,3), intent(inout) :: A

        !---- Local variables ----!
        type(rational)               :: det
        type(rational), dimension(3) :: row

        det = rational_determinant(A)

        if (det < (0//1)) then
            row(:) = A(:,1)
            A(:,1) = A(:,2)
            A(:,2) = row(:)
        end if

    end subroutine Set_Right_Handedness

    !!---- Subroutine Smallest_Integral_Vector(v)
    !!----     type(rational), dimension(:), intent(inout) :: v
    !!----
    !!---- Finds the smallest set of integers for the direction given by v.
    !!---- It assumes denominators in v have no prime numbers greater than 97.
    !!----
    !!---- Updated: January - 2019
    !!

    subroutine Smallest_Integral_Vector(v)

        !---- Arguments ----!
        type(rational), dimension(:), intent(inout) :: v

        !---- Local variables ----!
        integer                :: i
        integer, dimension(25) :: primos = [ 2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97 ]
        type(rational), dimension(size(v)) :: vAux

        do i = 1 , size(v)
            v = v(i)%Denominator * v
        end do

        do i = 1 , 25
            vAux = v / (primos(i)//1)
            do
                if (.not. IsInteger(vAux)) exit
                v    = vAux
                vAux = v / (primos(i)//1)
            end do
        end do

    end subroutine Smallest_Integral_Vector

end module CFML_Standard_Settings