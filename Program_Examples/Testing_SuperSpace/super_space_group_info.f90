!!----
!!---- Program: Super Space Group Info (SSPG_INFO)
!!----          Made after the Space Group Info program in  CFML.
!!----          Restricted to Magnetic SSGs in particular.
!!---- Author:
!!---- Revision:
!!




Module SGk_to_SSGs

use CFML_Math_General,              only: Determinant, Equal_Matrix, &
                                            Equal_Vector, Modulo_Lat, Zbelong
use CFML_Crystallographic_Symmetry, only: Space_Group_Type, set_SpaceGroup, &
                                             Write_SpaceGroup, Symmetry_Symbol, &
                                             Get_Generators_From_SpGSymbol,     &
                                             err_symm, err_symm_mess
use CFML_GlobalDeps,                only: cp
use CFML_Propagation_Vectors,       only: Group_k_Type, Set_Gk, K_Star, &
                                             Write_Group_k, k_EQUIV
use CFML_String_Utilities,          only: Get_Fraction_1Dig
use CFML_SuperSpaceGroups
use CFML_ssg_datafile
use CFML_Rational_Arithmetic_test
use Matrix_Mod

implicit none
contains

!>>!-----------------------------------------------------

!Subroutine Set_SSGs_from_Gkk_MYTEST(SpG,nk,kv)!,ssg,nss)
Subroutine Set_SSGs_from_Gkk_MYTEST(SpG,nk,kv,ssg,nss)
      type(Space_Group_Type),                                intent(in)  :: SpG
      integer,                                               intent(in)  :: nk
      real(kind=cp),dimension(:,:),                          intent(in)  :: kv
      type(SuperSpaceGroup_Type), dimension(:), allocatable, intent(out) :: ssg
      integer,                                               intent(out) :: nss
      !--- Local variables ---!
      character(len=132) :: line
      integer :: i,j,k,Dd
      type(Space_Group_Type),dimension(nk) :: Gkks !extdd Litt. Groups
      type(Space_Group_Type)               :: Gkk  !extdd Litt. Group (Intersec. of LGs for all nk)
      type(SuperSpaceGroup_Type)           :: trial_ssg
      Type (Group_k_Type)                  :: Gk
      real(kind=cp), dimension(3+nk)       :: tr
      integer,       dimension(3+nk,3+nk)  :: Mat
      !
      !--- Local variables (for subroutines)
      real(kind=cp), dimension(:,:), allocatable       :: vectors_H_R
      integer                                          :: max_order, size_vec_all_orders
      integer, dimension(:), allocatable               :: vec_all_orders
      !
      !--- Local variables (for subroutines)
      integer,           dimension(10)                 :: point_op
      integer,           dimension(10)                 :: gen_set
      character(len=50), dimension(10)                 :: gen
      character(len=70)                                :: op_symb
      integer                                          :: ngen
      type(SSym_Oper_Type), dimension(:), allocatable  :: setOpGen


      WRITE(*,*) " entering in Set_SSGs_from_Gkk_MYTEST"


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

     !>Derive possible superspace groups from Gkk
     call get_H_R_vectors(Gkk, kv, vectors_H_R) !(Gkk, vectors_H_R, kv) where kv=dim(:,:)
     call get_Max_Order(Gkk, max_order)
     call Allocate_Vector_Of_Orders(Gkk, vec_all_orders)
     call get_Group_Orders(Gkk,vec_all_orders)


     !--- select one
     !
     call Get_Generators_From_SpGSymbol(Gkk,gen,point_op,ngen)
     !call Get_Generators_From_SpGOpList(Gkk, point_op, ngen)
     !--------------

     call Allocate_SSG_SymmOps(1, 8*Gkk%multip, setOpGen)
     call genToSgen( 1, ngen, point_op, Gkk, vectors_H_R, setOpGen )
     call SSGs_maker(point_op, ngen, vec_all_orders, setOpGen,ssg,nss)


!>>
!>>      !First: colorless groups of type 1
!>>      !Second: black and white groups of type 3
!>>      !Third: black and white groups of type 4
!>>

End Subroutine Set_SSGs_from_Gkk_MYTEST

!---------------------------------------------------------------

Subroutine Get_Generators_From_SpGOpList(SpG, point_op, ngen)
      Type(Space_Group_Type),             intent (in) :: SpG
      integer, dimension(10),             intent(out) :: point_op
      integer,                            intent(out) :: ngen
      !--- Local variables ---!
      integer                 :: i,j,n,m,nOps
      integer, dimension(3,3) :: testMatrix

      !--- Algorithm description
      !    passing an SPG as argument to this subroutine, the generators
      !    will be the following:
      !    1st (the Identity) and 2nd always will be generators
      !    From the structure of the list of operators of the SPG,
      !    if the 4th /= (2nd)**3, then, the 4th is gnerator

      !if(allocated(gen_set)) deallocate(gen_set)
      !allocate(gen_set(SPGk_in%multip))  ! maximum size is the group's


      !> IMPORTANT: keep out the Id of the gen list (shift operator numbers)
      !
      point_op   =0
      !-->point_op(1)=1 !if considering the ID as part of the gen list
      !-->point_op(2)=2
      point_op(1)=2 !keep out the Id of the gen list
      point_op(2)=3
      ngen=2
      testMatrix=matmul(SpG%SymOp(3)%Rot,SpG%SymOp(3)%Rot) ! (the 2nd)**2
      testMatrix=matmul(SpG%SymOp(3)%Rot, testMatrix)      ! (the 2nd)**3

      if( .not. Equal_Matrix(SpG%SymOp(5)%Rot, testMatrix, 3 )) then
       !point_op(3)=5 ! non 'collapsed' format
       point_op(3)=5  ! yes 'collapsed' format
       ngen=3
      end if

End Subroutine Get_Generators_From_SpGOpList

!---------------------------------------------------------------

Subroutine get_H_R_vectors(SPGk_out, vec_k, vecs_H_R)
   Type(Space_Group_Type),                     intent (in)  :: SpGk_out
   real(kind=cp), dimension(3),                intent (in)  :: vec_k
   real(kind=cp), dimension(:,:), allocatable, intent (out) :: vecs_H_R
   !--- Local variables
   integer :: i,j,k,Dd
   real(kind=cp), dimension(3)                :: vec_k_input
   real(kind=cp), dimension(3)                :: vec_k_prime
   real(kind=cp), dimension(3)                :: vec_k_condition1
   real(kind=cp), dimension(3)                :: vec_k_condition2

      !> allocate H(R) vector with the number of operators (SPGk_out%multip)
      WRITE(*,*) "number of  R_S matrices = ", SPGk_out%multip
      if(allocated(vecs_H_R)) deallocate(vecs_H_R)
      allocate(vecs_H_R(4,SPGk_out%multip))
      vecs_H_R=-10

      do j=1, SPGk_out%multip+0
          vec_k_prime=matmul(vec_k, SPGk_out%SymOp(j)%Rot)
          vec_k_condition1=vec_k_prime-vec_k ! k1=k'-k  <=>  k' = kR = +k + H(R)
          vec_k_condition2=vec_k_prime+vec_k ! k2=k'+k  <=>  k' = kR = -k + H(R)

          ! if k1=H(R) belongs to reciprocal lattice, then ...
          if(        k_EQUIV(vec_k_prime, vec_k,  SPGk_out%SPG_lat)    ) then
               vecs_H_R(1:3,j)=vec_k_condition1
               vecs_H_R(4:4,j)=1
          else if(   k_EQUIV(vec_k_prime, -vec_k, SPGk_out%SPG_lat)    ) then
               vecs_H_R(1:3, j)=vec_k_condition2
               vecs_H_R(4:4, j)=-1
          else
               WRITE(*,*) "Error computing the H(R)!"
          end if
      end do

      WRITE(*,*) ""

      !> printing the H(R_j) vectors plus the 4th component RI
      do j=1, SPGk_out%multip
      WRITE(*,*) "H(R)|RI=", vecs_H_R(1:4,j)
      end do
End Subroutine get_H_R_vectors

!---------------------------------------------------------------

Subroutine get_Max_Order(SPG_in, maxorder)
! this function returns the maximum order
! (defined as the exponent to become the Identity)
! of a group of rotational operators

type(Space_Group_type), intent(in)    :: SPG_in
integer,                intent(out)   :: maxorder

real(kind=cp), dimension(3,3)     :: ID, CM, TM
integer                           :: cnt, i, j
!integer, dimension(SPG_in%multip) :: vec_all_orders

TM=SPG_in%Symop(1)%rot ! Current Matrix
CM=SPG_in%Symop(1)%rot ! Initiallized to the Identity
ID=SPG_in%Symop(1)%rot ! Identity matrix is always the first of the SG

maxorder=0

out:do i=1, SPG_in%multip
   CM=SPG_in%Symop(i)%rot ! Current Matrix
   TM=SPG_in%Symop(1)%rot ! Test Matrix (initiallized to the Id)
   cnt=0
   in:do j=2,6 !order is never greater than 6
       cnt=cnt+1
       TM=matmul(TM,CM)
       if(Equal_Matrix(TM,ID,3)) then
          !cnt is the order of the i-operator; load out
          !vec_all_orders(i)=cnt
          if(cnt > maxorder) then
             maxorder=cnt
          end if
          cycle out
       else
          continue
       end if

  end do in
end do out

end Subroutine get_Max_Order

!---------------------------------------------------------------

Subroutine Allocate_Vector_Of_Orders(SPG, vec_group_orders)

type(Space_Group_type),              intent(in)     :: SPG
integer, dimension(:), allocatable,  intent(in out) :: vec_group_orders

integer                                             :: size_vec_group_orders

size_vec_group_orders=SPG%multip
allocate(vec_group_orders(size_vec_group_orders))

end Subroutine Allocate_Vector_Of_Orders

!---------------------------------------------------------------

Subroutine  get_Group_Orders(SPG_in, vec_all_orders)
! this function returns the maximum order
! (defined as the exponent to become the Identity)
! of a group of rotational operators

type(Space_Group_type), intent(in)             :: SPG_in
integer, dimension(SPG_in%multip), intent(out) :: vec_all_orders

real(kind=cp), dimension(3,3)     :: ID, CM, TM
integer                           :: cnt, i, j, maxorder

TM=SPG_in%Symop(1)%rot ! Current Matrix
CM=SPG_in%Symop(1)%rot ! Initiallized to the Identity
ID=SPG_in%Symop(1)%rot ! Identity matrix is always the first of the SG

maxorder=0

out:do i=1, SPG_in%multip
   CM=SPG_in%Symop(i)%rot ! Current Matrix
   TM=SPG_in%Symop(1)%rot ! Test Matrix (initz. to Id)
   cnt=0
   in:do j=2,6 !order is never greater than 6
       cnt=cnt+1
       TM=matmul(TM,CM)
       if(Equal_Matrix(TM,ID,3)) then
          !cnt is the order of the i-operator; load out
          vec_all_orders(i)=cnt
          if(cnt > maxorder) then
             maxorder=cnt
          end if
          cycle out
       else
          continue
       end if

  end do in
end do out

end Subroutine  get_Group_Orders

!---------------------------------------------------------------

Subroutine genToSgen( nk, ngen, point_op, SPGk_out, vecs_H_R, setOpGen)
      integer,                                            intent (in) :: nk, ngen
      integer, dimension(10),                             intent (in) :: point_op
      Type(Space_Group_Type),                             intent (in) :: SPGk_out
      real(kind=cp), dimension(:,:), allocatable,         intent (in) :: vecs_H_R
      type(SSym_Oper_Type), dimension(:), allocatable, intent(in out) :: setOpGen


      !--- Local variables ---!
      integer :: i,j,k,Dd
      character(len=80),    dimension(5,5)             :: cMat1, cMat2

call identity_matrix( 5, setOpGen(1)%Mat(:,:) ) ! automate 5
do k=1, ngen
  j=point_op(k) ! j is the number of SG generator from the SG generator ptr list
  !
  setOpGen(k+1)%Mat(:,:)    = 0_ik//1_ik
  setOpGen(k+1)%Mat(1:3,1:3)=SPGk_out%SymOp(j)%Rot(:,:)//1 ! ID sumbatrix Rot 3*3
  setOpGen(k+1)%Mat(4,1:4)  =vecs_H_R(1:4,j)
  setOpGen(k+1)%Mat(1:3,5)  =SPGk_out%SymOp(j)%Tr(1:3) ! sic, at fifth column
  setOpGen(k+1)%Mat(4,5)    =0_ik//1_ik                ! location of t4
  setOpGen(k+1)%Mat(5,5)    =1_ik//1_ik                ! sic, to complete
enddo

WRITE(*,*) " "
WRITE(*,*) " Checking SSG Generators---------------------------------------------"
do k=1, ngen !size(setOpGen)
WRITE(*,*) " "
WRITE(*,*) " setOpGen", k
cMat1=print_rational(setOpGen(k)%Mat(:,:))
 do i=1,5
   write(*,"(5a16)") (trim( cMat1(i,j))//"   ",j=1,5)
 end do
enddo
WRITE(*,*) " End Checking SSG Generators------------------------------------------"
WRITE(*,*) " "
End Subroutine genToSgen

!---------------------------------------------------

Subroutine SSGs_maker(point_op,ngen,vec_all_orderss,setOpGen,ssg,nss)

      integer, dimension(10),                                intent(in)     :: point_op
      integer,                                               intent(in)     :: ngen
      integer, dimension(:), allocatable ,                   intent(in)     :: vec_all_orderss
      type(SSym_Oper_Type), dimension(:), allocatable,       intent(in out) :: setOpGen
      type(SuperSpaceGroup_Type), dimension(:), allocatable, intent(out)    :: ssg
      integer,                                               intent(out)    :: nss

      !---Local variables
      integer :: i,j,k,Dd, ii,jj,kk
      integer                                      :: finalnumb
      integer                                      :: initlnumb
      integer, dimension(:,:), allocatable         :: genGroupTable
      character(len=70)                            :: op_symb

      character(len=15)                            :: forma
      integer                                      :: orderG1, orderG2, orderG3
      type (rational),   dimension(:), allocatable :: ordVector

      type(SuperSpaceGroup_Type)                   :: SSG_out


  orderG1=vec_all_orderss(point_op(1))
  orderG2=vec_all_orderss(point_op(2))
  orderG3=vec_all_orderss(point_op(3))

  initlnumb=ngen+1

  select case (ngen)

     case (1)
     print*, "case: 1 generator,", ngen
     nss=orderG1
     if(allocated(SSG)) deallocate(SSG)
     allocate(SSG(nss))
     kk=0
     do i=1, orderG1
        call orderVector(orderG1, ordVector)
        setOpGen(2)%Mat(4,5)    = ordVector(i) ! location of t4

        finalnumb=0
        Call Gen_SSGroup(initlnumb, setOpGen, SSG_out ,"xyz", genGroupTable)
        if (Err_ssg) then
           WRITE(*,*)" ------", Err_ssg_mess
           cycle
        endif
        !
        kk=kk+1
        ssg(kk)=SSG_out
        !
        WRITE(*,*)"  "
     enddo



     case (2)
     print*, "case: 2 generators,", ngen
     nss=orderG1*orderG2
     if(allocated(SSG)) deallocate(SSG)
     allocate(SSG(nss))
     kk=0
     do i=1, orderG1
        call orderVector(orderG1, ordVector)
        setOpGen(2)%Mat(4,5)    = ordVector(i)
        do j=1, orderG2
           call orderVector(orderG2, ordVector)
           setOpGen(3)%Mat(4,5)    = ordVector(j)

           finalnumb=0
           Call Gen_SSGroup(initlnumb, setOpGen, SSG_out ,"xyz", genGroupTable)
           if (Err_ssg) then
              WRITE(*,*)" ------", Err_ssg_mess
              cycle
           endif
           !
           kk=kk+1
           ssg(kk)=SSG_out
           !
           WRITE(*,*)"  "
        enddo
     enddo


     case (3)
     print*, "case: 3 generators,", ngen
     nss=orderG1*orderG2*orderG3
     if(allocated(SSG)) deallocate(SSG)
     allocate(SSG(nss))
     kk=0
     do i=1, orderG1
        call orderVector(orderG1, ordVector)
        setOpGen(2)%Mat(4,5)    = ordVector(i)
        do j=1, orderG2
           call orderVector(orderG2, ordVector)
           setOpGen(3)%Mat(4,5)    = ordVector(j)
           do k=1, orderG3
              call orderVector(orderG3, ordVector)
              setOpGen(4)%Mat(4,5)    = ordVector(k)

              finalnumb=0
              Call Gen_SSGroup(initlnumb, setOpGen, SSG_out ,"xyz", genGroupTable)
              if (Err_ssg) then
                 WRITE(*,*)" ------", Err_ssg_mess
                 cycle
              endif
              !WRITE(*,*)" ------->  initial generators = ", ngen+1
              !WRITE(*,*)" ------->    final group size = ", finalnumb, "=",SSG_out%multip
              !WRITE(*,*)" ------->         num of SSGs = ",nss
              !
              kk=kk+1
              ssg(kk)=SSG_out
              !
              WRITE(*,*)"  "
           enddo
        enddo
     enddo
     if (kk .NE. nss) then
       WRITE(*,*) "bug here: dirty trick failed "
     else
       nss=kk ! subroutine (unnecessary) output
     endif



     case default
        print*, "Error! 0 generators"

  end select

End Subroutine SSGs_maker

!---------------------------------------------------


!---------------------------------------------------

Subroutine orderVector(ord, ordVec)
  ! identifies the vector of orders with its order (or size)
  integer, intent(in) :: ord
  type(rational), dimension(:), allocatable, intent(out) :: ordVec

   Type(rational),   dimension(2) :: o2
   Type(rational),   dimension(3) :: o3
   Type(rational),   dimension(4) :: o4
   Type(rational),   dimension(6) :: o6

   o2=[ 0//2, 1//2 ]
   o3=[ 0//3, 1//3, 2//3 ]
   o4=[ 0//4, 1//4, 2//4, 3//4 ]
   o6=[ 0//6, 1//6, 2//6, 3//6, 4//6, 5//6 ]

   select case (ord)

   case (2)
      if(allocated(ordVec)) deallocate(ordVec)
      allocate(ordVec(ord))
      ordVec=o2
   case (3)
      if(allocated(ordVec)) deallocate(ordVec)
      allocate(ordVec(ord))
      ordVec=o3
   case (4)
      if(allocated(ordVec)) deallocate(ordVec)
      allocate(ordVec(ord))
      ordVec=o4
   case (6)
      if(allocated(ordVec)) deallocate(ordVec)
      allocate(ordVec(ord))
      ordVec=o6
   end select


End Subroutine orderVector

!---------------------------------------------------
end module SGk_to_SSGs




!>>-------------------------------------------------
!>> MAIN
!>>-------------------------------------------------

Program SSPG_Info

   !---- Use Modules ----!
   use SGk_to_SSGs


   !---- Variables ----!
   implicit none


   !---variables for  STEP 1:
   !                  Prompts the user to type a SG and a  k-vector
   !                  Computes the SG's Small Group (or Group of k)
   !                  The single k has to be written as a general case
   !                  of several of them.
   !
   character(len=20)             :: spgr
   type(Space_Group_type)        :: grp_espacial, SPGk_in, SPGk_out
   type(Group_k_Type)            :: Gk_out
   character(len=1)              :: default_example ! len=1 ('y' or 'n')
   real(kind=cp), dimension(3)   :: vec_k
   integer                       :: i,j,k
   integer                       :: size_Gk


   !---variables for  STEP 2:
   !                  subroutine with all steps to create
   !                  all the compatible SSGs with the initial SG+k
   !
   integer nkvecs
   real, dimension(:,:), allocatable                      :: kvecs
   type(SuperSpaceGroup_Type), dimension(:), allocatable  :: SSGs_out
   integer                                                :: nSSGs



   !---- Procedure ----!
   do
      WRITE(*,*) " "
      WRITE(*,*) " "
      WRITE(*,*) " ========================================================================"
      WRITE(*,*) " "
      WRITE(*,*) " STEP 1: find the SMALL GROUP"
      WRITE(*,*) " "

      ! Overriding with a default input example
      !
      !write(unit=*,fmt="(a)") " => Enter a space group: "
      !write(unit=*,fmt="(a)",advance="no") " => Space Group (HM/Hall symbol or number): "
      !read(unit=*,fmt="(a)") spgr
      !overriding: default example
      !spgr="P n m a"
      !WRITE(*,*) "your input:", spgr

      !if (len_trim(spgr)==0) exit
      !write(unit=*,fmt="(a)",advance="no") " => Enter a a k-vector: "
      !read(unit=*,fmt=*) vec_k
      !overriding: default example
      !vec_k=(/ 0.321, 0.0, 0.5 /)

      spgr=""
      vec_k=0
      write(unit=*,fmt="(a)") " => Would you like the default Space Group and h-vector"
      write(unit=*,fmt="(a)") " => P n m a"
      write(unit=*,fmt="(a)") " => 0.321 0 0.5"
      write(unit=*,fmt="(a)") " => y/n?"
      read(unit=*,fmt=*) default_example
      if(default_example=="y") then
          write(unit=*,fmt="(a)") "Using default example"
          spgr="P n m a"
          vec_k=(/ 0.321, 0.0, 0.5 /)
      else
          !> New prompt
          write(unit=*,fmt="(a)") " => Enter a space group: "
          write(unit=*,fmt="(a)",advance="no") " => Space Group (HM/Hall symbol or number): "
          read(unit=*,fmt="(a)") spgr
          if (len_trim(spgr)==0) exit
          write(unit=*,fmt="(a)",advance="no") " => Enter a k-vector: "
          read(unit=*,fmt=*) vec_k
      end if

      !> Setting and writing the SG Information
      !> Casting the SG variable to decouple this part of the code from the following
      WRITE(*,*) " setting and writing grp_espacial:"
      call Set_SpaceGroup(spgr, grp_espacial)
      call Write_SpaceGroup(grp_espacial, full=.true.)
      size_Gk=grp_espacial%multip
      SPGk_in=grp_espacial ! casting


WRITE(*,*) " "
WRITE(*,*) " ========================================================================"
WRITE(*,*) " "
WRITE(*,*) " STEP 2:  using Set_SSGs_from_Gkk"
WRITE(*,*) " "



!--- the final goal is to write Set_SSGs_from_Gkk
!    we adapt here the new parameters: nk, kv
!
if(allocated(kvecs)) deallocate(kvecs)
allocate(kvecs(1,3))
nkvecs=1 ! num of k-vectors
kvecs(1,:)=vec_k


!---++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
call Set_SSGs_from_Gkk_MYTEST(SPGk_in, nkvecs, kvecs, SSGs_out, nSSGs)
!
!---++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
write(*,*) "final results:"
write(*,*) "num of SSGSs", nSSGs
do i=1, nSSGs
 call Write_SSG(SSGs_out(i),full=.true.)
enddo
!
!---++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



WRITE(*,*) " "
WRITE(*,*) " ========================================================================"
WRITE(*,*) " "




end do ! prompting loop



!!----------------------------------------------------------------------------

End Program SSPG_Info
