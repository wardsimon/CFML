!!----
!!---- Program: Super Space Group Info (SSPG_INFO)
!!----          Made after the Space Group Info program in  CFML.
!!----          Restricted to Magnetic SSGs in particular.
!!----
!!---- Author: 
!!---- Revision: 
!!

!-----------------------------------------------------------------

!-----------------------------------------------------------------
!---> Main
!-----------------------------------------------------------------

Module pepito

   use CFML_Crystallographic_Symmetry, only: Space_Group_Type, set_SpaceGroup, &
                                             Write_SpaceGroup, Symmetry_Symbol, &
                                             Get_Generators_From_SpGSymbol,     &
                                             err_symm, err_symm_mess
   use CFML_GlobalDeps,                only: cp
   use CFML_Propagation_Vectors,       only: Group_k_Type, Set_Gk, K_Star, &
                                             Write_Group_k, k_EQUIV
   use CFML_String_Utilities,          only: Get_Fraction_1Dig
   use CFML_Math_General,              only: Determinant, Equal_Matrix, & 
                                             Equal_Vector, Modulo_Lat, Zbelong
   use CFML_SuperSpaceGroups
   use CFML_ssg_datafile
   use CFML_Rational_Arithmetic_test
   use Matrix_Mod

implicit none ! 

contains

function GET_MAX_ORDER(SPG_in) result(maxorder)
! this function returns the maximum order 
! (defined as the exponent to become the Identity) 
! of a group of rotational operators

type(Space_Group_type), intent(in)    :: SPG_in
!logical,   optional,    intent(in)    :: full
integer                               :: maxorder

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

!if (present(full))  then
!  return vec_all_orders
!else
  return
!end if
end function GET_MAX_ORDER

!-----------------------------------------------------

function GET_POSSIBLE_t4_VALUES(maximum_order) result(vec_possible_t4_values)
! this function returns a vector with the possible t4 values
! corresponding to a maximum order (matrix potentiation)
!
integer, intent(in)             :: maximum_order
real, dimension(:), allocatable :: vec_possible_t4_values
!
integer :: init_size, inst_size, max_order, num, den, final_size
real    :: t4
real, dimension(:), allocatable  :: vectorA
real, dimension(:), allocatable  :: vectorB


init_size=0  ! start with a zero-sized vector
inst_size=init_size

max_order=maximum_order

!> using vectorB as copy buffer
allocate(vectorA(init_size))
allocate(vectorB(init_size))
vectorA=-11
vectorB=-11


do num=0, max_order-1
   do den=1, max_order
      t4=real(num)/real(den)
      !WRITE(*,*) "a possible t4 is: ", num, " / ",den, " = ", t4
      
      if( Zbelong( t4 ) .and. (t4/=0)) then
         !WRITE(*,*) "        exclude ", t4 
         !cycle
         continue
      else
         !WRITE(*,*) "        accept ", t4
         if ( ANY( vectorA==t4 ) ) then
            !WRITE(*,*) "        is present: do nothing "
            !cycle
            continue
         else
            !WRITE(*,*) "        is impresent: increase vec size and add "
             
            !> BEFORE: copy A -> B
            if(allocated(vectorB)) deallocate(vectorB)
            allocate(vectorB(inst_size))
            vectorB=vectorA
            vectorA(inst_size)=t4 
            
            !> DURING: increase size, reallocate A
            inst_size=inst_size+1
            if(allocated(vectorA)) deallocate(vectorA)
            allocate(vectorA(inst_size))
            
            !> AFTER: copy B -> A, add new to A 
            vectorA(1:inst_size-1)=vectorB ! sic, since different sizes
            vectorA(inst_size)=t4
            
            !WRITE(*,*) "initial vector------------------------------------------------------- a     ", vectorA
            !WRITE(*,*) "initial vector------------------------------------------------------- b     ", vectorB
            
         end if
      end if
      !        
   end do ! num
end do ! den

! casting the results
final_size=size(vectorA)
if(allocated(vec_possible_t4_values)) deallocate(vec_possible_t4_values)
allocate(vec_possible_t4_values(final_size))
vec_possible_t4_values=vectorA

return

end function GET_POSSIBLE_t4_VALUES

!-----------------------------------------------------
!-----------------------------------------------------
function GET_GROUP_ORDERS(SPG_in) result(vec_all_orders)
! this function returns the maximum order 
! (defined as the exponent to become the Identity) 
! of a group of rotational operators

type(Space_Group_type), intent(in)    :: SPG_in
!logical,   optional,    intent(in)    :: full
integer, dimension(SPG_in%multip)     :: vec_all_orders

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

return

end function GET_GROUP_ORDERS


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

    Subroutine Set_SSGs_from_Gkktt(SpG,nk,kv)!,ssg,nss)
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

      call get_H_R_vectors(Gkk, vecs_H_R)

      max_order              = GET_MAX_ORDER(Gkk) 
      vec_T4_possible_values = GET_POSSIBLE_t4_VALUES(max_order) 
      size_vec_all_orders=Gkk%multip
      !allocate(vec_all_orders(size_vec_all_orders))
      vec_all_orders=GET_GROUP_ORDERS(Gkk) 
      WRITE(*,*) " vector of SG orders = ", vec_all_orders

 
      !call Get_Generators_From_SpGSymbol(Gkk,gen,point_op,ngen)
      call Get_Generators_From_SpGOpList(Gkk, point_op,ngen) ! falla
      
      call Allocate_SSG_SymmOps(1, 8*Gkk%multip, setOpGen) !
      call genToSgen( 1, ngen, point_op, Gkk, setOpGen )
      call SSG_maker(point_op, ngen, setOpGen)



      !First: colorless groups of type 1
      !Second: black and white groups of type 3
      !Third: black and white groups of type 4

    End Subroutine Set_SSGs_from_Gkktt



    Subroutine Get_Generators_From_SpGOpList(SpG, point_op, ngen)
      Type(Space_Group_Type),             intent (in) :: SpG
      integer, dimension(10),             intent(out) :: point_op
      integer, intent(out)                            :: ngen
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

      write(*,*) "looking for generators in the list "
      point_op   =0
      point_op(1)=1 
      point_op(2)=2 
      ngen=2
      testMatrix=matmul(SpG%SymOp(2)%Rot,SpG%SymOp(2)%Rot) ! (the 2nd)**2
      testMatrix=matmul(SpG%SymOp(2)%Rot, testMatrix)      ! (the 2nd)**3

      if( .not. Equal_Matrix(SpG%SymOp(4)%Rot, testMatrix, 3 )) then
       !point_op(3)=4 ! non 'collapsed' format
       point_op(3)=4  !     'collapsed' format
       ngen=3
      end if

    End Subroutine Get_Generators_From_SpGOpList 


Subroutine genToSgen( nk, ngen, point_op, SPGk_out, setOpGen)
      integer,                                        intent (in) :: nk, ngen
      integer, dimension(10),                         intent (in) :: point_op
      Type(Space_Group_Type),                         intent (in) :: SPGk_out
      type(SSym_Oper_Type), dimension(:), allocatable, intent(in out):: setOpGen 


      !--- Local variables ---!
 
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


Subroutine SSG_maker(point_op,ngen,setOpGen)

integer,           dimension(10),                intent(in)     :: point_op
integer,                                         intent(in)     :: ngen
type(SSym_Oper_Type), dimension(:), allocatable, intent(in out) :: setOpGen 

!---Locals
integer                                 :: finalnumb
integer, dimension(:,:), allocatable    :: genGroupTable
character(len=70)                       :: op_symb

!call SSG_maker(setOpGen, setOpGen_size)

  orderG1=vec_all_orders(point_op(1)) 
  orderG2=vec_all_orders(point_op(2)) 
  orderG3=vec_all_orders(point_op(3)) 


  select case (ngen)
  
     case (1) 
     print*, "case: 1 generator" 
     do i=1, orderG1
        call orderVector(orderG1, ordVector)
        setOpGen(2)%Mat(4,5)    = ordVector(i) ! location of t4
        
        finalnumb=0  
        Call Gen_Group(3, setOpGen, finalnumb, genGroupTable) 
        !Call Gen_SSGroup(ngen,gen,SSG,x1x2x3_type,table)
        if (Err_ssg) then
           WRITE(*,*)" ------", Err_ssg_mess 
           cycle
        endif
        WRITE(*,*)" ------->  initial generators = ", ngen+1
        WRITE(*,*)" ------->    final group size = ", finalnumb
        !check generated
        do ii=1, finalnumb
           call Get_SSymSymb_from_Mat(setOpGen(ii)%Mat, op_symb,"xyz")
           write(*,*) ii, trim(op_symb)
        enddo
        WRITE(*,*)"  "
     enddo



     case (2)
     print*, "case: 2 generators" 
     do i=1, orderG1 
        call orderVector(orderG1, ordVector)
        setOpGen(2)%Mat(4,5)    = ordVector(i)    
        do j=1, orderG2  
           call orderVector(orderG2, ordVector)
           setOpGen(3)%Mat(4,5)    = ordVector(j)      
           
           finalnumb=0  
           Call Gen_Group(3, setOpGen, finalnumb, genGroupTable) 
           if (Err_ssg) then
              WRITE(*,*)" ------", Err_ssg_mess 
              cycle
           endif
           WRITE(*,*)" ------->  initial generators = ", ngen+1
           WRITE(*,*)" ------->    final group size = ", finalnumb
           !check 
           do ii=1, finalnumb
              call Get_SSymSymb_from_Mat(setOpGen(ii)%Mat, op_symb,"xyz")
              write(*,*) ii, trim(op_symb)
           enddo
           WRITE(*,*)"  "
        enddo
     enddo

  
     case (3) 
     print*, "case: 3 generators"
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
              Call Gen_Group(4, setOpGen, finalnumb, genGroupTable) 
              if (Err_ssg) then
                 WRITE(*,*)" ------", Err_ssg_mess 
                 cycle
              endif
              WRITE(*,*)" ------->  initial generators = ", ngen+1
              WRITE(*,*)" ------->    final group size = ", finalnumb
              !check 
              do ii=1, finalnumb
                 call Get_SSymSymb_from_Mat(setOpGen(ii)%Mat, op_symb,"xyz")
                 write(*,*) ii, trim(op_symb)
              enddo
              WRITE(*,*)"  "
           enddo
        enddo
     enddo

     case default
        print*, "Error! 0 generators" 
     
  end select
 

End Subroutine SSG_maker


Subroutine get_H_R_vectors(SPGk_out, vecs_H_R)
   Type(Space_Group_Type),                     intent (in)  :: SpGk_out
   real(kind=cp), dimension(:,:), allocatable, intent (out) :: vecs_H_R

   !--- Locals
   real(kind=cp), dimension(3)                :: vec_k_input
   real(kind=cp), dimension(3)                :: vec_k_prime
   real(kind=cp), dimension(3)                :: vec_k_condition1
   real(kind=cp), dimension(3)                :: vec_k_condition2
   !Type (Space_Group_Type)                    :: SPGk_input


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


end module pepito

Program SSPG_Info

   !---- Use Modules ----!
   
   use CFML_Crystallographic_Symmetry, only: Space_Group_Type, set_SpaceGroup, &
                                             Write_SpaceGroup, Symmetry_Symbol, &
                                             Get_Generators_From_SpGSymbol,     &
                                             err_symm, err_symm_mess
   use CFML_GlobalDeps,                only: cp
   use CFML_Propagation_Vectors,       only: Group_k_Type, Set_Gk, K_Star, &
                                             Write_Group_k, k_EQUIV
   use CFML_String_Utilities,          only: Get_Fraction_1Dig
   use CFML_Math_General,              only: Determinant, Equal_Matrix, & 
                                             Equal_Vector, Modulo_Lat, Zbelong
   use CFML_SuperSpaceGroups
   use CFML_ssg_datafile
   use CFML_Rational_Arithmetic_test
   use Matrix_Mod
   use pepito

   !use CFML_Crystallographic_Symmetry 

   !---- Variables ----!
   implicit none

   !>
   !> final result must be written as a SuperSpaceGroup_Type
   !> use notation similar to SG write()
   !> type(SuperSpaceGroup_Type),dimension(:),allocatable :: SSG_out 
   !> use write_ssg()
   !>


   !---variables for  STEP 1: 
   !                  Prompts the user to type a SG and a  k-vector
   !                  Computes the SG's Small Group (or Group of k)
   !
   character(len=20)             :: spgr
   type(Space_Group_type)        :: grp_espacial, SPGk_in, SPGk_out
   type(Group_k_Type)            :: Gk_out
   character(len=1)              :: default_example ! len=1 ('y' or 'n')  
   real(kind=cp), dimension(3)   :: vec_k
   integer                       :: size_Gk
   !
   !                  For the retrieval of the SG generators
   integer,           dimension(10) :: point_op
   integer,           dimension(10) :: gen_set

   character(len=50), dimension(10) :: gen
   character(len=70)                :: op_symb
   integer                          :: ngen



   !---variables for  STEP 2: 
   !                  Compute the SG's Small Group Extended for the k-vector.
   !                  Cast the result to a vector of type(SSym_Oper_Type)
   !
   integer :: size_Gkmink
   integer :: i,j,k


   !---variables for  STEP 3: 
   !                  Compute the reciprocal vectors {H(R)}
   !                  that fulfill conditions 1 and  2, namely the small group 
   !                  (G_{k}) and extended small group (G_{k,-k}), respectively
   !
   real(kind=cp), dimension(:,:), allocatable :: vecs_H_R


   !---variables for  STEP 4:
   !                  Compute the possible values t4 and store in vector
   !
   integer                             :: m,n, max_order, size_vec_all_orders
   real, dimension(:), allocatable     :: vec_T4_possible_values
   integer, dimension(:), allocatable  :: vec_all_orders



   !---variables for  STEP 5:
   !                  Compound the (3+d+1) matrices as a type(SSym_Oper_Type)
   !                  leaving just t4(=44) indetermined. 
   !                  Note: include the extra operator 'closeness of SSG'
   !                  Note: this step has been skipped and merged with 6
   !
   type(SSym_Oper_Type), dimension(:), allocatable  :: setOpGen 
   character(len=80),    dimension(5,5)             :: cMat1, cMat2



   !---variables for  STEP 6-7:
   !                  work with just the generators of the group plus 1'
   !
   integer nGenerators, nmax, ncmb, ngen_ssg, num, den, ii
   character(len=15)                            :: forma
   integer                                      :: orderG1, orderG2, orderG3
   type (rational),   dimension(:), allocatable :: ordVector


   !---variables for  STEP 8:
   !                  final subroutine
   !
   integer nk
   real, dimension(:,:), allocatable :: kv
   !Type (SuperSpaceGroup_Type)                     :: SSG_out


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
      
      call K_Star(vec_k, SPGk_in, Gk_out, ext=.true.)  ! ext=.true. gives G_{k,-k} 
      call Write_Group_k(Gk_out)
      call Set_Gk(Gk_out, SPGk_out, ext=.true.)
      call Write_SpaceGroup(SPGk_out, full=.true.) ! SPGk_out is now the EG

 

!--- using the final goal: Set_SSGs_from_Gkk
!  
if(allocated(kv)) deallocate(kv) 
allocate(kv(1,3)) 
nk=1 ! num of k-vectors
kv(1,:)=vec_k 




!---++++++++++++++++++++++++++++++++++++++++++++++++++
!
call Set_SSGs_from_Gkktt(SPGk_out, nk, kv)!,ssg,nss)
!
!---++++++++++++++++++++++++++++++++++++++++++++++++++
              





WRITE(*,*) " "
WRITE(*,*) " ========================================================================"
WRITE(*,*) " "
WRITE(*,*) " STEP 9:"
WRITE(*,*) " "


!> writing final results in SSG info format
!>
!WRITE(*,*) ""
!WRITE(*,*) ""
!SSG_out%SymOp=setOpGen
!call Write_SSG(SSG_out) 
!WRITE(*,*) ""
!WRITE(*,*) ""
!

end do ! prompting loop



!--- FUNCTIONS--------------------------------------------

!contains


!!----------------------------------------------------------------------------

End Program SSPG_Info
