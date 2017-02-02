!!----
!!---- Program: SPG_INFO
!!----          Example of simple program using CFML
!!----
!!---- Author: Juan Rodriguez-Carvajal
!!---- Revision: Nov-2008
!!

Module STEP_1_compute_small_group
!---modules

implicit none
!---variables

contains
!---functions

End Module STEP_1_compute_small_group


Program SPG_Info

   !---- Use Modules ----!
   use STEP_1_compute_small_group
   use CFML_Crystallographic_Symmetry, only: Space_Group_Type, set_SpaceGroup, &
                                             Write_SpaceGroup
   use CFML_GlobalDeps,                only: cp
   use CFML_Propagation_Vectors,       only: Group_k_Type, Set_Gk, K_Star, &
                                             Write_Group_k, k_EQUIV
   use CFML_String_Utilities,          only: Get_Fraction_1Dig
   use CFML_Math_General,              only: Determinant, Equal_Matrix, & 
                                             Equal_Vector, Modulo_Lat, Zbelong


   !---- Variables ----!
   implicit none

   !---variables for  STEP 1: 
   !                  Prompting the user to enter a space group and
   !                  computing the small group (Group k)
   !
   character(len=20)             :: spgr
   type(Space_Group_type)        :: grp_espacial, SPGk_in, SPGk_out
   type(Group_k_Type)            :: Gk_out
   character(len=1)              :: default_example ! len=1, just 'y' or 'n'
   real(kind=cp), dimension(3)   :: vec_k
   integer                       :: size_Gk

   !---variables for  STEP 2: 
   !                  Compute the small group extended
   !
   integer :: size_Gkmink 

   !---variables for  STEP 3: 
   !                  Compute the reciprocal vectors {H(R)}
   !                  that fulfill condition 1 and condition 2, namely the
   !                  small group (Gk) and extended small group (Gkmink)
   !                  respectively 
   !
   real(kind=cp), dimension(3)   :: vec_k_input
   real(kind=cp), dimension(3)   :: vec_k_prime
   real(kind=cp), dimension(3)   :: vec_k_condition1
   real(kind=cp), dimension(3)   :: vec_k_condition2
   real(kind=cp), dimension(:,:), allocatable :: vecs_H_R
   Type (Space_Group_Type)       :: SPGk_input
   integer                       :: i,j

   !---variables for  STEP 4:
   !                  Compute the possible values t4
   !real, dimension(6)   :: m_values
   character(len=512)   :: aux_string
   integer              :: max_order
   integer              :: m,n
   real                 :: t4

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
          !write(unit=*,fmt="(a)") "Hacking the code to put the right prompt to the user"
          write(unit=*,fmt="(a)") " => Enter a space group: "
          write(unit=*,fmt="(a)",advance="no") " => Space Group (HM/Hall symbol or number): "
          read(unit=*,fmt="(a)") spgr
          if (len_trim(spgr)==0) exit
          write(unit=*,fmt="(a)",advance="no") " => Enter a a k-vector: "
          read(unit=*,fmt=*) vec_k
      end if

      !> Setting the Space Group Information and
      !> writing the SpaceGroup Information.
      WRITE(*,*) " setting and writing grp_espacial:"
      call set_spacegroup(spgr,grp_espacial)
      call Write_SpaceGroup(grp_espacial, full=.true.)
      size_Gk=grp_espacial%multip  ! size of Gk, correct? needed for the vector of indices



      WRITE(*,*) " "
      WRITE(*,*) " ============= END OF ORIGINAL EXAMPLE. ========"
      WRITE(*,*) " "

      !> Recasting the SpaceGroup variable, to decouple this new part of the code
      !> from the original 'space_group_info' example.
      !> So, we just rename the space group input by the user as 'SPGk_in':
      SPGk_in=grp_espacial


      WRITE(*,*) " ========================================================================"
      WRITE(*,*) " "
      WRITE(*,*) " STEP 2: find the SMALL GROUP EXTENDED "
      WRITE(*,*) " "

      !> Computing the K-Star (i.e. the small group or Gk)
      !> NOTE: if using ext=.true. we get the extended group Gk-k
      !> The output will be loaded into Gk_out so we have to write it out:
      WRITE(*,*) " "
      !WRITE(*,*) " +++++++++++++++++++++++++++++++1"
      WRITE(*,*) " "
      call K_Star(vec_k, SPGk_in, Gk_out, ext=.true.)  !with ext=true computes Gk-k, the extended group
      call Write_Group_k(Gk_out)
      call Set_Gk(Gk_out, SPGk_out, ext=.true.)
      call Write_SpaceGroup(SPGk_out, full=.true.)
      size_Gkmink=SPGk_out%multip
      !where SPGk_out es Gk-k is the extended group

      WRITE(*,*) " "
      WRITE(*,*) " ========================================================================"
      WRITE(*,*) " "
      WRITE(*,*) " STEP 3: find the H(R) vectors and R_S matrices "
      WRITE(*,*) " "

      !call here the subrutine that makes the R_S matrices (simplified)
      WRITE(*,*) "number of  R_S matrices = ", SPGk_out%multip

      !passing parameters (put this as subroutine)
      !vec_k_input=vec_k
      !SPGk_input=SPGk_out

      if(allocated(vecs_H_R)) deallocate(vecs_H_R)
      allocate(vecs_H_R(4,SPGk_out%multip)) 
      vecs_H_R=-10

      do j=1, SPGk_out%multip
          vec_k_prime=matmul(vec_k, SPGk_out%SymOp(j)%Rot) ! ese SPGk_out
          vec_k_condition1=vec_k_prime-vec_k ! k1=k'-k  <=>  k' = kR = +k + H(R)   
          vec_k_condition2=vec_k_prime+vec_k ! k2=k'+k  <=>  k' = kR = -k + H(R)

          ! if k1=H(R) belongs to reciprocal lattice
          if(        k_EQUIV(vec_k_prime, vec_k,  SPGk_out%SPG_lat)    ) then  
               vecs_H_R(1:3,j)=vec_k_condition1
               vecs_H_R(4:4,j)=1 
          else if(   k_EQUIV(vec_k_prime, -vec_k, SPGk_out%SPG_lat)    ) then
               vecs_H_R(1:3, j)=vec_k_condition2
               vecs_H_R(4:4, j)=-1     
          else
               WRITE(*,*) "Error computing the H(R)!!"
          end if
      end do

      WRITE(*,*) ""

      !printing the H(R_j) vectors plus the 4th component RI
      do j=1, SPGk_out%multip
      WRITE(*,*) "H(R)|RI=", vecs_H_R(1:4,j)
      end do


      WRITE(*,*) " "
      WRITE(*,*) " ========================================================================"
      WRITE(*,*) " "
      WRITE(*,*) " STEP 4: find the t4(=m/n) values"
      WRITE(*,*) " "
      !NOTE: n = maximum order of the set of matrices
      WRITE(*,*) " "

      !m_values  = (/0,1,2,3,4,6/) !not needed
      max_order = GET_MAX_ORDER(SPGk_out)
 
      WRITE(*,*) " max_order",max_order
      WRITE(*,*) " "

      !NOTE: there's no need to compute the possibilities
      !      we just loop the m=0,...,max_order
      !      and the from     n=1,...,max_order 
      t4=-1
      do m=0,max_order
        do n=1,max_order
        t4=real(m)/real(n)
        WRITE(*,*) "a possible t4 is: ",m,n,t4
        end do
     end do

     


! la solucion es: los posibles grupos tal son:
! las Rs 4*4 y las Tr , con la t4 siendo los posibles t4 

!intento de solution:
do i=1, SPGk_out%multip
WRITE(*,*) " "
WRITE(*,*) " SOL.",i,"     Rs=",SPGk_out%SymOp(i)%Rot
WRITE(*,*) " SOL.",i,"   H(R)=",vecs_H_R(1:4,i)
WRITE(*,*) " SOL.",i," T(1:3)=",SPGk_out%SymOp(i)%Tr
WRITE(*,*) " SOL.",i," T(4)  = 0, 0.5, 1, 2"  ! NOO
WRITE(*,*) " "
end do

! cierre, operador inv. magnetica




      WRITE(*,*) " "
      WRITE(*,*) " ========================================================================"
      WRITE(*,*) " "
      WRITE(*,*) " STEP 5:"
      WRITE(*,*) " "



   end do ! prompting loop

!--- FUNCTIONS--------------------------------------------

contains

function GET_MAX_ORDER(SPG_in) result(maxorder)
! this function returns the maximum order 
! (defined as the exponent to become the Identity) 
! of a group of rotational operators

type(Space_Group_type), intent(in)    :: SPG_in
integer                               :: maxorder

real(kind=cp), dimension(3,3) :: ID, CM, TM 
integer                       :: cnt, i, j

TM=SPG_in%Symop(1)%rot ! Current Matrix
CM=SPG_in%Symop(1)%rot ! Initiallized to the Identity
ID=SPG_in%Symop(1)%rot ! Identity matrix is always the first of the SG

maxorder=0

out:do i=1, SPG_in%multip
   CM=SPG_in%Symop(i)%rot ! Current Matrix
   TM=SPG_in%Symop(1)%rot ! Test Matrix (initz. to Id)
   cnt=0 
   in:do j=2,6 !order is never gerater then 6
       cnt=cnt+1
       TM=matmul(TM,CM)
       if(Equal_Matrix(TM,ID,3)) then
          if(cnt > maxorder) then
             maxorder=cnt
          end if
          cycle out
       else
          continue
       end if
  
  end do in
end do out

!WRITE(*,*) " exiting GET_MAX_ORDER"

return
end function GET_MAX_ORDER

!-----------------------------------------------------
!function GET_MAX_ORDER(maximum_order_of_matrix_set) result(maxorder)

!return
!end function GET_MAX_ORDER

!-----------------------------------------------------

!-----------------------------------------------------

End Program SPG_Info
