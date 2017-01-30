!!----
!!---- Program: SPG_INFO
!!----          Example of simple program using CFML
!!----          to make the superspace magnetic groups
!!----
!!----
!!---- Author:
!!---- Revision:
!!
!----- Default input:
!
! => Enter a space group:
! => Space Group (HM/Hall symbol or number): P n m a
! => Enter a a k-vector: 0.321 0 0.5





Program SSPG_Info
   !---- Use Modules ----!
   use CFML_Crystallographic_Symmetry, only: Space_Group_Type, set_SpaceGroup, &
                                             Write_SpaceGroup
   use CFML_String_Utilities,          only: Get_Fraction_1Dig

   use CFML_Propagation_Vectors,       only: Group_k_Type, Set_Gk, K_Star, &
                                             Write_Group_k, k_EQUIV

   use CFML_GlobalDeps,                only: cp

   use CFML_Math_General, only: Determinant, Equal_Matrix, Equal_Vector, Modulo_Lat, Zbelong

   !---- Variables ----!

   implicit none

   integer :: i, j, m, n, sizegg, iter, cnt1, cnt2, cnt3, nmax, k, k1, k2, t1, t2, possibleCombinations, seclimit, l
   integer :: nmatrices, nmatricesRS, npossibleTs, size_Gk, size_Gkmink, size_m, size_mn, theTrueSize, dim_buff, sizeResults

   real                       :: anImprobableNumber, buff1, buff2

   !--- Variables of the original example:
   !
   !--- To prompt the user: "Enter a space group"
   !
   character(len=20)      :: spgr
   character(len=512)     :: aux_string
   type(Space_Group_type) :: grp_espacial
   character(len=1) :: default_example ! len=1 since choice is just 'y' or 'n'

   !--- Variables for the new part pf the code:
   !
   !--- A vector to compute the K_Star() of entered space group.
   !
   real(kind=cp), dimension(3)   :: vec_k
   real(kind=cp), dimension(3)   :: vec_k_prime
   real(kind=cp), dimension(3)   :: vec_k_condition1
   real(kind=cp), dimension(3)   :: vec_k_condition2
   real(kind=cp), dimension(3)   :: vec_H_test

   real, dimension(:), allocatable :: m_values        ! m_values = (/0,1,2,3,4,6/) ! 5 is excluded
   real, dimension(:), allocatable :: n_values_Gkmink ! matrix indexes of only Gk-k
   real, dimension(:), allocatable :: n_values_Gk     ! matrix indexes of only Gk
   real, dimension(:), allocatable :: vec_T4_possible_values             ! possible values m/n
   real, dimension(:), allocatable :: vec_T4_possible_values_definitive  ! possible values m/n, trimmed
   real(kind=cp), dimension(4)     :: vec_H_R                 ! a buffer vec, there are many, one for operator
   real(kind=cp), dimension(4)     :: vec_T_hat, vec_T_hat_0, vec_T_hat_1, vec_T_hat_2 , PV  ! a buffer vec, is the %tr part plus a possible t value. PV=Prueba Vector
   real(kind=cp), dimension(3)     :: vec_buff0_0, vec_buff0_1, vec_buff0_2  ! a buffer vec, is the %tr part plus a possible t value
   real(kind=cp), dimension(7)     :: vec_buff0_3,  vec_buff0_11,  vec_buff0_22  ! store 3+4 final result
   real(kind=cp), dimension(5)     :: vec_buff0_5  ! store 1+4 purged final result

   real, dimension(:,:,:),   allocatable :: matrices_R_S   ! The {H(R)} que seran (4:4:NUMERO OPERADORES)
   real, dimension(:,:),     allocatable :: FINAL_RESULT
   real, dimension(:,:),     allocatable :: FINAL_RESULT5
   real, dimension(:,:),     allocatable :: FINAL_RESULT_PURGED
   !real, dimension(:,:),     allocatable :: FINAL_RESULT_PURGED_DEFINITIVE
   real, dimension(:,:),     allocatable :: vectors_T, vectors_T_hat_0, vectors_T_hat_1, vectors_T_hat_2     ! The translations
   real, dimension(:,:),     allocatable :: PP, PM  ! a test matrix for multiplication
   integer, dimension(:,:),  allocatable :: multiplitable  !
   real, dimension(128) :: movn
   integer              :: n_movn

   logical :: equal

   !--- Declaring variables for Space Group and Group k
   !
   Type (Space_Group_Type) :: SPGk_in, SPGk_out
   Type (Group_k_Type)     :: Gk_out







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




      !> Recast the SpaceGroup variable, to decouple this new part of the code
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

      !---> MOD(A,P) computes the remainder of the division of A by P.
      !---> vector v starts in 1: v(1), v(2), v(3)

      WRITE(*,*) " "
      WRITE(*,*) " ========================================================================"
      WRITE(*,*) " "
      WRITE(*,*) " STEP 3: find the H(R) vectors and R_S matrices "
      !WRITE(*,*) " now with the SPGk_out: STEP 3 find the H(R)"
      WRITE(*,*) " "

      ! Pointer to the extended little group:  ",Gk%p(1:2*Gk%ngk)
      !do j=1, size(Gk_out%p(1:2*Gk_out%ngk))

      ! allocate here al the {R_S} matrices and set all to zero. Also the 3*3 T vectors
      !se van a expandir las matrixec_R_S, hay que realocatar para expandir. Usamos un nmax gordo
    nmax=51
    !
    if(allocated(matrices_R_S)) deallocate(matrices_R_S)
       allocate(matrices_R_S(4,4,SPGk_out%multip+nmax)) ! 4,4,num operators; put a variabel
     matrices_R_S=0
     !write(*,*) "check-------------------", matrices_R_S

   if(allocated(vectors_T)) deallocate(vectors_T)
       allocate(vectors_T(4, SPGk_out%multip)) ! 3rd is the %tr, whwreas the 4th is the possible t value
     vectors_T=-10 ! an improbable initialization
     !write(*,*) "check-------------------", vectors_T
   if(allocated(vectors_T_hat_0)) deallocate(vectors_T_hat_0)
       allocate(vectors_T_hat_0(SPGk_out%multip, 4)) ! 3rd is the %tr, whwreas the 4th is the possible t value
     vectors_T_hat_0=-100.0 ! an improbable initialization
     !write(*,*) "check-------------------", vectors_T
   if(allocated(vectors_T_hat_1)) deallocate(vectors_T_hat_1)
       allocate(vectors_T_hat_1(SPGk_out%multip, 4)) ! 3rd is the %tr, whwreas the 4th is the possible t value
     vectors_T_hat_1=-101.0 ! an improbable initialization
     !write(*,*) "check-------------------", vectors_T
   if(allocated(vectors_T_hat_2)) deallocate(vectors_T_hat_2)
       allocate(vectors_T_hat_2(SPGk_out%multip, 4)) ! 3rd is the %tr, whwreas the 4th is the possible t value
     vectors_T_hat_2=-102.0 ! an improbable initialization
     !write(*,*) "check-------------------", vectors_T



     !call here the subrutine that makes the R_S matrices
     WRITE(*,*) "    number of  R_S matrices = ", SPGk_out%multip
     WRITE(*,*) "    making the R_S matrices:"
     call  Get_Extended_Group_Matrices(SPGk_out, vec_k, nmax, matrices_R_S)

     WRITE(*,*) " "
     WRITE(*,*) " ========================================================================"
     WRITE(*,*) " "
     WRITE(*,*) " STEP 4: find possible t (=m/n) values consistent with the SSG law"
     WRITE(*,*) " "


     nmatrices=GET_NMATRICES(matrices_R_S)
     !WRITE(*,*) " testing index"
     !WRITE(*,*) " test index-------------- ", INDEX_MATRIX  (matrices_R_S(:, :, 2),4)    ! 4=dim
     !WRITE(*,*) " test indexes------------ ", INDEX_MATRICES(matrices_R_S, 4, nmatrices) ! 4=dim;
     !WRITE(*,*) " R_S indexes:         ", INDEX_MATRICES(matrices_R_S, 4, nmatrices) ! 4=dim;


     if(allocated(m_values)) deallocate(m_values)
     allocate(m_values(6))  !m_values = (/0,1,2,3,4,6/) ! fixed values; 5 is excluded.
     if(allocated(n_values_Gkmink)) deallocate(n_values_Gkmink)
      allocate(n_values_Gkmink(nmatrices))        ! nmatrices=8 in example, size of Gk-k
     if(allocated(n_values_Gk)) deallocate(n_values_Gk)
     allocate(n_values_Gk(SPGk_out%multip))   ! is 4 only, size of Gk

     m_values=0


     m_values = (/0,1,2,3,4,6/) ! fixed values; 5 is excluded.
     n_values_Gkmink = INDEX_MATRICES(matrices_R_S, 4, nmatrices)
     n_values_Gk     = INDEX_MATRICES(matrices_R_S, 4, nmatrices)
     size_m=SIZE(m_values)
     size_mn=size_m*SIZE(n_values_Gk)
     WRITE(*,"(a,6i4)") " m values:            ", nint(m_values)        !," ----------", size_m ! is fixed
     WRITE(*,"(a,128i4)") " n indexes:           ", nint(n_values_Gkmink) !," ----------", size_Gkmink
     !WRITE(*,*) " n values of Gk:      ", n_values_Gk     !," ----------", size_Gk  ! beware!!! is 4 or 8??
     !WRITE(*,*) " possible m/n values:   "

     if(allocated(vec_T4_possible_values)) deallocate(vec_T4_possible_values)
     allocate(vec_T4_possible_values(size_mn))
     anImprobableNumber=0.0051
     vec_T4_possible_values=anImprobableNumber ! initializ



     vec_T4_possible_values=GET_POSSIBLE_T4_VALUES(n_values_Gk, size_mn)
     theTrueSize=true_size(vec_T4_possible_values)
     !WRITE(*,*) ""
     !WRITE(*,*) ""
     !WRITE(*,*) "vector initial ", vec_T4_possible_values
     !WRITE(*,*) "after trimming ", TRIMMER(vec_T4_possible_values,theTrueSize)
     !WRITE(*,*) " possible m/n values: "
     !WRITE(*,*) " possible m/n values---> ", vec_T4_possible_values, theTrueSize
    ! WRITE(*,*) " possible t=(m/n) values: ", TRIMMER(vec_T4_possible_values,theTrueSize)

     !write the result to a clean vector of size 'theTrueSize'
     if(allocated(vec_T4_possible_values_definitive)) deallocate(vec_T4_possible_values_definitive)
     allocate(vec_T4_possible_values_definitive(theTrueSize))
     vec_T4_possible_values_definitive=-1 ! just in case, to check if anything fails
     vec_T4_possible_values_definitive=TRIMMER(vec_T4_possible_values,theTrueSize)
     WRITE(*,*) " possible t values, deft: ", vec_T4_possible_values_definitive
     aux_string=" "
     l=0
     call get_m_over_n(m_values,n_values_Gkmink,movn,n_movn)
     do i=1,n_movn
       call Get_Fraction_1Dig(movn(i),aux_string(l+2:))
       l=len_trim(aux_string)
     end do
     write(*,"(a)") trim(aux_string)

     ! futura sol. elegante
     !WRITE(*,*) " n values of Gk-k", INDEX_MATRICES(matrices_R_S, 4, nmatrices, ext=1)
     !WRITE(*,*) " n values of Gk",   INDEX_MATRICES(matrices_R_S, 4, nmatrices, ext=0) ! we need this one

     !WRITE(*,*) " testing print matrices:"
     !WRITE(*,*) "", PRINT_MATRICES(matrices_R_S)



     WRITE(*,*) " "
     WRITE(*,*) " ========================================================================"
     WRITE(*,*) " "
     WRITE(*,*) " STEP 5: expansion"
     WRITE(*,*) " "
     ! here we play with the expansion-group function.
     ! she is able to say if matrices form group or not and if they do, expand and gives truth table
     ! let's use here the RS matrices as input
     !WRITE(*,*) "", PRINT_MATRICES(matrices_R_S)

     !WRITE(*,*) " starting expansion (of the 8 matrices_R_S)"
     !WRITE(*,*) "", PRINT_MATRICES(matrices_R_S)
     !WRITE(*,*) "", EXPAND_MATRICES(matrices_R_S) ! weird bug, memory corruptio, but works; move on


!       WRITE(*,*) " consistence SSG law"
!       WRITE(*,*) " multiplicidades", SPGk_out%multip
!       do i=1, SPGk_out%multip
!       !     matrices_R_S = SPGk_out%Symop(i)%rot   !ya la tenemos
!            vectors_T(:,i)= SPGk_out%Symop(i)%tr
!            !WRITE(*,*) " Vectors  TRANS:"  , SPGk_out%Symop(i)%tr   !just to check
!       end do
!       WRITE(*,*) " Vectors  TRANS:"  , vectors_T
!     !la 4a componente de vectors_T es:
!     !vec_T4_possible_values_definitive

       !comprobar partmatriz es grupal, con isGroup
       !comprobar partvectiz es grupal, con moduloLat


!       do i=1, size(vec_T4_possible_values_definitive)
!       do j=1, size(vec_T4_possible_values_definitive) !beware, size of matrices =! that of vectors
!       !WRITE(*,*) " Vectors yudayvcs", vec_T_hat , GET_NMATRICES(matrices_R_S) !size(matrices_R_S)
!       !PP=matmul(matrices_R_S(:,:,i), matrices_R_S(:,:,j))
!       !PT=matmul(matrices_R_S(:,:,i),vectors_T(:,i))+vectors_T(:,j) ! salvo lattice
!       end do
!       end do


!       nmatricesRS=GET_NMATRICES(matrices_R_S) !repeating, just in case
!       npossibleTs=size(vec_T4_possible_values_definitive)
!       do i=1, nmatricesRS!size(vec_T4_possible_values_definitive)
!        do j=1, nmatricesRS!size(vec_T4_possible_values_definitive) !beware, size of matrices =! that of vectors
!
!          PM=matmul(matrices_R_S(:,:,i), matrices_R_S(:,:,j))  !--> comprobar grupal
!          do k1=1,npossibleTs
!           do k2=1,npossibleTs
!             !vectors_T(k1,i)=vec_T4_possible_values_definitive(k1) filling 4th value in 2 steps!
!             !vectors_T(k2,j)=vec_T4_possible_values_definitive(k2)
!             buff1=vec_T4_possible_values_definitive(k1)
!             buff2=vec_T4_possible_values_definitive(k2)
!             vectors_T(k1,i)=buff1
!             vectors_T(k2,j)=buff2
!             WRITE(*,*) " Vec", vectors_T(k1,i), vectors_T(k2,j)
!             !PV=matmul(matrices_R_S(:,:,i), vectors_T(:,j)) + vectors_T(:,i) !pero vectors_T(k1,j) k1=valores t y k2
!           end do
!          end do
!        end do
!       end do


     !PLAN: we retrieve the table from the EXPAND_MATRICES function and we use it to find the PV's
     !WRITE(*,*) " starting expansion (of the 8 matrices_R_S)"
     !WRITE(*,*) "", PRINT_MATRICES(matrices_R_S)
     WRITE(*,*) "", EXPAND_MATRICES(matrices_R_S, vec_T4_possible_values_definitive)

     nmatricesRS=GET_NMATRICES(matrices_R_S)
     if(allocated(multiplitable)) deallocate(multiplitable)
     allocate(multiplitable(nmatricesRS,nmatricesRS)) ! n*n table of already determined size of group
     multiplitable=EXPAND_MATRICES(matrices_R_S, vec_T4_possible_values_definitive)


     WRITE(*,*) "   checking the group table"
     do i=1,nmatricesRS
        WRITE(*,*),(multiplitable(i,j),j=1,nmatricesRS)
     end do

   vec_T_hat_1=-12 !initialliz
   vec_T_hat_2=-12

   !check
   !WRITE(*,*) " Vectors  TRANS before:"  , vectors_T    !32 vectors, ok
   !WRITE(*,*) " Vectors  TRANS before:"  , vec_T_hat_1  !inits, ok
   !WRITE(*,*) " Vectors  TRANS before:"  , vec_T_hat_2
   !WRITE(*,*) " Vectors  TRANS before:"  , vec_T4_possible_values_definitive !ok
   !WRITE(*,*) " vec_bu",vec_buff0_1


   ! for some reason this has to be computed in an independent loop, otherwise it fails
   ! we are creating a single vectors_T_hat_0 with the nmatricesRS size, even though
   ! the translations are repeated and only 4 are independent (but later we a re going to prune)
   ! we set the fouth component as 0. It will be later replaced by the possible values of t=m/n
   WRITE(*,*) ""
   do i=1,nmatricesRS
      !WRITE(*,*) "traslations --- ",i, " ",SPGk_out%Symop(i)%tr  !ok, beware: multiplicity is 4
      vec_buff0_0=0
      vec_buff0_0(1:3) = SPGk_out%Symop(i)%tr(1:3)
      vectors_T_hat_0(i, 1:3) = SPGk_out%Symop(i)%tr
      vectors_T_hat_0(i, 4  ) = 0 !set 4th component to zero
   end do
   do i=1, nmatricesRS
      !WRITE(*,*) "Visualization: ",(vectors_T_hat_0(i,j),j=1,4)
   end do
   do i=1, nmatricesRS
      !vectors_T_hat_1(i)=vectors_T_hat_0(i)
   end do
   WRITE(*,*) ""





   ! to store results, just k index and vector possibles
    if(allocated(FINAL_RESULT)) deallocate(FINAL_RESULT)  !useless, rewrite
       allocate(FINAL_RESULT(4096, 7))!  4096=8*8*8*8 just number k (1*1) and the possible vectors+
    FINAL_RESULT=-11 ! 7 cause 1st is the ijk index and the rest is the possible vector PV

    if(allocated(FINAL_RESULT5)) deallocate(FINAL_RESULT5)
       allocate(FINAL_RESULT5(4096, 5))!  4096=8*8*8*8 just number k (1*1) and the possible vectors+
    FINAL_RESULT5=-131 ! 5 cause 1st is the ijk index and the rest is the possible vector PV


    possibleCombinations=0 !here yes but change name
    npossibleTs=size(vec_T4_possible_values_definitive)
    cnt1=0 !the position of the unrepeated results

iter=0!size of results
sizeResults=0

    do i=1, nmatricesRS ! for each 4*4 matrix forming group (rotational SSG)
     do j=1, nmatricesRS ! for each 4*4 matrix forming group (rotational SSG)

          ! now we create the two vectors using the T_hat_0 template
          vectors_T_hat_1(i,:)=vectors_T_hat_0(i,:)
          vectors_T_hat_2(j,:)=vectors_T_hat_0(j,:)

          k=multiplitable(i,j)   ! the index of the matrix product
          PM=matrices_R_S(:,:,k) ! R1*R2
          !WRITE(*,*) "Prueba Matrix (PM=R1*R2)   = ",i,j,k,"   ",PRINT_MATRIX(PM)

          !vec_buff0_3=0
          vec_buff0_5=0

          outer:do k1=1,npossibleTs
           inner:do k2=1,npossibleTs

             possibleCombinations = possibleCombinations + 1   !counter

             buff1=vec_T4_possible_values_definitive(k1)
             buff2=vec_T4_possible_values_definitive(k2)
             vectors_T_hat_1(i, 4)=buff1 ! filling the 4th component to make the SSG Translational part
             vectors_T_hat_2(j, 4)=buff2 !
             !WRITE(*,*) " T hat 1 full", vectors_T_hat_1
             !WRITE(*,*) " T hat 2 full", vectors_T_hat_2
             !WRITE(*,*) " T values", possibleCombinations, " ", buff1, buff2

             !>to check just this case
             !if(i==3 .and. j==6) then !--------------------------------------------
               PV=matmul(matrices_R_S(:,:,i), vectors_T_hat_2(j,:)) + vectors_T_hat_1(i,:) ! the Proof Vector
               !WRITE(*,*) " PV before Modulo-Lat", PV ! checked ok, so far so good
               !WRITE(*,*) " Proof Vector (PV=R1*T_hat2 + T_hat1) =", PV
               !WRITE(*,*) "hat2",vectors_T_hat_2(j,:)
               !WRITE(*,*) "hat1",vectors_T_hat_1(i,:)
               PV=Modulo_Lat(PV)
               !WRITE(*,*) " PV after Modulo-Lat", PV  ! some are vector zero; store the non-zero vectors with its table-matrix
               if( SUM(PV) /= 0.0 ) then
                  !note that the existence of zeroes in the 3 first positions only depends on case ij, and not on t4
                  !and, probably, no ij-case causes a full 3 zero vector. This condition on the code may be, then, redundant.
                  !WRITE(*,*) " PV after Modulo-Lat non zero                    ", PV
                  !final result is a (1+4)-vector:
           !       vec_buff0_5(1)=k    ! 1st component is the k-index from multiplication table
           !       vec_buff0_5(2:5)=PV ! the rest is the PV (Proof Vector)
                  !1WRITE(*,*) " PV after Modulo-Lat non zero bf", vec_buff0_5
                  !WRITE(*,*) " "

                  !check on-the-fly if vector already exists
                  !and store unique ones
                  !cnt1=0 !size of final results
                  !WRITE(*,*) " hello",FINAL_RESULT5(44,:)


                  !find a good way to implement a dynamic array
                  !by simply filling it (no need to know the final size)
                  !we would not need to fix the limit as 4096


seclimit=35 ! automate this,
superdoext:do m=1,seclimit!23

     iter=sizeResults
     !vec_buff0_5=FINAL_RESULT5(i,:)
     vec_buff0_5(1)=k    ! 1st component is the k-index from multiplication table
     vec_buff0_5(2:5)=PV ! the rest is the PV (Proof Vector)

     !WRITE(*,*) "Taking element-----------       ",i, " ", vec_buff0_5
     equal=.false.
     checker:do n=1,sizeResults+1!5!sizeResults
          !WRITE(*,*), "  checker current size is ", sizeResults
          !WRITE(*,*), "  comparing ",  vec_buff0_5
          !WRITE(*,*), "  comparing ",  FINAL_RESULT5(n,:)


          if(all( vec_buff0_5  .EQ. FINAL_RESULT5(n,:)  )) then
             !WRITE(*,*), "    ok equals, cycling "!, i , n!, "  ", FINAL_RESULT5(n,:)
             equal=.true.
             cycle superdoext !YESS!! keep looking
          else
             !WRITE(*,*), "    ko nequals"!, i , n!, "  ", FINAL_RESULT5(n,:)
             equal=.false.
             !sizeResults=sizeResults+1
             !FINAL_RESULT5(sizeResults,:)=vec_buff0_5
             !WRITE(*,*), "    now: FINAL_RESULT5 at pos", sizeResults
             !WRITE(*,*), "     is:                     ", vec_buff0_5
             !cycle checker  !keep looking
          end if
     end do checker

          if(equal) then
               !WRITE(*,*), "    ok equals, cycling "
               !NADA cycle superdoext
          else
             sizeResults=sizeResults+1
             FINAL_RESULT5(sizeResults,:)=vec_buff0_5
             !WRITE(*,*), "    now: FINAL_RESULT5 at pos", sizeResults, "added"
             !WRITE(*,*), "    instantaneous results are:"
             !do k1=1,sizeResults !more indexes!
              !   WRITE(*,*) "    instant res--> ", (FINAL_RESULT5(k1,k2),k2=1,5)
             !end do
             cycle superdoext
          end if

end do superdoext




               end if !non-zero vectors
             !end if ! check particular case i==3 .and. j==6  !--------------------------------------



           end do inner!k2
          end do outer!k1

       end do ! j=1, nmatricesRS
       end do ! i=1, nmatricesRS



WRITE(*,*) ""
WRITE(*,*) ""


WRITE(*,*) "The result is the set {R_s|t=(t1,t2,t3,t4)} where the 3*3 submatrix of R_s is R_g."
WRITE(*,*) "   A natural index 'n' checks the number of total results."
WRITE(*,*) "   The R_g matrix, defined as the matrix ID of the group."
WRITE(*,*) "   The 4 components of the 't' vector: t=(t1,t2,t3,t4)"
WRITE(*,*) "      NOTE: the n=35 limit is a forced declaration     (to be upgraded as a dynamic array)"
WRITE(*,*) "      NOTE: the -131.00000 is a forced initialization  (to be upgraded as a dynamic array)"
WRITE(*,*) ""
WRITE(*,*) "         |  n    |   matrix R_g ID |  t=(t1,t2,t3,t4)    "
 WRITE(*,*) ""
do k1=1,seclimit !53! 4096!cnt1
   WRITE(*,"(a,i4,a,5f10.3)") " ",k1, "    ",(FINAL_RESULT5(k1,k2),k2=1,5) !sort-ear
end do

!check, get tr vectors and matrices and get results
!WRITE(*,*) "   che", SPGk_out%Symop(3)%tr(1:3)
!WRITE(*,*) "   che", SPGk_out%Symop(4)%tr(1:3)


!check: es normal que de los vectores final result cambie la penultima celda?
!check: por que el caso de matriz 1 o 5 hay solo multiplicidad 3 y no 4?
!       es acaso por la forma de las matrices que son sin elementos fuera de la diagonal?
!possible bug?: entering Pnma with vector (alfa, 0,0) crashes for some values of alpha (with alpha=0.1 works)
!              in any other case it crashes (allocate) is the bug in this code or in the interaction
!              with other modules of the CrysFML library?





      WRITE(*,*) " "
      WRITE(*,*) " "
      WRITE(*,*) " ================== END. ======================"
      WRITE(*,*) " "





   end do ! of the whole program






!--- FUNCTIONS--------------------------------------------

contains

  subroutine get_m_over_n(m_vector,n_vector,movn,nind)
    real, dimension(:), intent(in) :: m_vector
    real, dimension(:), intent(in) :: n_vector
    real, dimension(:), intent(out) :: movn
    integer,            intent(out) :: nind

    !Local variables
    integer :: lm,ln,i,j,k
    real    :: coc
    lm=size(m_vector)
    ln=size(n_vector)
    nind=1
    movn(1)=0
    do i=2,lm
      do_j:do j=1,ln
        if( m_vector(i) > n_vector(j)) cycle
        coc= m_vector(i)/n_vector(j)
        if(Zbelong(coc)) cycle
        do k=nind,1,-1
          if(abs(coc-movn(k)) < 0.00001) then
            cycle do_j
          end if
        end do
        nind=nind+1
        movn(nind) = coc
      end do do_j
    end do
  end subroutine get_m_over_n

!----------------------------------------------------------------------------
! TO DO: built in functions to get nrows, ncols and nmatrices in a set (:,:,:)
!        we better define a new type: set of matrices. However this is practically done

!----------------------------------------------------------------------------
!--------------------------------------------------
function true_size(vector) result (trueSize)
! Returns the 'true size' of a vector. It is supposed that the vector was allocated
! with a maximum size, then initiallized to an improbable number
! (because it can be filled with the value '0')
! and we want now to trim it. For instance, in this case,
! the allocation was to 7 but the true size is 3:
!
! (0.0,  1.0,  4.0,  5.1E-03,  5.1E-03,  5.1E-03,  5.1E-03)

real, dimension(:), intent(in) ::  vector
integer :: trueSize

integer :: i, cnt

cnt=0

do i=1, size(vector)
  if (vector(i) /= anImprobableNumber) then
    cnt=cnt+1
  end if
end do
!WRITE(*,*) " val-------------------->>>>>>>>>>>", anImprobableNumber


trueSize=cnt
return
end function true_size
!--------------------------------------------------

!extra complexity, uses allocatable dummy arguments with a module.
function GET_POSSIBLE_T4_VALUES(vec_N_values, sizeMN)  result(vec_P_values)
!use arrayMod
real(kind=cp), dimension(:),  intent(in)    :: vec_N_values
integer,                      intent(in)    :: sizeMN
real(kind=cp), dimension(sizeMN)            :: vec_P_values ! possible values.
!if(allocated(vec_P_values)) then

real(kind=cp), dimension(6)         :: vec_M_values
real(kind=cp), dimension(sizeMN)    :: vec_P_values_buffer

integer :: i,j, sizeN, sizeM,  cnt, right_size
real :: val

vec_M_values = (/0,1,2,3,4,6/)

sizeN=SIZE(vec_N_values)
sizeM=SIZE(vec_M_values)

!allocate(vec_P_values(sizeN*sizeM))


vec_P_values_buffer=anImprobableNumber !initiallize to improbable number, bug was here

!seeding
cnt=1
vec_P_values_buffer(cnt)=0  ! starts with a zero, NOT FULL INITIALIZATED

val=anImprobableNumber

!WRITE(*,*) " valm--------------------", vec_P_values_buffer
!WRITE(*,*) " valm--------------------"


do i=1, sizeN
 do j=1, sizeM

  val = vec_M_values(j)/vec_N_values(i)
  !WRITE(*,*) " val-------------------", val

  if ( ANY( vec_P_values_buffer == val ) ) then  !remove duplicates
    cycle
  else
  cnt=cnt+1
  vec_P_values_buffer(cnt)=val
  end if

 end do
end do

!endif
!right_size=SIZE(vec_P_values_buffer )
right_size=cnt

vec_P_values = TRIMMER(vec_P_values_buffer, right_size)
return
end function GET_POSSIBLE_T4_VALUES

!----------------------------------------------------------------------------

function TRIMMER(vector, trueSize) result (vector_trimmed)
  real, dimension(:), intent(in)  ::  vector
  integer,            intent(in)  ::  trueSize
  real, dimension(trueSize) ::  vector_trimmed

  integer :: i, cnt

  !write(*,*), "trimming"
  do i=1, trueSize
  vector_trimmed(i)=vector(i)
  end do


return
end function TRIMMER

!----------------------------------------------------------------------------

function GET_ROWS_MATRIX(A) result(nrows)
real(kind=cp), dimension(:,:),  intent(in)    :: A  ! 2D
integer :: nrows
nrows=UBOUND(A,DIM=1)
return
end function GET_ROWS_MATRIX

!----------------------------------------------------------------------------

function GET_COLS_MATRIX(A) result(ncols)
real(kind=cp), dimension(:,:),  intent(in)    :: A ! 2D
integer :: ncols
ncols=UBOUND(A,DIM=2)
return
end function GET_COLS_MATRIX

!----------------------------------------------------------------------------

function GET_NMATRICES(A) result(nmatrices)
!returns the number of matrices of a set
!needs the global variable nmax

real(kind=cp), dimension(:,:,:),  intent(in)    :: A ! 3D
integer :: nmatrices
nmatrices=UBOUND(A,DIM=3)-nmax
return
end function GET_NMATRICES

!----------------------------------------------------------------------------

function IS_MATRIX_ZERO(A, dim) result(isMatrixZero)

use CFML_Math_General, only: Equal_Matrix

real(kind=cp), dimension(:,:),  intent(in)    :: A ! a matrix
integer,  intent(in)                            :: dim
logical :: isMatrixZero, isZero

integer :: cnt, i !unused
real(kind=cp), dimension(dim,dim) :: MatrixZero

!making the 0-Matrix
do i=1,dim
   MatrixZero(i,i)=0
end do

isMatrixZero=.false.


if(Equal_Matrix( A(:,:), MatrixZero(:,:), dim) ) then
   isMatrixZero=.true.
end if


return
end function IS_MATRIX_ZERO

!----------------------------------------------------------------------------

function GET_NMATRICES_NOTZERO(A, dim) result(nmatrices_out)
!returns the number of matrices of a set
!does not count matrices equal to zeroes

use CFML_Math_General, only: Equal_Matrix

real(kind=cp), dimension(:,:,:),  intent(in)    :: A ! a set of matrices
integer,  intent(in)                            :: dim
integer :: nmatrices_out

integer :: cnt, i, allocated_matrices, nmatrices_in
real(kind=cp), dimension(dim,dim)          :: MatrixZero


  nmatrices_in =0
  nmatrices_out=0

  !making the 0-Matrix
  do i=1,dim
     MatrixZero(i,i)=0
  end do

  !A=0 !prohibido por ser intention in tlw
  allocated_matrices=UBOUND(A,DIM=3)! number of initial matrices (allocated ones, filled or empty)
  nmatrices_in=allocated_matrices   ! just for clarifying the notation

  cnt=0
  do i=1, nmatrices_in
  if(.not. Equal_Matrix( A(:,:, i), MatrixZero(:,:), dim) ) then
  !check how many of the A matrices are different
  ! from an unusual matrix, (quasi)zero
     cnt=cnt+1
  end if
  end do
  nmatrices_out=cnt

  WRITE(*,*) " GET_NMATRICES_NOTZERO: num of input  (allocated) matrices------", nmatrices_in
  WRITE(*,*) " GET_NMATRICES_NOTZERO: num of output (not zero)  matrices------", nmatrices_out


return
end function GET_NMATRICES_NOTZERO

!----------------------------------------------------------------------------

function PRINT_MATRICES(matrices_A) result(printed_OK)
!---parameters
real(kind=cp), dimension(:,:,:),  intent(in)    :: matrices_A
!integer, optional,              intent(inout) :: order
logical :: printed_OK

!---local variables
integer :: i,j,m,n, nmatrices, nmatrices_total, nrows, ncols

printed_OK=.false. ! default
nmatrices_total=GET_NMATRICES(matrices_A)
nmatrices=nmatrices_total-nmax ! nmax is global, beware
nrows=GET_ROWS_MATRIX(matrices_A(:,:,1)) ! taking the first
ncols=GET_COLS_MATRIX(matrices_A(:,:,1)) ! taking the first

!WRITE(*,*) " printing nmatrices_total, nmatrices, nrows, ncols", nmatrices_total, " ",nmatrices, " ",nrows, " ", ncols

!WRITE(*,*) " printing a set of matrices:"
           !do j=1,nmatrices
           do j=1,nmatrices_total
              !WRITE(*,*) " "
              WRITE(*,*) " "
              do n=1,ncols
                 write(*, '(1000F14.7)')( matrices_A(m,n, j), m=1,nrows)
              end do
           end do

   !WRITE(*,*) " nmax, nmatrices_total,  nmatrices:", nmax,  nmatrices_total,  nmatrices
printed_OK=.true.
return
end function PRINT_MATRICES
!----------------------------------------------------------

function PRINT_MATRIX(A) result(printed_OK)
!function PRINT_MATRIX(A, order) result(printed_OK)

! if it is fed a set of matrices, it must compute the order and number, or deduce, or be put as optargs

real(kind=cp), dimension(:,:),  intent(in)    :: A
!integer, optional,              intent(inout) :: order
logical :: printed_OK
integer :: m, n, order   ! assumed square

if(GET_ROWS_MATRIX(A) /= GET_COLS_MATRIX(A)) then
   WRITE(*,*) "Error! matrix is not square. "
   WRITE(*,*) "  matrix info: ------------------------", UBOUND(A(:,:))
   WRITE(*,*) "  rows: -------------------------------", GET_ROWS_MATRIX(A) !4 it was bugged here!
   WRITE(*,*) "  cols: -------------------------------", GET_COLS_MATRIX(A) !3
end if


order=INT(SQRT(REAL(SIZE(A))))
!NOTE: in a 3*4 matrix, the dimension is 3 times 4
!NOTE: in a 4*4 matrix, the order is 4

!if( .not. present(order) )  order=INT(SQRT(REAL(SIZE(A))))
!WRITE(*,*) " order of matrix=", order
!! Beware: we are assuming squared matrices.
!! We should protect the code here in case a n*m matrix is inserted

printed_OK=.false.

!do j=1, total_matrices
  do n=1,order
     write(*, '(1000F14.7)')( A(m,n), m=1,order)
  end do
!end do

printed_OK=.true.

  !each column (DIM=1) or row (DIM=2),
  !WRITE(*,*) " mmmmm-------------------------------", UBOUND(A,DIM=1),UBOUND(A,DIM=2)


return
end function PRINT_MATRIX

!----------------------------------------------------------

function SUB_MATRIX(A)  result(subA)
!extract a submatrix from matrix A. Default is 3, but put options.

real(kind=cp), dimension(:,:),    intent(in)  :: A
real(kind=cp), dimension(3,3)                 :: subA

subA = A(1:3, 1:3)
WRITE(*,*) "", PRINT_MATRIX(subA)

return
end function SUB_MATRIX

!----------------------------------------------------------
function INDEX_MATRICES(matrices_A, dimOfMatrices, numOfMatrices)  result(vector_indexes)
implicit none
!---parameters
real(kind=cp), dimension(:,:,:),  intent(in)  :: matrices_A
integer(kind=cp),                 intent(in)  :: dimOfMatrices  ! dimension of input matrices (assuming all sme)
integer(kind=cp),                 intent(in)  :: numOfMatrices  ! fortran does not let do it in any other way!!
real(kind=cp), dimension(numOfMatrices)       :: vector_indexes


!---local
integer  :: l_bound, u_bound, i
!l_bound=lbound(matrices_A,dim=3)
!u_bound=ubound(matrices_A,dim=3)-nmax


!initiallization to zeroes
vector_indexes=0



  do i=1,numOfMatrices
     vector_indexes(i)=INDEX_MATRIX(matrices_A(:,:,i), dimOfMatrices)
     !WRITE(*,*) "mono", PRINT_MATRIX(matrices_A(:,:,i))
  end do





return
end function INDEX_MATRICES

!----------------------------------------------------------

function INDEX_MATRIX(A, dim)  result(index)
  !computes the index of a matrix, i.e. the number of times it has to be
  !multiplied by itself in order to reach the Identity matrix.
  !By default computes only a 3*3 submatrix (upper-left: the rotational part)
  !of an indetermined (:,:) matrix.

  ! if it is fed a set of matrices, it must store the indexes in a vector

  use CFML_Math_General, only: Equal_Matrix

  implicit none

  !WRITE(*,*) " -----------------------------------"


  real(kind=cp), dimension(:,:),    intent(in)  :: A
  integer(kind=cp),                 intent(in)  :: dim  ! dimension of input matrix (redundant, just to check)
  integer                                       :: index, counter, limit

  real(kind=cp), dimension(3,3)              :: subA ! a submatrix 3*3 (rotational part) of A
                                                     !    (assuming order of A not smaller than 3)
  real(kind=cp), dimension(dim,dim)          :: P
  real(kind=cp), dimension(dim,dim)          :: Identity

  !make an Identity matrix of dimension dim
  !if(allocated(A)) deallocate(A)
  !allocate(A(dim,dim))

  !--- creating Id matrix
  !
  Identity=0
  do i=1,dim
     identity(i,i)=1
  end do

  counter=0
  limit=1000
  index=0

  P=Identity ! sic  P=A


           do
             P(:,:) = matmul(A(:,:),P(:,:))
             counter=counter+1
             if(.not. Equal_Matrix( P(:,:), Identity(:,:), dim) ) then
                 !WRITE(*,*) " not found. Cycling. ", counter
                 if(counter .le. limit) then
                    cycle
                 else
                    WRITE(*,*) " not found below security limit. Exiting ", counter, "", limit
                    exit
                 end if
             else
                 index=counter
                 !WRITE(*,*) " yes found, index= ", index
                 exit
             end if

           end do



  return  ! beware: returns 0 if index not found
end function INDEX_MATRIX

!----------------------------------------------------------------------------

! the extended group matrices as a set. Better as a subroutine
!Function Get_Extended_Group_Matrices(SPGk_input, vec_k_input)  result(matrices_R_S_output)
Subroutine Get_Extended_Group_Matrices(SPGk_input, vec_k_input, nmax, matrices_R_S_output)

use CFML_Math_General, only: Equal_Matrix

!implicit none !only for functions


!---parameters
Type (Space_Group_Type)               :: SPGk_input
real(kind=cp), dimension(3)           :: vec_k_input
real, dimension(:,:,:),   allocatable :: matrices_R_S_output  ! The {H(R)} que seran (4:4:NUMERO OPERADORES)



!---local variables
integer                       :: nmax, i, j, m, n
real(kind=cp), dimension(3)   :: vec_k_prime
real(kind=cp), dimension(3)   :: vec_k_condition1
real(kind=cp), dimension(3)   :: vec_k_condition2
real(kind=cp), dimension(3)   :: vec_H_test

real(kind=cp), dimension(4)        :: vec_H_R  ! a buffer vec, therea aer many, one for operator, allocat
real, dimension(:,:),  allocatable :: vectors_T      ! The translations
real, dimension(:,:),  allocatable :: P ! a test matrix for multiplication



!nmax=100
if(allocated(matrices_R_S_output)) deallocate(matrices_R_S_output)
allocate(matrices_R_S_output(4,4,SPGk_input%multip+nmax)) ! 4,4,num operators; put a variable
matrices_R_S=0
if(allocated(vectors_T)) deallocate(vectors_T)
allocate(vectors_T(3, SPGk_out%multip)) ! 3,3,nomp
vectors_T=0



      do j=1, SPGk_out%multip


        !WRITE(*,*) "  " !   limite=",SPGk_out%nk
        vec_k_prime=matmul(vec_k_input, SPGk_input%SymOp(j)%Rot) ! ese SPGk_out
        !WRITE(*,*) " k_prime = ", vec_k_prime

        vec_k_condition1=vec_k_prime-vec_k ! k1=k'-k  <=>  k' = kR = +k + H(R)   ver si H es de red reciproca
        vec_k_condition2=vec_k_prime+vec_k ! k2=k'+k  <=>  k' = kR = -k + H(R)

        !para ver si es de la red reciproca usese k_Equiv(H,k) metiendo k=(000) fijo

        !WRITE(*,*) " k1 =      ", vec_k_condition1
        !WRITE(*,*) " k2 =      ", vec_k_condition2

        vec_H_test=(/ 0, 0, 0 /)
        !WRITE(*,*) "vec_H_test  =      ", vec_H_test

        if(   k_EQUIV(vec_k_condition1, vec_H_test, SPGk_input%SPG_lat)   ) then  ! si k1=H(R) es de la red reciproca
             !WRITE(*,*) "     found  cond 1 . expanding and storing   "
             vec_H_R(1:3)=vec_k_condition1
             vec_H_R(4:4)=1 ! for cond 1

             !join H(R(3*3)) and vec_H_R(1*4)  as a the j-esim matrices_R_S

             !WRITE(*,*) " R 3*3 :"
               do  n=1,3
                 !WRITE(*,"(a,3i4)") "      ",     SPGk_input%SymOp(j)%Rot(n,:) ! la matiz R
                 matrices_R_S_output(n,:,j)= SPGk_input%SymOp(j)%Rot(n,:)
               end do

            !WRITE(*,*) "     H(R)=", vec_H_R
            matrices_R_S_output(:,4,j) = vec_H_R  ! add thevector

            !WRITE(*,*) " next "

            !show R_S j-esi matrix
            !WRITE(*,*) " R_S 4*4 no ordered    ",  matrices_R_S_output(:,:,j)
            !WRITE(*,*) " R_S 4*4 ordered:"
            !do n=1,4
            !   write(*, '(1000F14.7)')( matrices_R_S_output(m,n, j) ,m=1,4)
            !end do



       else
             !WRITE(*,*) "     found cond 2  . (si no es la 1 es la 2) expanding and storing   "
             vec_H_R(1:3)=vec_k_condition2
             vec_H_R(4)=-1 ! for cond 2

             !join H(R(3*3)) and vec_H_R(1*4)  as a the j-esim matrices_R_S

             !WRITE(*,*) " R 3*3 :"
               do  n=1,3
             !    WRITE(*,"(a,3i4)") "      ",     SPGk_input%SymOp(j)%Rot(n,:) ! la matiz R
                 matrices_R_S_output(n,:,j)= SPGk_input%SymOp(j)%Rot(n,:)
               end do

            !WRITE(*,*) "     H(R)=", vec_H_R
            matrices_R_S_output(:,4,j) = vec_H_R  ! add thevector

            !WRITE(*,*) " next "

            !show R_S j-esi matrix
            !show R_S j-esi matrix
            !WRITE(*,*) " R_S 4*4 no ordenada    ",  matrices_R_S(:,:,j)
            !WRITE(*,*) " R_S 4*4 ordered:"
            !do n=1,4
            !   write(*, '(1000F14.7)')( matrices_R_S_output(m,n, j) ,m=1,4)
            !end do

       end if


       ! hacer lo del vector t fases
       ! hacer lo del vector t fases    SPGk_input%SymOp(j)%Rot(n,:) pero con Trans
       !
       !WRITE(*,*) " multiplicidad", SPGk_out%multip
       !do i=1, SPGk_out%multip
       ! !   !matrices_R_S = SPGk_input%Symop(i)%rot
       !     vectors_T(:,i)= SPGk_input%Symop(i)%tr
       !     WRITE(*,*) " Vectores de T (individuales):"  , SPGk_input%Symop(i)%tr
       !end do
       !WRITE(*,*) " Vectores de T (todos):"  , vectors_T    ! PUT IN ORDER, repeated!



        ! make the SSG type (tupe superspace group)

      end do ! del j


      !WRITE(*,*) " MATRICES exported, out of the loop where they were created (exported so t osys)"
      !     do j=1,SPGk_input%multip
      !        WRITE(*,*) ""
      !        do n=1,4
      !           write(*, '(1000F14.7)')( matrices_R_S_output(m,n, j) ,m=1,4)
      !        end do
      !     end do
      !
      !WRITE(*,*) " printing matrices:"
      WRITE(*,*) " ", PRINT_MATRICES(matrices_R_S_output)
      !WRITE(*,*) " printing matrices done"
return
!end function Get_Extended_Group_Matrices
end subroutine Get_Extended_Group_Matrices

!----------------------------------------------------------------------------

!function EXPAND_MATRICES(matrices_A, vec_T4_possible_values_definitively) result(isGroup)
function EXPAND_MATRICES(matrices_A, vec_T4_possible_values_definitively) result(table)


use CFML_Math_General, only: Determinant, Equal_Matrix
implicit none

!---parameters
real(kind=cp), dimension(:,:,:),   intent(in)    :: matrices_A
real(kind=cp), dimension(4),       intent(in)    :: vec_T4_possible_values_definitively
integer,       dimension(:,:),   allocatable     ::  table
logical :: isGroup

!---local variables
   integer, dimension(:,:),   allocatable  :: R
   integer, dimension(:,:),   allocatable  :: A,P,identity,minus_identity
   integer, dimension(:,:,:), allocatable  :: G, SS, SG
   integer :: ngroup, nmaxx , i, j, k, kk1, kk2, m, n, nrows, nrowsread, nrowsread2, ppp, q, iter, iter2, sizegg, &
              trials, idet, min, max, ssIndex, sizess, sizesg, zcnt, nonzeros, cnt, cntr, newelem, &
              sizeggg, sizesss, binSizeSS, binSizeGG, binSizeSG, dimread, nmatrices_before, nmatrices_after
   real    :: det, rnd
   logical :: eq
   character(len=80) :: fmtt, message
   integer, dimension(:),     allocatable :: ssetIndexes ! array to store the indexes of extracted subset
   integer, dimension(:),     allocatable :: ssetIndexesShuffled ! indexes of the subset after shuffling

   integer, dimension(:),     allocatable :: indexesGG ! array to store the indexes of the Global Group
   integer, dimension(:),     allocatable :: indexesSS ! array to store the indexes of the starting SubSet
   integer, dimension(:),     allocatable :: indexesSG ! array to store the indexes of the local SubGroup

   integer, dimension(:),     allocatable :: indexesGGminusSG ! remaining indexes of the logical subtraction SG-GG
   integer, dimension(:),     allocatable :: indexesGGminusSGnoZeros ! previous, having removed zeroes

   integer, dimension(:),     allocatable :: arraytest0,arraytest1,arraytest2,arraytest3,arraytest4,arraytest5
   integer, dimension(:),     allocatable :: binGG, binSZ, binSI, binSS, binSG, binCT, binPR


   !IO features
   real, dimension(:,:),   allocatable :: readtable   ! to read the exported-to-file table
   real, dimension(:,:,:), allocatable :: readmatrices   ! to read the exported-to-file matrices

   real, dimension(4) :: vec_TRANS_buff, vec_buff_t1, vec_buff_t2


!nmax: is a security limit, number of matrices in the expanding group.
!      it is set to a fixed, high number
!   n: is the size of the matrix. It is fixed to 4.
!nmaxx=100 !is global!!

isGroup=.false. ! by default, input matrices do not form a group
n=GET_ROWS_MATRIX(matrices_A(:,:,1)) ! taking the first element (must be the dimension: n=4)

nmatrices_before=0 !initiallization to avoid cycle-persitence
nmatrices_after=0

      !--- Allocating. Allocate  G as a set of nmax matrices of n*n
      !
      !if(allocated(matrices_A)) deallocate(matrices_A) ! set of generators (all the input A matrices)
      !
      if(allocated(G)) deallocate(G) ! set of generators (all the input A matrices)
      if(allocated(SS)) deallocate(SS) ! a subset
      if(allocated(SG)) deallocate(SG) ! subgroup
      if(allocated(R)) deallocate(R) ! real
      if(allocated(P)) deallocate(P) ! prueba
      if(allocated(identity)) deallocate(identity)
      if(allocated(minus_identity)) deallocate(minus_identity)
      allocate(G(n,n,nmax),SG(n,n,nmax),R(n,n), P(n,n),identity(n,n),minus_identity(n,n))
      !allocate(matrices_A(n,n,nmax))
      !
      !--- A is a proof to check with the identity
      if(allocated(A)) deallocate(A)
      allocate(A(n,n))
      !

      !--- creating I and -I
      !
      identity=0
      do i=1,n
         identity(i,i)=1
      end do
      minus_identity=-identity

   !--- Now starts the expansion of the generator set (G) to get a group ()
   !
   ! sizegg is the instantaneous size of the expanding set (a better name would be sizeEXPSET)
   ! it is initialized to 3 which is the initial number of generators
   !G=matrices_A
   !G(:,:,2)=0                    ! trick to see if group computes. Beware: there is no zero-purger and size may increase
   !G(:,:,1)=matrices_A(:,:,1)=0  ! trick to see if group computes.
   !

   !preunbug, trick n befor in order to make small matrices
   !G=0
   !G(:,:,1) = reshape( (/ 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4 /), (/ 4, 4/))
   !G(:,:,0) = reshape( (/ 1,0,0, 0,1,0,  0,0,1 /),    (/ 3, 3/))
   !G(:,:,2) = reshape( (/ 1,0,0, 0,1,0,  0,0,-1 /),   (/ 3, 3/))
   !G(:,:,3) = reshape( (/ -1,0,0, 0,-1,0,  0,0,-1 /), (/ 3, 3/))

   !WRITE(*,*) " priny dummy", PRINT_MATRICES(real(G)) !doen not work inside another function


   ! G is the set to be expanded. We set manually to some input matrices in order to test the code
   !G=0 ! always initiallize
   !G(:,:,1)=matrices_A(:,:,1) ! beware, answer depends on initial sset
   !G(:,:,2)=matrices_A(:,:,2)
   !G(:,:,3)=matrices_A(:,:,3)
   !G(:,:,4)=matrices_A(:,:,4)
   !G(:,:,5)=matrices_A(:,:,5)
   !G(:,:,6)=matrices_A(:,:,6)
   !G(:,:,7)=matrices_A(:,:,7)
   !G(:,:,8)=matrices_A(:,:,8)
   !!G(:,:,2)=0

   ! function called from anothere function does not work, why?
   ! it returns the T7F value but not the text messages
   !WRITE(*,*) " printing G -------------------", PRINT_MATRICES(real(G))


   !WRITE(*,*) " "
   !WRITE(*,*) " test ismatrixZero---------------------------------====="
   !dim_buff=GET_ROWS_MATRIX(real(G(:,:,1)))
   !WRITE(*,*) " ", IS_MATRIX_ZERO(real(G(:,:,2)), dim_buff)
   !WRITE(*,*) " test ismatrixZero---------------------------------====="
   !WRITE(*,*) " "
   !WRITE(*,*) " "
   !WRITE(*,*) " in the beginning, matrices are: ", G !PRINT_MATRICES(real(G))
   !WRITE(*,*) " "

   WRITE(*,*) " "
   WRITE(*,*) " "
   !WRITE(*,*) " "

   !!uncomment this when everything is debugged
   G=0 ! always initiallize
   !!copy matrices_A to G's
   dim_buff=GET_ROWS_MATRIX(real(G(:,:,1)))!
   !!WRITE(*,*) " print======", dim_buff, nmax !nmax=51
   do i=1,nmax ! AUTOMATE THIS LIMIT!!
      if(.not. IS_MATRIX_ZERO(matrices_A(:,:,i), dim_buff)) then! if matrix is not zero
         G(:,:,i)=matrices_A(:,:,i)
      end if
   end do






   !WRITE(*,*) " print===================================== all Gs", PRINT_MATRICES(real(G)) !ok  but ... prints just number, not matrices. Fix the sizes vs true sizes


   ! SELECT HERE tricjek ones, G's or A's
   !dim_buff=GET_ROWS_MATRIX(real(matrices_A(:,:,1)))!
   !nmatrices_before=GET_NMATRICES_NOTZERO(matrices_A, dim_buff)
   !!!
   dim_buff=GET_ROWS_MATRIX(real(G(:,:,1)))! for instance
   nmatrices_before=GET_NMATRICES_NOTZERO(real(G), dim_buff) ! stable, bug fixed
   WRITE(*,*) " print====== cycle-unstable--------------", dim_buff, nmatrices_before, UBOUND(G,DIM=3)
   !WRITE(*,*) " print====== cycle-unstable?", nmatrices_before

   !sizegg=GET_NMATRICES_NOTZERO(real(G), dim_buff)
   sizegg=nmatrices_before
   WRITE(*,*) ""
   WRITE(*,*) " sizegg init -------------------", sizegg
   WRITE(*,*) " nmax   init -------------------", nmax
   !WRITE(*,*) "if(sizegg >= nmax) then procedure fails!"
   WRITE(*,*) ""

   iter=0
   fmtt=" "
   write(unit=fmtt,fmt="(a,i5,a)") "(a,i3,a,i3,a,",n*n,"i3,a,i4)"
   do_ext:do   ! an infinite loop here to use the iter (as variable index)
     iter=sizegg
     do_i:do i=1,sizegg
       do_j: do j=1,sizegg ! multiply the matrices of the G set each other
           P(:,:)= matmul(G(:,:,i),G(:,:,j))! assign
           eq=.false. ! set default guess: they are not equal
           do k=1,sizegg ! comprobar si esta en el grupo generador
                if(Equal_Matrix (P(:,:), G(:,:,k), n)) then
                   eq=.true. ! if repeated it does not expand the set, keep multiplying
                   exit ! exist current loop of k-index
                     ! for the next program we need a kind of isGroup function
                !else ! if product of G(i)*G(j) is not in the set means is not yet a group:
                   !write(unit=*,fmt="(a)") "   current set is not a group (yet) "
                endif
           end do
           if(eq) cycle do_j
           !else, include it in the generator set
           !
           sizegg=sizegg+1
           if(sizegg > nmax) then
                write(unit=*,fmt="(a,i3,a,1000i3)") " =>  Procedure fails! ... expansion above expected limit"
                exit do_ext
           endif
           G(:,:,sizegg)=P  ! include it to the set (identified as i*j) and display its current size
           write(unit=*,fmt=fmtt) " =>  Including the matrix: G(",i,")xG(",j,")=", P(:,:), "   Current size of Set=",sizegg
           cycle do_ext
       end do do_j
     end do do_i
     if (iter == sizegg) then
         !table here? no
         exit
     endif
   end do do_ext !main loop Procedure

   !if(sizegg >= nmax) then
   !    write(unit=*,fmt="(a,i3,a,1000i3)") " =>  Procedure fails! ... expansion above maximum order of the exptd group"
       !cycle
   !end if
   !
   !write(unit=*,fmt="(a)") "    "
   !
   !do i=1,sizegg
      !write(unit=*,fmt="(a,i3,a,1000i3)") " =>  Matrix(", i,")=", G(:,:,i)
      !write(*,*) " =>  Matrix(", i,")=", G(:,:,i)
   !end do
   !
   !write(*,*) " final sizegg-------------------------------------- ", sizegg
   !write(unit=*,fmt="(a)") "  "
   !write(unit=*,fmt="(a)") "  "

   if(sizegg > nmax) then
     write(*,*) "   no table for this case (expansion exceeds sec. limit)"
   else

   !--- Multiplication table allocation
   !
   write(unit=*,fmt="(a)") " Multiplication table: "
   !
   if(allocated(table)) deallocate(table)
   allocate(table(sizegg,sizegg)) ! n*n table of already determined size of group

   !Multiplication table fill
   table=0
   do i=1,sizegg
     do j=1,sizegg
       P=matmul(G(:,:,i),G(:,:,j))
       eq=.false.
       do k=1,sizegg
          if(Equal_Matrix (P(:,:), G(:,:,k), n)) then
            eq=.true.
            exit
          end if
       end do
       if(.not. eq) then
         write(unit=*,fmt="(2(a,i3),a)") " => Warning! Matrices ",i," and ",j," give a result not included in the group!" ! to debug previous code, just in case
       else
         Table(i,j)=k  ! the shadow of a thought
       end if
     end do
   end do
   fmtt=" "
   write(unit=fmtt,fmt="(a,i5,a)") "(",sizegg+1,"i4)"
   write(unit=*,fmt=fmtt) 0,(j,j=1,sizegg)
   do i=1,sizegg
      write(unit=*,fmt=fmtt) i,(table(i,j),j=1,sizegg)
   end do

   endif
   !==================================================================================
   ! END OF MATRIX EXPANSION AND MULTIPLICATION TABLE
   ! initial group of G matrices expanded to a group

!TODO: purge in case of a 0-matrix (trick used to test the code, but increases the size artificially)
!      a kind of trimmer of zeroes should be implemented

WRITE(*,*) ""
dim_buff= size(G(:,:,1),dim=1) !GET_ROWS_MATRIX(real(G(:,:,1)))!
nmatrices_after=GET_NMATRICES_NOTZERO(real(G), dim_buff) ! sic


WRITE(*,*) ""
WRITE(*,*) "summary----------------------------------------------------------------------------"
WRITE(*,*) "   nmatrices_before=", nmatrices_before
WRITE(*,*) "   nmatrices_after =", nmatrices_after

if (nmatrices_before == nmatrices_after) then
   isGroup=.true. !input matrices formed group
   WRITE(*,*) "   input matrices formed group", nmatrices_before," " ,nmatrices_after
else if (nmatrices_before < nmatrices_after) then
   isGroup=.false. !input matrices do not fomed group and were expanded
   WRITE(*,*) "   input matrices did not form group; were expanded from:", nmatrices_before," to:" ,nmatrices_after
else if(nmatrices_before > nmatrices_after) then
   isGroup=.false. !input matrices do not fomed group and were expanded
   WRITE(*,*) "   ERROR: you have empty matrices artif. increasing the size! SOL: prune", nmatrices_before," " ,nmatrices_after
else
   WRITE(*,*) "   no more options!"
end if
WRITE(*,*) "------------------------------------------------------------------------------------"
WRITE(*,*) ""



! here we can use the table in order to find the t values that are consistent with SSG law
!WRITE(*,*) "   the table"
!   do i=1,sizegg
!      write(unit=*,fmt=fmtt) ,(table(i,j),j=1,sizegg)
!   end do

!computar todos los vectores, guardarlos, luego, discriminarlos
!vec_buff_t2= (/1,2,3,4/)!
!vec_buff_t1= (/2,2,4,4/)! just export the table
!do i=1, sizegg ! just run among all matrices what we found as group
!  do kk1=1, size(vec_T4_possible_values_definitively)
!  do kk2=1, size(vec_T4_possible_values_definitively)
!     vec_TRANS_buff=(matmul(G(:,:,i),vec_buff_t2) + vec_buff_t1) ! add 4th component, and discriminate it
!                                                                 ! according to lattice module
!     WRITE(*,*) "   thething---- ",vec_TRANS_buff
!  end do
!  end do
!end do





nmatrices_before=0 !fin-initiallization to avoid cycle-persitence
nmatrices_after=0


return
end function EXPAND_MATRICES


!-----------------------------
End Program SSPG_Info
