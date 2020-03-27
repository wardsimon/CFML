!!----
!!---- Program: Super Space Group Info (SSPG_INFO)
!!----          Made after the Space Group Info program in  CFML.
!!----          Restricted to Magnetic SSGs in particular.
!!---- Author:
!!---- Revision:
!!



!>>-------------------------------------------------
!>> MAIN
!>>-------------------------------------------------

Program SSPG_Info

  !---- Use Modules ----!

  use CFML_Math_General,              only: Determinant, Equal_Matrix, &
                                            Equal_Vector, Modulo_Lat, Zbelong
  use CFML_Crystallographic_Symmetry, only: Space_Group_Type, set_SpaceGroup, &
                                            Write_SpaceGroup, Symmetry_Symbol, &
                                            Get_Generators_From_SpGSymbol,     &
                                            err_symm, err_symm_mess
  use CFML_GlobalDeps,                only: cp
  use CFML_Propagation_Vectors,       only: Group_k_Type, Set_Gk, K_Star, &
                                            Write_Group_k, k_EQUIV

  !use CFML_SuperSpaceGroups
  use CFML_SuperSpaceGroups_LOCAL ! module name, not file

  use CFML_Rational_Arithmetic_test
  use Matrix_Mod


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

   type(SuperSpaceGroup_Type), dimension(8)  :: SSGs_out_test !only case nss=8
   integer                                   :: lun
   !character(len=80),    dimension(6,6)                  :: cMat1
   character(len=80),    dimension(5,5)                  :: cMat1
   character(len=80),    dimension(:,:),allocatable      :: cMat2


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
!
!--- CASE nk=1
!write(*,*) "testing with 1 k-vec"
!nkvecs=1
!if(allocated(kvecs)) deallocate(kvecs)
!allocate(kvecs(3,nkvecs))
!kvecs=0 ! init
!kvecs(:,1)=vec_k


!--- CASE nk>1
write(*,*) "testing with 2 k-vecs"
nkvecs=2
if(allocated(kvecs)) deallocate(kvecs)
allocate(kvecs(3,nkvecs))
kvecs=0 ! init
kvecs(:,1)=vec_k
kvecs(:,2)=vec_k + [0.201, 0., 0.]


!---++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
call Set_SSGs_from_Gkk_MYTEST(SPGk_in, nkvecs, kvecs, SSGs_out, nSSGs)
!
!---++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!!---++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!
WRITE(*,*) "nk=",nkvecs
write(*,*) "final results:"
write(*,*) "num of SSGSs", nSSGs
do i=1, nSSGs
 call Write_SSG(SSGs_out(i),full=.true.)
enddo
!!
!!---++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


end do ! prompting loop

write(*,*) ""
write(*,*) "end testing with 2 k-vecs"
write(*,*) ""

!!----------------------------------------------------------------------------

End Program SSPG_Info
