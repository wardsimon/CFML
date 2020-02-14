SubModule (CFML_gSpaceGroups) Init_Proc
   Contains

   !!----
   !!---- SET_IDENTITY_MATRIX
   !!----
   !!---- 19/04/2019
   !!
   Module Subroutine Set_Identity_Matrix(D)
      !---- Arguments ----!
      integer, intent(in) :: D    ! Dimension of the Matrix

      !---- Local Variables ----!
      integer :: i

      !> Init
      if (allocated(identity_matrix)) deallocate(identity_matrix)
      allocate(identity_matrix(d,d))

      call Rational_Identity_Matrix(identity_matrix)

   End Subroutine Set_Identity_Matrix

   !!----
   !!---- Init_SpaceGroup
   !!----
   !!----    Initializes the components that are not allocatable arrays
   !!----
   !!---- 19/04/2019
   !!
   Module Subroutine Init_SpaceGroup(Grp)
      !---- Arguments ----!
      class(Group_type),  intent(in out) :: Grp

      !> Init
      Grp%multip      = 0
      Grp%d           = 0

      Select Type(mGrp => Grp)

         type is (Spg_Type)
            mGrp%standard_setting=.false.
            mGrp%numspg      = 0
            mGrp%numshu      = 0
            mGrp%Numops      = 0
            mGrp%centred     = 1
            mGrp%anticentred = 0
            mGrp%mag_type    = 0
            mGrp%num_lat     = 0
            mGrp%num_alat    = 0
            mGrp%Parent_num  = 0       ! Number of the parent Group
            mGrp%spg_lat     = " "
            mGrp%shu_lat(1)  = " "
            mGrp%shu_lat(2)  = " "
            mGrp%Parent_spg  = "       "
            mGrp%tfrom_parent= "       "
            mGrp%setting     = "       "
            mGrp%spg_symb    = "       "
            mGrp%BNS_symb    = "       "
            mGrp%BNS_num     = "       "
            mGrp%OG_symb     = "       "
            mGrp%OG_num      = "       "
            mGrp%Centre      = "       " ! Alphanumeric information about the center of symmetry
            mGrp%crystalsys  = "       "
            mGrp%pg          = "       " !Point Group
            mGrp%mag_pg      = "       " !Magnetic Point Group
            mGrp%laue        = "       " !laue Group
            mGrp%mat2std     = "       "
            mGrp%mat2std_shu = "       "
            mGrp%generators_list = "      "
            mGrp%Hall       ="          "
            mGrp%SSG_symb   ="          "
            mGrp%SSG_Bravais="          "
            mGrp%SSG_nlabel ="          "
            mGrp%Bravais_num=0      ! Number of the Bravais class

         type is (SuperSpaceGroup_Type)
            mGrp%nk=0               !  number of k-vectors
            mGrp%nq=0               !  number of q-coefficients

         class default
            Err_CFML%Ierr=1
            Err_CFML%Msg="The class passed to 'Init_SpaceGroup' is not implemented"
      End Select

   End Subroutine Init_SpaceGroup

   !!----
   !!---- SET_CONDITIONS_GROUP
   !!----
   !!---- 19/04/2019
   !!
   Module Subroutine Set_Conditions_NumOP_EPS(maxop,epsg)
      !---- Arguments ----!
      integer,       optional, intent(in) :: maxop
      real(kind=cp), optional, intent(in) :: epsg

      !> Set a  maximum number of operators
      if (present(maxop)) maxnum_op=maxop

      !> Set a new EPS
      if (present(epsg)) then
         call Set_Eps_Math(epsg)
      else
         call Set_Eps_Math(0.001_cp)
      end if

   End Subroutine Set_Conditions_NumOp_EPS

End SubModule Init_Proc