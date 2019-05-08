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
   !!---- INIT_SPACEG
   !!----
   !!----    Initializes the components that are not allocatable arrays
   !!----
   !!---- 19/04/2019
   !!
   Module Subroutine Init_SpaceG(Grp)
      !---- Arguments ----!
      class(Group_type),  intent(in out) :: Grp
      
      !> Init
      Grp%multip=0
      Grp%     d=0
      
      Select Type(Grp)
         type is (Spg_Type)
            Grp%numspg      = 0
            Grp%numshu      = 0
            Grp%Numops      = 0
            Grp%centred     = 1
            Grp%anticentred = 1
            Grp%mag_type    = 0
            Grp%num_lat     = 0
            Grp%num_alat    = 0
            Grp%spg_lat     = " "
            Grp%shu_lat(1)  = " "
            Grp%shu_lat(2)  = " "
            Grp%spg_symb    = "    "
            Grp%shu_symb    = "    "
            Grp%crystalsys  = "    "
            Grp%pg          = "    "
            Grp%laue        = "    "
            Grp%mat2std     = "    "
            Grp%mat2std_shu = "    "
            Grp%generators_list = "    "
      End Select
      
   End Subroutine Init_SpaceG
   
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