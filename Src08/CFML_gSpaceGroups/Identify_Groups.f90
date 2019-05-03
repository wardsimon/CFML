SubModule (CFML_SpaceG) SPG_045
   Contains
   
   !!---- 
   !!---- Identify_Cryst_SpG
   !!----
   !!---- For a given crystallographic group G in an arbitrary
   !!---- setting, identifies the space group and computes the
   !!---- transformation matrix to the standard setting.
   !!----
   !!---- It follows Acta Cryst. A55 383-395 (1999). P,M,A and
   !!---- C matrices in the paper correspond to Pinv, Minv,
   !!---- Ainv and Cinv in this subroutine.
   !!----
   !!---- 22/04/2019 
   !!
   Module Subroutine Identify_Cryst_SpG(G)
      !---- Arguments ----!
      type(spg_type),    intent(in out) :: G

      !---- Local variables ---!
      integer                          :: n
      character                        :: lattyp
      type(rational), dimension(3,3)   :: P,Mp,Mc,M
      type(rational), dimension(3,3,6) :: A
      logical :: pout

      !>===== DEBUG =====
      pout=.false.
      pout=(pout .or. CFML_DEBUG)
      !>=================

      P=Get_P_Matrix(G,nospin=.true.)
      if (Err_CFML%Ierr /= 0) return

      Mp=Get_Mp_Matrix(G,P)
      if (Err_CFML%Ierr /= 0) return

      Mc=Get_Mc_Matrix(G%laue,Mp)
      if (Err_CFML%Ierr /= 0) return

      M = matmul(Mp,Mc)
      lattyp=Get_Lattice_Type(M)
      if (Err_CFML%Ierr /= 0) return
      
      call Get_A_Matrix_Crys(G%laue,A,n)
      !call Match_Crystallographic_Space_Group(G,P,M,A(:,:,1:n),n)

   end subroutine Identify_Cryst_SpG
   
   !!---- 
   !!---- Identify_Group
   !!----
   !!---- Initialize the identification of the group by calling
   !!---- the appropiate subroutine according to the nature  of
   !!---- the group -crystallographic, magnetic, superspace-.
   !!----
   !!---- 22/04/2019 
   !!
   Module Subroutine Identify_Group(G)
      !---- Arguments ----!
      type(spg_type),    intent(in out) :: G

      !---- Local Variables ----!
      logical :: pout

      !>===== DEBUG =====
      pout=.false.
      pout=(pout .or. CFML_DEBUG)
      !>=================

      !> Init
      call Clear_Error()

      if (G%d < 4) then
         Err_CFML%Ierr = 1
         Err_CFML%Msg = "Identify_Group@SPACEG: Group dimension < 4."
         return

      else if (G%d == 4) then ! Shubnikov groups
         call Identify_Shubnikov_Group(G)

      else ! Superspace groups
         Err_CFML%Ierr = 1
         Err_CFML%Msg = "Identify_Group@SPACEG: Superspace groups not implemented yet"
         return
      end if
      
   End Subroutine Identify_Group
     
   !!---- 
   !!---- Identify_Laue_Class
   !!----
   !!---- Sets the Laue class from the crystallographic point group
   !!----
   !!---- 22/04/2019 
   !!
   Module Subroutine Identify_Laue_Class(G)
      !---- Arguments ----!
      type(spg_type), intent(inout) :: G

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
            Err_CFML%Ierr = 1
            Err_CFML%Msg ="Identify_Laue_Class@SPACEG: Inconsistent crystallographic point group."
      end select
   End Subroutine Identify_Laue_Class
     
   !!----
   !!---- IDENTIFY_MAGNETIC_SPACE_GROUP
   !!----
   !!---- Identifies the Shubnikov group of group G. Before
   !!---- using this subroutine, the subroutine Identify_
   !!---- _Crystallographic_Space_Group must be called
   !!----
   !!---- 22/04/2019 
   !!----
   Module Subroutine Identify_Shubnikov_Group(G)
      !---- Arguments ----!
      type(spg_type),    intent(in out) :: G

      !---- Local variables ---!
      character                        :: lattyp
      type(rational), dimension(3,3)   :: P,Mp,Mc,M
      logical :: pout

      !>===== DEBUG =====
      pout=.false.
      pout=(pout .or. CFML_DEBUG)
      !>=================

      if (.not. Magnetic_DBase_allocated) then
         call Read_Magnetic_Data()
         if (Err_CFML%IErr /=0) return
      end if

      call Identify_Crystallographic_PG(G)
      if (Err_CFML%Ierr /= 0) return

      call Identify_Laue_Class(G)
      if (Err_CFML%Ierr /= 0) return

      P=Get_P_Matrix(G)
      if (Err_CFML%Ierr /= 0) return

      Mp=Get_Mp_Matrix(G,P)
      if (Err_CFML%Ierr /= 0) return

      Mc=Get_Mc_Matrix(G%laue,Mp)
      if (Err_CFML%Ierr /= 0) return

      M = matmul(Mp,Mc)
      lattyp=Get_Lattice_Type(M)
      if (Err_CFML%Ierr /= 0) return

      call Match_Shubnikov_Group(G,P,M)

      return
   End Subroutine Identify_Shubnikov_Group
    
End Submodule SPG_045    