SubModule (CFML_gSpaceGroups) SPG_Identify_Group
   implicit none
   Contains

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
      class(spg_type),    intent(in out) :: G

      !---- Local Variables ----!
      character(len=5) :: car
      integer          :: n
      logical          :: pout

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
         if (G%Numspg == 0) then
            if (G%Numshu > 0) then
               car=adjustl(nlabel_bns(G%numshu))
               n=index(car,'.')
               car=car(:n-1)
               read(unit=car,fmt='(i3)') G%numspg
            end if
         end if
      else ! Superspace groups
         Err_CFML%Ierr = 1
         Err_CFML%Msg = "Identify_Group@SPACEG: Superspace groups not implemented yet"
         return
      end if

   End Subroutine Identify_Group

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
      class(spg_type),    intent(in out) :: G

      !---- Local variables ---!
      type(rational), dimension(3,3)   :: P,Mp,Mc,M
      logical :: pout

      !>===== DEBUG =====
      pout=.false.
      pout=(pout .or. CFML_DEBUG)
      !>=================

      if(.not. Magnetic_DBase_allocated) then
        call Read_Magnetic_Data()
        if (Err_CFML%IErr /=0) return
      end if

      call Identify_PointGroup(G)
      if (Err_CFML%Ierr /= 0) return

      call Identify_LaueClass(G)
      if (Err_CFML%Ierr /= 0) return

      call Identify_Crystal_System(G)
      if (Err_CFML%Ierr /= 0) return

      P=Get_P_Matrix(G)
      if (Err_CFML%Ierr /= 0) return

      Mp=Get_Mp_Matrix(G,P)
      if (Err_CFML%Ierr /= 0) return

      Mc=Get_Mc_Matrix(G%laue,Mp)
      if (Err_CFML%Ierr /= 0) return

      M = matmul(Mp,Mc)
      G%spg_lat=Get_Lattice_Type(M)
      if (Err_CFML%Ierr /= 0) return

      call Match_Shubnikov_Group(G,P,M)
      if (Err_CFML%Ierr /= 0) then
         write(unit=*,fmt="(a)") " => "//trim(err_CFML%msg)
      end if

      return
   End Subroutine Identify_Shubnikov_Group

End Submodule SPG_Identify_Group