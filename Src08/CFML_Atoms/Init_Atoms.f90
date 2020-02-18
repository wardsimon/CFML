!!----
!!----
!!----
SubModule (CFML_Atoms) Atm_002
   Contains

   !!----
   !!---- INIT_ATOM_TYPE
   !!----    Initialize Atom_Type
   !!----
   !!---- 12/06/2019
   !!
   Module Subroutine Init_Atom_Type(Atm,d)
      !---- Arguments ----!
      class(Atm_Type), intent(in out)   :: Atm
      integer,         intent(in)       :: d !Number of k-vectors

      Atm%Lab      =" "
      Atm%ChemSymb =" "
      Atm%Z        =0
      Atm%Mult     =0
      Atm%Charge   =0
      Atm%X        =0.0_cp
      Atm%Occ      =1.0_cp

      !> Thermal model
      Atm%UType    ="B"
      Atm%ThType   ="iso"
      Atm%U_iso    = 0.0_cp
      Atm%U        = 0.0_cp

      !> Magnetic
      Atm%Magnetic =.false.
      Atm%Moment   = 0.0_cp

      !> SFac
      Atm%SfacSymb =" "
      Atm%Ind_ff   = 0

      !> Info
      Atm%AtmInfo  = " "
      Atm%wyck      = " "

      select type (Atm)
         type is (atm_std_type)
            Atm%X_Std     = 0.0_cp
            Atm%Occ_Std   = 0.0_cp
            Atm%U_iso_Std = 0.0_cp
            Atm%U_Std     = 0.0_cp
            Atm%Moment_Std= 0.0_cp

         type is (matm_std_type)
            Atm%X_Std     = 0.0_cp
            Atm%Occ_Std   = 0.0_cp
            Atm%U_iso_Std = 0.0_cp
            Atm%U_Std     = 0.0_cp
            Atm%Moment_Std= 0.0_cp

            Atm%n_oc      = 0
            Atm%n_mc      = 0
            Atm%n_dc      = 0
            Atm%n_uc      = 0
            Atm%poc_q     = 0
            Atm%pmc_q     = 0
            Atm%pdc_q     = 0
            Atm%puc_q     = 0
            Atm%Ocs       = 0.0_cp
            Atm%Ocs_std   = 0.0_cp
            Atm%Mcs       = 0.0_cp
            Atm%Mcs_std   = 0.0_cp
            Atm%Dcs       = 0.0_cp
            Atm%Dcs_std   = 0.0_cp
            Atm%Ucs       = 0.0_cp
            Atm%Ucs_std   = 0.0_cp
            if(allocated (Atm%Xs)) deallocate(Atm%Xs)
            if(allocated (Atm%Us)) deallocate(Atm%Us)
            if(allocated (Atm%Moms)) deallocate(Atm%Moms)
            allocate(Atm%Xs(3+d),Atm%Us(3+d,3+d),Atm%Moms(3+d))
            Atm%Xs=0.0; Atm%Us=0.0; Atm%Moms=0.0

         type is (atm_ref_type)
            Atm%X_Std     = 0.0_cp
            Atm%Occ_Std   = 0.0_cp
            Atm%U_iso_Std = 0.0_cp
            Atm%U_Std     = 0.0_cp
            Atm%Moment_Std= 0.0_cp

            Atm%L_X       =0
            Atm%L_Occ     =0
            Atm%L_U_iso   =0
            Atm%L_U       =0
            Atm%L_moment  =0
            Atm%M_X       =0.0_cp
            Atm%M_Occ     =0.0_cp
            Atm%M_U_iso   =0.0_cp
            Atm%M_U       =0.0_cp
            Atm%M_moment  =0

         type is (Matm_ref_type)
            Atm%X_Std     = 0.0_cp
            Atm%Occ_Std   = 0.0_cp
            Atm%U_iso_Std = 0.0_cp
            Atm%U_Std     = 0.0_cp
            Atm%Moment_Std= 0.0_cp

            Atm%n_oc      = 0
            Atm%n_mc      = 0
            Atm%n_dc      = 0
            Atm%n_uc      = 0
            Atm%poc_q     = 0
            Atm%pmc_q     = 0
            Atm%pdc_q     = 0
            Atm%puc_q     = 0
            Atm%Ocs       = 0.0_cp
            Atm%Ocs_std   = 0.0_cp
            Atm%Mcs       = 0.0_cp
            Atm%Mcs_std   = 0.0_cp
            Atm%Dcs       = 0.0_cp
            Atm%Dcs_std   = 0.0_cp
            Atm%Ucs       = 0.0_cp
            Atm%Ucs_std   = 0.0_cp

            Atm%L_X       =0
            Atm%L_Occ     =0
            Atm%L_U_iso   =0
            Atm%L_U       =0
            Atm%M_X       =0.0_cp
            Atm%M_Occ     =0.0_cp
            Atm%M_U_iso   =0.0_cp
            Atm%M_U       =0.0_cp

            Atm%L_Ocs       = 0.0_cp
            Atm%L_Mcs       = 0.0_cp
            Atm%L_Dcs       = 0.0_cp
            Atm%L_Ucs       = 0.0_cp
            Atm%M_Ocs       = 0.0_cp
            Atm%M_Mcs       = 0.0_cp
            Atm%M_Dcs       = 0.0_cp
            Atm%M_Ucs       = 0.0_cp
            if(allocated (Atm%Xs)) deallocate(Atm%Xs)
            if(allocated (Atm%Us)) deallocate(Atm%Us)
            if(allocated (Atm%Moms)) deallocate(Atm%Moms)
            allocate(Atm%Xs(3+d),Atm%Us(3+d,3+d),Atm%Moms(3+d))
            Atm%Xs=0.0; Atm%Us=0.0; Atm%Moms=0.0
      end select

   End Subroutine Init_Atom_Type

End SubModule Atm_002