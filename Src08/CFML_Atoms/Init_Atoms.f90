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
   Module Subroutine Init_Atom_Type(Atm)
      !---- Arguments ----!
      class(Atm_Type), intent(in out)   :: Atm

      Atm%Lab      =" "
      Atm%ChemSymb =" "
      Atm%Z        =0
      Atm%Mult     =0
      Atm%Charge   =0
      Atm%X        =0.0_cp
      Atm%Occ      =1.0_cp
      
      !> Thermal model
      Atm%UType    ="U"      
      Atm%ThType   ="iso"      
      Atm%U_iso    = 0.0_cp
      Atm%U        = 0.0_cp
      
      !> Magnetic
      Atm%Magnetic =.false.
      Atm%Moment   = 0.0_cp
      
      !> SFac
      Atm%SfacSymb =" "
      Atm%Ind_ff   = 0
      
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
            
            Atm%wyck      = " "
            Atm%n_mc      = 0       
            Atm%n_dc      = 0
            Atm%Mcs       = 0.0_cp  
            Atm%Mcs_std   = 0.0_cp  
            Atm%Dcs       = 0.0_cp  
            Atm%Dcs_std   = 0.0_cp  
            
         type is (atm_ref_type) 
            Atm%X_Std     = 0.0_cp 
            Atm%Occ_Std   = 0.0_cp     
            Atm%U_iso_Std = 0.0_cp     
            Atm%U_Std     = 0.0_cp
            Atm%Moment_Std= 0.0_cp 
            
            Atm%LX       =0      
            Atm%LOcc     =0
            Atm%LU_iso   =0 
            Atm%LU       =0
            Atm%MX       =0.0_cp 
            Atm%MOcc     =0.0_cp
            Atm%MU_iso   =0.0_cp
            Atm%MU       =0.0_cp
            
         type is (Matm_ref_type)  
            Atm%X_Std     = 0.0_cp 
            Atm%Occ_Std   = 0.0_cp     
            Atm%U_iso_Std = 0.0_cp     
            Atm%U_Std     = 0.0_cp
            Atm%Moment_Std= 0.0_cp 
            
            Atm%wyck      = " "
            Atm%n_mc      = 0       
            Atm%n_dc      = 0
            Atm%Mcs       = 0.0_cp  
            Atm%Mcs_std   = 0.0_cp  
            Atm%Dcs       = 0.0_cp  
            Atm%Dcs_std   = 0.0_cp 
            
            Atm%LX       =0      
            Atm%LOcc     =0
            Atm%LU_iso   =0 
            Atm%LU       =0
            Atm%MX       =0.0_cp 
            Atm%MOcc     =0.0_cp
            Atm%MU_iso   =0.0_cp
            Atm%MU       =0.0_cp
      end select
      
   End Subroutine Init_Atom_Type
   
End SubModule Atm_002   