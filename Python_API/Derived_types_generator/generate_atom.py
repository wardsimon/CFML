from DerivedType import DerivedType

d = DerivedType("Atom", "Atom_Type", "atom_typedef", "atom")


d.addParam("label","Lab","string","label")         
d.addParam("chemical_symbol","ChemSymb","string","Chemical Symbol")    
d.addParam("sfac_symbol","SfacSymb","string","Symbol for Scattering Factor")
d.addParam("wyckoff_letter","wyck","string","Wyckoff letter")

d.addParam("active","Active","bool","Flag to control is the atom is included or excluded")

d.addParam("atomic_number","Z","object","Atomic number")
d.addParam("site_multiplicity","Mult","object","Multiplicity of the site")

d.addParam("xyz","X","ndarray","Fractional coordinates")
d.addParam("xyz_std_dev","X_Std","ndarray","Standard deviation of the fractional coordinates")
d.addParam("multiplier_xyz","MX","ndarray","Multipliers for xyz")
d.addParam("lsq_xyz","LX","ndarray","Labels in the LSQ list (refined parameters) for xyz")

d.addParam("occ","Occ","object","Occupation factor")
d.addParam("occ_std_dev","Occ_Std","object","Standard deviation of the occupation factor")
d.addParam("multiplier_occ","MOcc","object","Multiplier for occ")
d.addParam("lsq_occ","LOcc","object","Label in the LSQ list for occ")

d.addParam("biso","Biso","object","Isotropic B-factor")
d.addParam("biso_std_dev","Biso_std","object","Standard deviation for biso")
d.addParam("multiplier_biso","MBiso","object","Multiplier for biso")
d.addParam("lsq_biso","LBiso","object","Label in the LSQ list for biso")



# Betas are defined by the following expression of the temperature factor:
# Taniso= exp( -(beta11 h^2 + beta22 k^2 + beta33 l^2 + 2 (beta12 h k + beta13 h l + beta23 k l)) )
# Taniso= exp( -(bet(1) h^2 + bet(2) k^2 + bet(3) l^2 + 2 (bet(4) h k + bet(5) h l + bet(6) k l)) )

# Us are defined by the following expression of the temperature factor:
# Taniso= exp( -2pi^2 (h^2 (a*)^2 U11+ k^2 (b*)^2 U22+ l^2 (c*)^2 U33+
#                2 (h k (a*) (b*) U12+ h l (a*) (c*) U13+  k l (b*) (c*) U23)) )

# Us and Bs are related by a scaling factor B = 8*Pi^2 U

d.addParam("u_type","Utype","string","Type of anisotropic thermal parameters: u_ij, b_ij, beta, none")
d.addParam("thermal_type","ThType","string","Type of behavior wrt to the temperature: isotr, aniso or other")

d.addParam("u_array","U","ndarray","Array of anisotropic thermal coefficients U11, U22, U33, U12, U13, U23 or B11, B22, B33, B12, B13, B23 depending on the value of Utype")
d.addParam("u_array_std_dev","U_std","ndarray","Array of standard deviation of the anisotropic thermal coefficients")
d.addParam("multiplier_u_array","MU","ndarray","Array of multipliers for the anisotropic thermal coefficients")
d.addParam("lsq_u_qrrqy","LU","ndarray","Array of labels in the LSQ list for the anisotropic thermal coefficients")

d.addParam("u_equiv","Ueq","object","Equivalent isotropic atomic displacement parameter, U(equiv), in angstroms squared, calculated from anisotropic atomic displacement parameters. \n\
U(equiv) = (1/3) sum~i~[sum~j~(U^ij^ a*~i~ a*~j~ a~i~ a~j~)]")

d.addParam("formal_charge","Charge","object","The net integer charge assigned to this atom. This is the formal charge assignment normally found in chemical diagrams.")
d.addParam("moment","Moment","object","Theoretical value of the magnetic moment. i.e., number of unpaired electrons for transition metal oxides")

d.addParam("index_sfac","Ind","ndarray","Index for different purposes (e.g. 1:pointer to the scattering factor, 2:pointer to the magnetic scattering factor")

d.addParam("nvarF","NVar","object","Dimension of VarF")
d.addParam("varF","VarF","ndarray","Free variables used for different purposes (1,2,3 reserved for occupations, not refinable)")
d.addParam("multiplier_varF","MVarF","ndarray","Multipliers of the free variables varF")
d.addParam("lsq_varF","LVarF","ndarray","Labels in the LSQ list of the free variables varF")

d.addParam("info","AtmInfo","string","Information String")


d.addParam("m_xyz","m_xyz","ndarray","Magnetic moment along x,y,z")
d.addParam("m_xyz_std_dev","sm_xyz","ndarray","Standard deviation of the magnetic moment along x,y,z")
d.addParam("multiplier_m_xyz","Mm_xyz","ndarray","Multipliers of the magnetic moment along x,y,")
d.addParam("lsq_m_xyz","Lm_xyz","ndarray","Label in the LSQ list of the magnetic moment along x,y,")

d.generate_template("file_atom")