from DerivedType import DerivedType

d = DerivedType("SpaceGroup", "Space_Group_Type", "crystallographic_symmetry", "spacegroup")

# Special functions
d.addParam("wyckoff_info","wyckoff","object","Wyckoff Information")
d.addParam("symmetry_operators","SymOp","object","Symmetry operators (192)")
d.addParam("symmetry_operators_as_text","SymopSymb","object","String form of symmetry operators")

# Casual getters
d.addParam("lattice_translation", "latt_trans", "ndarray", "Lattice translation matrix")
d.addParam("hexa", "hexa", "bool", "Is space group hexagonal")
d.addParam("space_group_number", "numspg", "object", "Number of the space group")
d.addParam("space_group_symbol", "spg_symb", "string", "Hermann-Mauguin Symbol (str)")
d.addParam("hall_symbol","hall","string","Hall symbol")
d.addParam("generalized_hall_symbol", "ghall","string","Generalised Hall symbol")
d.addParam("crystal_system","CrystalSys","string","Crystal System")
d.addParam("laue_class","Laue","string","Laue class")
d.addParam("point_group","PG","string","Point group")
d.addParam("extra_information","Info","string","Extra Information")
d.addParam("space_group_setting","SG_setting","string","Information about the SG setting (IT,KO,ML,ZA,Table,Standard,UnConventional)")
d.addParam("space_group_lattice_type","SPG_lat","string","Lattice type")
d.addParam("space_group_lattice_type_symbol","SPG_latsy","string","Lattice type Symbol")
d.addParam("number_of_lattice_points","NumLat","object","Number of lattice points in a cell")
d.addParam("bravais","bravais","string","String with Bravais symbol + translations")
d.addParam("centre","centre","string","Alphanumeric information about the center of symmetry")
d.addParam("centric","centred","object","Centric or Acentric [ =0 Centric(-1 no at origin),=1 Acentric,=2 Centric(-1 at origin)]")
d.addParam("inversion_centre","centre_coord","ndarray","Fractional coordinates of the inversion centre")
d.addParam("number_of_symops","numops","object","Number of reduced set of Symmetry Operators")
d.addParam("multiplicity","multip","object","Multiplicity of the general position")
d.addParam("operators_minimum_number","num_gen","object","Minimum number of operators to generate the Group")
d.addParam("asymetric_unit","R_asym_unit","ndarray","Asymmetric unit in real space")

d.generate_template("file_spg")
