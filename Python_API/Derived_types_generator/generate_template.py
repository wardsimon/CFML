from DerivedType import DerivedType

d = DerivedType("SpaceGroup", "Space_Group_Type", "crystallographic_symmetry")
d.addParam("lattice_translation", "latt_trans", "ndarray", "Lattice translation matrix")
d.addParam("hexa", "hexa", "bool", "Is space group hexagonal")
d.addParam("number_of_space_group", "numspg", "object", "Number of the space group")
d.addParam("space_group_symbol", "spg_symb", "string", "Space group symbol as text")
d.generate_template("file1")

d = DerivedType("Cell", "Cell_Type", "crystal_metrics")
d.addParam("lattpar", "cell", "ndarray", "")
d.addParam("lattangle", "angl", "ndarray", "")
d.generate_template("file2")
    
    
    #generate_template("hall_symbol
    #generate_template("generalized_hall_symbolp
    #generate_template("crystal_system
    #generate_template("laue_class
    #generate_template("point_group
    #generate_template("extra_information
    #generate_template("space_group_setting
    #generate_template("space_group_lattice_type 
    #generate_template("space_group_lattice_type_symbol
    #generate_template("number_of_lattice_points
    #generate_template("lattice_translation
    #generate_template("bravais
    #generate_template("centre
    #generate_template("centric
    #generate_template("inversion_centre
    #generate_template("number_of_reduced_set
    #generate_template("multiplicity
    #generate_template("operators_minimum_number
    #generate_template("symmetry_operators
    #generate_template("symmetry_operators_as_text
    #generate_template("wyckoff_type
    #generate_template("asymetric_unit