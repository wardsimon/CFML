from DerivedType import DerivedType

d = DerivedType("Cell", "Cell_Type", "crystal_metrics", "cell")

d.addParam("lattpar", "cell", "ndarray", "")
d.addParam("lattangle", "angl", "ndarray", "")

d.generate_template("file_cell")