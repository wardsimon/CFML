from DerivedType import DerivedType

d = DerivedType("DiffractionPattern", "Diffraction_Pattern_Type", "diffraction_patterns", "diffractionpattern")

# Casual getters
d.addParam("title",     "title",     "string",  "Identification of the pattern"                             )
d.addParam("diff_kind", "diff_kind", "string",  "type of radiation"                                         )
d.addParam("scat_var",  "scat_var",  "string",  "x-space: 2theta, TOF, Q, s, d-spacing, SinT/L, etc"        )
d.addParam("xax_text",  "xax_text",  "string",  "x-axis legend, eg. Lambda (Angstroms)"                     )
d.addParam("yax_text",  "yax_text",  "string",  "y-axis legend, eg. Intensity (arb. units)"                 )
d.addParam("instr",     "instr",     "string",  "file type"                                                 )
d.addParam("filename",  "filename",  "string",  "file name"                                                 )
d.addParam("filepath",  "filepath",  "string",  "file path"                                                 )
d.addParam("xmin",      "xmin",      "int",     ""                                                          )
d.addParam("xmax",      "xmax",      "int",     ""                                                          )
d.addParam("ymin",      "ymin",      "int",     ""                                                          )
d.addParam("ymax",      "ymax",      "int",     ""                                                          )
d.addParam("scal",      "scal",      "int",     ""                                                          )
d.addParam("monitor",   "monitor",   "int",     ""                                                          )
d.addParam("norm_mon",  "norm_mon",  "int",     ""                                                          )
d.addParam("col_time",  "col_time",  "int",     ""                                                          )
d.addParam("step",      "step",      "int",     ""                                                          )
d.addParam("zerop",     "zerop",     "int",     ""                                                          )
d.addParam("Tsamp",     "Tsamp",     "int",     "Sample Temperature"                                        )
d.addParam("Tset",      "Tset",      "int",     "Setting Temperature (wished temperature)"                  )
d.addParam("npts",      "npts",      "int",     "Number of points"                                          )
d.addParam("ct_step",   "ct_step",   "bool",    "Constant step"                                             )
d.addParam("gy",        "gy",        "bool",    "logicals for graphics"                                     )
d.addParam("gycalc",    "gycalc",    "bool",    "logicals for graphics"                                     )
d.addParam("gbgr",      "gbgr",      "bool",    "logicals for graphics"                                     )
d.addParam("gsigma",    "gsigma",    "bool",    "logicals for graphics"                                     )
d.addParam("sig_var",   "sig_var",   "bool",    "If .true. the content of sigma is in fact the variance"    )
d.addParam("al_x",      "al_x",      "bool",    "logicals for allocation"                                   )
d.addParam("al_y",      "al_y",      "bool",    "logicals for allocation"                                   )
d.addParam("al_ycalc",  "al_ycalc",  "bool",    "logicals for allocation"                                   )
d.addParam("al_bgr",    "al_bgr",    "bool",    "logicals for allocation"                                   )
d.addParam("al_sigma",  "al_sigma",  "bool",    "logicals for allocation"                                   )
d.addParam("al_istat",  "al_istat",  "bool",    "logicals for allocation"                                   )
d.addParam("conv",      "conv",      "ndarray", "Wavelengths or Dtt1, Dtt2 for converting to Q,d, etc"      )
d.addParam("x",         "x",         "ndarray", "Scattering variable (2theta...)"                           )
d.addParam("y",         "y",         "ndarray", " Experimental intensity"                                   )
d.addParam("sigma",     "sigma",     "ndarray", "observations sigma or variance (depending on sig_var)"     )
d.addParam("istat",     "istat",     "ndarray", "Information about the point 'i'"                           )
d.addParam("ycalc",     "ycalc",     "ndarray", "Calculated intensity"                                      )
d.addParam("bgr",       "bgr",       "ndarray", "Background"                                                )
d.addParam("nd",         "nd",       "ndarray", "Number of detectors contributing to the point 'i'"         )

d.generate_template("file_diffpatt")