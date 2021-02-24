import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))

import matplotlib.pyplot as plt

import CFML_api

filename = sys.argv[1]

cif_file = CFML_api.CIFFile(filename)

cell = cif_file.cell
space_group = cif_file.space_group
atom_list = cif_file.atom_list
job_info = cif_file.job_info

print(job_info.pattern_types)

cell.print_description()
space_group.print_description()
atom_list.print_description()

powder_pattern_simulation_conditions = CFML_api.PowderPatternSimulationConditions()

reflection_list = CFML_api.ReflectionList(
    cell, space_group, True, 0, powder_pattern_simulation_conditions.getSinThetaOverLambdaMax())
reflection_list.compute_structure_factors(space_group, atom_list)

diffraction_pattern = CFML_api.DiffractionPattern(powder_pattern_simulation_conditions, reflection_list, cell.reciprocal_cell_vol)

plt.plot(diffraction_pattern.x, diffraction_pattern.ycalc)
plt.show()
