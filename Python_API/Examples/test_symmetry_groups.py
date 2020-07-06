import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))

import CFML_api

a=CFML_api.SymmetryGroups(5)
b=CFML_api.SymmetryGroups(10)
print("======")
a.printDescription()
print("======")
b.printDescription()
print("======")
print(b.lattice_translation)
print("======")
print(a.lattice_translation)


#print("======")
#a.lattice_translation = [1,2,3]