import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))

import CFML_api
import numpy as np 

cellv = np.asarray([1,2,3], dtype='float32')
angl = np.asarray([1,2,3], dtype='float32')

print("======")

a=CFML_api.Cell(cellv, angl)

print("======")
a.printDescription()
