import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))

import CFML_api

while True:
    a=CFML_api.SpaceGroup(5)
    b=CFML_api.SpaceGroup(10)
