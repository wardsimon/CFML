import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))

import CFML_api
import matplotlib.pyplot as plt
plt.figure(1)

powder_pattern = CFML_api.PowderPatternSimulator()
powder_pattern.compute("Data/Si_Laue")
plt.plot(powder_pattern.x, powder_pattern.y, label="Si_Laue")

powder_pattern = CFML_api.PowderPatternSimulator()
powder_pattern.compute("Data/ponsin")
plt.plot(powder_pattern.x, powder_pattern.y, label="ponsin")

plt.legend()

plt.show()