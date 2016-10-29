import numpy as np
import matplotlib.pyplot as plt


deltaP = np.linspace(10, 100, 10)
Q = np.array([0.45, 0.76, 1.01, 1.15, 1.41, 1.57, 1.66, 1.85, 2.05, 2.25])

plt.figure()
plt.plot(deltaP,Q)
plt.xlabel('Change in Pressure (kPa)')
plt.ylabel('Volumetric Flow Rate (L/min)')
plt.title('Q as a function of Change in Pressure')
plt.show()