import matplotlib.pyplot as plt
import numpy as np

D = np.array([6, 10, 3, 1.0/16.0*25.4, 5.0/64*25.4, 3.0/32.0*25.4, 1.0/8.0*25.4, 3.0/16.0*25.4, 8.0, 20.0, 15.0])/1000.0

Bd = np.array([3.5, 8, 3.5, 2.875, 3.0, 2.75, 3.75, 5, 5.5, 12.5, 12.0])*25.4/1000.0

g = 9.81  # m/s^2
h = 0.17  # m
rho = 999.0  # kg/m^3
sigma = 7.28*1.00**(-2.0)  # N/m^2
U = (2.0*g*h)**2.0

We = rho*U*U*D/sigma

plt.plot(D, We)
plt.plot(D, Bd/D, '*')
plt.xlabel('Hole Diameter (m)')
plt.ylabel('Break-up Distance (m)')
plt.title('Break-up Distance of Water Jet to Hole Diameter')
plt.xlim([0,21])
plt.show()
