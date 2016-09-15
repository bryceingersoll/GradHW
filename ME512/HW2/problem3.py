import numpy as np
import matplotlib.pyplot as plt

h = np.linspace(0,0.4,100)

g = 9.81
v1 = 3
A = 0.01

r = ((v1*A)/(np.pi*(2*g*h-v1**2.0)**2.0))**0.5

plt.plot(r,h)
plt.xlabel('r (m)')
plt.ylabel('h (m)')
plt.title('r(h)')
plt.show()
