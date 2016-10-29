import numpy as np
import matplotlib.pyplot as plt

C1 = 1.0
C2 = 1.0

r = 1.0

R = np.linspace(r,2*r,50)

u_theta = R*C1 + C2/R

rho = 1000.0

P = rho*(C1*np.log(R)-1.0/2.0/R**2.0*C2)

plt.figure()

plt.plot(R, u_theta,'g')
plt.plot()
plt.ylabel('Velocity Profile in Theta Direction')
plt.xlabel('Radial Position')
plt.title('Velocity Profile')
plt.show()

plt.figure()

plt.plot(R,P,'g')
plt.plot()
plt.ylabel('Pressure Profile in Theta Direction')
plt.xlabel('Radial Position')
plt.title('Pressure Profile')
plt.show()

omega = 1.0

u_p = 1.0*1.0*1.0*R

plt.plot(R, u_p,'g')
plt.plot()
plt.ylabel('Velocity Profile in Theta Direction')
plt.xlabel('Radial Position')
plt.title('Velocity Profile for Potential Flow')
plt.show()

gamma = 2*np.pi*1.0*R*R

P_p = -1/(8.0*np.pi*np.pi)*rho*gamma*gamma/(R*R)

plt.plot(R,P,'g')
plt.plot()
plt.ylabel('Pressure Profile in Theta Direction')
plt.xlabel('Radial Position')
plt.title('Pressure Profile for Potential Flow')
plt.show()