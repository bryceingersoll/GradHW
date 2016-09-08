import matplotlib.pyplot as plt
import numpy as np

Q = np.array([0, 100, 150, 200, 250, 300, 325, 350]) # m^3/hr
P = np.array([361, 349, 328, 293, 230, 145, 114, 59]) # kPa
P = P*1000.0 # Pa

D = 0.1 #m

Area = D**2.0*np.pi/4.0

rho = 999 # kg/m^3

omega = 750 # rpm
omega = omega/60.0 # rps

pi1_750 = P*D**4.0/(rho*Q**(-2.0))
pi2_750 = omega*D**(3.0)/(Q)

plt.plot(pi2_750, pi1_750, '*')
plt.xlabel('$\Pi_2$')
plt.ylabel('$\Pi_1$')
plt.title('$\Pi_1$ as a function of $\Pi_2$')
plt.show()

omega = 450 # rpm
omega = omega/60.0 # rps

pi1_450 = P*D**4.0/(rho*Q**(-2.0))
pi2_450 = omega*D**(3.0)/(Q)

line450, = plt.plot(pi2_450, pi1_450, '--x', label='$\omega$ = 450')
line750, = plt.plot(pi2_750, pi1_750, '--*', label='$\omega$ = 750')

legend = plt.legend(handles=[line450, line750],loc=1)

plt.xlabel('$\Pi_2$')
plt.ylabel('$\Pi_1$')
plt.title('$\Pi_1$ as a function of $\Pi_2$')
plt.show()