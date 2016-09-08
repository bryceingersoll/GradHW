import matplotlib.pyplot as plt
import numpy as np

Q = np.array([0, 100, 150, 200, 250, 300, 325, 350]) # m^3/hr
P = np.array([361, 349, 328, 293, 230, 145, 114, 59]) # kPa

# without first entry
Q = np.array([100, 150, 200, 250, 300, 325, 350]) # m^3/hr
Q = Q/3600.0 # m^3/s
P = np.array([349, 328, 293, 230, 145, 114, 59]) # kPa
P = P*1000.0 # Pa

D = 0.1 #m

Area = D**2.0*np.pi/4.0

rho = 999 # kg/m^3

omega = 750 # rpm
omega = omega/60.0 # rps

pi1_750 = P*Q**(-2.0/3.0)*omega**(-4.0/3.0)*rho**(-1.0)
pi2_750 = D*Q**(-1.0/3.0)*omega**(1.0/3.0)

# first plot

plt.plot(pi2_750, pi1_750, '*')
plt.xlabel('$\Pi_2$')
plt.ylabel('$\Pi_1$')
plt.xlim([0.5, 0.8])
plt.ylim([0, 140])
plt.title('$\Pi_1$ as a function of $\Pi_2$')
plt.show()

# line trend analysis

p1 = np.polyfit(x = pi2_750, y = pi1_750, deg=1)
p2 = np.polyfit(x = pi2_750, y = pi1_750, deg=2)


p2_linspace = np.linspace(start=0.5, stop=0.8, num=50)
p2_fit = p2[0]*p2_linspace**(2.0) + p2[1]*p2_linspace + p2[2]

linfit, = plt.plot([0.5, 0.8], [p1[0]*0.5+p1[1], p1[0]*0.8+p1[1]], label='Linear Fit')
quadfit, = plt.plot(p2_linspace, p2_fit,'--', label='Quadratic Fit')
data, = plt.plot(pi2_750, pi1_750, '*', label= 'Given Data')
plt.xlim([0.5, 0.8])
plt.ylim([0, 140])
plt.xlabel('$\Pi_2$')
plt.ylabel('$\Pi_1$')
plt.title('$\Pi_1$ as a function of $\Pi_2$')
legend = plt.legend(handles=[data, linfit, quadfit],loc=4)
plt.show()

# P - Q curves

omega = 700 # rpm
omega = omega/60.0 # rps

P_700 = Q**(2.0/3.0)*omega**(4.0/3.0)*rho*(p2[0]*(D*Q**(-1.0/3.0)*omega**(1.0/3.0))**(2.0)+p2[1]*(D*Q**(-1.0/3.0)*omega**(1.0/3.0))
            +p2[2])

omega = 750 # rpm
omega = omega/60.0 # rps

P_750 = Q**(2.0/3.0)*omega**(4.0/3.0)*rho*(p2[0]*(D*Q**(-1.0/3.0)*omega**(1.0/3.0))**(2.0)+p2[1]*(D*Q**(-1.0/3.0)*omega**(1.0/3.0))
            +p2[2])

omega = 800 # rpm
omega = omega/60.0 # rps

P_800 = Q**(2.0/3.0)*omega**(4.0/3.0)*rho*(p2[0]*(D*Q**(-1.0/3.0)*omega**(1.0/3.0))**(2.0)+p2[1]*(D*Q**(-1.0/3.0)*omega**(1.0/3.0))
            +p2[2])

omega = 850 # rpm
omega = omega/60.0 # rps

P_850 = Q**(2.0/3.0)*omega**(4.0/3.0)*rho*(p2[0]*(D*Q**(-1.0/3.0)*omega**(1.0/3.0))**(2.0)+p2[1]*(D*Q**(-1.0/3.0)*omega**(1.0/3.0))
            +p2[2])

omega = 900 # rpm
omega = omega/60.0 # rps

P_900 = Q**(2.0/3.0)*omega**(4.0/3.0)*rho*(p2[0]*(D*Q**(-1.0/3.0)*omega**(1.0/3.0))**(2.0)+p2[1]*(D*Q**(-1.0/3.0)*omega**(1.0/3.0))
            +p2[2])

Q = Q*3600.0 # m^3/hr

p700, = plt.plot(Q,P_700/1000, '--*', label='$\omega$ = 700 rpm')
p750, = plt.plot(Q,P_750/1000, '--o', label='$\omega$ = 750 rpm')
p800, = plt.plot(Q,P_800/1000, '--x', label='$\omega$ = 800 rpm')
p850, = plt.plot(Q,P_850/1000, '--D', label='$\omega$ = 850 rpm')
p900, = plt.plot(Q,P_900/1000, '--v', label='$\omega$ = 900 rpm')
legend = plt.legend(handles=[p700, p750, p800, p850, p900],loc=3)
plt.ylabel('Pressure Head (kPa)')
plt.xlabel('Q (cubic meters per hour)')
plt.title('Pressure as a function of flow rate, impeller speed')
plt.show()

# # different impeller speeds
#
# # 450
#
# omega = 450 # rpm
# omega = omega/60.0 # rps
#
# pi1_450 = P*Q**(-2.0/3.0)*omega**(-4.0/3.0)*rho**(-1.0)
# pi2_450 = D*Q**(-1.0/3.0)*omega**(1.0/3.0)
#
# # 600
#
# omega = 600 # rpm
# omega = omega/60.0 # rps
#
# pi1_600 = P*Q**(-2.0/3.0)*omega**(-4.0/3.0)*rho**(-1.0)
# pi2_600 = D*Q**(-1.0/3.0)*omega**(1.0/3.0)
#
# # 900
#
# omega = 900 # rpm
# omega = omega/60.0 # rps
#
# pi1_900 = P*Q**(-2.0/3.0)*omega**(-4.0/3.0)*rho**(-1.0)
# pi2_900 = D*Q**(-1.0/3.0)*omega**(1.0/3.0)
#
# line450, = plt.plot(pi2_450, pi1_450, '--x', label='$\omega$ = 450 rpm')
# line600, = plt.plot(pi2_600, pi1_600, '--*', label='$\omega$ = 600 rpm')
# line750, = plt.plot(pi2_750, pi1_750, '--o', label='$\omega$ = 750 rpm')
# line900, = plt.plot(pi2_900, pi1_900, '--v', label='$\omega$ = 900 rpm')
#
#
#
# legend = plt.legend(handles=[line450, line600, line750, line900],loc=1)
#
# plt.xlabel('$\Pi_2$')
# plt.ylabel('$\Pi_1$')
# plt.title('$\Pi_1$ as a function of $\Pi_2$')
# plt.show()