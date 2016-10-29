import numpy as np
import matplotlib.pyplot as plt

h = 2.5*0.001

mu_l = 0.5
mu_u = 5.0

k = -1000.0

C1 = k*h/(2.0*mu_u)*((mu_u-mu_l)/(mu_u+mu_l))
C3 = k*h/(2.0*mu_l)*((mu_u-mu_l)/(mu_u+mu_l))

C2 = k*h*h/(mu_u+mu_l)
C4 = C2

y_u = np.linspace(0,h,5000)
y_l = np.linspace(-h,0,5000)

u_u = k/(2*mu_u)*y_u*y_u+C1*y_u+C2

u_l = k/(2*mu_u)*y_l*y_l+C3*y_l+C4

mu_1 = 0.5
mu_2 = 5.0
y_1 = y_l
y_2 = y_u

u1 = k/(2.0*mu_1)*(y_1**(2.0)+y_1*h*(mu_2-mu_1)/(mu_2+mu_1))-k*h**(2.0)/(mu_2+mu_1)
u2 = k/(2.0*mu_2)*(y_2**(2.0)+y_2*h*(mu_2-mu_1)/(mu_2+mu_1))-k*h**(2.0)/(mu_2+mu_1)

plt.figure()

plt.plot(u1,y_1,'g')
plt.plot(u2,y_2,'--')
plt.plot()
plt.xlabel('Velocity (m/s)')
plt.ylabel('Height (m)')
plt.title('Velocity Profile')
plt.show()