import numpy as np
import math
import matplotlib.pyplot as plt

g = 9.81 # m/s^2
h = 0.2*10**(-3.0) #m
m = 5.0 #kg
mu = 239.39*10**(-3.0) # Pa*s
A = 0.2**(2.0) # m^2

t = np.linspace(0,0.5,100)
t = np.array(t)
v = np.zeros([len(t),1])
a = np.zeros([len(t),1])
term = np.zeros([len(t),1])

for i in range(0, len(t)):
    v[i] = g*m*h/(2.0*mu*A)*(1.0-math.exp(-A*mu*t[i]/(m*h)))
    a[i] = g/2.0*math.exp(-A*mu*t[i]/(m*h))


# first plot

plt.plot(t, v)
plt.xlabel('Time (s)')
plt.ylabel('Velocity (m/s)')
plt.title('Block Velocity as a Function of Time')
#plt.xlim([0.5, 0.8])
#plt.ylim([0, 140])
#plt.title('$\Pi_1$ as a function of $\Pi_2$')
plt.show()