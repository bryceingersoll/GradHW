import numpy as np
import matplotlib.pyplot as plt

p_atm = 101.325
p0 = 700.0 + p_atm

k = -0.001242

t = np.linspace(0,30,30)

p_saying = np.ones([len(t),1])*p0

for i in range(0, len(t)):
    p_saying[i] = p_saying[i] - i*6.89476/2.0

p_actual = p_atm/(1-(p0-p_atm)/p0*np.exp(k*t))

actual, = plt.plot(t,p_actual, label='Actual')
saying, = plt.plot(t,p_saying, label='Saying')
plt.xlabel('time (days)')
plt.ylabel('total pressure (kPa)')
plt.title('Total Pressure vs. Time')
legend = plt.legend(handles=[actual, saying],loc=1)
plt.show()