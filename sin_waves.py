import numpy as np
import matplotlib.pyplot as plt

freq = 10*2*np.pi
period = 1/freq


t_array = np.linspace(0,period*10,1000)

plt.plot(t_array,np.sin(t_array*freq))

plt.show()