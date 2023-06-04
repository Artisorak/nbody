import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('../out/energy.dat')
ekin = data[:,0]/data[0,0]
epot = data[:,1]/data[0,1]
etot = (data[:,0] + data[:,1]) / (data[0,0] + data[0,1])
steps = np.linspace(0, len(data), len(data))

plt.plot(steps*10, ekin, label='Kinetic Energy')
plt.plot(steps*10, epot, label='Potential Energy')
plt.plot(steps*10, etot, label='Total Energy')

plt.xlabel('Steps')
plt.ylabel(r'Energy $E/E_0$')

plt.legend()

plt.savefig('../plots/energyplot.png', dpi=300)
plt.show()