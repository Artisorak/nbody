import numpy as np
import matplotlib.pyplot as plt

forces = np.loadtxt('../out/forces.txt')
print("forces: ", forces.shape)

fx = forces[:,0]
fy = forces[:,1]
fz = forces[:,2]
fr = np.sqrt(fx**2 + fy**2 + fz**2)

data = np.loadtxt('../data/data.txt')
print("data: ", data.shape)

x = data[:,2]
y = data[:,3]
z = data[:,4]
r = np.sqrt(x**2 + y**2 + z**2)

m = data[:,1]

# plt.plot(np.log(r), fr, '.', label='Force Profile from calculation', markersize=1)
plt.plot(r, fr, '.', label='Force Profile from calculation', markersize=1)

plt.xlabel(r'Radius $log(r)$')
plt.xlabel(r'Radius $r$')
plt.ylabel(r'Force $F(r)$')
plt.xlim(right=400)

plt.legend()
plt.savefig('../plots/forceprofile.png', dpi=200)
plt.show()