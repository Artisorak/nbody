import numpy as np
from matplotlib import pyplot as plt

data = np.loadtxt('data0.txt')
x = data[:,2]
y = data[:,3]

plt.scatter(x,y, s=1)
plt.axis('equal')
plt.show()