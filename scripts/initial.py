import numpy as np
from matplotlib import pyplot as plt

data = np.loadtxt("../data/twosheets.ascii")
print(data.shape)

N = int(data.shape[0]/9)
print("using %d particles"%N)

x = data[1+1*N:2*N]
y = data[1+2*N:3*N]
z = data[1+3*N:4*N]

fig = plt.figure()
ax = plt.axes(projection="3d")

ax.scatter3D(x, y, z, s=1, edgecolor="none")
ax.set_xlim(-1, 1)
ax.set_ylim(-1, 1)
ax.set_zlim(-1, 1)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
ax.equal_aspect = True

plt.show()