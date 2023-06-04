import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation

data = np.loadtxt("../out/positions.dat")
print(data.shape)

skip = data[0,1]

N = int(data[0,0]/skip)
print("using %d of %d particles"%(N, N*skip))

nSteps = int(data.shape[0]/N)
print(nSteps)

dim = data.shape[1]

x = data[1:, 0]
y = data[1:, 1]
z = data[1:, 2]

fig = plt.figure()
ax = plt.axes(projection="3d")

size = 0.5

def animate(i):
    xi = x[i*N:(i+1)*N]
    yi = y[i*N:(i+1)*N]
    zi = z[i*N:(i+1)*N]
    ax.cla()
    ax.scatter3D(xi, yi, zi, s=1, edgecolor="none")
    ax.set_xlim(-size, size)
    ax.set_ylim(-size, size)
    ax.set_zlim(-size, size)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.equal_aspect = True
    ax.set_title("step: %d"%(i*10))

anim = animation.FuncAnimation(fig, animate, frames=nSteps)
mp4writer = animation.FFMpegWriter(fps=24)
anim.save("../plots/3Danimation.mp4", writer=mp4writer, dpi=300)

plt.show()