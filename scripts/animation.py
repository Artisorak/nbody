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

fig = plt.figure()

size = 1

def animate(i):
    xi = x[i*N:(i+1)*N]
    yi = y[i*N:(i+1)*N]
    plt.cla()
    plt.scatter(xi, yi, s=1, edgecolor="none")
    # plt.scatter(np.sum(xi)/N, np.sum(yi)/N, s=1, c="red")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.axis("equal")
    plt.xlim(-size, size)
    plt.ylim(-size, size)
    plt.title("step: %d"%(i*10))

anim = animation.FuncAnimation(fig, animate, frames=nSteps)
mp4writer = animation.FFMpegWriter(fps=24)
anim.save("../plots/animation.mp4", writer=mp4writer, dpi=300)