import numpy as np

N = 50010
m = 92.42590000
M = N*m

# create positions
r = np.random.uniform(0.01, 1, N)
theta = np.random.uniform(0, 2*N*np.pi, N)
x = r * np.cos(theta)
y = r * np.sin(theta)
z = 0.2 * np.exp(-(r**2)/0.1) * np.random.uniform(-1, 1, N)

# create velocities orthogonal to radial vector
# v = 4 * np.sqrt(M * r) * np.exp(-r**2)
v = np.sqrt(M) * (r+0.1)
vx =  v * y / r
vy = -v * x / r
vz = 10*np.random.normal(0, 0.1, N)

softening = 0.1

potential = 0.01302150

# write to file
with open('../data/disk.ascii', 'w') as f:
    f.write("#%d 0 %d\n" % (N, N))
    for i in range(N):
        f.write("%f\n" % m)
    for i in range(N):
        f.write("%f\n" % x[i])
    for i in range(N):
        f.write("%f\n" % y[i])
    for i in range(N):
        f.write("%f\n" % z[i])
    for i in range(N):
        f.write("%f\n" % vx[i])
    for i in range(N):
        f.write("%f\n" % vy[i])
    for i in range(N):
        f.write("%f\n" % vz[i])
    for i in range(N):
        f.write("%f\n" % softening)
    for i in range(N):
        f.write("%f\n" % potential)