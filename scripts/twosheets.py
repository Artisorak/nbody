import numpy as np

N = 50010
m = 92.42590000
M = N*m

Lx = int(np.sqrt(N/2))
Ly = int(N/2/Lx)
N = Lx*Ly*2
D = 0.1
A = 0.5

x = np.array([])
x1 = np.linspace(D, D+A, Lx)
for i in range(Lx):
    x = np.concatenate((x, x1))

y = np.array([])
y1 = np.linspace(-A/2, A/2, Ly)
for i in range(Ly):
    y = np.concatenate((y, [y1[i]]*Lx))

z1 = np.array([0.1]*int(N/2))

x2 = np.linspace(-D, -D-A, Lx)
for i in range(Lx):
    x = np.concatenate((x, x2))

y2 = np.linspace(-A/2, A/2, Ly)
for i in range(Ly):
    y = np.concatenate((y, [y1[i]]*Lx))

z2 = np.array([-0.1]*int(N/2))

z = np.concatenate((z1, z2))

v = 3000
vx1 = np.array([-v]*int(N/2))
vy1 = np.array([0]*int(N/2))
vz1 = np.array([0]*int(N/2))

vx2 = np.array([v]*int(N/2))
vy2 = np.array([0]*int(N/2))
vz2 = np.array([0]*int(N/2))

vx = np.concatenate((vx1, vx2))
vy = np.concatenate((vy1, vy2))
vz = np.concatenate((vz1, vz2))

softening = 0.1
potential = 0.01302150

# write to file
with open('../data/twosheets.ascii', 'w') as f:
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