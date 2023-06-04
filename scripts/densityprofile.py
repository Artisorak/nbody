import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('../data/data.txt')

x = data[:,2]
y = data[:,3]
z = data[:,4]
r = np.sqrt(x**2 + y**2 + z**2)

m = data[:,1]
M = np.sum(m)

print("Number of particles:",len(m))

halfmassradius = 0
tmpmass = 0
rsorted = np.sort(r)
for i in range(len(rsorted)):
    tmpmass += m[i]
    if tmpmass > M/2:
        halfmassradius = rsorted[i]
        break
print("Half mass radius:",halfmassradius)

rmax = 128
dr = 0.075
nbins = int(rmax/dr)

hist, bins = np.histogram(r, bins=nbins, range=(0,rmax), weights=m)

for i in range(len(hist)):
    volume = 4/3*np.pi*(bins[i+1]**3 - bins[i]**3)
    hist[i] /= volume

a = halfmassradius/(1+np.sqrt(2))
rhoHernquist = M / (2*np.pi) * a / (bins[1:]-dr/2) / (bins[1:]-dr/2 + a)**3

print("Number of bins:",nbins)
print("Bin size:",dr)

# for i in range(len(hist)):
#     if hist[i] == 0:
#         hist[i] = np.nan

# plot poissonian error bars
plt.errorbar(bins[1:], rhoHernquist, yerr=np.sqrt(rhoHernquist), color='grey', fmt='.', label='Hernquist with Poisson error', markersize=0)
# plot density profile of data 
plt.loglog(bins[1:]-dr/2, hist, 'b-', label='Density Profile from data', linewidth=1)
# plot Hernquist profile
plt.loglog(bins[1:]-dr/2, rhoHernquist, 'r-', label='Hernquist Density Profile', linewidth=1)

plt.xlabel(r'Radius $r$')
plt.ylabel(r'Density $\rho(r)$')

plt.legend()
plt.savefig('../plots/densityprofile.png', dpi=200)
plt.show()