import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('../data/data.txt')

x = data[:,2]
y = data[:,3]
z = data[:,4]

m = data[:,1]
M = np.sum(m)

print("Number of particles:",len(m))

maxr = 16
dr = maxr/100

a = dr

r = np.arange(dr, maxr, dr)
rhoHernquist = M / (2*np.pi) * a / r / (r + a)**4
rho = np.zeros(len(r))

for i in range(len(x)):
    ri2 = x[i]**2 + y[i]**2 + z[i]**2
    for j in range(len(r)):
        if (ri2 > r[j]**2 and ri2 < (r[j] + dr)**2):
            rho[j] += m[i]
            break


for i in range(len(r)-1):
    volume = 4*np.pi*(r[i+1]**3 - r[i]**3)/3
    rho[i] /= volume

#rho /= 4*np.pi*r**2*dr

# Normalize
# rhoHernquist /= np.sum(rhoHernquist)
# rho /= np.sum(rho)

print("Number of bins:",len(rho))
print("Bin size:",dr)

# plot poissonian error bars
plt.errorbar(r, rhoHernquist, yerr=np.sqrt(rhoHernquist), color='grey', fmt='.', label='Hernquist with Poisson error', markersize=0)
# plot Hernquist profile
plt.loglog(r, rhoHernquist, 'r.-', label='Hernquist Density Profile')
# plot density profile of data 
plt.loglog(r, rho, 'b.-', label='Density Profile from data')

plt.xlabel('Radius r')
plt.ylabel('Density rho(r)')

plt.legend()
plt.savefig('../plots/densityprofile.png', dpi=200)
plt.show()