import numpy as np
import matplotlib.pyplot as plt

# Read the data from the file
data = np.loadtxt('../out/forcedeviation.txt')

index = data[:,0]
cf = data[:,1]
af = data[:,2]

print('direct force: ', cf)
print('tree force: ', af)

# calculate histogram for both forces
hist1, bin_edges1 = np.histogram(cf, bins=100)
hist2, bin_edges2 = np.histogram(af, bins=100)

line = np.linspace(np.max([np.min(cf), np.min(af)]), np.min([np.max(cf), np.max(af)]), 2)

# plot
# plt.plot(bin_edges1[:-1], hist1, label='correct f')
# plt.plot(bin_edges2[:-1], hist2, label='actual f')
# plt.loglog(bin_edges1[:-1], hist1, label='f')
# plt.loglog(bin_edges2[:-1], hist2, label='actual f')
colors = plt.cm.viridis(np.linspace(0, 1, len(index)))
plt.scatter(cf, af, label='particles', s=1, edgecolor='none', c=colors)
plt.plot(line, line, label='ideal', linewidth=0.5, color='red')
plt.loglog()
plt.xlabel('direct force')
plt.ylabel('tree force')
plt.axis('equal')
plt.legend()
plt.savefig('../plots/forcecomparison.png', dpi=300)
plt.show()