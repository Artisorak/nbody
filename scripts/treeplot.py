import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

tree = np.loadtxt("../out/tree.txt")
pos = np.loadtxt("../data/data.txt")

x = pos[:,2]
y = pos[:,3]

maxCoordinate = np.max(np.abs(pos[:,2:4]))
size = 1.5

fig, ax = plt.subplots()

# ax.scatter(tree[:,0], tree[:,1], color='r', s=0.1)
ax.scatter(x, y, color='b', edgecolor="none", s=0.33)
ax.set_xlim(-4/3*size, 4/3*size)
ax.set_ylim(-size, size)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_aspect('equal')

for i in range(tree.shape[0]):
    size = tree[i,3]
    ax.add_patch(Rectangle((tree[i,0] - size, tree[i,1] - size), 2*size, 2*size, linewidth=0.67, edgecolor='k', facecolor='none', alpha=0.5))
    # ax.text(tree[i,0], tree[i,1], str(tree[i,4]), fontsize=8, horizontalalignment='center', verticalalignment='center')

plt.savefig("../plots/treeplot.png", dpi=300)
plt.show()