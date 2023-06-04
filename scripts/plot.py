import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('../out/data.txt')
steps = np.linspace(0, len(data), len(data))

plt.plot(steps, data)
plt.xlabel('Steps')
plt.savefig('../plots/plot.png', dpi=300)