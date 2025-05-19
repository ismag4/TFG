import numpy as np
from getdist import MCSamples
from getdist import plots
import matplotlib.pyplot as plt

h = np.loadtxt('lcdm_cuv.txt', skiprows = 1, usecols = 0, unpack = False)
om = np.loadtxt('lcdm_cuv.txt', skiprows = 1, usecols = 1, unpack = False)
ok = np.loadtxt('lcdm_cuv.txt', skiprows = 1, usecols = 3, unpack = False)
ol = np.loadtxt('lcdm_cuv.txt', skiprows = 1, usecols = 4, unpack = False)

samples = np.array([h, om, ok, ol]).T
names = ['h', 'omega_m', 'omega_k', 'omega_lambda']
labels = ['h', r'\Omega_m', r'\Omega_K', r'\Omega_\Lambda']

mc_samples = MCSamples(samples=samples, names=names, labels=labels)

print(mc_samples.getTable(limit=1).tableTex())

g = plots.get_subplot_plotter()
g.triangle_plot(mc_samples, filled=True)
plt.show()
