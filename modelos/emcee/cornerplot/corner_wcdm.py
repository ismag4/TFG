import numpy as np
from getdist import MCSamples
from getdist import plots
import matplotlib.pyplot as plt

h = np.loadtxt('wcdm_param.txt', skiprows = 1, usecols = 0, unpack = False)
om = np.loadtxt('wcdm_param.txt', skiprows = 1, usecols = 1, unpack = False)
w = np.loadtxt('wcdm_param.txt', skiprows = 1, usecols = 2, unpack = False)
ol = np.loadtxt('wcdm_param.txt', skiprows = 1, usecols = 3, unpack = False)

samples = np.array([h, om, w, ol]).T
names = ['h', 'omega_m','w', 'omega_lambda']
labels = ['h', r'\Omega_m','w', r'\Omega_\Lambda']

mc_samples = MCSamples(samples=samples, names=names, labels=labels)

print(mc_samples.getTable(limit=1).tableTex())

g = plots.get_subplot_plotter()
g.triangle_plot(mc_samples, filled=True)
plt.show()
