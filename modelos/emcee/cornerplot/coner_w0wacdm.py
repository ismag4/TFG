import numpy as np
from getdist import MCSamples
from getdist import plots
import matplotlib.pyplot as plt

h = np.loadtxt('w0wacdm_param.txt', skiprows = 1, usecols = 0, unpack = False)
om = np.loadtxt('w0wacdm_param.txt', skiprows = 1, usecols = 1, unpack = False)
w0 = np.loadtxt('w0wacdm_param.txt', skiprows = 1, usecols = 3, unpack = False)
wa = np.loadtxt('w0wacdm_param.txt', skiprows = 1, usecols = 4, unpack = False)
ol = np.loadtxt('w0wacdm_param.txt', skiprows = 1, usecols = 5, unpack = False)

samples = np.array([h, om, w0, wa, ol]).T
names = ['h', 'omega_m','w0', 'wa', 'omega_lambda']
labels = ['h', r'\Omega_m','w_0', 'w_a', r'\Omega_\Lambda']

mc_samples = MCSamples(samples=samples, names=names, labels=labels)

print(mc_samples.getTable(limit=1).tableTex())

g = plots.get_subplot_plotter()
g.triangle_plot(mc_samples, filled=True)
plt.show()
