#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# from matplotlib import rc, rcParams
# from matplotlib.image import NonUniformImage
# from matplotlib import cm, ticker
# from matplotlib.font_manager import FontProperties
# import matplotlib
# config = {
#     "font.family": ["serif", "Times New Roman"],
#     "font.size": 20,
#     "mathtext.fontset": 'stix',
#     "font.serif": ['Computer Modern'],
#     "text.latex.preamble": r"\usepackage{amsmath}"
# }
# rcParams.update(config)

chan = 2

from scipy.interpolate import interp1d
data50 = np.load("Data/number50.npy")
data100 = np.load("Data/number100.npy")
data150 = np.load("Data/number150.npy")
bkg50 = 7444.89
bkg100 = 2131.04
bkg150 = 533.085

if chan == 1:
    data = data50
    bkg = bkg50
elif chan == 2:
    data = data100
    bkg = bkg100
elif chan == 3:
    data = data150
    bkg = bkg150

sgn = interp1d(np.linspace(0., 1., data.shape[0]), data[:, 1], kind="cubic")
# sgn = interp1d(np.linspace(0., 1., data50.shape[0]), data50[:, 1], kind="quadratic")

xx = np.linspace(0., 1., 1000)
yy1 = sgn(xx)

fig = plt.figure(figsize=(5,5))
ax = fig.add_axes([0.16, 0.16, 0.82, 0.82])
# ax.plot(np.linspace(0., 1., data50.shape[0]), data50[:, 1], "o")
# ax.plot(xx, yy1, "-")
# ax.plot(xx, yy2, "-")


mashgrid = pd.DataFrame(
    index=np.linspace(0., 1., 30),
    columns=np.linspace(0., 1., 30)
).unstack().reset_index().rename(columns={'level_0': 'yy', 'level_1': 'xx', 0: 'z'})
mashgrid['z'] = (sgn(mashgrid["xx"]) * mashgrid['yy'] * mashgrid['yy']) / (sgn(mashgrid["xx"]) * mashgrid['yy'] * mashgrid['yy'] + bkg)**0.5


from matplotlib.tri import Triangulation, TriAnalyzer, UniformTriRefiner
plottri = Triangulation(mashgrid['xx'], mashgrid['yy'])
refiner = UniformTriRefiner(plottri)
tri_refine_zz, zz_refine = refiner.refine_field(mashgrid['z'], subdiv=1)

ax.tricontour(tri_refine_zz, zz_refine, levels=[1.6, 2.0, 2.4], transform=ax.transAxes)
# ax.tricontourf(tri_refine_zz, zz_refine, 100, transform=ax.transAxes)
ax.set_xlim(data[0, 0], data[-1, 0])
ax.set_ylim(0., 1.)

plt.show()
# print(mashgrid)

