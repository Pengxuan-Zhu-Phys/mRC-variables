#!/usr/bin/env python3

import cmath
from distutils.log import Log
from turtle import color, left, right
from matplotlib.transforms import Transform
import numpy as np
from matplotlib.ticker import FixedLocator, AutoMinorLocator
import pandas as pd
import os
import sys
from shutil import which
import matplotlib.pyplot as plt
import json
import math 
import matplotlib
from matplotlib import rc, rcParams
from matplotlib.image import NonUniformImage
from matplotlib import cm, ticker
from matplotlib.font_manager import FontProperties
from yaml import DirectiveToken
config = {
    "font.family":["serif", "Times New Roman"],
    "font.size": 20,
    "mathtext.fontset":'stix',
    "font.serif": ['Computer Modern'],
    "text.latex.preamble": r"\usepackage{amsmath}"
}
rcParams.update(config)

fig = plt.figure(figsize=(10, 8))
ax = fig.add_axes([ 0.13, 0.13, 0.715, 0.83])
axc = fig.add_axes([0.85, 0.15, 0.02, 0.79])

bkg = pd.read_csv("bkg.csv")
# dat = pd.read_csv("../Grid_Smuon11/GridData_pure.csv")
dat01 = pd.read_csv("../Grid_Smuon11/GridData_pure.csv")
dat02 = pd.read_csv("../Grid_Smuon12/GridData_pure.csv")
dat03 = pd.read_csv("../Grid_Smuon13/GridData_pure.csv")
dat = pd.concat([dat01, dat02, dat03], ignore_index=True)
# dat = pd.concat([dat01, dat02], ignore_index=True)

SRs = ['SRH01', "SRH02", "SRM01", "SRM02", "SRM03", "SRM04", "SRM05", "SRL01", "SRL02", "SRL03", "SRL04", "SRL05"]
SRc = ['SRH03', "SRM06", "SRL06"]
# print(dat.shape, bkg.shape)
xsect = pd.DataFrame({
    "x":    dat['msmuR'],
    'y':    dat['mn1'],
    "z":    dat['XSectSmuon'] * 1.e3,
    "SRH01":    dat['SRH01'],
    "SRH02":    dat['SRH02'],
    "SRH03":    dat['SRH03'],
    "SRM01":    dat['SRM01'],
    "SRM02":    dat['SRM02'],
    "SRM03":    dat['SRM03'],
    "SRM04":    dat['SRM04'],
    "SRM05":    dat['SRM05'],
    "SRM06":    dat['SRM06'],
    "SRL01":    dat['SRL01'],
    "SRL02":    dat['SRL02'],
    "SRL03":    dat['SRL03'],
    "SRL04":    dat['SRL04'],
    "SRL05":    dat['SRL05'],
    "SRL06":    dat['SRL06'],
})
# print(xsect)
cmap = plt.get_cmap("PuBu")
cmax = 1000
cmin = 1.e-3

for idx, row in xsect.iterrows():
    # print(dict(pd.Series(row)))
    if row['x'] > 119.7:
        row['z'] = 1.8 * row['z']
    if row['z'] < cmin:
        cc = cmin
    else:
        cc = (math.log10(row["z"]) - math.log10(cmin))  / 6
    fc = cmap(cc)[0:3]
    if row['x'] < 110:
        dxm = 1.0 
        dxp = 1.0
    elif row['x'] ==  110:
        dxm = 1.0
        dxp = 0.5
    elif row['x'] > 110 and row['x'] < 115:
        dxm = 0.5
        dxp = 0.5
    elif row['x'] == 115:
        dxm = 0.5
        dxp = 0.25
    elif row['x'] > 115 and row['x'] < 130:
        dxm = 0.25
        dxp = 0.25
    elif row['x'] == 130:
        dxm = 0.5
        dxp = 1.0
    else:
        dxm = 1.0
        dxp = 1.0
    ax.fill(
        [row['x'] - dxm, row['x'] + dxp, row['x'] + dxp, row['x'] - dxm, row['x'] - dxm],
        [row['y'] - 1.0, row['y'] - 1.0, row['y'] + 1.0, row['y'] + 1.0, row['y'] - 1.0],
        facecolor=fc,
        edgecolor=None
    )
from matplotlib.colors import LogNorm
sc = ax.scatter(xsect['x'], xsect['y'], c=xsect['z'], s=1, cmap=cmap, norm=LogNorm(vmin=cmin, vmax=cmax), zorder=0)
# print(xsect)
for sr in SRs:
    xsect.loc[(xsect.x > 119.7), sr] = 1.8 * xsect[sr]
    # print(xsect[sr])
    xsect[sr] = xsect[sr] / (xsect[sr] + bkg[sr].sum())**0.5
    # xsect[sr] = xsect[sr] / (bkg[sr].sum())**0.5
    # print(bkg[sr].sum())
# print(xsect[SRs].max(axis=1))

xsect['ss'] = xsect[SRs].max(axis=1)
# df.loc[(df.Event == 'Dance'),'Event']='Hip-Hop'

xsect.loc[(xsect.ss > 10.), 'ss'] = 10
print(xsect['ss'].min(), xsect['ss'].max())
from matplotlib.tri import Triangulation, TriAnalyzer, UniformTriRefiner

plottri = Triangulation(xsect['x'], xsect['y'])
refiner = UniformTriRefiner(plottri)
tri_refine_zz, zz_refine = refiner.refine_field(xsect['ss'], subdiv=1)

cs = ax.tricontour(tri_refine_zz, zz_refine, levels=[2, 5], colors='white', linewidths=5, zorder=20)
cs = ax.tricontour(tri_refine_zz, zz_refine, levels=[2, 5], colors=['red', 'black'], linewidths=3, zorder=20)
# cs = ax.tricontourf(tri_refine_zz, zz_refine, levels=100, zorder=20)


for sr in SRc:
    # print(xsect[sr])
    # xsect.loc[(xsect.x > 119.7), sr] = 1.8 * xsect[sr]

    xsect[sr] = xsect[sr] / (xsect[sr] + bkg[sr].sum())**0.5

xsect['sc'] = xsect[SRc].max(axis=1)
xsect.loc[(xsect.sc > 10.), 'sc'] = 10
tri_refine_zc, zc_refine = refiner.refine_field(xsect['sc'], subdiv=1)

# cs = ax.tricontour(tri_refine_zc, zc_refine, levels=[2, 5], colors='white', linewidths=5, zorder=20)
cs = ax.tricontour(tri_refine_zc, zc_refine, levels=[2, 5], colors=['magenta', 'lime'], linewidths=3, zorder=10)

# ddt = cs.allsegs[0][1]
# ax.plot(ddt[:,0], ddt[:,1])
ax.set_xlim(100, 140)
ax.set_ylim(0., 140)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(labelsize=18, direction="in", which='both', bottom=True, top=True, left=True, right=True)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4)
plt.colorbar(sc, axc, orientation='vertical')
axc.tick_params(labelsize=18, direction="in", which='both', bottom=False, top=False, left=False, right=True)
axc.tick_params(which='major', length=7)
axc.tick_params(which='minor', length=4)

ax.plot([100, 120, 140], [99.9, 119.1, 119.9], ":", c='gray', linewidth=5)
ax.plot([120, 120, 140], [0., 119.1, 99.9], ":", c='gray', linewidth=5)

ax.text(0, 1.015, r"$e^+ e^- \to {\tilde{\mu}}_{L,R}^+ \tilde{\mu}_{L, R}^-, \tilde{\mu} \to \mu \tilde{\chi}_1^0 $", 
            fontsize=18, color='black',transform=ax.transAxes)

ax.text(101, 130, r"$\sqrt{s} =  240~{\rm GeV}, 5.05~{\rm ab}^{-1}$", fontsize=24, fontdict={"fontfamily": "DejaVu Sans", "weight": 'bold', "fontstyle": "normal"})
ax.text(105, 109, r"$m_{\tilde{\chi}_1^0} = m_{\tilde{\mu}} - m_{\mu}$", c='gray', fontsize=20, rotation=15)
ax.text(125, 123, r"$m_{\tilde{\chi}_1^0} = \left. \sqrt{s} \right/ 2 - m_{\mu} $", c='gray', fontsize=20, rotation=0)

ax.set_xlabel(r"$m_{\tilde{\mu}}~[{\rm GeV}]$", fontsize=30, loc='right')
ax.set_ylabel(r"$m_{\tilde{\chi}_1^0}~[{\rm GeV}]$", fontsize=30, loc='top')
axc.set_ylabel(r"$\sigma(e^+ e^- \to {\tilde{\mu}} \tilde{\mu}^{(\star)})~[{\rm fb}]$", fontsize=30, loc='top')

# s2 = np.loadtxt("2sigma.csv", delimiter=',')
# ax.plot(s2[0:4, 0], s2[0:4, 1], '-', c='magenta', linewidth=3, zorder=10)
# ax.plot(s2[4:, 0], s2[4:, 1], '-', c='magenta', linewidth=3, zorder=10)
# s5 = np.loadtxt("5sigma.csv", delimiter=',')
# ax.plot(s5[0:25, 0], s5[0:25, 1], '-', c='lime', linewidth=3, zorder=10)
# ax.plot(s5[25:, 0], s5[25:, 1], '-', c='lime', linewidth=3, zorder=10)

ax.text(125, 85, "This work:", fontsize=14, fontdict={"fontfamily": "DejaVu Sans", "weight": 'bold', "fontstyle": "normal"}, zorder=20)
ax.plot([128, 132], [80, 80], '-', c='white', linewidth=5, zorder=20)
ax.plot([128, 132], [80, 80], '-', c='red', linewidth=3, zorder=20)
ax.text(133, 80, r"$2\sigma$", fontsize=16, zorder=20, va='center', ha='left')
ax.plot([128, 132], [73, 73], '-', c='white', linewidth=5, zorder=20)
ax.plot([128, 132], [73, 73], '-', c='black', linewidth=3, zorder=20)
ax.text(133, 73, r"$5\sigma$", fontsize=16, zorder=20, va='center', ha='left')

ax.text(125, 55, "arXiv:2203.10580", fontsize=14, fontdict={"fontfamily": "DejaVu Sans", "weight": 'bold', "fontstyle": "normal"}, zorder=20)
ax.plot([128, 132], [50, 50], '-', c='magenta', linewidth=3, zorder=20)
# ax.plot([128, 132], [60, 60], '-', c='red', linewidth=3, zorder=20)
ax.text(133, 50, r"$2\sigma$", fontsize=16, zorder=20, va='center', ha='left')
# ax.plot([128, 132], [53, 53], '-', c='white', linewidth=5, zorder=20)
ax.plot([128, 132], [43, 43], '-', c='lime', linewidth=3, zorder=20)
ax.text(133, 43, r"$5\sigma$", fontsize=16, zorder=20, va='center', ha='left')

# plt.show()
plt.savefig("excl.pdf")