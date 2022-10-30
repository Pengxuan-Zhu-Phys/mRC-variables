#!/usr/bin/env python3

import numpy as np
from matplotlib.ticker import FixedLocator, AutoMinorLocator
import pandas as pd
import os
import sys
from shutil import which
import matplotlib.pyplot as plt
import json


from matplotlib import rc, rcParams
from matplotlib.image import NonUniformImage
from matplotlib import cm, ticker
from matplotlib.font_manager import FontProperties
import matplotlib
config = {
    "font.family": ["serif", "Times New Roman"],
    "font.size": 20,
    "mathtext.fontset": 'stix',
    "font.serif": ['Computer Modern'],
    "text.latex.preamble": r"\usepackage{amsmath}"
}
rcParams.update(config)

pwd = os.path.abspath(os.path.dirname(__file__))

fig = plt.figure(figsize=(5, 4))
ax = fig.add_axes([0.11, 0.18, 0.865, 0.8])


xbin = 250
xrange = (0, 125)

chan = "mRecoil"


ax.set_ylabel(r"A.U.", fontsize=20, loc='center')

ww = pd.read_csv(f"{pwd}/Data/WW_reportDF.csv")
zz = pd.read_csv(f"{pwd}/Data/ZZ_reportDF.csv")
smu50_0 = pd.read_csv(f"{pwd}/Data/smu50_chi0_reportDF.csv")
smu115_0 = pd.read_csv(f"{pwd}/Data/smu115_chi0_reportDF.csv")
smu115_40 = pd.read_csv(f"{pwd}/Data/smu115_chi40_reportDF_2.csv")
smu115_80 = pd.read_csv(f"{pwd}/Data/smu115_chi80_reportDF_2.csv")
smu115_110 = pd.read_csv(f"{pwd}/Data/smu115_chi110_reportDF_2.csv")

# wwhist = np.histogram(
#     ww[chan], xbin, range=xrange
# )
# wwxx = np.array([np.array(wwhist[1][:-1]),
#                  np.array(wwhist[1][1:])]).ravel(order="F")
# wwyy = np.array([np.array(wwhist[0]), np.array(wwhist[0])]).ravel(order="F")
# wwyy = wwyy / np.max(wwyy)

# zzhist = np.histogram(
#     zz[chan], xbin, range=xrange
# )
# zzxx = np.array([np.array(zzhist[1][:-1]),
#                  np.array(zzhist[1][1:])]).ravel(order="F")
# zzyy = np.array([np.array(zzhist[0]), np.array(zzhist[0])]).ravel(order="F")
# zzyy = zzyy / np.max(zzyy)

# smu1hist = np.histogram(
#     smu50_0[chan], xbin, range=xrange
# )
# smuxx1 = np.array([np.array(smu1hist[1][:-1]),
#                    np.array(smu1hist[1][1:])]).ravel(order="F")
# smuyy1 = np.array(
#     [np.array(smu1hist[0]), np.array(smu1hist[0])]).ravel(order="F")
# smuyy1 = smuyy1 / np.max(smuyy1)

smu2hist = np.histogram(
    smu115_0["mLSPmax"], xbin, range=xrange
)
smuxx2 = np.array([np.array(smu2hist[1][:-1]),
                   np.array(smu2hist[1][1:])]).ravel(order="F")
smuyy2 = np.array(
    [np.array(smu2hist[0]), np.array(smu2hist[0])]).ravel(order="F")
smuyy2 = smuyy2 / np.max(smuyy2)

smu3hist0 = np.histogram(
    smu115_40['mLSPmax'], xbin, range=xrange
)
smuxx3 = np.array([np.array(smu3hist0[1][:-1]),
                   np.array(smu3hist0[1][1:])]).ravel(order="F")
smuyy3 = np.array(
    [np.array(smu3hist0[0]), np.array(smu3hist0[0])]).ravel(order="F")
smuyy3 = smuyy3 / np.max(smuyy3)

smu4hist = np.histogram(
    smu115_80['mLSPmax'], xbin, range=xrange
)
smuxx4 = np.array([np.array(smu4hist[1][:-1]),
                   np.array(smu4hist[1][1:])]).ravel(order="F")
smuyy4 = np.array(
    [np.array(smu4hist[0]), np.array(smu4hist[0])]).ravel(order="F")
smuyy4 = smuyy4 / np.max(smuyy4)

smu5hist = np.histogram(
    smu115_110['mLSPmax'], xbin, range=xrange
)
smuxx5 = np.array([np.array(smu5hist[1][:-1]),
                   np.array(smu5hist[1][1:])]).ravel(order="F")
smuyy5 = np.array(
    [np.array(smu5hist[0]), np.array(smu5hist[0])]).ravel(order="F")
smuyy5 =  smuyy5 / np.max(smuyy5)



lstyle = {"linewidth": 1.8, "solid_joinstyle": "miter"}


# ax.plot(wwxx, wwyy, '-', c="#c1272d", **lstyle)
# ax.plot(zzxx, zzyy, '-', c="#0000a7", **lstyle)
# ax.plot(smuxx1, smuyy1, '-', c="#00b000", **lstyle)
ax.plot(smuxx2, smuyy2, '-', c="#008176", **lstyle)
ax.text(30, 1.1, r"$\tilde{\mu}\tilde{\mu}$", fontsize=16, c="#008176",
        verticalalignment='bottom', horizontalalignment="center")
ax.text(30, 1.01, r"(115, 0)", c="#008176", fontsize=16,
        verticalalignment='bottom', horizontalalignment="center")

ax.plot(smuxx3, smuyy3, c="#0000a7", linestyle="solid", **lstyle)
ax.text(52, 1.1, r"$\tilde{\mu}\tilde{\mu}$", fontsize=16, c="#0000a7",
        verticalalignment='bottom', horizontalalignment="center")
ax.text(52, 1.01, r"(115, 40)", c="#0000a7", fontsize=16,
        verticalalignment='bottom', horizontalalignment="center")

ax.plot(smuxx4, smuyy4, c="#00b000", linestyle="solid", **lstyle)
ax.text(80, 1.1, r"$\tilde{\mu}\tilde{\mu}$", fontsize=16, c="#00b000",
        verticalalignment='bottom', horizontalalignment="center")
ax.text(80, 1.01, r"(115, 80)", c="#00b000", fontsize=16,
        verticalalignment='bottom', horizontalalignment="center")

ax.plot(smuxx5, smuyy5, c="#c1272d", linestyle="solid", **lstyle)
ax.text(110, 1.1, r"$\tilde{\mu}\tilde{\mu}$", fontsize=16, c="#c1272d",
        verticalalignment='bottom', horizontalalignment="center")
ax.text(110, 1.01, r"(115, 110)", c="#c1272d", fontsize=16,
        verticalalignment='bottom', horizontalalignment="center")


ax.set_xlabel(r"$m_{\rm LSP}^{\rm max}~[{\rm GeV}]$", fontsize=20, loc='right')



ax.text(2, 1.2, r"$\sqrt{\bf s}$=240 GeV", fontsize=14, fontdict={"fontfamily": "DejaVu Sans",
                                                              "weight": "bold",
                                                              "fontstyle": "normal", "fontvariant": "small-caps",
                                                              "fontstretch": "extra-condensed"})

ax.tick_params(
    which='both',
    direction="in",
    left=True,
    right=True,
    bottom=True,
    top=True
)
ax.tick_params(which="major", length=6, width=1.2)
ax.tick_params(which="minor", length=3.6, width=1.2)


ax.xaxis.set_major_locator(FixedLocator(np.linspace(0, 120, 7)))
ax.xaxis.set_minor_locator(FixedLocator(np.linspace(0, 125, 26)))

ax.yaxis.set_minor_locator(FixedLocator(np.linspace(0, 1, 26)))
ax.yaxis.set_major_locator(FixedLocator(np.linspace(0, 1, 6)))
ax.set_yticklabels([])


ax.set_xlim(xrange)
ax.set_ylim(0, 1.3)

plt.savefig("Figure/mLSPmaxs.pdf")
plt.show()
