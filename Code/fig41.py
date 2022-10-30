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
ax = fig.add_axes([0.11, 0.08, 0.865, 0.9])

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

# smu2hist = np.histogram(
#     smu115_0[chan], xbin, range=xrange
# )
# smuxx2 = np.array([np.array(smu2hist[1][:-1]),
#                    np.array(smu2hist[1][1:])]).ravel(order="F")
# smuyy2 = np.array(
#     [np.array(smu2hist[0]), np.array(smu2hist[0])]).ravel(order="F")
# smuyy2 = smuyy2 / np.max(smuyy2)

smu3hist0 = np.histogram(
    smu115_40['mRCmin'], xbin, range=xrange
)
smuxx3 = np.array([np.array(smu3hist0[1][:-1]),
                   np.array(smu3hist0[1][1:])]).ravel(order="F")
smuyy3 = np.array(
    [np.array(smu3hist0[0]), np.array(smu3hist0[0])]).ravel(order="F")
smuyy3 = 0.6 * smuyy3 / np.max(smuyy3)

smu3hist40 = np.histogram(
    smu115_40['mRCmin40'], xbin, range=xrange
)
smuxx340 = np.array([np.array(smu3hist40[1][:-1]),
                   np.array(smu3hist40[1][1:])]).ravel(order="F")
smuyy340 = np.array(
    [np.array(smu3hist40[0]), np.array(smu3hist40[0])]).ravel(order="F")
smuyy340 = 0.6 *  smuyy340 / np.max(smuyy340)

smu4hist = np.histogram(
    smu115_80['mRCmin'], xbin, range=xrange
)
smuxx4 = np.array([np.array(smu4hist[1][:-1]),
                   np.array(smu4hist[1][1:])]).ravel(order="F")
smuyy4 = np.array(
    [np.array(smu4hist[0]), np.array(smu4hist[0])]).ravel(order="F")
smuyy4 = 0.8* smuyy4 / np.max(smuyy4)

smu4hist40 = np.histogram(
    smu115_80['mRCmin40'], xbin, range=xrange
)
smuxx440 = np.array([np.array(smu4hist40[1][:-1]),
                   np.array(smu4hist40[1][1:])]).ravel(order="F")
smuyy440 = np.array(
    [np.array(smu4hist40[0]), np.array(smu4hist40[0])]).ravel(order="F")
smuyy440 = 0.8*  smuyy440 / np.max(smuyy440)

smu4hist80 = np.histogram(
    smu115_80['mRCmin80'], xbin, range=xrange
)
smuxx480 = np.array([np.array(smu4hist80[1][:-1]),
                   np.array(smu4hist80[1][1:])]).ravel(order="F")
smuyy480 = np.array(
    [np.array(smu4hist80[0]), np.array(smu4hist80[0])]).ravel(order="F")
smuyy480 = 0.8* smuyy480 / np.max(smuyy480)

smu5hist = np.histogram(
    smu115_110['mRCmin'], xbin, range=xrange
)
smuxx5 = np.array([np.array(smu5hist[1][:-1]),
                   np.array(smu5hist[1][1:])]).ravel(order="F")
smuyy5 = np.array(
    [np.array(smu5hist[0]), np.array(smu5hist[0])]).ravel(order="F")
smuyy5 = 1.0 * smuyy5 / np.max(smuyy5)

smu5hist40 = np.histogram(
    smu115_110['mRCmin40'], xbin, range=xrange
)
smuxx540 = np.array([np.array(smu5hist40[1][:-1]),
                   np.array(smu5hist40[1][1:])]).ravel(order="F")
smuyy540 = np.array(
    [np.array(smu5hist40[0]), np.array(smu5hist40[0])]).ravel(order="F")
smuyy540 = 1.0 * smuyy540 / np.max(smuyy540)

smu5hist80 = np.histogram(
    smu115_110['mRCmin80'], xbin, range=xrange
)
smuxx580 = np.array([np.array(smu5hist80[1][:-1]),
                   np.array(smu5hist80[1][1:])]).ravel(order="F")
smuyy580 = np.array(
    [np.array(smu5hist80[0]), np.array(smu5hist80[0])]).ravel(order="F")
smuyy580 = 1.0 * smuyy580 / np.max(smuyy580)

smu5hist110 = np.histogram(
    smu115_110['mRCmin110'], xbin, range=xrange
)
smuxx5110 = np.array([np.array(smu5hist110[1][:-1]),
                   np.array(smu5hist110[1][1:])]).ravel(order="F")
smuyy5110 = np.array(
    [np.array(smu5hist110[0]), np.array(smu5hist110[0])]).ravel(order="F")
smuyy5110 = 1.0 * smuyy5110 / np.max(smuyy5110)


lstyle = {"linewidth": 1.8, "solid_joinstyle": "miter"}


# ax.plot(wwxx, wwyy, '-', c="#c1272d", **lstyle)
# ax.plot(zzxx, zzyy, '-', c="#0000a7", **lstyle)
# ax.plot(smuxx1, smuyy1, '-', c="#00b000", **lstyle)
# ax.plot(smuxx2, smuyy2, '-', c="#008176", **lstyle)

ax.plot(smuxx3, smuyy3, c="#0000a7", linestyle="solid", **lstyle)
ax.text(10, 0.1, r"$m_{\rm RC}^{\rm min}$", c="#0000a7", fontsize=14, fontfamily="stix", rotation=30, 
        verticalalignment='bottom', horizontalalignment="left")
ax.plot(smuxx340, smuyy340, c="#0000a7", linestyle="dotted", **lstyle)
ax.text(60, 0.1, r"$m_{\rm RC}^{\rm min}\left(40\right)$", c="#0000a7", fontsize=14, fontfamily="stix", rotation=12, 
        verticalalignment='bottom', horizontalalignment="left")

ax.plot(smuxx4, smuyy4, c="#00b000", linestyle="solid", **lstyle)
ax.text(30, 0.76, r"$m_{\rm RC}^{\rm min}$", c="#00b000", fontsize=14, fontfamily="stix", rotation=50, 
        verticalalignment='top', horizontalalignment="left")
ax.plot(smuxx440, smuyy440, c="#00b000", linestyle="dotted", **lstyle)
ax.text(60, 0.80, r"$m_{\rm RC}^{\rm min}\left(40\right)$", c="#00b000", fontsize=14, fontfamily="stix", rotation=0, 
        verticalalignment='bottom', horizontalalignment="left")
ax.plot(smuxx480, smuyy480, c="#00b000", linestyle="dashed", **lstyle)
ax.text(116, 0.70, r"$m_{\rm RC}^{\rm min}\left(80\right)$", c="#00b000", fontsize=14, fontfamily="stix", rotation=90, 
        verticalalignment='top', horizontalalignment="left")


ax.plot(smuxx5, smuyy5, c="#c1272d", linestyle="solid", **lstyle)
ax.text(24, 1.01, r"$m_{\rm RC}^{\rm min}$", c="#c1272d", fontsize=14, fontfamily="stix", rotation=-70, 
        verticalalignment='top', horizontalalignment="left")
ax.plot(smuxx540, smuyy540,  c="#c1272d",linestyle="dotted", **lstyle)
ax.text(43, 1.01, r"$m_{\rm RC}^{\rm min}\left(40\right)$", c="#c1272d", fontsize=14, fontfamily="stix", rotation=-80, 
        verticalalignment='top', horizontalalignment="left")
ax.plot(smuxx580, smuyy580, c="#c1272d",linestyle="dashed", **lstyle)
ax.text(84, 1.01, r"$m_{\rm RC}^{\rm min}\left(80\right)$", c="#c1272d", fontsize=14, fontfamily="stix", rotation=-87, 
        verticalalignment='top', horizontalalignment="left")
ax.plot(smuxx5110, smuyy5110, c="#c1272d",linestyle="dashdot", **lstyle)
ax.text(107, 1.01, r"$m_{\rm RC}^{\rm min}\left(100\right)$", c="#c1272d", fontsize=14, fontfamily="stix", rotation=-90, 
        verticalalignment='top', horizontalalignment="left")

ax.text(15, 1.1, r"$\tilde{\mu}\tilde{\mu}$", fontsize=16, c="#0000a7",
        verticalalignment='bottom', horizontalalignment="left")
ax.text(15, 1.01, r"(115, 40)", c="#0000a7", fontsize=16,
        verticalalignment='bottom', horizontalalignment="left")
ax.plot([5, 12], [1.13, 1.13], c="#0000a7", linestyle="solid", **lstyle)
ax.plot([5, 12], [1.08, 1.08], c="#0000a7", linestyle="dotted", **lstyle)

ax.text(55, 1.1, r"$\tilde{\mu}\tilde{\mu}$", fontsize=16, c="#00b000",
        verticalalignment='bottom', horizontalalignment="left")
ax.text(55, 1.01, r"(115, 80)", c="#00b000", fontsize=16,
        verticalalignment='bottom', horizontalalignment="left")
ax.plot([45, 52], [1.155, 1.155], c="#00b000", linestyle="solid", **lstyle)
ax.plot([45, 52], [1.105, 1.105], c="#00b000", linestyle="dotted", **lstyle)
ax.plot([45, 52], [1.055, 1.055], c="#00b000", linestyle="dashed", **lstyle)

ax.text(95, 1.1, r"$\tilde{\mu}\tilde{\mu}$", fontsize=16, c="#c1272d",
        verticalalignment='bottom', horizontalalignment="left")
ax.text(95, 1.01, r"(115, 110)", c="#c1272d", fontsize=16,
        verticalalignment='bottom', horizontalalignment="left")
ax.plot([85, 92], [1.18, 1.18],  c="#c1272d", linestyle="solid", **lstyle)
ax.plot([85, 92], [1.13, 1.13], c="#c1272d",linestyle="dotted", **lstyle)
ax.plot([85, 92], [1.08, 1.08],  c="#c1272d",linestyle="dashed", **lstyle)
ax.plot([85, 92], [1.03, 1.03], c="#c1272d",linestyle="dashdot", **lstyle)


# ax.set_xlabel(r"$m_{\rm RC}^{\rm min}~[{\rm GeV}]$", fontsize=20, loc='right')



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

plt.savefig("Figure/mRCmins.pdf")
plt.show()
