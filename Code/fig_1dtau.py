#!/usr/bin/env python3

from cProfile import label
import numpy as np
from matplotlib.ticker import LinearLocator, FixedLocator, AutoMinorLocator, MaxNLocator
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
plt.rcParams['axes.formatter.min_exponent'] = 2
from matplotlib.ticker import LogFormatterSciNotation

class CustomTicker(LogFormatterSciNotation):
    def __call__(self, x, pos=None):
        if x not in [0.1,1,10]:
            return LogFormatterSciNotation.__call__(self,x, pos=None)
        else:
            return "{x:g}".format(x=x)

pwd = os.path.abspath(os.path.dirname(__file__))


from matplotlib.colors import LinearSegmentedColormap
cmap_name = "self"
collist = ["#02004f", "#000b57", "#00165e", "#002166", "#012c6e", "#003875", '#00437d', "#004e85", "#00598c", "#026494", "#006f9c", "#007aa3", "#0085ab", "#0090b3", "#009cba", "#27a699", "#4bb178", "#70bc56", "#97c637", "#d2c70d", "#e8bd08", "#ffb300", "#ffbf2a", "#fecc55", "#ffd980", "#ffe6aa", "#ffe6aa"]
cmap = LinearSegmentedColormap.from_list(cmap_name, collist, N=256)


fig = plt.figure(figsize=(10, 9))
ax  = fig.add_axes([0.198, 0.21, 0.657, 0.73])

data = pd.read_csv(f"{pwd}/Data/Tau_reportDF.csv")
xsect = 1886.47 * 0.1739 * 0.1739
ax.text(0.02, 1.001, r"$\tau^+\tau^- \to \mu^+\mu^- \nu\nu\nu\nu $", fontsize=26, ha='left', va='bottom', transform=ax.transAxes)
nmc = 113317


ax.hist(data['mRCmin'], 120, range=(0, 5), color='#2e66ff', histtype="stepfilled", alpha=0.5, orientation='vertical', zorder=10)
ax.hist(data['mRCmin'], 120, range=(0, 5), color='#231aa5', histtype="step", alpha=1, orientation='vertical', zorder=10)

ax.hist(data['mRCmax'], 120, range=(0, 5), color='#012c6e', histtype="stepfilled", alpha=0.5, orientation='vertical', zorder=10)
ax.hist(data['mRCmax'], 120, range=(0, 5), color='#012c6e', histtype="step", alpha=1, orientation='vertical', zorder=10)

ax.set_xlim(0, 5)

plt.show()
