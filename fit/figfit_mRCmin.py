#!/usr/bin/env python3

from matplotlib.colors import LinearSegmentedColormap
from matplotlib.legend_handler import HandlerTuple
from matplotlib.ticker import LogFormatterSciNotation
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
    "font.family": ["serif", "stix", "Times New Roman"],
    "font.size": 20,
    "mathtext.fontset": 'stix',
    "font.serif": ['Computer Modern'],
    "text.latex.preamble": r"\usepackage{amsmath}"
}
rcParams.update(config)
plt.rcParams['axes.formatter.min_exponent'] = 3


class CustomTicker(LogFormatterSciNotation):
    def __call__(self, x, pos=None):
        if x not in [0.1, 1, 10]:
            return LogFormatterSciNotation.__call__(self, x, pos=None)
        else:
            return "{x:g}".format(x=x)

def cords(x1, x2):
    return np.array([x1, x2]).ravel(order="F")

cmap_name = "self"
collist = ["#02004f", "#000b57", "#00165e", "#002166", "#012c6e", "#003875", '#00437d', "#004e85", "#00598c", "#026494", "#006f9c", "#007aa3", "#0085ab",
           "#0090b3", "#009cba", "#27a699", "#4bb178", "#70bc56", "#97c637", "#d2c70d", "#e8bd08", "#ffb300", "#ffbf2a", "#fecc55", "#ffd980", "#ffe6aa", "#ffe6aa"]
cmap = LinearSegmentedColormap.from_list(cmap_name, collist, N=256)

pwd = os.path.abspath(os.path.dirname(__file__))

fig = plt.figure(figsize=(10, 10))
ax = fig.add_axes([0.17, 0.45, 0.81, 0.53])
axr = fig.add_axes([0.17, 0.15, 0.81, 0.3])


base = 80390 
vars = [80360, 80420]

dat_base = pd.read_csv("./Hist/{}_mRCmin.csv".format(base))
def draw(mW, c):
    dat_var = pd.read_csv("./Hist/{}_mRCmin.csv".format(mW))

    xx = cords(dat_var['xlow'], dat_var['xhigh'])
    ax.plot(xx, cords(dat_var['val'], dat_var['val']), '-', lw=2, solid_joinstyle='miter', c=c)
    ax.fill_between(
        xx, 
        cords(dat_var['val'] - dat_var['errminus'], dat_var['val'] - dat_var['errminus']), 
        cords(dat_var['val'] + dat_var['errplus'], dat_var['val'] + dat_var['errplus']), 
        facecolor=c, edgecolor=None, alpha=0.3
    )
    axr.plot(
        xx, 
        cords(dat_var['val'] /dat_base['val'], dat_var['val'] /dat_base['val']),
        solid_joinstyle='miter',
        lw=3,
        color=c,
    )
    axr.plot(
        xx, 
        cords((dat_var['val'] - dat_var['errminus'])/dat_base['val'], (dat_var['val'] - dat_var['errminus'])/dat_base['val']),
        solid_joinstyle='miter',
        color=c,
        alpha=0.7,
        lw=0.1
    )
    axr.plot(
        xx, 
        cords((dat_var['val'] + dat_var['errplus'])/dat_base['val'], (dat_var['val'] + dat_var['errplus'])/dat_base['val']),
        solid_joinstyle='miter',
        color=c,
        alpha=0.7,
        lw=0.1
    )
    axr.fill_between(
        xx, 
        cords((dat_var['val'] - dat_var['errminus'])/dat_base['val'], (dat_var['val'] - dat_var['errminus'])/dat_base['val']),
        cords((dat_var['val'] + dat_var['errplus'])/dat_base['val'], (dat_var['val'] + dat_var['errplus'])/dat_base['val']),
        facecolor=c, edgecolor=None,
        alpha=0.3
    )
draw(base, 'red')
draw(80350, 'blue')
draw(80430, 'green')

ax.set_xlim([dat_base['xlow'].min(), dat_base['xhigh'].max()])
ax.set_ylim(0., 50000)
axr.set_xlim([dat_base['xlow'].min(), dat_base['xhigh'].max()])
axr.set_ylim(0.95, 1.05)

# ax.set_xticks(np.linspace(76, 84, 5))
ax.set_xticks([])
axr.set_xticks(np.linspace(76, 84, 5))
axr.set_yticks(np.linspace(0.96, 1.04, 5))
# ax.xaxis.set_minor_locator(FixedLocator(np.linspace(75, 85, 51)))
axr.xaxis.set_minor_locator(FixedLocator(np.linspace(75, 85, 51)))
ax.yaxis.set_minor_locator(AutoMinorLocator())
axr.yaxis.set_minor_locator(AutoMinorLocator())


ax.tick_params(
    which='both',
    direction="in",
    labelsize=26,
    left=True,
    right=True,
    bottom=True,
    top=True
)
ax.tick_params(which="major", length=10, width=1.2)
ax.tick_params(which="minor", length=4, width=1.2)
axr.tick_params(
    which='both',
    direction="in",
    labelsize=26,
    left=True,
    right=True,
    bottom=True,
    top=True
)
axr.tick_params(which="major", length=10, width=1.2)
axr.tick_params(which="minor", length=4, width=1.2)

axr.set_xlabel(r"$m_{\rm RC}^{\rm min}~[{\rm GeV}]$", fontsize=40, loc='right')
ax.set_ylabel(r"${\rm Events~\left/~0.2~GeV \right.}$", fontsize=40, loc='top')
axr.set_ylabel(r"$\rm Ratio$", fontsize=40, loc='center')

from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
hdl = [
    (
        Rectangle((0, 0), 1, 1, facecolor="red", edgecolor=None, alpha=0.3),
        Line2D([0], [0], color='red', lw=4) 
    ),
    (
        Rectangle((0, 0), 1, 1, facecolor='blue', edgecolor=None, alpha=0.3),
        Line2D([0], [0], color='blue', lw=4) 
    ),
    (
        Rectangle((0, 0), 1, 1, facecolor='green', edgecolor=None, alpha=0.3),
        Line2D([0], [0], color='green', lw=4) 
    )
]
labs = [
    r"$m_W = 80390~{\rm MeV}$",
    r"$m_W = 80360~{\rm MeV}$",
    r"$m_W = 80420~{\rm MeV}$",
]
ax.legend(hdl, labs, framealpha=0, loc='best', fontsize=30)

# plt.show()
plt.savefig("mRCmin_PRD.pdf")

