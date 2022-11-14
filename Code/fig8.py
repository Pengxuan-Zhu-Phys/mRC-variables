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
    "font.family": ["serif", "Times New Roman"],
    "font.size": 20,
    "mathtext.fontset": 'stix',
    "font.serif": ['Computer Modern'],
    "text.latex.preamble": r"\usepackage{amsmath}"
}
rcParams.update(config)
plt.rcParams['axes.formatter.min_exponent'] = 2


class CustomTicker(LogFormatterSciNotation):
    def __call__(self, x, pos=None):
        if x not in [0.1, 1, 10]:
            return LogFormatterSciNotation.__call__(self, x, pos=None)
        else:
            return "{x:g}".format(x=x)


cmap_name = "self"
collist = ["#02004f", "#000b57", "#00165e", "#002166", "#012c6e", "#003875", '#00437d', "#004e85", "#00598c", "#026494", "#006f9c", "#007aa3", "#0085ab",
           "#0090b3", "#009cba", "#27a699", "#4bb178", "#70bc56", "#97c637", "#d2c70d", "#e8bd08", "#ffb300", "#ffbf2a", "#fecc55", "#ffd980", "#ffe6aa", "#ffe6aa"]
cmap = LinearSegmentedColormap.from_list(cmap_name, collist, N=256)

pwd = os.path.abspath(os.path.dirname(__file__))

fig = plt.figure(figsize=(10, 4))
ax  = fig.add_axes([0.096, 0.16, 0.89, 0.8])

SRs = ['SRH01', "SRH02", "SRM01", "SRM02", "SRM03", "SRM04", "SRM05", "SRL01", "SRL02", "SRL03", "SRL04", "SRL05"]
# Plot the SRs for signals
sgn = pd.read_csv(f"Data/sgn.csv")
mchi10 = [0, 40, 80, 110]
lincols = ['#009051', '#f44336', '#1A237E', '#795548']
for idx, row in sgn.iterrows():
    ax.plot(
        np.array([np.linspace(-0.5, 10.5, 12), np.linspace(0.5, 11.5, 12)]).ravel(order="F"),
        np.array([row[SRs], row[SRs]]).ravel(order="F"),
        '-',
        c='w',
        linewidth=6,
        solid_joinstyle="miter"
    )
    ax.plot(
        np.array([np.linspace(-0.5, 10.5, 12), np.linspace(0.5, 11.5, 12)]).ravel(order="F"),
        np.array([row[SRs], row[SRs]]).ravel(order="F"),
        '--',
        c=lincols[idx],
        linewidth=4.5,
        solid_joinstyle="miter"
    )

# print(sgn.loc[sgn['sgn'] == "150_0"])

# for ii in range(len(mchi10)):
#     chi = mchi10[ii]
#     ax.plot(
#         xx, yy, '-',
#         c='w',
#         linewidth=6,
#         solid_joinstyle="miter"
#     )
# for ii in range(len(mchi10)):
#     chi = mchi10[ii]
#     sgn = pd.read_csv(f"Data/Draw/{mode}/115_{chi}.csv")
#     xx = np.array([sgn['xlow'], sgn['xhigh']]).ravel(order='F')
#     yy = np.array([sgn['val'], sgn['val']]).ravel(order='F')
#     ax.plot(
#         xx, yy, '--',
#         c=lincols[ii],
#         linewidth=4.5,
#         solid_joinstyle="miter"
#     )


bkg = pd.read_csv(f"Data/bkg.csv")
# print(bkg)


colors = ['#FFC078', '#A9E34B', '#22B8CF', '#E599F7', '#3F51B5', '#ffc107', '#B197FC']
yy0 = np.zeros(12)
bkgs = []
for idx, row in bkg.iterrows():
    bkgs.append(row['bkg'])
    for ii in range(len(SRs)):
        ax.fill(
            [ii-0.5, ii+0.5, ii+0.5, ii-0.5],
            [yy0[ii], yy0[ii], yy0[ii]+row[SRs[ii]], yy0[ii]+row[SRs[ii]]],
            c=colors[idx]
        )
        yy0[ii] += row[SRs[ii]]
# print(bkgs)
eventlabels={
    'zh': r'$ZH$', 
    'wmv': r'$W\mu\nu$', 
    'zz': r'$ZZ$', 
    '4l': r"$\mu\mu \nu \nu$", 
    'zll': r"$Z\ell\ell$", 
    'tata': r"$\tau\tau$", 
    'ww': r"$WW$"
}
ax.set_xticks(
    np.linspace(0.3, 11.3, 12),
    labels=['SRH-01', "SRH-02", "SRM-01", "SRM-02", "SRM-03", "SRM-04", "SRM-05", "SRL-01", "SRL-02", "SRL-03", "SRL-04", "SRL-05"],
    rotation=26,
    ha='right',
    va='top'
)
ax.tick_params(
    which='both',
    direction="in",
    labelsize=16,
    left=True,
    right=True,
    bottom=False,
    top=False
)
ax.tick_params(which="major", length=9, width=1.2)
ax.tick_params(which="minor", length=4, width=1.2)
# Plot the SM total in SRs
ax.plot(
    np.array([np.linspace(-0.5, 10.5, 12), np.linspace(0.5, 11.5, 12)]).ravel(order="F"),
    np.array([yy0, yy0]).ravel(order='F'),
    "-",
    c='black',
    linewidth=2
)
# print(yy0)
ax.plot([1.5, 1.5], [0., yy0[1]], ':', linewidth=2, c='gray', alpha=0.7)
ax.plot([6.5, 6.5], [0., yy0[7]], ':', linewidth=2, c='gray', alpha=0.7)

for ii in range(len(SRs)): 
    ax.plot([ii+0.5, ii+0.5], [0., 1.e7], '-', linewidth=0.5, c='gray', alpha=0.2)
    ax.plot([ii/12., ii/12], [0., 0.02], '-', linewidth=1.5, c='black', transform=ax.transAxes)
    ax.plot([ii/12., ii/12], [0.98, 1.0], '-', linewidth=1.5, c='black', transform=ax.transAxes)
    err = yy0[ii] ** 0.5 
    ax.fill(
        [ii-0.5, ii+0.5, ii+0.5, ii-0.5],
        [yy0[ii]-err, yy0[ii]-err, yy0[ii]+err, yy0[ii]+err],
        c='gray',
        alpha=0.3,
        hatch='/'
    )

ax.set_xlim(-0.5, 11.5)
ax.set_ylim(5, 5.e7)
ax.set_yscale("log")
ax.set_ylabel(r"Events", loc="top", fontsize=30)
ax.yaxis.set_minor_locator(matplotlib.ticker.LogLocator(numticks=999, subs="auto"))
ax.yaxis.set_major_locator(matplotlib.ticker.LogLocator(numticks=999))

from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
handles = []
for cc in lincols:
    handles.append((
        Line2D([0], [0], color="w", lw=6),
        Line2D([0], [0], color=cc, linestyle='--', lw=4.5)
    ))

hdl = [(
    Rectangle((0,0),1,1,facecolor=None,ec='black', hatch='///', fill=False),  
    Rectangle((0,0),1,1,color="gray",ec=None, alpha=0.3),  
    Line2D([0], [0], color="black", lw=3)
    )]
labs = [
    r"$m_{\tilde{\mu}} = 115~{\rm GeV}, m_{\tilde{\chi}_1^0} = 0~{\rm GeV}$", 
    r"$m_{\tilde{\mu}} = 115~{\rm GeV}, m_{\tilde{\chi}_1^0} = 40~{\rm GeV}$", 
    r"$m_{\tilde{\mu}} = 115~{\rm GeV}, m_{\tilde{\chi}_1^0} = 80~{\rm GeV}$", 
    r"$m_{\tilde{\mu}} = 115~{\rm GeV}, m_{\tilde{\chi}_1^0} = 110~{\rm GeV}$", 
    'SM Total'
    ]
labs += list(eventlabels.values())
handles  += hdl

# ax.legend(handler[::-1], label[::-1], loc='upper right', ncol=2, framealpha=0)
handles += [Rectangle((0,0),1,1,color=c,ec=None) for c in colors]
ax.legend(handles, labs, framealpha=0, ncol=3, loc="best", fontsize='xx-small')
# ax.legend()


# plt.show()
plt.savefig("Figure/SRs.pdf")