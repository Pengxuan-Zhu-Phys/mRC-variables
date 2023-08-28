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

fig = plt.figure(figsize=(10, 8))
ax = fig.add_axes([0.15, 0.16, 0.818, 0.8])

# mode = "dMRC"
# ax.set_xlabel(
#     r"$m_{\rm RC}^{\rm max} - m_{\rm LSP}^{\rm max}~[{\rm GeV}]$", fontsize=48, loc='right')

# mode = "mRCmax"
# ax.set_xlabel(
#     r"$m_{\rm RC}^{\rm max}~[{\rm GeV}]$", fontsize=48, loc='right')

# mode = "mRCmin"
# ax.set_xlabel(
#     r"$m_{\rm RC}^{\rm min}~[{\rm GeV}]$", fontsize=48, loc='right')

mode = "mLSPmax"
ax.set_xlabel(
    r"$m_{\rm LSP}^{\rm max}~[{\rm GeV}]$", fontsize=48, loc='right')

bkgname = ['Zh', 'Wmuv', 'ZZ', 'mumuvv', 'Zll', 'TauTau', 'WW']
eventlabels={
    'zh': r'$ZH$', 
    'wmv': r'$W\mu\nu$', 
    'zz': r'$ZZ$', 
    '4l': r"$\mu\mu \nu \nu$", 
    'zll': r"$Z\ell\ell$", 
    'tata': r"$\tau\tau$", 
    'ww': r"$WW$"
}
mchi10 = [0, 40, 80, 110]
lincols = ['#009051', '#f44336', '#1A237E', '#795548']
for ii in range(len(mchi10)):
    chi = mchi10[ii]
    sgn = pd.read_csv(f"Data/Draw/{mode}/115_{chi}.csv")
    xx = np.array([sgn['xlow'], sgn['xhigh']]).ravel(order='F')
    yy = np.array([sgn['val'], sgn['val']]).ravel(order='F')
    ax.plot(
        xx, yy, '-',
        c='w',
        linewidth=6,
        solid_joinstyle="miter"
    )

for ii in range(len(mchi10)):
    chi = mchi10[ii]
    sgn = pd.read_csv(f"Data/Draw/{mode}/115_{chi}.csv")
    xx = np.array([sgn['xlow'], sgn['xhigh']]).ravel(order='F')
    yy = np.array([sgn['val'], sgn['val']]).ravel(order='F')
    ax.plot(
        xx, yy, linestyle=(0, (5, 1)),
        c=lincols[ii],
        linewidth=4,
        solid_joinstyle="miter"
    )



# colors = ['#FFC078', '#A9E34B', '#22B8CF', '#E599F7', '#B197FC', '#EF9A9A', '#00bcd4']
colors = ['#FFC078', '#A9E34B', '#22B8CF', '#E599F7', '#3F51B5', '#ffc107', '#B197FC']


yy0 = np.zeros(24)
sumerr2 = np.zeros(24)




for ii in range(len(bkgname)):
    bkg = bkgname[ii]
    dat = pd.read_csv(f"Data/Draw/{mode}/{bkg}.csv")
    xx = np.array([dat['xlow'], dat['xhigh']]).ravel(order='F')
    yy1 = yy0 + dat['val']
    yyl = np.array([yy0, yy0]).ravel(order='F')
    yyh = np.array([yy1, yy1]).ravel(order='F')
    # ax.fill_between(xx,yyl,yyh, fc="#00BBCC", edgecolor=None)
    ax.fill_between(xx, yyl, yyh, label=bkg, fc=colors[ii])
    yy0 = yy1
    sumerr2 += dat['errminus']**2

sumerr2 += yy0/(5)
err = (sumerr2**0.5)/2


# errl = np.array([errlow,errlow]).ravel(order='F')
# errh = np.array([errhigh,errhigh]).ravel(order='F')

ax.fill_between(
    xx, np.array([yy0-err, yy0-err]).ravel(order="F"), np.array([yy0 +
                                                                 err, yy0+err]).ravel(order="F"),
    step="mid", edgecolor="black", facecolor='w',
    label="SM", hatch='///', zorder=10
)
ax.fill_between(
    xx, np.array([yy0-err, yy0-err]).ravel(order="F"), np.array([yy0 +
                                                                 err, yy0+err]).ravel(order="F"),
    step="mid", edgecolor="gray",
    alpha=0.3, label="SM", zorder=11
)
ax.plot(xx, np.array([yy0, yy0]).ravel(
    order="F"), c='black', linewidth=0.8, label="SM", zorder=11)
ax.set_xlim(0, 120)
# ax.set_ylim(0, 5.e4)

ax.set_ylim(20, 5.e5)
ax.set_yscale("log")

ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.xaxis.set_minor_locator(AutoMinorLocator())
# ax.yaxis.set_major_locator(MaxNLocator(6))
ax.yaxis.set_minor_locator(matplotlib.ticker.LogLocator(numticks=999, subs="auto"))
ax.yaxis.set_major_locator(matplotlib.ticker.LogLocator(numticks=999))
ax.yaxis.set_major_formatter(CustomTicker())

ax.tick_params(
    which='both',
    direction="in",
    labelsize=32,
    left=True,
    right=True,
    bottom=True,
    top=True
)
ax.tick_params(which="major", length=14, width=2.0)
ax.tick_params(which="minor", length=8, width=1.2)

ax.set_ylabel(r"Events $\left/ 5~{\rm GeV}\right.$", fontsize=48, loc='top')

handler, label = ax.get_legend_handles_labels()
print(handler, label)
smtt = (handler[7], handler[8])
handler.pop(-1)
handler.pop(-1)
handler.append(smtt)
label.pop(-1)    


# ax1.legend(handler1, label1, loc='upper left', framealpha=0.001)
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

handles = []
for cc in lincols:
    handles.append((
        Line2D([0], [0], color="w", lw=3),
        Line2D([0], [0], color=cc, linestyle=(0, (5, 1)), lw=2)
    ))

hdl = [(
    Rectangle((0,0),1,1,facecolor=None,ec='black', hatch='///', fill=False),  
    Rectangle((0,0),1,1,color="gray",ec=None, alpha=0.3),  
    Line2D([0], [0], color="black", lw=1)
    )]
labs = [
    r"$m_{\tilde{\chi}_1^0} = {\bf 0~GeV}$", 
    r"$m_{\tilde{\chi}_1^0} = {\bf 40~GeV}$", 
    r"$m_{\tilde{\chi}_1^0} = {\bf 80~GeV}$", 
    r"$m_{\tilde{\chi}_1^0} = {\bf 110~GeV}$", 
    'SM Total'
    ]
labs += list(eventlabels.values())
handles  += hdl

# ax.legend(handler[::-1], label[::-1], loc='upper right', ncol=2, framealpha=0)
handles += [Rectangle((0,0),1,1,color=c,ec=None) for c in colors]
# ax.legend(handles, labs, framealpha=0, ncol=3, loc="best", fontsize='small')
ax.legend(handles, labs, framealpha=0, ncol=3, loc="upper right", fontsize=18)


# plt.legend([(ll1, ll2)], ['SM total'])
# plt.plot(xx,yy0,'-')
# plt.show()
plt.savefig(f'Figure/cut{mode}.pdf')
