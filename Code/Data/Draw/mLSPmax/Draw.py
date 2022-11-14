#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

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

namelist = ['Zh', 'Wmuv', 'ZZ', 'mumuvv', 'Zll', 'TauTau', 'WW']
bkg0 = pd.read_csv(namelist[0]+'.csv')
xx = np.array([bkg0['xlow'],bkg0['xhigh']]).ravel(order='F')
yy0 = np.zeros(bkg0.shape[0])
sumerr2 = np.zeros(bkg0.shape[0])
fig = plt.figure(figsize=(10, 8))
ax = fig.add_axes([0.16, 0.16, 0.82, 0.82])
for cname in namelist :
    bkg = pd.read_csv(cname+'.csv')
    yy1 = yy0 + bkg['val']
    yyl = np.array([yy0,yy0]).ravel(order='F')
    yyh = np.array([yy1,yy1]).ravel(order='F')
    # ax.fill_between(xx,yyl,yyh, fc="#00BBCC", edgecolor=None)
    ax.fill_between(xx,yyl,yyh)
    yy0 = yy1
    sumerr2 += bkg['errminus']**2
sumerr2 += yy0/(bkg0['xhigh']-bkg0['xlow'])
halferr = (sumerr2**0.5)/2
errlow = yy0 - halferr
errhigh = yy0 + halferr
errl = np.array([errlow,errlow]).ravel(order='F')
errh = np.array([errhigh,errhigh]).ravel(order='F')
ax.fill_between(xx,errl,errh)
ax.set_xlim(0, 120)
ax.set_ylim(1, 10**5)
ax.set_yscale("log")
ax.set_xlabel(r"$m_{\rm LSP}^{\rm max}~[{\rm GeV}]$", fontsize=30, loc='right')
ax.set_ylabel(r"Events $\left/ 5~{\rm GeV}\right.$", fontsize=30, loc='top')
# plt.plot(xx,yy0,'-')
# plt.show()
plt.savefig('mLSPmax.pdf')
