#!/usr/bin/env python3 

import os, sys 

import scipy.stats
import math 

CI = scipy.stats.chi2.cdf(1., 1.)

chi2 = scipy.stats.chi2.ppf(CI, 15)
print(chi2)
# pwd = os.path.abspath(os.path.dirname(__file__))

# for ii in range(13):
#     wmass = 80330 + 10*ii
#     yodapath = os.path.abspath(os.path.join(pwd, "../{}_tot.yoda".format(wmass)))
#     cmd = "rivet-mkhtml --errs {}".format(yodapath)
#     os.system(cmd)
#     mvcmd = "cp -r ./rivet-plots/CEPC_MASS {}".format(wmass)
#     os.system(mvcmd)
    # print(cmd, "\n", mvcmd)