#!/usr/bin/env python3 

import os, sys 
from operator import mod
import pyhf 
import json, jsonschema
import numpy as np 
import scipy

jpath = os.path.abspath("~/Jarvis-HEP/src".replace("~", os.path.expanduser("~")))
pwd = os.path.abspath(os.path.dirname(__file__))

sys.path.append(jpath)
from IOs import YodaFile

XSect = 395.000
Lint  = 5050.



def load_wmass(wmass):
    wmass = "{}".format(wmass)

    yodamin = YodaFile()
    yodamin.file = os.path.join(pwd, "{}/Nhist_Wmass_mRCmin.dat".format(wmass))
    # yodamin.file = os.path.join(pwd, "{}/Nhist_mRCmin.dat".format(wmass))
    yodamin.read()

    yodamax = YodaFile()
    yodamax.file = os.path.join(pwd, "{}/Nhist_Wmass_mRCmax.dat".format(wmass))
    yodamax.read()
    data = {}

    for kk, vv in yodamin.hist1ddata.items():
        if "mRCmin" in kk:
            data['mRCmin'] = vv['data']
    for kk, vv in yodamax.hist1ddata.items():
        if "mRCmax" in kk:
            data['mRCmax'] = vv['data']
    
    return data

truedata = load_wmass(80390)
# print(truedata)
idmin = 0
idmax = -1

signal = list(truedata['mRCmin']['val'] * XSect * Lint)[idmin:idmax]
background = list(truedata['mRCmin']['val'] * XSect * Lint)[idmin:idmax]
bkg_err = list(truedata['mRCmin']['errminus'] * XSect * Lint)[idmin:idmax]
# print(signal[0]**0.5, bkg_err[0])

dof = len(signal)
CI = scipy.stats.chi2.cdf(1., 1.)
chi2 = scipy.stats.chi2.ppf(CI, dof)
# signal = list(truedata['mRCmin']['val'])
# background = list(truedata['mRCmin']['val'])
# bkg_err = list(truedata['mRCmin']['errminus'])

modelmin = pyhf.simplemodels.uncorrelated_background(
    signal, background, bkg_err
)

loglike = []
for ii in range(13):
    wma = 80330 + 10 * ii
    test = load_wmass(wma)
    observations = list(test['mRCmin']['val'] * XSect * Lint)[idmin:idmax]
    data = pyhf.tensorlib.astensor(observations + modelmin.config.auxdata)
    pars_bonly, twice_nll = pyhf.infer.mle.fixed_poi_fit(0., data, modelmin, return_fitted_val=True)
    # print(data)
    likelihood = -2.0 * modelmin.logpdf(pars_bonly, data)
    loglike.append(likelihood[0])

loglike -= min(loglike)
print(loglike)


idmin = 20
idmax = 30
signal = list(truedata['mRCmax']['val']* XSect * Lint)[idmin:idmax]
background = list(truedata['mRCmax']['val']* XSect * Lint)[idmin:idmax]
bkg_err = list(truedata['mRCmax']['errminus']* XSect * Lint)[idmin:idmax]
# bkg_err = list(truedata['mRCmin']['errminus']/10.0 )

dof2 = len(signal)
CI = scipy.stats.chi2.cdf(1., 1.)
chi22 = scipy.stats.chi2.ppf(CI, dof2)

modelmax = pyhf.simplemodels.uncorrelated_background(
    signal, background, bkg_err
)
loglike2 = []

for ii in range(13):
    wma = 80330 + 10 * ii
    test = load_wmass(wma)
    observations = list(test['mRCmax']['val']* XSect * Lint)[idmin:idmax]
    data = pyhf.tensorlib.astensor(observations + modelmax.config.auxdata)
    pars_bonly, twice_nll = pyhf.infer.mle.fixed_poi_fit(0., data, modelmax, return_fitted_val=True)
    # print(data)
    likelihood = -2.0 * modelmax.logpdf(pars_bonly, data)
    loglike2.append(likelihood[0])

loglike2 -= min(loglike2)
print(loglike2)
# print(loglike2 / chi22 )

ll = loglike2 + loglike
print(ll)
# print(chi2, chi22)
CI = scipy.stats.chi2.cdf(1., 1.)
chi2t = scipy.stats.chi2.ppf(CI, dof+dof2)
# print(chi2t)
# print(ll/chi2t)

result = {
    "wmass":    np.linspace(80330., 80450., 13),
    "loglik_min":   loglike,
    "sigma_min":    loglike / chi2,
    "loglik_max":   loglike2,
    "sigma_max":    loglike2/ chi22,
    "loglik_tot":   ll, 
    "sigma_tot":    ll/chi2t
}
from pandas import DataFrame
df = DataFrame(result)
df.to_csv("FitWmass.csv")



