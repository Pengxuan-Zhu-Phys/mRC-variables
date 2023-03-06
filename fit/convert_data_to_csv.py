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

for ii in range(13):
    wma = 80330 + 10 * ii
    test = load_wmass(wma)
    # print(test)
    for col in ['val', 'errminus', 'errplus']:
        test['mRCmin'][col] = test['mRCmin'][col] * XSect * Lint
        test['mRCmin'].to_csv("./Hist/{}_mRCmin.csv".format(wma))

        test['mRCmax'][col] = test['mRCmax'][col] * XSect * Lint
        test['mRCmax'].to_csv("./Hist/{}_mRCmax.csv".format(wma))
