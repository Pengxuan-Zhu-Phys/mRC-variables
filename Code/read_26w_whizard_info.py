#!/usr/bin/env python3 

import os, sys 
import pandas as pd 
pwd = os.path.abspath(os.path.dirname(__file__))
import json 

wmss = {
    "01":   "80362",
    "02":   "80367",
    "03":   "80372",
    "04":   "80378",
    "05":   "80382",
    "06":   "80387",
    "07":   "80392"
}

def read_mass(nn):
    pos = "corr. to luminosity [fb-1] ="
    lumis = {}
    for ii in range(26):
        logf = os.path.abspath(os.path.join(pwd, "../Whizard3/26M{}/whizard_{:0>3}.log".format(nn, ii)))
        print(logf)
        with open(logf, 'r') as f1:
            logg = f1.readlines()
            for line in logg:
                if pos in line:
                    print(line)
                    lumi = float(line.split()[-1])
                    lab  = "{}_{:0>2}".format(wmss[nn], ii)
                    lumis[lab] = lumi
    with open("{}.json".format(wmss[nn]), "w") as f1:
        json.dump(lumis, f1)


if __name__ == "__main__":
    read_mass("07")