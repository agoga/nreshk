from scipy import interpolate
import numpy as np
from calc_shk import calc_shk
from astropy.time import Time
import plotting as plot
import helpers as h
import numpy as np
import scipy.constants as sc
import astropy.io.fits
import numpy as np
import os
import glob

import helpers as h
import logging

from calc_shk import calc_targOlapf
from calc_shk import calc_shk
import pickle
import pprint
import pickle
import pprint
import os
import astropy.io.fits

import pipeline as pipe
import plotting as plot

from matplotlib import pyplot as plt

from astropy.time import Time


def make_random_data(n):
    r = np.random.rand(n)
    b = []
    for p in r:
        d = 58890 +20*p
        b.append(h.analyzedData(d))
        b.append(h.analyzedData(d+.1,))

    return b


def test_daily_data_sum():
    r = make_random_data(10)

    pipe.sum_daily_data(r,'',[])


#divide out by the original alpha then multiply by new
def adjust_sites(data,alphas):
    for d in data:
        
        unShk = d.shk/d.alpha
        d.alpha = alphas[d.site]
        d.shk = unShk*d.alpha

        if d.average is True:
            print("MJD: "+str(d.mjd)+" SHK: " + str(d.shk))
    return data

