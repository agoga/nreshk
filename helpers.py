# helpers.py 
# Adam Goga
# Holds the two data structures built for spectra data(raw/analyzed) and other
# helper functions. Most of the constants for the pipeline which deal with 
# spectral analysis are defined here as well

#File Functions
#
#def print_header
#class rawData
#class analyzedData
#def mkdir_p
#def bad_spec_detection_v2
#def mjd_from_hdu
#def find_nearest_mjd
#def closestKey
#def is_folder_star
#def get_immediate_subdirectories
#def find_closest

import os
import copy

import numpy as np
import scipy as sc

import astropy.io.fits
import astropy.io.fits.header

from scipy import signal
from os import makedirs,path
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.time import Time
from errno import EEXIST



debug = False#DONT TRUST THIS, BAD ADAM

c=2.99792458e5  #speed of light (km/s)

cahLam=396.847# 396.967 #   #Ca II H line wavelength (nm, vacuum)
cakLam=393.366#393.485 #    #Ca II K line wavelength (nm, vacuum)


#center of blue continuum band (nm, vacuum)
#offset by .5 nm for new 68 order data

lamB=390.8#-.116 subtraction is the offset from our lab values     
#center of red continuum band (nm, vacuum)
#TODO TAKE FROM VAUGHAN 1978 subtraction is the offset from our lab values 
lamR=399.4#-.116subtraction is the offset from our lab values  

conWid=1#2.0         #wavelength width of continuum bands (nm, vacuum)
lineWid=.109       #FWHM of line core window functions (nm, vacuum)
    
sigToFWHM = 2.355#used to display a real FWHM value on plots where smoothing occurs

lowGOrd = 63
highGOrd = 66

oldScale = 1#6/5

siteColors = {  'lsc':["b","Cerro Tololo Interamerican Obs\'"],
                'cpt':["g","South African Astro Obs\'"],
                'elp':["r","McDonald Obs\'"],
                'tlv':["k","Wise Obs\'"]}


# factor to make shk into equivalent width (a guess currently!)
siteAlpha ={     'lsc':44,#The alpha value in the Ca HK SHK calculation is telescope dependent
                'cpt':10,
                'elp':39,
                'tlv':39}

#TODO NEW STARS NEED TEFF AND CONFIRM MOST OF THESE TEFF
tEffLookup = {"1835":5837,
              "12235":6097,
              "20630":5742,
              "26913":5631,#picked one from 15 at http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%40700565&Name=HD++26913&submit=display+all+measurements#lab_meas
              "37394":5249,
              "43587":5892,
              "75332":6258,
              "78366":6014,
              "82443":5287,
              "82885":5499,
              "88737":6345,#http://simbad.u-strasbg.fr/simbad/sim-basic?Ident=HD+88737
              "98230B":6005,#https://en.wikipedia.org/wiki/Xi_Ursae_Majoris
              "100180":6002,
              "114710":6075,
              "115383":6234,
              "115404":4976,#https://en.wikipedia.org/wiki/HD_115404
              "126053":5714,
              "136202":6052,
              "149661":5277,
              "152391":5479,
              "154417":5989,
              "165341A":5419,#http://simbad.u-strasbg.fr/simbad/sim-basic?Ident=HD+165341
              "176051":6000,#https://en.wikipedia.org/wiki/HD_176051
              "182101":6427,#http://simbad.u-strasbg.fr/simbad/sim-basic?Ident=HD+182101
              "187691":6139,
              "194012":6250,#ONE OF THREE http://simbad.u-strasbg.fr/simbad/sim-id?mescat.distance=on&mescat.velocities=on&Ident=%402756705&Name=HD+194012&submit=display+all+measurements#lab_meas
              "206860":5974,
              "22049":5065,
              "17051":6200,#could not find in references, google's value
              "30495":5834,
              "49933":6674,
              "76151":5761,#could not find in references, google's value
              "120136":6420,#could not find in references, google's value
              "190406":5940
            }

def print_header(header):
    for h in header.keys():
        print(h, header[h])    
        
#hopefully not SHK dependant, move to pipeline        
class rawData:
    header:astropy.io.fits.header.Header#
    spec:list
    fitsFile:str
    format:bool#change to int for more than 2 formats
    waveGrid:list
    mjd:float#
    site:str#header['SITEID']
    date:str#header['DATE-OBS']
    star:str
    nx:int
    nOrd:int
    alpha:float

    #def __init__(self, target=None):
    ##    if orig is None:
      #      self._constructor()
      #  else:
      #      self._copy_constructor(target)

    def __init__(self=None,waveGrid=None,spec=None,header=None,fileName=None,format=None,copy=None):
        self.mjd = None
        
        if copy is None:
            self.header = header
            #TODO check the header as you go and throw error if missing any
            self.mjd = header['MJD-OBS']
            self.site = header['SITEID']
            self.alpha = siteAlpha[self.site]
            self.date = header['DATE-OBS']
            self.star = header['OBJ1']
            self.nOrd = header['NORD']
            self.nx = header['NX']
            self.spec=spec
            self.header=header
            self.fitsFile=fileName
            self.format=format 
            self.waveGrid = waveGrid   
        else:
            self.header = copy.header
            self.star = copy.star
            self.mjd = copy.mjd
            self.site = copy.site
            self.date = copy.date
            self.spec=copy.spec
            self.fitsFile=copy.fitsFile
            self.format=copy.format 
            self.waveGrid = copy.waveGrid
            self.nOrd = copy.nOrd
            self.nx = copy.nx
            self.alpha = copy.alpha

    


#TODO maybe make rawData which specData and printData inherit from.
#not a NRES SHK specific function            
class analyzedData(rawData):
    flat:[]

    
    shk:float#
    window:[]#window order: Ca H, Ca K, R band, B band
    offset:[]#list of offsets for each window


    day:int
    hour:int
    
    shk:float#
    decimalYr:Time
    label:str

    bad:bool
    lamGrid:list#
    targOlapf:list#

    average:bool#

    def __init__(self, raw=None,lamGrid=None,flat=None,targOlapf=None, shk=None,offset=None,average=None,bad=None,window=None):
        if raw is not None:
            super().__init__(copy=raw)


        if hasattr(self,'mjd') and self.mjd is not None:
            mjd = self.mjd
            self.day = int(np.floor(mjd))
            self.hour = int(str(mjd).split('.')[1])
            self.decimalYr = Time(mjd,format='mjd')
            self.decimalYr.format = 'decimalyear'
        else:
            self.mjd = None

        
        self.bad=bad
        self.flat = flat
        self.lamGrid = lamGrid
        self.targOlapf = targOlapf
        self.offset = offset
        self.shk = shk
        self.average = average
        self.window = window

    def raw(self):#for 
        return self

    
    #if this is an average then we put it in the day's folder as day_combined_data and save it slightly diff
    #location that the data should go to
    def data_path(self):
        if self.average is False:
            return self.outputDir()+str(self.hour)+"_data"
        else:
            return self.outputDir()+str(self.day)+"_combined_data"

    #location pdf reports go to
    def report_path(self):
        if self.average is False:
            return  self.outputDir()+str(self.hour)+"_"+str(self.nOrd)+"_report.pdf"
        else:
            return  self.outputDir()+str(self.day)+"_"+str(self.nOrd)+"_combined_report.pdf"

    def label(self):
        return 'MJD: ' + str(self.mjd) + ' and decYr ' + str(self.decimalYr) + ' with '+str(self.nOrd)+' orders. SHK: ' + str(self.shk) + ' and offsets:' + str(self.offset)

    def starDir(self):
        return  "output/"+self.star+"/"

    def outputDir(self):
        return "output/"+self.star+"/"+str(self.day)+"/"

    def pdfTitle(self):
        t = ''
        if self.bad:
            t = 'bad '
        t += 'NRES spectra, ' + self.site +', '+self.date+' ('+ '{:.6}'.format(self.decimalYr.value) +', ' +str(self.mjd)+'), S='+'{:.4}'.format(self.shk)
        return t
            
        

#https://stackoverflow.com/questions/11373610/save-matplotlib-file-to-a-directory
def mkdir_p(mypath):
    '''Creates a directory. equivalent to using mkdir -p on the command line'''
    try:
        makedirs(mypath)
    except OSError as exc: # Python >2.5
        if exc.errno == EEXIST and path.isdir(mypath):
            pass
        else: raise
            

        
#the idea behind this function is to auto-correlate a spectra.
#if the spectra is not complete noise then the correlation function will have a small peek at 0
#if the spectra is noise then the peek should be much more massive since no where besides at no
#offset will the noise overlap.
def bad_spec_detection_v2(lamGrid, targ):
    


    correlation = signal.fftconvolve(targ,targ[::-1],'full')
    #correlation = np.correlate(targ,targ,'full')
    
    #TODO maybe take the maximum out of the mean for a different ratio.

    
    mean = np.mean(correlation[39900:40100])
    maxi = max(correlation[39900:40100])
    #print("bad correlation: " + str(maxi/mean)) 
    
    #plt.figure()
    #plt.plot(range(len(correlation)),correlation)
    #plt.xlim([39900,40100])
    #plt.show()
    #plt.close()
    
    
    #TODO for future projects, fiddling with this number will remove or allow different spectra
    #1.7 might be a better go
    if maxi/mean > 1.7:
        return True
    else:
        return False
    
#deprecated but potentially useful
#we would like to find the 2 lowest local min in this region
#that should be the H and K lines. If they are farther apart than we expect them to be
#then we need to throw out this spectra

#ran into issues with it finding other peeks and the noisey data not being bad enough
#Also the code is buggy I think..opps
def bad_spec_detection(lamGrid, targ):
    
    
    #1.5 angstroms allowed diff from the desired spectral features
    allowedError = .24
    
    lowI = np.abs(lamGrid - 392).argmin()
    highI =np.abs(lamGrid - 405).argmin()
    bad = True
    width = 1000
    
    dLam = lamGrid[1] - lamGrid[0]
    
    
    gausdTarg = sc.ndimage.filters.gaussian_filter(targ,.05/dLam)
    newLook = gausdTarg[lowI:highI]
    adjLam = lamGrid[lowI:highI]

    first = int(np.nanargmin(newLook))
    
    low = first-width
    high = first+width
    if low < 0:
        low = 0
    if high > len(newLook) - 1:
        high = len(newLook)
    newLook[low:high] = np.nan
    
    second = int(np.nanargmin(newLook))
    
    #print('first: ' + str(first) + ' second: ' + str(second))

    
    #fig = plt.figure()
    #plt.axvline(x=adjLam[0],color='r')
    #plt.axvline(x=adjLam[-1],color='r')
   # 
    #plt.axvline(x=adjLam[first],color='b')
    #plt.axvline(x=adjLam[second],color='b')
    #plt.plot(adjLam, newLook, 'k-')
    #plt.show()
    #plt.close()
    trueDist = abs(cahLam - cakLam)
    realDist = abs(adjLam[second] - adjLam[first])
    
    
    if abs(trueDist - realDist) <= allowedError:
        bad = False
        print('spec with diff: ' + str(abs(trueDist - realDist)))
    else:
        print('bad spec with diff: ' + str(abs(trueDist - realDist)) + ' first: ' + str(first) + ' second: ' + str(second))
        #gausdTarg = sc.ndimage.filters.gaussian_filter(targ,.15/dLam)
        #newLook = gausdTarg[lowI:highI]
        #fig = plt.figure()
        #plt.axvline(x=adjLam[0],color='r')
        #plt.axvline(x=adjLam[-1],color='r')
        #plt.axvline(x=adjLam[first],color='b')
        #plt.axvline(x=adjLam[second],color='b')
        #plt.plot(adjLam, newLook, 'k-')
        #plt.show()
        #plt.close()
    return bad

    #todo clean and place into healpers file

#do you want to crash? name something just 'print'...
def dprint(input):
    if debug:
        print(input)

def mjd_from_hdu(hdu):
    return hdu[0].header['MJD-OBS']

#https://stackoverflow.com/questions/47725773/finding-an-integer-key-of-a-python-dict-that-is-closest-to-a-given-integer
def find_nearest_mjd(dd,mjd):
    low = max([d for d in dd if d<= mjd])
    high = min([d for d in dd if d>= mjd])
    nearkey = low if mjd - low <= high - mjd else high
    return nearkey

def closestKey(dic, key):
    diff = {k:abs(k - key) for k in dic}
    return min(diff, key=diff.get)

def is_folder_star(folder):
    try:
        #if first letter is number we can run pipeline on it
        int(folder[0])
        return True
    except ValueError:
        return False
    
def get_immediate_subdirectories(a_dir):
    return [name for name in os.listdir(a_dir)
            if os.path.isdir(os.path.join(a_dir, name))]

#https://stackoverflow.com/questions/21388026/find-closest-float-in-array-for-all-floats-in-another-array
def find_closest(A, target):
    #A must be sorted
    idx = A.searchsorted(target)
    idx = np.clip(idx, 1, len(A)-1)
    left = A[idx-1]
    right = A[idx]
    idx -= target - left < right - target
    return idx