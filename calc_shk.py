# calc_shk.py 
# Adam Goga
# This file holds the functions which go from raw data and flat to 
# Ca II H&K activity indes shk. These functions assume the spectra 
# has been aligned already.

#File Functions (not in order)
#
#def calc_targOlapf - Take raw data from NRES into lamda & target arrays 
#def hk_windows - main function called for V&R-bands and H&K windows
#def create_window - Used in window alignment to find a window around a lamda val
#def calc_shk - the main function which takes a target calls hk_windows and calc shk val
#def smart_hk_windows - old verify remove is ok
#def smart_calc_shk - old verify remove is ok
#def multi_window_calc_shk - for summing multiple observations


#the original calc_shk function has been split into calc_targOlapf and calc_shk to better 
#debug intermediate values 

# This routine manages the computation of the Ca II H&K activity indes shk.
# Inputs are:
#  mjd = MJD-OBS for the input spectrum.
#  siteid = site from which observation comes.
#  teff = Estimated Teff of target star.
#  extrct(nx,nord) = extracted but not flat-fielded spectra of target. (ADU)
#  lam(nx,nord) = wavelength solution corresp to extrct. (nm)
# Outputs are:
#   shk = (fH + fK)/(fB + fR) where fH, fK are the counts in the cores of
#       the H and K lines, resp., and fB, fR are fluxes in nearby blue and
#       red continuum bands.  Note that shk requires some considerable work
#       to be converted to the more physically useful index R'_HK.
#   eshk = an estimate of the photon+read noise applicable to shk.




import os
import sys

import numpy as np
import scipy.io as sc
import astropy.io.fits

from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.convolution import convolve, Box1DKernel
from astropy.modeling.models import BlackBody
from astropy import units as u
from astropy.time import Time

import helpers as h


#  extrct(nord,nx) = extracted but not flat-fielded spectra of target. (ADU)
#  lam(nord,nx) = wavelength solution corresp to extrct. (nm)
def calc_targOlapf(raw,lamGrid, flatOlap):

    #number of good orders
    highOrd = raw.nOrd#virtually always we want to go to highest(lowest wavelength) order

    gOrd = np.arange(h.lowGOrd,highOrd)

    nx=raw.nx
    lam=raw.waveGrid
    extrct=raw.spec
    
    #old T.Brown constants
    #rdnoi=7.*np.sqrt(5.*5.)           # read noise per resolution element in e-
    #resel=.0015                       #resolution element (nm)


    #T.Brown
    #find the flat field that will be applied to the order-overlapped spectrum.
    #flat field found in above code TODO 
    nLamg=len(lamGrid)

    #T.Brown
    #make the order-overlapped stellar spectrum, divide it by the flat
    nLam = len(lamGrid)
    targOlap=np.zeros(int(nLam),dtype=np.float64)


    from scipy import interpolate
    #https://stackoverflow.com/questions/18326714/idl-interpol-equivalent-for-python
    for i in range(len(gOrd)) :
        lamC = lam[gOrd[i],:]
        dLamdx = np.gradient(lamC)

        scale = dLamdx/dLamdx[int(nx/2.)]

        #original T.T.Brown Code
        #flux=interpol(extrct(*,gord(i))*scale,lam(*,gord(i)),lamgrid)
        interpfunc = interpolate.interp1d(lam[gOrd[i],:],extrct[gOrd[i],:]*scale, kind='linear', fill_value='extrapolate')
        flux=interpfunc(lamGrid)


        #deprecated verify and delete soon
        #sl=where(lamgrid le min(lam(*,gord(i))),nsl)
        #sh=where(lamgrid ge max(lam(*,gord(i))),nsh)
        #if(nsl gt 0) then flux(sl)=0.d0
        #if(nsh gt 0) then flux(sh)=0.d0
        
        #ADAM CODE OF ABOVE - POTENTIALLY BUGGY
        mini = np.min(lam[gOrd[i],:])
        maxi = np.max(lam[gOrd[i],:])

        sl=(lamGrid<=mini).nonzero()
        sh=(lamGrid>=maxi).nonzero()

        if len(sl)>0 :
            for i in sl:
                flux[i]=0.
        if len(sh)>0 :
            for i in sh:
                flux[i]=0.

        #E ADAM CODE
        targOlap=targOlap+flux

    #wherever sg is greater than .01 set targOlapf value
    sg=(flatOlap>=.01).nonzero()
    targOlapf = np.zeros(int(nLamg),dtype=np.float64)
    #same code as below? TODO TEST
    #targOlapf[sg] = targOlap[i]/flatOlap[i]
    for i in sg :
        targOlapf[i] = targOlap[i]/flatOlap[i]


    #fig = plt.figure()
    #print(flatOlap)

    #plt.plot(lamGrid, targOlapf, 'k-')
    #plt.axvline(x=389.9,color='blue')
    #plt.xlim(389,408)
    #plt.xlabel('wavelength [nm]')
    #plt.ylabel('tragOlapf')
    ##plt.show()
    #plt.close()
    
    return targOlapf
    

#could easily pass a function for the pass we'd like but now it's either triangle or flat
def create_window(lamGrid,loc,width,triangle=True,rvcc=None):
    #T.Brown
    #make output array
    nLam = len(lamGrid)
    window = np.zeros(nLam,dtype=np.float32)
    z = 1
    if rvcc is not None:
        z = 1. +rvcc/h.c


    #emmission features
    #creates a triangle filter around CA H and K
    # H 
    if triangle:
        d = abs(lamGrid-loc*z)/(width)
        s = (d<=1.0).nonzero()

        if len(s) > 0:#         |
            window[s]=1.-d[s]#  |
    else:#why do we need *2?    V
        d = (abs(lamGrid-loc*z)*2.)/(width)
        s = (d<=1.0).nonzero()
        if len(s) > 0:
            window[s]=1.

    return window



# rvcc = redshift of target star relative to lab. (km/s)
def hk_windows(lamGrid,rvcc=None,scaleWidth=1):
    #T.Brown
    #make output array
    nLam = len(lamGrid)
    windows = np.zeros((4,nLam),dtype=np.float32)

    z = 1
    if rvcc is not None:
        z = 1. +rvcc/h.c

    
    #emmission features
    #creates a triangle filter around CA H and K
    # H 
    d0 = abs(lamGrid-h.cahLam*z)/(h.lineWid*scaleWidth)
    s = (d0<=1.0).nonzero()
    if len(s) > 0:
        #print('g')
        #print(1.-d0[s])
        windows[0,s]=1.-d0[s] 
    # K
    d1 = abs(lamGrid-h.cakLam*z)/(h.lineWid*scaleWidth)
    s = (d1<=1.0).nonzero()
    if len(s) > 0:
        windows[1,s]=1.-d1[s] 
    
    #continum bands
    #no filter applied to continum bands
    #r-band
    d2 = (abs(lamGrid-h.lamR*z)*2.)/(h.conWid*scaleWidth)
    s = (d2<=1.0).nonzero()
    if len(s) > 0:
        windows[2,s]=1.
    #v-band
    d3 = (abs(lamGrid-h.lamB*z)*2.)/(h.conWid*scaleWidth)
    s = (d3<=1.0).nonzero()
    if len(s) > 0:
        windows[3,s]=1.
        
    return windows

#This is the actual function to calculate the SHK
#Takes the standard lamda grid, and correlated/noise free target data(targOlapf)
#also takes teff for use in the planck function and a radial velocity
#which is not currently used. I have left it in case we decide later to use it
#over cross-correlation
#rvcc is optional parameter if we want radial velocity shift

#OLD COMMENT ABOUT NRES ORDERS
# For NRES, fB is not easily accessible, since the relevant echelle order is
# not extracted.  Therefore, fB is approximated as FR*K(Teff), where the 
# function K is the ratio of Planck functions in the red and blue bandpasses.
#END T.Brown COMMENTS
def calc_shk(lamGrid, targOlapf, raw, windows, rvcc=None):
    starName=raw.star.strip('/')
    
    
    gain=3.4           # e-/ADU

    
    #atm we're not going to apply the radial velocity and just use the adjusted spectra
    rv = 0#rv/10000 Need RV dictionary
    tempEff = 0 #base value if not found BUT BAD
    tempEff = h.tEffLookup[starName]
    if tempEff is 0:
        print('Update temperature for HD ' + starName)
        tempEff = 6200
    



    if raw.nOrd == 67:
        #this function gives us the regions of our arrays which hold the information we need to sum
        #windows = hk_windows(lamGrid-offsets[0],rv)
        raw.fh=(targOlapf*windows[0]).sum()

        #windows = hk_windows(lamGrid-offsets[1],rv)
        raw.fk=(targOlapf*windows[1]).sum()

        #windows = hk_windows(lamGrid-offsets[2],rv)
        raw.fr=(targOlapf*windows[2]).sum()
    elif raw.nOrd == 68:
        #this function gives us the regions of our arrays which hold the information we need to sum
        #windows = smart_hk_windows(lamGrid,rv)

        raw.fh=(targOlapf*windows[0]).sum()
        raw.fk=(targOlapf*windows[1]).sum()
        raw.fr=(targOlapf*windows[2]).sum()
        raw.fb=(targOlapf*windows[3]).sum()

    
    #windows = hk_windows(lamGrid,rv)
    if h.debug:
        plt.figure()
        cur=windows[0]
        plt.plot(lamGrid[cur!=0],targOlapf[cur!=0])#*cur[cur!=0],'k-')
        plt.show()
        plt.close()
        plt.figure()
        cur=windows[1]
        plt.plot(lamGrid[cur!=0],targOlapf[cur!=0])#*cur[cur!=0],'k-')
        plt.show()
        plt.close()
        plt.figure()
        cur=windows[2]
        plt.plot(lamGrid[cur!=0],targOlapf[cur!=0])#*cur[cur!=0],'k-')
        plt.show()
        plt.close()
    
    if raw.nOrd == 67:
        #the SHK calculation with pseudo V-Band
        #https://docs.astropy.org/en/v4.1/modeling/blackbody_deprecated.html
        #plFactor = blackbody_lambda(h.lamB*u.nm,tempEff*u.K).value/blackbody_lambda(h.lamR*u.nm,tempEff*u.K).value
        #original above
        bb=BlackBody(tempEff*u.K, scale=1.0 * u.erg / (u.cm ** 2 * u.AA * u.s * u.sr))
        plFactor=bb(h.lamB*u.nm).value/bb(h.lamR*u.nm).value
        raw.fb = raw.fr*plFactor

    num = (raw.fh+raw.fk)*gain
    den = (raw.fr+raw.fb)*gain
    alpha = h.siteAlpha[raw.site]

    if hasattr(raw,'format') and raw.format == True:
        alpha = alpha * h.oldScale

    shk = alpha*(raw.fh+raw.fk)/(raw.fr+raw.fb)
    #print("shk: "+ str(shk))

    return shk, windows, raw.fr/raw.fb

def multi_window_calc_shk(lamGrid, targOlapf, raw, windows, rvcc=None):
    starName=raw.star.strip('/')
    
    
    gain=3.4           # e-/ADU

    
    #atm we're not going to apply the radial velocity and just use the adjusted spectra
    rv = 0#rv/10000 Need RV dictionary
    tempEff = 0 #base value if not found BUT BAD
    tempEff = h.tEffLookup[starName]
    if tempEff is 0:
        print('Update temperature for HD ' + starName)
        tempEff = 6200
    
    fh=0
    fk=0
    fr=0
    fb=0


    if raw.nOrd == 67:

        for i in range(len(targOlapf)):
            #this function gives us the regions of our arrays which hold the information we need to sum
            #windows = hk_windows(lamGrid-offsets[0],rv)
            fh+=(targOlapf[i]*windows[i][0]).sum()

            #windows = hk_windows(lamGrid-offsets[1],rv)
            fk+=(targOlapf[i]*windows[i][1]).sum()

            #windows = hk_windows(lamGrid-offsets[2],rv)
            fr+=(targOlapf[i]*windows[i][2]).sum()
    elif raw.nOrd == 68:
        for i in range(len(targOlapf)):
            #this function gives us the regions of our arrays which hold the information we need to sum
            #windows = smart_hk_windows(lamGrid,rv)

            fh+=(targOlapf[i]*windows[i][0]).sum()
            fk+=(targOlapf[i]*windows[i][1]).sum()
            fr+=(targOlapf[i]*windows[i][2]).sum()
            fb+=(targOlapf[i]*windows[i][3]).sum()

    
    #windows = hk_windows(lamGrid,rv)
    if h.debug:
        plt.figure()
        cur=windows[0]
        plt.plot(lamGrid[cur!=0],targOlapf[cur!=0])#*cur[cur!=0],'k-')
        plt.show()
        plt.close()
        plt.figure()
        cur=windows[1]
        plt.plot(lamGrid[cur!=0],targOlapf[cur!=0])#*cur[cur!=0],'k-')
        plt.show()
        plt.close()
        plt.figure()
        cur=windows[2]
        plt.plot(lamGrid[cur!=0],targOlapf[cur!=0])#*cur[cur!=0],'k-')
        plt.show()
        plt.close()
    
    if raw.nOrd == 67:
        #the SHK calculation with pseudo V-Band
        #plFactor = blackbody_lambda(h.lamB*u.nm,tempEff*u.K).value/blackbody_lambda(h.lamR*u.nm,tempEff*u.K).value
        #https://docs.astropy.org/en/v4.1/modeling/blackbody_deprecated.html
        #plFactor = blackbody_lambda(h.lamB*u.nm,tempEff*u.K).value/blackbody_lambda(h.lamR*u.nm,tempEff*u.K).value
        #original above
        bb=BlackBody(tempEff*u.K, scale=1.0 * u.erg / (u.cm ** 2 * u.AA * u.s * u.sr))
        plFactor=bb(h.lamB*u.nm).value/bb(h.lamR*u.nm).value
        fb = fr*plFactor

    num = (fh+fk)*gain
    den = (fr+fb)*gain
    alpha = h.siteAlpha[raw.site]

    if hasattr(raw,'format') and raw.format == True:
        alpha = alpha * h.oldScale

    shk = alpha*(fh+fk)/(fr+fb)
    #print("shk: "+ str(shk))

    return shk, windows, fr/fb

#smarts specific hk windows with V-Band included
def old_hk_windows(lamGrid,rvcc=0):
    #make output array
    nLam = len(lamGrid)
    windows = np.zeros((nLam,3),dtype=np.float32)
    z = 1. +rvcc/h.c
    
    
    #emmission features
    #creates a triangle filter around CA H and K
    # H 
    d0 = abs(lamGrid-h.cahLam*z)/h.lineWid
    s = (d0<=1.0).nonzero()
    print('old ' + str(len(s)))
    if len(s) > 0:
        windows[s,0]=1.-d0[s] 
    # K
    d1 = abs(lamGrid-h.cakLam*z)/h.lineWid
    s = (d1<=1.0).nonzero()
    if len(s) > 0:
        windows[s,1]=1.-d1[s] 
    
    #continum bands
    #no filter applied to continum bands
    #r-band
    d2 = abs(lamGrid-h.lamR*z)*2./h.conWid
    s = (d2<=1.0).nonzero()
    if len(s) > 0:
        windows[s,2]=1.
    #v-band
    #d3 = abs(lamGrid-h.lamB*z)*2./h.conWid
    #s = (d3<=1.0).nonzero()
    #if len(s) > 0:
    #    windows[s,3]=1.
        
    return windows




#DEPRECATED
#smarts specific calc_shk which just calls the special windows function
#to find the v-band window and uses the real v-band data
def old_calc_shk(lamGrid, targOlapf, raw, rvcc=0):
    starName=raw.star.strip('/')
    
    
    gain=3.4           # e-/ADU

    
    #atm we're not going to apply the radial velocity and just use the adjusted spectra
    rv = 0#rv/10000 Need RV dictionary
    tempEff = 0 #base value if not found BUT BAD
    tempEff = h.tEffLookup[starName]
    if tempEff is 0:
        print('Update temperature for HD ' + starName)
        tempEff = 6200
    


    windows = hk_windows(lamGrid,rv)
    if raw.nOrd == 67:
        #this function gives us the regions of our arrays which hold the information we need to sum
        fh=(targOlapf*windows[:,0]).sum()
        fk=(targOlapf*windows[:,1]).sum()
        fr=(targOlapf*windows[:,2]).sum()


    
    windows = hk_windows(lamGrid,rv)
    
    #plt.figure()
    #cur=windows[:,0]
    #plt.plot(lamGrid[cur!=0],targOlapf[cur!=0])#*cur[cur!=0],'k-')
    #plt.show()
    #plt.close()    

    #the SHK calculation with pseudo V-Band
    #plFactor = blackbody_lambda(h.lamB*u.nm,tempEff*u.K).value/blackbody_lambda(h.lamR*u.nm,tempEff*u.K).value
    #https://docs.astropy.org/en/v4.1/modeling/blackbody_deprecated.html
    #plFactor = blackbody_lambda(h.lamB*u.nm,tempEff*u.K).value/blackbody_lambda(h.lamR*u.nm,tempEff*u.K).value
    #original above
    bb=BlackBody(tempEff*u.K, scale=1.0 * u.erg / (u.cm ** 2 * u.AA * u.s * u.sr))
    plFactor=bb(h.lamB*u.nm).value/bb(h.lamR*u.nm).value
    fb = fr*plFactor

    num = (fh+fk)*gain
    den = (fr+fb)*gain
    alpha = h.siteAlpha[raw.site]

    if hasattr(raw,'format') and raw.format == True:
        alpha = alpha * h.oldScale

    shk = alpha*(fh+fk)/(fr+fb)

    return shk, windows, fr/fb