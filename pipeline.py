# pipeline.py 
# Adam Goga
# 
# Intended to be the home of general functions that most spectra analysis will need.
#
#def import_aligning_spectra - load from file a spectra which aids in aligning our observations
#def create_integration_windows - Creates the specific integration windows for each interesting location. 
# for NRES that's V&R bands and CA H&K locations
#
#def find_window_offsets - uses calc_del_lam to align each window against the alignment spectra
#def calc_del_lam - First interpolates the align spectra onto the obsevations lamda grid then use 
# convolution of the obs and align and fit a parabola around the center to discover accurate offset of the window
#
#def sum_daily_data - goes through entire list of analyzed data and sums the observations of each night to provide
# an averaged SHK value



import os
import os.path
import re
import scipy
import scipy.signal
from scipy import signal
import scipy.ndimage.filters

import numpy as np
import matplotlib as mpl
import scipy as sc
import numpy.polynomial.polynomial as poly

from matplotlib import pyplot as plt
from astropy.convolution import convolve, Box1DKernel
from scipy import interpolate
from astropy.time import Time

from calc_shk import calc_shk
from calc_shk import hk_windows
from calc_shk import create_window
from calc_shk import multi_window_calc_shk
import plotting as plot
import helpers as h

#remove once calc shk hack gone
from astropy.modeling.blackbody import blackbody_lambda
from astropy import units as u
#E remove
#this file holds pipeline functions which are not telescope dependant

#imports from a lot of data the lab frame spectra so we may cross-correlate onto this and
#find better locations of the ca hk lines
#taken from Ricky Egeland's code, probably could be optimized since he
#did more with this code than we need

#TODO determine if this needs change for other uses 
def import_aligning_spectra(fluxdir, minwl=None, maxwl=None, resolution=1, residual=False):   
    wavelength = []
    res_flux = []
    irradiance = []
    
    for f in os.listdir(fluxdir):
        if not re.match('^lm', f): continue
        fpath = os.path.join(fluxdir, f)
        d = np.genfromtxt(fpath, unpack=True)
        wavelength = np.append(wavelength, d[0]) # in nm
        res_flux = np.append(res_flux, d[1]) # from a normalized spectra
        irradiance = np.append(irradiance, d[2]) # from a regular spectra; mu-W / cm^2 / nm
        
    wavelength = np.array(wavelength)
    res_flux = np.array(res_flux)
    irradiance = np.array(irradiance)
    sort = np.argsort(wavelength)
    wavelength = wavelength[sort]
    res_flux = res_flux[sort]
    irradiance = irradiance[sort]

    angstroms = wavelength * 10 # in Angstroms
    
    if minwl is None:
        minwl = angstroms[0]
    if maxwl is None:
        maxwl = angstroms[-1]
    sel = (angstroms >= minwl) & (angstroms <= maxwl)
    dw = angstroms[1] - angstroms[0]
    if not residual:
        series = irradiance / 10. # nm^-1 => Ang.^-1
    else:
        series = res_flux

    if resolution > 0:
        series = scipy.ndimage.filters.gaussian_filter(series, resolution/dw)

    #print('degrading source: ' + str(resolution/dw))
    
    return angstroms[sel], series[sel]
    
    


def create_integration_windows(tarGrid, targ, offsets=None):
    if offsets is None:
        offsets = np.zeros(4,dtype=np.float32)
    
    windows = []

    windows.append(create_window(tarGrid-offsets[0],h.cahLam,h.lineWid,True))
    windows.append(create_window(tarGrid-offsets[1],h.cakLam,h.lineWid,True))
    windows.append(create_window(tarGrid-offsets[2],h.lamR,h.conWid,True))
    windows.append(create_window(tarGrid-offsets[3],h.lamB,h.conWid,True))

    return windows
    



def find_window_offsets(labGrid, lab, tarGrid, targ, scale, smooth):
    llen = len(labGrid)
    tlen = len(tarGrid)

    windows = [[],[]]
    #Ca H
    windows[0].append(create_window(labGrid, h.cahLam, scale*h.lineWid, True))
    windows[1].append(create_window(tarGrid, h.cahLam, scale*h.lineWid, True))
    #Ca k
    windows[0].append(create_window(labGrid, h.cakLam, scale*h.lineWid, True))
    windows[1].append(create_window(tarGrid, h.cakLam, scale*h.lineWid, True))
    #R band
    windows[0].append(create_window(labGrid, h.lamR, scale*h.conWid, False))
    windows[1].append(create_window(tarGrid, h.lamR, scale*h.conWid, False))
    #B Band
    windows[0].append(create_window(labGrid, h.lamB, scale*h.conWid, False))
    windows[1].append(create_window(tarGrid, h.lamB, scale*h.conWid, False))

    #windows[0].append(create_window(labGrid, h.lamR, scale*h.conWid, False))
    #windows[1].append(create_window(tarGrid, h.lamB, scale*h.conWid, False))

    #labWin = [lhWindow,lkWindow,lrWindow]
    #tarWin = [thWindow,tkWindow,trWindow]
    if h.debug:
        plt.figure()
        cur=windows[1][0]
        plt.plot(tarGrid[cur!=0],targ[cur!=0])#*cur[cur!=0],'k-')
        plt.show()
        plt.close()
        plt.figure()
        cur=windows[1][1]
        plt.plot(tarGrid[cur!=0],targ[cur!=0])#*cur[cur!=0],'k-')
        plt.show()
        plt.close()
        plt.figure()
        cur=windows[1][2]
        plt.plot(tarGrid[cur!=0],targ[cur!=0])#*cur[cur!=0],'k-')
        plt.show()
        plt.close()
    

    offsets = []

    for i in range(4):
        cLab=windows[0][i]
        cTar=windows[1][i]

        l1 = labGrid[cLab!=0]
        l2 = lab[cLab!=0]
        l3 = tarGrid[cTar!=0]
        l4 = targ[cTar!=0]
        if h.debug:
            plt.figure()
            cur=windows[0][i]
            plt.plot(labGrid[cur!=0],lab[cur!=0])#*cur[cur!=0],'k-')
            plt.show()
            plt.close()

            plt.figure()
            cur=windows[1][i]
            plt.plot(tarGrid[cur!=0],targ[cur!=0])#*cur[cur!=0],'k-')
            plt.show()
            plt.close()

        out = calc_del_lam(l1,l2,l3,l4,.005,i)
        offsets.append(out[0])
        labInterp = out[1]
    
    #offsets.append(0)#temp blue hack
    return offsets#, labInterp, tarW

    

#############################################
#This is the big function we calculate our delta lambe value to place this current
#spectra in the pipeline in the lab reference frame.
#First we will interpolate to get these on the same grid
#Then we do a cross correlation
#With the cross correlation function we will fit a parabola right around the peak
#We do this because we want to have more precision for offset than our grid will give us
#the peak of this new parabola is the real offset
#lots of plotting function commented out for debugging since this function does a lot and
#could certainly be improved on.
#TODO not called if rvcc provided for star
def calc_del_lam(labGrid, lab, tarGrid, targ, smooth,tmpDebug=0) :
    tmpGridScale = 1
    dLam = tarGrid[1] - tarGrid[0]

    h.dprint('dlam: ' + str(dLam))
    #print(np.shape(targ))
    gausdTarg = sc.ndimage.filters.gaussian_filter(targ,(smooth/dLam)/2)#TODO hmmm
    


    


    dLabLam = labGrid[1] - labGrid[0]
    
    #TODO SHOULD 2.55 be h.sigToFWHM?
    gausedLab = lab#sc.ndimage.filters.gaussian_filter(lab,(dLam/dLabLam)/(2.55*2)

    #get the lab spectrum into angs div by 10 on angstrom grid for our purposes
    interpfunc = interpolate.interp1d(labGrid, gausedLab,fill_value='extrapolate') # kind='linear')#9/1/2020 swapped 'ext..' w 'lin..''
    labInterp=interpfunc(tarGrid)


    #ZERO OUT THE EDGES OF OUR LAB SPECTRA
    #TODO do this with strict values to remove edge errors which may affect correlation
    #0 out the lab grid edges based on first and last nonzero element of targolap
    #h.dprint(labInterp)
    labInterp[:targ.nonzero()[0][0]]=0
    #from last nonzero element to end
    labInterp[targ.nonzero()[0][-1]:]=0
    #h.dprint(labInterp)
    
    #do not consider 0's while taking mean
    labInterp[labInterp==0]=np.nan
    targ[targ==0]=np.nan



    #plt.figure(figsize=(12,6))
    #plt.plot(1000*labInterp[targ!=0]-rmsx,'k-')
    #plt.plot(targ[targ!=0]-rmsy,'g-')
    #plt.show()
    #plt.close()
    
    rmsx = np.nanmean(labInterp)
    rmsy =  np.nanmean(targ)

    #place the 0's back
    labInterp[np.isnan(labInterp)]=0
    targ[np.isnan(targ)]=0




    #plt.figure(figsize=(12,6))
    #plt.plot(1000*labInterp[targ!=0]-rmsx,'k-')
    #plt.plot(targ[targ!=0]-rmsy,'g-')
    #plt.show()
    #plt.close()
    
    #THIS IS THE CROSS CORRELATION SECTION
    ##
    ##correlate must have same sized arrays input
    correlation = signal.fftconvolve(targ[targ!=0]-rmsy,(labInterp[targ!=0]-rmsx)[::-1],'full')
    #correlation = np.correlate(targ[targ!=0]-rmsy,(labInterp[targ!=0]-rmsx),'full')
    #length of the correlation array is length input array times 2 plus 1
    #if the two arrays are already aligned then the peak of correlation function should be middle
    middle = int((len(correlation)-1)/2)
    
    
    #width is used for local maximum finding
    #1 is a number that worked here for smarts and NRES data. MAY NEED TO BE ADJUSTED for new data
    width = int(1/dLam)
    
    #max value is the index of the maximum value(local max around middle of array if width used)
    mval = middle-width+np.argmax(correlation[middle-width:middle+width])
    
    
    #want to make a quadratic around mval to be more precise with 'peak' of correlation
    #if we just use max val, precision is limited to grid spacing.
    #need to do poly only in certain range around center because wings will take over the fit
    fitWidth = 5
    
    xRange = mval + np.arange(2*fitWidth)- fitWidth
    polyFunc = np.polyfit(xRange, correlation[mval-fitWidth:mval+fitWidth],2)

    #corr = correlation[mval-fitWidth:mval+fitWidth]
    #h.dprint('xRange: ' + str(len(xRange)))
    #h.dprint('corr: ' + str(len(corr))) 
    #polyFunc = np.polyfit(xRange, corr,2)
    
    #take the derivative of the function and set equal to 0 for the peak.
    #f=a*x^2 +b*x+c
    #f' = 2*a*x+b = 0 -> x = -b/2a 
    xVal = -polyFunc[1]/(2*polyFunc[0])
    
    if h.debug:
        print(xVal)
        #ffit = np.polyval(polyFunc,xRange)
        plt.figure()
        #plt.xlim(mval-fitWidth*3,mval+fitWidth*3)
        #plt.plot(xRange, ffit,'g-')
        plt.plot(range(len(correlation)), correlation, 'k-')
        plt.show()
        plt.close()
    #mval= np.argmax(ffit)
    
    
    #the actual lamda offset is how far from middle we are in pixel space times the
    #pixel to grid ratio
    offset = (xVal-middle)*(tarGrid[1]-tarGrid[0])
    h.dprint('offset '+str(offset))

    #tmpMax = np.argmax(out)
    #print('index of maximum: ' + str(tmpMax) + ' and adjusted delLam: ' + str(tmpMax/len(out)))
    #print(out)
    if h.debug:
        
        #used for printing/debuging
        #scale = np.mean(gausdTarg)/np.mean(labInterp)
        #print('offset: ' + str(offset) + ' and scale: ' + str(scale))

        #fig, ax = plt.subplots(figsize=(25,5))
        #ax.ticklabel_format(useOffset=False)
        #plt.title("Unadjusted stellar spectra over reference spectra")
        #plt.xlabel("Wavelength (nm)")
        #plt.ylabel("Scaled irradiance")
        #plt.xlim([385,405])
        #plt.ylim([0,2800])
        #plt.plot(tarGrid, gausdTarg, 'g-')
        #plt.plot(tarGrid,labInterp*scale, 'k-')
        #plt.show()
        #plt.close()

        fig, ax = plt.subplots(figsize=(25,5))
        ax.ticklabel_format(useOffset=False)
        plt.title("Unadjusted stellar spectra over reference spectra")
        plt.xlabel("Wavelength (nm)")
        plt.ylabel("Scaled irradiance")
        #plt.xlim([385,405])
        #plt.ylim([0,2800])
        plt.plot(tarGrid-offset, targ, 'g-')
        if tmpDebug is 1:
            plt.axvline(x=h.cakLam, color='blue')
        if tmpDebug is 0:
            plt.axvline(x=h.cahLam, color='blue')
        #SCALE JUST FOR VIEWING
        #plt.plot(tarGrid,labInterp*scale, 'k-')
        plt.show()
        plt.close()


    
    
    
    #fig, ax = plt.subplots(figsize=(25,5))
    #ax.ticklabel_format(useOffset=False)
    #plt.title("Correlated stellar spectra over reference spectra")
    #plt.xlabel("Wavelength (nm)")
    #plt.ylabel("Scaled irradiance")
    
    #plt.xlim([392,398])
    #plt.ylim([0,2800])
    
    #plt.plot(tarGrid-offset, gausdTarg, 'g-')
    #plt.axvline(x=393.369, color='blue')
    #plt.axvline(x=396.85, color='blue')

    #plt.plot(tarGrid-offset, gausdTarg, 'k-',color='green')
    #SCALE JUST FOR VIEWING
    #plt.plot(tarGrid,labInterp*scale, 'k-')
    #plt.show()
    #plt.close()
    
    return offset,labInterp


#Takes in the data array with reference spectra and set name, reference spectra and setname used 
#for printing.
#This function goes through the data array for each observation and finds any other obs on that day
#It will combined these observations into one targOLap and run the rest of the pipeline on the 
#averaged data to give a better signal to noise for these observations.
#the output is a longer data array to be used later for printing the final 
#This is how the data structure is created and maintained:
#E
#correlation[0] is the delta lamda to place spectra in lab frame and lamGrid is the lab frame
def sum_daily_data(inData,starName,labSpec):  
    rv = 0
    #array to hold which MJD are done
    done = []

    #make a clean copy
    allData = np.asarray(inData)

    datesList = []
    sortedList = []


    for curD in allData:
        #find the closest MJD's to the current

        sameD = [curD]#initialize the same day list 

        curMjd = curD.mjd #dont double up if we've alrady DONE this day
        print(str(curMjd)+' has ' +str(curD.nOrd))
        if(curMjd in done):
            continue

        #grab the same day's data(will happen for multiple obs)
        for t in allData:
            if t.hour != curD.hour and t.day == curD.day:
                if t.nOrd is not curD.nOrd:
                    continue
                sameD.append(t)
            

        #if there is no other obs on this day or if cur
        if len(sameD) is 1: 
            done.append(curMjd)#must be the MJD not class obj
            continue
        


        header = curD.header#for faking the data
        mjd = curD.mjd

        if curD.header is not None:
            site = curD.site 
        else:
            for d in sameD:
                done.append(d.mjd)
            print('sun_daily_data missing header very bad')
            continue#go past this day something very wrong
        
        targList = []
        lamGrid = curD.lamGrid#use the same lamGrid for all of the observations to average
        #combinedWin = curD.window
        windowList =[]

        
        #combining targs
        #this is the big deal, increases signal to noise for each days spectra
        for day in sameD:
            #offset = day.offset
            ##dWindow = day.window
            #dLam = day.lamGrid

            targList.append(day.targOlapf)
            windowList.append(day.window)
#adam add the current day to this too
#############CAUTION here be dragons
            #Going to add to the target overlap FOR ONLY the windows around each line and band
            #for i in range(len(dWindow)):

                
                #for the combined window we will take the min value of each of the windows
                #so for a triangle window it may be a bit smaller than normal
                # potentially max is better here? Ricky Q TODO 
                #combinedWin[i] = [max(value) for value in zip(combinedWin[i], dWindow[i])]

                #need to interpolate because all the observations have slightly different lamda grids
                #put them all onto a blank spectra
                #interpfunc = interpolate.interp1d(dLam-offset[i], day.targOlapf, kind='linear',fill_value='extrapolate')    
                #combinedTarg = combinedTarg + interpfunc(lamGrid)
                        
            
            #these two are now done
            done.append(day.mjd)

        shkRet = multi_window_calc_shk(lamGrid,targList,curD.raw(),windowList)
        #end combining loop
        #wins = hk_windows(lamGrid,0,1)
        #print(np.shape(lamGrid))
        #print(np.shape(inData[0].lamGrid))
        #print(np.shape(combinedTarg))
        ##print(np.shape(inData[0].targOlapf))
        #print(max(combinedWin[0]))
        #print(max(inData[0].window[0]))
        
        #find SHK with new offset to lamda grid
        #shkRet = calc_shk(lamGrid, combinedTarg, curD, wins)

        shk = shkRet[0]
        windows = curD.window

        combData = h.analyzedData(curD.raw(),lamGrid,curD.flat,curD.targOlapf,shk,curD.offset,True,False,windows)


        outputDir = combData.outputDir()
        decimalYr = combData.decimalYr

        print(combData.label())
        

        
        #OUTPUT FOR FURTHER PRINTING
        if h.pdfMode == 0:
            plot.pdf_from_intermediate_data(curD.lamGrid, labSpec,combData,.3)

        #h.mkdir_p(outputDir)
        inData.append(combData)
    return inData