#this file holds pipeline functions which are not telescope dependant

#imports from a lot of data the lab frame spectra so we may cross-correlate onto this and
#find better locations of the ca hk lines
#taken from Ricky Egeland's code, probably could be optimized since he
#did more with this code than we need

#TODO determine if this needs change for other uses 
def import_aligning_spectra(fluxdir, minwl=None, maxwl=None, resolution=1, residual=False):
    import numpy as np
    import matplotlib as mpl
    from matplotlib import pyplot as plt
    import os
    import os.path
    import re
    import scipy
    import scipy.signal
    import scipy.ndimage.filters
    import helpers as h#for constants
    
    
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
def calc_del_lam(labGrid, lab, tarGrid, targ, smooth) :
    import scipy as sc
    #import astropy.io.fits
    import numpy as np
    import helpers as h#for constants
    import numpy.polynomial.polynomial as poly
    #import os
    #from calc_shk import calc_shk
    #from calc_shk import calc_targOlapf
    #from mk_flatolap import mk_flatolap
    from matplotlib import pyplot as plt
    from astropy.convolution import convolve, Box1DKernel
    from scipy import interpolate
    
    tmpGridScale = 1
    dLam = tarGrid[1] - tarGrid[0]
    #print(np.shape(targ))
    gausdTarg = sc.ndimage.filters.gaussian_filter(targ,smooth/dLam/2)
    
    dLabLam = labGrid[1] - labGrid[0]
    
    #TODO SHOULD 2.55 be h.sigToFWHM?
    gausedLab = sc.ndimage.filters.gaussian_filter(lab,(dLam/dLabLam)/2.55)

    #get the lab spectrum into angs div by 10 on angstrom grid for our purposes
    interpfunc = interpolate.interp1d(labGrid, gausedLab, kind='linear')#,fill_value='extrapolate')
    labInterp=interpfunc(tarGrid)



    #ZERO OUT THE EDGES OF OUR LAB SPECTRA
    #TODO do this with strict values to remove edge errors which may affect correlation
    #from 0 to first nonzero element of targolap
    labInterp[:targ.nonzero()[0][0]]=0
    #from last nonzero element to end
    labInterp[targ.nonzero()[0][-1]:]=0

    
    #do not consider 0's while taking mean
    labInterp[labInterp==0]=np.nan
    targ[targ==0]=np.nan
    
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
    correlation = np.correlate(targ[targ!=0]-rmsy,(labInterp[targ!=0]-rmsx),'full')
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
    
    #take the derivative of the function and set equal to 0 for the peak.
    #f=a*x^2 +b*x+c
    #f' = 2*a*x+b = 0 -> x = -b/2a 
    xVal = -polyFunc[1]/(2*polyFunc[0])
    
    #print(xVal)
    #ffit = np.polyval(polyFunc,xRange)
    #plt.figure()
    #plt.xlim(mval-fitWidth*3,mval+fitWidth*3)
    #plt.plot(xRange, ffit,'g-')
    #plt.plot(range(len(correlation)), correlation, 'k-')
    #plt.show()
    #plt.close()
    #mval= np.argmax(ffit)
    
    
    #the actual lamda offset is how far from middle we are in pixel space times the
    #pixel to grid ratio
    offset = (xVal-middle)*(tarGrid[1]-tarGrid[0])
    #print('offset: ' + str(offset))

    
    
    #used for printing/debuging
    #scale = np.mean(gausdTarg)/np.mean(labInterp)

    #tmpMax = np.argmax(out)
    #print('index of maximum: ' + str(tmpMax) + ' and adjusted delLam: ' + str(tmpMax/len(out)))
    #print(out)
    #fig, ax = plt.subplots(figsize=(25,5))
    #ax.ticklabel_format(useOffset=False)
    #plt.title("Unadjusted stellar spectra over reference spectra")
    #plt.xlabel("Wavelength (nm)")
    #plt.ylabel("Scaled irradiance")
    #plt.xlim([392,398])
    #plt.ylim([0,2800])
    #plt.plot(tarGrid, gausdTarg, 'g-')
    #plt.axvline(x=393.369, color='blue')
    #plt.axvline(x=396.85, color='blue')
    #SCALE JUST FOR VIEWING
    #plt.plot(tarGrid,labInterp*scale, 'k-')
    #plt.show()
    #plt.close()
    
    
    
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
    
    return offset,targ,labInterp,gausdTarg


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
    from scipy import interpolate
    import numpy as np
    from calc_shk import calc_shk
    from astropy.time import Time
    import plotting as plot
    import helpers as h
    import numpy as np
    
    rv = 0
    done = []#array to hold which MJD are done

    #make a clean copy
    allData = np.asarray(inData)

    datesList = []
    sortedList = []


    for curD in allData:
        #find the closest MJD's to the current

        
        #sameD = None
        sameD = [curD]#initialize the same day list 
        curMjd = curD.mjd #dont double up if we've alrady DONE this day

        if(curMjd in done):
            continue

        #grab the same day's data(will happen for multiple obs)
        for t in allData:
            if t.hour != curD.hour and t.day == curD.day:
                #sameD=t
                sameD.append(t)
                #break

        #if there is no other obs on this day or if cur
        #if sameD is None:
        if len(sameD) is 1: 
            done.append(curMjd)
            continue
        


        header = curD.header#for faking the data
        mjd = curD.mjd

        if curD.header is not None:
            site = curD.site 
        else:
            for d in sameD:
                done.append(d)
            print('sun_daily_data missing header very bad')
            continue#go past this day something very wrong
        
        combinedTarg = np.zeros(len(curD.targOlapf))
        lamGrid = curD.lamGrid#use the same lamGrid for all of the observations to average
        #sameMjd = sameD.mjd
        for day in sameD:
        
            #need to interpolate because all the observations have slightly different lamda grids
            #put them all onto the referece spectra's grid held in curLamGrid
            interpfunc = interpolate.interp1d(day.lamGrid-day.offset, day.targOlapf, kind='linear',fill_value='extrapolate')
            combinedTarg+=interpfunc(lamGrid)
            

            #combinedTarg = curTarg        

            #The issue is the two targetOlaps are on different grids. 
            #Interp them on the same grid then add
            #interpfunc = interpolate.interp1d(sameD.lamGrid-sameD.offset, sameD.targOlapf, kind='linear',fill_value='extrapolate')
            #sameTarg=interpfunc(curD.lamGrid)

            #this is the big deal, increases signal to noise for each days spectra
            #combinedTarg = combinedTarg+sameTarg

                    
            #these two are now done
            done.append(day)
            #done.append(curMjd)
            #done.append(sameMjd)

        #end combining loop

        #find SHK with new offset to lamda grid
        shkRet = calc_shk(lamGrid, combinedTarg,curD)

        shk = shkRet[0]
        windows = shkRet[1]

        combData = h.analyzedData(curD.raw(),lamGrid,[],combinedTarg,shk,0,True,False,windows)


        outputDir = combData.outputDir()
        decimalYr = combData.decimalYr

        print(combData.label())
        

        
        #OUTPUT FOR FURTHER PRINTING
        if h.pdfMode == 0:
            plot.pdf_from_intermediate_data(curD.lamGrid, labSpec,combData,.3)

        #h.mkdir_p(outputDir)
        inData.append(combData)
    return inData