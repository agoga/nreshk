# NRES_SHK_Pipeline.py 
# Adam Goga
# 
# Intended to be NRES SHK specific pipeline functions though that is not likely the case currently. 
# Order of pipeline and folder structure taylored for nres_shk
#
#
#def NRES_SHK_MkFlat - wrapper function for eventual command line use
#def create_flat_dict_file - reads all the flats and creates a disctionary for pipeline
#def old_multifile_NRES_to_data - Reads old NRES specific observations in .fits format 
#def load_obs_for_pipeline - Reads NRES specific observations in .fits format(will call old_mltifile if detected)
#def mk_flatolap - combines the proper orders of flat file to make flat overlap
#def NRES_SHK_Pipeline - the meat and potatoes which handles interactions from raw folders to analyzed data


import sys
import os
import glob
import logging
import pickle
import pprint
import astropy.io.fits
import time

import scipy.constants as sc
import numpy as np
from scipy import interpolate
from astropy.time import Time
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
from astropy.convolution import convolve, Box1DKernel
            

import helpers as h
import pipeline as pipe
import plotting as plot
from calc_shk import calc_targOlapf
from calc_shk import calc_shk
from calc_shk import old_calc_shk
from calc_shk import old_hk_windows
#This is the flat field dictionary creator.
#it must be run before the pipeline is run since we need the flatDict var
#do not run create_flat_dict_file unless you have added a new flat field to the flats folder
def NRES_SHK_MkFlat(flatPath,flatPickle):
    #UNCOMMENT TO RECREATE THE FLAT DICTIONARY
    create_flat_dict_file(flatPath,flatPickle)

    f = open(flatPath+flatPickle,"rb")
    flatDict = pickle.load(f)
    f.close()

    return flatDict

#This function should not be run every time the pipeline is run.
#It only needs to be run when you have added new flat files to the data
#A two dimensional dictionary is created for quick lookup of flat file 
#names by their site and MJD's and stored in a pickle file to be loaded
#in every time the pipeline is run
def create_flat_dict_file(flatDir,fileName):

    #find every flat file and put into dictionary
    flatDict = dict() 
    flatDict.update({'cpt':dict()})
    flatDict.update({'tlv':dict()})
    flatDict.update({'lsc':dict()})
    flatDict.update({'elp':dict()})


    #this loops through all the flat_....fits files in the flats directory
    #and creates the 2 dimensional dictionary  
    for path, subdirs, files in os.walk(flatDir):
        for name in files:
            if name.startswith("flat_") and name.endswith(".fits"):
                flHDu1 = astropy.io.fits.open(os.path.join(path,name))
                header = flHDu1[0].header
                
                mjd = float(header['MJD-OBS'])
                site = header['SITEID']
                
                siteDict = flatDict[site]
                siteDict.update({mjd:os.path.join(path,name)})
                
    f = open(flatDir+fileName,"wb")
    pickle.dump(flatDict,f)
    f.close()
    

def old_multifile_NRES_to_data(obsFileName):
    #print('old format')
    #load the data in the old format
    
    data=None
    splt = os.path.splitext(obsFileName)
    noFlat = splt[0]+'-noflat'+splt[1]
    wave  = splt[0]+'-wave'+splt[1]


    fileName = obsFileName
    specHDu = astropy.io.fits.open(noFlat)
    wvHDu = astropy.io.fits.open(wave)

    spec = specHDu[0].data[0]
    header = specHDu[0].header
    
    
    #there are two different names for this data, it's mostly in the first name but
    #sometimes in the second
    #need the length here because sometimes the wave grid is in the first element
    #of this array but most of the time it's not 
    if len(wvHDu)>1 and type(wvHDu[1]) == astropy.io.fits.hdu.table.BinTableHDU:
        try:
            waveGrid = wvHDu[1].data['WavelenStar'][0]
        except:
            waveGrid=wvHDu[1].data['Wavelength'][0]
            
    elif len(wvHDu)>1:
        waveGrid = (wvHDu[1].data)[0]
    else:
        #print(wvHDu[0].data.dtype)
        waveGrid = (wvHDu[0].data)[0]

    data = h.rawData(waveGrid,spec,header,obsFileName,True)  
    
    specHDu.close()
    wvHDu.close()
    return data
    
#This is the NRES specific function goes through the current folder(which should be one observation)
#and loads, from fits files, the wavelength grid, spectra, header information
#This function is only needed when there is both old and new data from NRES and we don't know which
#is which. When all the data is in the new format then half this function can be removed.
def load_obs_for_pipeline(obsFileName):    
    waveGrid = []
    spec = []
    header = []
    fileName = ''#spec file name
    oldFormat = False

    data=None
    splitF = os.path.splitext(obsFileName)
    #create the no flat file to check it's existence
    noFlat = splitF[0]+'-noflat'+splitF[1]
    #see if the folder contains a -noflat or -wave file
    #There may be a faster way to determine if this folder is old or new but I don't see a safer way

    if os.path.exists(noFlat):
        oldFormat = True
        

    #we've now looped through all the files if any are old then we try old format
    if oldFormat:
        return old_multifile_NRES_to_data(obsFileName)
    else: #ELSE IT'S NEW FORMAT and super easy
        
        oHDu = astropy.io.fits.open(obsFileName)
        waveGrid = oHDu[7].data
        spec = oHDu[1].data
        header = oHDu[0].header
        data = h.rawData(waveGrid,spec,header,obsFileName,oldFormat)  
        oHDu.close()
   #end big if              
    #if len(waveGrid) == 0 or len(spec) == 0 or len(header) == 0:
    #   print('BIG PROBLEM IN load_folder_for_pipeline')
         
    return data


# This routine makes a wavelength grid spanning lamran[2], with constant
# wavelength interval dlam (nm).  It reads an extracted flat-field file
# filin from the reduced/flat directory appropriate to input site, and
# the corresponding wavelength scale from the reduced/thar directory.
# It trims and smooths the 4 bluemost orders of the flat field, and then
# interpolates and sums the 4 flats onto the wavelength grid, taking account
# of the effect of varying dlambda/dx in the original spectra.
def mk_flatolap(raw, flat):
    #WORKING MK_FLATOLAP
    ##PORT OF DR. TIM T.Brown'S NRES HK CODE 
    ##comments marked T.Brown are Dr. T.Brown's

    lam = raw.waveGrid
    highOrd = raw.nOrd#virtually always we want to go to highest(lowest wavelength) order

    gOrd = np.arange(h.lowGOrd,highOrd)


    #TODO Q for Ricky, this isnt right correct? We should be finding the actual lamda range from flat and using?
    lamRan=[380.,420.]
    dLam =0.001
    nGord= len(gOrd)
    nx=raw.nx

    #bounds provides the cut off for each order, anthing below the low and above the high 
    #index of each order in the flat file will be ignored
    #bounds=[[615,3803],[670,3770],[733,3740],[750,3660]]

    bounds =[]
    for i in range(nGord):
        bounds.append([600,3800])
    bad = False
    #read filesfits

    
    #T.Brown
    #make wavelength grid
    nLam = ((lamRan[1]-lamRan[0])/dLam)+1

    lamGrid = lamRan[0]+dLam*np.arange(nLam,dtype=np.float64)

    #T.Brown
    #make wavelength derivative, scale function to assure integral flux*dlambda
    #is preserved in transformation to constant wavelength bins.

    #TODO verify data types https://www.harrisgeospatial.com/docs/PythonDataConvert.html
    dLambx = np.zeros((nGord,nx),dtype=np.float64)
    scale = np.zeros((nGord,nx),dtype=np.float64)



    for y in range(len(gOrd)) :
        ##Potential BUG first and last in each order are not accurate w idl code and precision is less
        #UPDATE Not likely an issue
        #
        #dlamb only accurate to 3 sigfig
        #scale is accurate against idl to about 3 decimal places
        dLambx[y,:] = np.gradient(lam[y,:]) 
        scale[y,:] =  dLambx[y,:]/dLambx[y,int(nx/2)]


    #T.Brown 
    #isolate the desired orders, set contents to zero outside boundaries.
    gFlat = flat[0,gOrd,:]
    sgFlat = np.zeros((nGord,nx),dtype=np.float64)


    stuff = []
    for i in range(len(gOrd)) :
        gFlat[i,0:bounds[i][0]]=0
        gFlat[i,bounds[i][1]:]=0
        #BOX CAR SMOOTHING INSTEAD OF IDL SMOOTH()
        #https://joseph-long.com/writing/AstroPy-boxcar/
        sgFlat[i,:]=convolve(convolve(convolve(gFlat[i,:]*scale[i,:], Box1DKernel(25)), Box1DKernel(25)), Box1DKernel(25))      

    flatOlap = np.zeros(int(nLam),dtype=np.float64)    



    #T.Brown
    #interpolate onto flatolap

    #https://stackoverflow.com/questions/18326714/idl-interpol-equivalent-for-python
    
    #need extrapolate likely due to edge cases of lamda grid not being perfectly aligned. 
    #since linear spacing should be fine?
    for i in range(len(gOrd)) :
        interpfunc = interpolate.interp1d(lam[gOrd[i],:],sgFlat[i,:], kind='linear', fill_value='extrapolate')
        flatOlap = flatOlap+interpfunc(lamGrid)

    #T.Brown
    #make output data array -- [lambda,flatolap].  Write it out.

    #there should not be negatives so remove them
    flatOlap = flatOlap.clip(min=0)

    #this returns to the pipeline that we have a bad flat    
    if max(flatOlap) > 1.5:
        bad = True

    ##OUTPUTS
    #lamGrid-the x val of all these plots
    #flatOlap-the overlapped flat on these ranges
    #bad-if the flat is crap
    return lamGrid, flatOlap, bad



def find_and_make_flat(obs):
    return
    




#this is the BIG wrapper. This handles all the folder structure stuff and manages the
#pipeline. The general outline of pipeline is;
#
# import data to arrays
# create flat field
# remove instrumental error by dividing out flat field
# align stellar spectrum to reference by cross correlation
# find ca HK emission reversal and calculate s-index(shk)
# detect bad spectra and erronous data and remove from final product
# Sum spectra from each night(2-3 obs per night) to increase signal/noise
# Calculate S-index for each nights combined data
# print intermediate pdfs(for single and combined obs) for debugging into output folder
# plot final time series of star with all good data included
#
#loop through each star as desired

#Inputs
def NRES_SHK_Pipeline(dataPath,outputPath,flatDict,lab,skip,forceRun,manualAdj,only=None):
    """Main NRES SHK Function with various options to create time-series for input stars.\n
    dataPath- location of data folder\n
    outputPath- location of output folder\n
    flatDict- dictionary of flat files created with create_flat_dict_file\n
    lab- reference spectra to align against\n
    skip- List of MJD's to skip due to bad data or testing\n
    forceRun- stars to run even if they have an output folder already\n
    manualAdj- list of MJD's coupled with 4 offsets(analyzedData.offset) to manually adjust bad shifting\n
    """
    starName = ''

    #resolution for printing
    res = .01
    #label =''
    stars = h.get_immediate_subdirectories(dataPath)
    starData = {}#dictionary for returning from this function
    
    #we're going to loop through all folders in a stars data and determine if each is 
    #old or new format once thats determined we get the wave/flat/spec data differently
    #then we should not care about the format anymore and will do the pipeline as normal
    for s in stars:
        #need blank arrays to append to since we don't know length due to bad spectra
        obsFiles=[]
        analyzed=[]
        badD=[]

        if not h.is_folder_star(s):
            print('bad folder \'' + s + '\'')
            continue
        #if the star has an output folder and isnt forced to run we won't run again
        if os.path.exists(outputPath+s) and s not in forceRun:
            print(s + " already has output")
            continue
        
        print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        print("Running pipeline on HD " + s)
        print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        starName = s

        starPath = os.path.join(dataPath,s)
         
        #Find all the files that may be correct observations which we will cull during loading
        obsUn = glob.glob(starPath+'\\**\\*-e91.fits', recursive=True)
        obsFiles = [r for r in obsUn if not "output" in r]
        print('Number of observations: ' + str(len(obsFiles)))
        #star loop
        for obsFileName in obsFiles:
            
            obsRaw = load_obs_for_pipeline(obsFileName)
            if obsRaw is None:
                continue
            
            if only is not None:
                if obsRaw.mjd not in only:
                    continue
            
            if obsRaw.mjd in skip:
                continue
            
            print('-------------------------------------------------------------------------------')
            obsRaw.star=starName#sometimes the fits file will have a different star than the file name

            #if you have data from a different system you may just want to simulate the header
            #dictionary.
    
            #BEGIN FORMAT INDEPENDANT CODE(HOPEFULLY)

            #try to find a flat otherwise fail
            #try:
            if obsRaw.nOrd is 68:
                continue#TODO REMOVE THIS WHEN YOU WANNA DO NEW SPEC

            fD = flatDict[obsRaw.site]
            fK = h.closestKey(fD,float(obsRaw.mjd))

            flatFilePath = fD[fK]

            ##8/21/2020 was fHDu = astropy.io.fits.open(os.path.join(path,flatFilePath))
            fHDu = astropy.io.fits.open(os.path.join(starPath,flatFilePath))
            flat = fHDu[0].data
            fHDu.close()

            #if the flat and obs are many days apart. Should probably add this check back in
            #if abs(mjd-fK) > 25:
            #    print('BAD BAD closest for: '+str(mjd) +' is '+ str(abs(mjd-fK)))
            #elif abs(mjd-fK) > 2:
            #    print('closest for: '+str(mjd) +' is '+ str(abs(mjd-fK)))


            #give multiple arrays of flats whose lam values are stored in multiple wave grid arrays
            flatRet = mk_flatolap(obsRaw, flat)

            #return one flat array with lambda grid
            flatOlap = flatRet[1]
            lamGrid = flatRet[0]

            #get the target data minus the flat
            targOlapf = calc_targOlapf(obsRaw, lamGrid, flatOlap)


            #wins = old_hk_windows(lamGrid)
                                            #porque esto?
            #bigger windows to align with
            #aligningWins = pipe.create_shk_windows(lab[0]/10,lab[1], lamGrid, targOlapf, scale=10)
            offsets = pipe.find_window_offsets(lab[0]/10,lab[1], lamGrid, targOlapf, scale=12, smooth=res)
            
            integrationWins = pipe.create_integration_windows(lamGrid, targOlapf,offsets)

            #cross correlation returns a lambda offset(dlam) and save the lab spectra
            #quick fix TODO FIX
            labout = pipe.calc_del_lam(lab[0]/10,lab[1], lamGrid, targOlapf,res)
            labSpec = labout[1]
            #print(multiCorr)

            
            
            
            #ok now if we've been given specific offsets replace them now
            #if the MJD's are not unique(multiple could have same mjd) then this is borked
            try:
                i = manualAdj[0].index(obsRaw.mjd)
                offs = manualAdj[1][i]

                print('manual offset for mjd ' + str(obsRaw.mjd))
                #manual
                for o in range(len(offs)):
                    if offs is not None:
                        offsets[o]=offs[0]
            except:#we dont need this s
                pass #all gggggg
 
                

            #radial velocity calculation
            #lamRef = 396.85
            #delta lamda / ref lamda * speed of light
            # rv = dLam/ lamRef * sc.c 
            #rv from meters to km/s as desired by hk_windows
            #find SHK with new offset to lamda grid
            shkRet = calc_shk(lamGrid, targOlapf, obsRaw, integrationWins)

            shk = shkRet[0]
            windows = shkRet[1]

            
            # Time to toss bad spec so they won't be summed
            badSpec = False
            badReason =' Bad: '

            #not included but may need to be, check if any values of targolapf
            #are negative. Makes sense to me that those spectra should be tossed
            badSpec = h.bad_spec_detection_v2(lamGrid,targOlapf)          

            
            if badSpec:
                badSpec = True
                badD.append(obsRaw.mjd)         #fugly
                badReason+=(' bad spec detector.')
            #among other things?!
            elif shk < 0 or shk > 1:
                badD.append(obsRaw.mjd)
                badSpec = True  
                badReason+=('shk above 1 or below 0. ')
            elif flatRet[2] == True:
                badD.append(obsRaw.mjd)
                badSpec = True
                badReason+=(' bad flat. ')



            oData = h.analyzedData(obsRaw,lamGrid,flatOlap,targOlapf,shk,offsets,False,badSpec,integrationWins)

            outputDir = oData.starDir()
            decimalYr = oData.decimalYr
            
            print(oData.label())
            

            #OUTPUT FOR FURTHER PRINTING
            if h.pdfMode == 0:
                plot.pdf_from_intermediate_data(lamGrid, labSpec,oData,.3)


            if(badSpec):
                print(str(obsRaw.mjd)+badReason)
                continue

            analyzed.append(oData)
        #end of loop - for obsFileName in obsFiles:
        
        #daily averaging
        analyzed = pipe.sum_daily_data(analyzed,starName,labSpec)

        #save all the stars data to file for reloading 
        #just in case it hasnt been made
        h.mkdir_p(oData.starDir())

        #deprecated?! who knows
        t = [[],[]]
        for d in analyzed:
            t[0].append(d.mjd)
            t[1].append(d.shk)

        f = open(oData.starDir() + starName+'.pkl',"wb")
        pickle.dump(analyzed,f)
        f.close()

        #output bad data, should go into file which is checked
        print(badD)
        plot.plot_timeseries(analyzed,badD)
        
        
        starData.update({starName: analyzed})
    #end of loop - for s in stars: 
    return analyzed#pass back analyzed data           

                