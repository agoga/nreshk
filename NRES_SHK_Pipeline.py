#NRES SHK Specific pipeline. 
# Order of pipeline and folder structure taylored for nres_shk

import scipy.constants as sc
import astropy.io.fits
import numpy as np
import os
import glob
import helpers as h
import logging

from calc_shk import calc_targOlapf
from calc_shk import calc_shk
from calc_shk import smart_calc_shk
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
    

def old_multifile_NRES_to_data(filePath):
    #print('old format')
    #load the data in the old format
    import os
    import astropy.io.fits

    if file.endswith("-noflat.fits"):
        fileName = filePath
        specHDu = astropy.io.fits.open(os.path.join(path,dirs,file))
    elif file.endswith("-wave.fits"):
        wvHDu = astropy.io.fits.open(os.path.join(path,dirs,file))

    spec = specHDu[0].data[0]
    header = specHDu[0].header
    
    
    #need the length here because sometimes the wave grid is in the first element
    #of this array but most of the time it's not 
    if len(wvHDu)>1 and type(wvHDu[1]) == astropy.io.fits.hdu.table.BinTableHDU:
    #print('wavehdu')   
    #print(wvHDu.info())
    #print(wvHDu[1].data.dtype)
    #print(type(wvHDu[1].data))
    #print('done')
    #print('stuff ' + str(i))
    #print((wvHDu[1].data)[0])
    
        #there are two different names for this data, it's mostly in the first name but
        #sometimes in the second
        try:
            waveGrid = wvHDu[1].data['WavelenStar'][0]
        except:
            waveGrid=wvHDu[1].data['Wavelength'][0]
            
    elif len(wvHDu)>1:
        waveGrid = (wvHDu[1].data)[0]
    else:
        #print(wvHDu[0].data.dtype)
        waveGrid = (wvHDu[0].data)[0]
        
    specHDu.close()
    wvHDu.close()
    
#This is the NRES specific function goes through the current folder(which should be one observation)
#and loads, from fits files, the wavelength grid, spectra, header information
#This function is only needed when there is both old and new data from NRES and we don't know which
#is which. When all the data is in the new format then half this function can be removed.
def load_obs_for_pipeline(obsFileName):
    import os
    import astropy.io.fits
    
    waveGrid = []
    spec = []
    header = []
    fileName = ''#spec file name
    oldFormat = False


    splitF = os.path.splitext(obsFileName)
    #create the no flat file to check it's existence
    noFlat = splitF[0]+'-noflat'+splitF[1]
#see if the folder contains a -noflat or -wave file
    #There may be a faster way to determine if this folder is old or new but I don't see a safer way

    if os.path.exists(noFlat):
        print('wowz')
        oldFormat = True
        
    

    #we've now looped through all the files if any are old then we try old format
    if oldFormat:
        return
        old_multifile_NRES_to_data(obsFileName)
    else: #ELSE IT'S NEW FORMAT and super easy
        oHDu = astropy.io.fits.open(obsFileName)
        waveGrid = oHDu[7].data
        spec = oHDu[1].data
        header = oHDu[0].header
        oHDu.close()
   #end big if              
    #if len(waveGrid) == 0 or len(spec) == 0 or len(header) == 0:
    #   print('BIG PROBLEM IN load_folder_for_pipeline')
                    
    return waveGrid, spec, header, obsFileName, oldFormat


# This routine makes a wavelength grid spanning lamran[2], with constant
# wavelength interval dlam (nm).  It reads an extracted flat-field file
# filin from the reduced/flat directory appropriate to input site, and
# the corresponding wavelength scale from the reduced/thar directory.
# It trims and smooths the 4 bluemost orders of the flat field, and then
# interpolates and sums the 4 flats onto the wavelength grid, taking account
# of the effect of varying dlambda/dx in the original spectra.  Results are
# written to a FITS file named in outfile (full pathname), and a new line
# describing this output is added to $STAGE2ROOT/flatolap.csv.
def mk_flatolap(lam, flat, idl=''):
    #WORKING MK_FLATOLAP
    ##PORT OF DR. TIM BROWN'S NRES HK CODE 
    ##comments marked brown are Dr. Brown's
    
    import numpy as np
    import scipy.io as sc
    import sys
    from matplotlib.backends.backend_pdf import PdfPages

    from matplotlib import pyplot as plt

    from astropy.convolution import convolve, Box1DKernel



    ##intializations and hardcoded inputs TODO fix hardcoded?
    ##mk_flatolap
    lamRan=[380.,420.]
    dLam =0.001
    gOrd=[63,64,65,66]
    nGord= len(gOrd)
    nx=4096
    bounds=[[615,3803],[670,3770],[733,3740],[750,3660]]

    bad = False
    #read filesfits
    #TRANSPOSE THE STRAIGHT DATA BC IDL CODE HAS REVERSED DIMENSIONS
    #TODO FIX THIS AND ALL ARRAY SHAPES TO FOLLOW PYTHON LOADING
    lam=np.transpose(lam)

    #brown
    #make wavelength grid
    nLam = ((lamRan[1]-lamRan[0])/dLam)+1
    lamGrid = lamRan[0]+dLam*np.arange(nLam,dtype=np.float64)

    #brown
    #make wavelength derivative, scale function to assure integral flux*dlambda
    #is preserved in transformation to constant wavelength bins.

    #verify data types https://www.harrisgeospatial.com/docs/PythonDataConvert.html
    dLambx = np.zeros((nx,nGord),dtype=np.float64)
    scale = np.zeros((nx,nGord),dtype=np.float64)



    #BUG potentially the IDL code is finding the flat for the first 4 orders 
    #when we want to use the last 4
    for y in range(len(gOrd)) :
        #TODO first and last in each order are not accurate w idl code and precision is less
        #dlamb only accurate to 3 sigfig
        #scale is accurate against idl to about 3 decimal places
        dLambx[:,y] = np.gradient(lam[:,y]) 
        scale[:,y] =  dLambx[:,y]/dLambx[int(nx/2),y]


    #brown 
    #isolate the desired orders, set contents to zero outside boundaries.
    

    #TODO mnad;lsad;lsadksadjkl transpose hack blargggggggg
    gFlat = np.transpose(flat[0,gOrd,:])
    sgFlat = np.zeros((nx,nGord),dtype=np.float64)
    #print(np.shape(gFlat))
    #print(np.shape(sgFlat))

    stuff = []
    for i in range(len(gOrd)) :
        #bug cant do bounds[i,0] for some reason
        gFlat[0:bounds[i][0],i]=0
        gFlat[bounds[i][1]:,i]=0
        #BOX CAR SMOOTHING INSTEAD OF IDL SMOOTH()
        #https://joseph-long.com/writing/AstroPy-boxcar/
        sgFlat[:,i]=convolve(convolve(convolve(gFlat[:,i]*scale[:,i], Box1DKernel(25)), Box1DKernel(25)), Box1DKernel(25))
        #fig = plt.figure()
        #plot python data
        #plt.plot(gFlat[:,i]*scale[:,i], 'k-')
        #plt.ylabel('sgFlats')
        #plt.show()
        #plt.close()
    
    

    flatOlap = np.zeros(int(nLam),dtype=np.float64)    



    #brown
    #interpolate onto flatolap

    #https://stackoverflow.com/questions/18326714/idl-interpol-equivalent-for-python
    from scipy import interpolate

    #need extrapolate likely due to edge cases of lamda grid not being perfectly aligned. 
    #since linear spacing should be fine?
    for i in range(len(gOrd)) :
        interpfunc = interpolate.interp1d(lam[:,gOrd[i]],sgFlat[:,i], kind='linear', fill_value='extrapolate')
        flatOlap = flatOlap+interpfunc(lamGrid)

    #brown
    #make output data array -- [lambda,flatolap].  Write it out.
    #output = np.zeros((int(nLam),2),dtype=np.float64)
    #print(np.shape(output))
    #output[:,0] = lamGrid    
    #output[:,1] = flatOlap  

    #plt.figure(figsize=(30,5))
    #plt.style.use('classic')
    #plot python data
    #plt.plot(lamGrid, flatOlap, 'k-')
    #plt.title('Orders 63-67 of flat field added togeather')
    #plt.xlim([391.5,408])
    ##plt.xlabel('wavelength [nm]')
    #plt.ylabel('flatOlap')
    #plt.show()
    #plt.close()



    #idlHdu = astropy.io.fits.open(idlfilepath+idlFlatFile)
    if idl != '':
        idlData = idl['flatolap']#idlHdu[0].data[1]#
        plt.plot(lamGrid, idlData, 'k-', color='blue')
        plt.plot(lamGrid, abs(flatOlap-idlData), color='green')
        
    #this returns to the pipeline that we have a bad flat    
    if max(flatOlap) > 1.5:
        print("bad flat detected mk_flatolap.py")
        bad = True
        #return [],[]
        
    ##OUTPUTS
    #lamGrid-the x val of all these plots
    #flatOlap-the overlapped flat on these ranges
    return lamGrid, flatOlap, bad

    #this is the BIG wrapper. This handles all the folder structure stuff and manages the
    #pipeline. The general outline of pipeline is(search tag to find general location);
    #
    #TAG_0_ import data to arrays
    #TAG_1_ create flat field
    #TAG_2_ remove instrumental error by dividing out flat field
    #TAG_3_ align stellar spectrum to reference by cross correlation
    #TAG_4_ find ca HK emission reversal and calculate s-index(shk)
    #TAG_5_ detect bad spectra and erronous data and remove from final product
    #TAG_6_ Sum spectra from each night(2-3 obs per night) to increase signal/noise
    #TAG_7_ Calculate S-index for each nights combined data
    #TAG_8_ print intermediate pdfs(for single and combined obs) for debugging into output folder
    #TAG_9_ plot final time series of star with all good data included
    #
    #loop through each star as desired
def NRES_SHK_Pipeline(dataPath,outputPath,flatDict,lab,badD,forceRun):
    starName = ''

    #resolution for printing
    res = .01
    label =''

    stars = h.get_immediate_subdirectories(dataPath)
    
    #loop through every folder in the data files path. 
    #these folders hold all the data for each star
    for s in stars:

        if not h.is_folder_star(s):
            print('bad folder \'' + s + '\'')
            continue
        #if the star has an output folder and isnt forced to run we won't run again
        if os.path.exists(outputPath+s+'/') and s not in forceRun:
            print("HD " + s + " already has output")
            continue
        
        print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        print("Running pipeline on HD " + s)
        starName = s
        #need blank arrays to append to since we don't know length due to bad spectra
        #data = [[],[],[],[]]
        tryData =[]

        dailyData = dict()

        starPath = os.path.join(dataPath,s)
        obsFiles = []   

        #Find all the files that may be correct
        obsUn = glob.glob(starPath+'\\**\\*-e91.fits', recursive=True)
        obsFiles = [r for r in obsUn if not "output" in r]

        #star loop
        for obsFileName in obsFiles:
            print('-------------------------------------------------------------------------------')
            #curFiles = os.listdir(os.path.join(path,dirs))
            #print('loading ' + dirs)
            #print('current files')
            #print(curFiles)

            #TAG_0_
            folderRet = load_obs_for_pipeline(obsFileName)
            if folderRet is None:
                continue
            waveGrid = folderRet[0]#key for pipeline
            spec = folderRet[1]#key for pipeline
            header = folderRet[2]
            obsFileName = folderRet[3]
            old = folderRet[4]#is old format or not, for debugging
        


            #if you have data from a different system you may just want to simulate the header
            #dictionary.
            mjd = header['MJD-OBS']
            site = header['SITEID']
            date = header['DATE-OBS']

            
            label=str(mjd).replace('.','/')

            tempEff = h.tEffLookup[starName.strip('/')]

            #title for pdfs
            title = ''
            #BEGIN FORMAT INDEPENDANT CODE(HOPEFULLY)



            #try to find a flat otherwise fail
            #try:
            fD = flatDict[site]
            #flatDates = list(flatDict.keys())
            #flatDates.sort(key=lambda x: abs(x - mjd))
            fK = h.closestKey(fD,float(mjd))

            flatFilePath = fD[fK]
            
            ##8/21/2020 was fHDu = astropy.io.fits.open(os.path.join(path,flatFilePath))
            fHDu = astropy.io.fits.open(os.path.join(starPath,flatFilePath))

            flat = fHDu[0].data#key for pipeline

            #if abs(mjd-fK) > 25:
            #    print('BAD BAD closest for: '+str(mjd) +' is '+ str(abs(mjd-fK)))
            #elif abs(mjd-fK) > 2:
            #    print('closest for: '+str(mjd) +' is '+ str(abs(mjd-fK)))
            
            #TAG_1_
            #give multiple arrays of flats whose lam values are stored in multiple wave grid arrays
            flatRet = mk_flatolap(waveGrid, flat)

            #except:
            # print('bad flat find with site: ' + site + ' and date ' + str(obsDate))
                #continue


            #return one flat array with lambda grid
            flatOlap = flatRet[1]
            lamGrid = flatRet[0]

            #TAG_2_
            #get the target data minus the flat
            targOlapf = calc_targOlapf(lamGrid, waveGrid, spec, flatOlap, label)

            #TAG_3_
            #cross correlation returns a lambda offset(dlam) and save the lab spectra
            correlation = pipe.calc_del_lam(lab[0]/10,lab[1], lamGrid, targOlapf,res)
            dLam = correlation[0]
            labSpec = correlation[2]


            #radial velocity calculation
            #lamRef = 396.85
            #delta lamda / ref lamda * speed of light
            # rv = dLam/ lamRef * sc.c 
            #rv from meters to km/s as desired by hk_windows

            #atm we're not going to apply the radial velocity and just use the adjusted spectra
            rv = 0#rv/10000 


            
            #TAG_4_
            #find SHK with new offset to lamda grid
            shkRet = calc_shk(lamGrid-dLam, targOlapf, rv, teff=tempEff)
            shk = shkRet[0]
            windows = shkRet[1]
            


           
            #if old:
            #   shk=shk*2
            #print(title)
            
            #TAG_5_
            #time to toss bad spec so they won't be summed
            badSpec = False
            #not included but may need to be, check if any values of targolapf
            #are negative. Makes sense to me that those spectra should be tossed
            badSpec = h.bad_spec_detection_v2(lamGrid-dLam,targOlapf)          
            
            #among other things?!
            if shk < 0 or shk > 1 or flatRet[2] == True:
                badSpec = True          

            if badSpec:
                badD.append(mjd)

            

            title = ''

            #create the directories for pdf plotting and save every intermediate data array
            #first = label.split("/")[0]
            #second = label.split("/")[1]

            #this observations data and printer
            obsP = h.obsPrinter(mjd,shk,dLam,header,site,date,starName,False,badSpec,windows)
            obsD = h.specData(mjd,header,lamGrid,targOlapf,dLam,shk,False)
            outputDir = obsP.outputDir()
 
            

            h.mkdir_p(outputDir)

            
            decimalYr = obsP.decimalYr
            dataPath=outputDir+ obsP.hour+"_data"


            #OUTPUT FOR FURTHER PRINTING
            np.savez(dataPath, targOlapf=targOlapf,flatOlap=flatOlap, lamGrid=lamGrid, adjLamGrid=lamGrid-dLam,windows=windows)
            plot.pdf_from_intermediate_data(lamGrid, labSpec,lamGrid-dLam, targOlapf, windows,obsP.pdfTitle(), outputDir ,obsP.hour,flatOlap,.3)

            if(badSpec):
                print('bad: ' + str(mjd))
                continue

            tryData.append(obsD)
            
        #end of star loop 
        
        #major portions of the pipeline are kept in these two functions for readability.
        #We could likely functionize from individual observation folder-> data output as well
        #but if people are interested in putting data through this pipeline from other telescopes
        #many parts of the above code must be changed. The below would not need to be changed.
        #TAG_6_

        #tryData = pipe.sum_daily_data(tryData,starName,labSpec)
        #data = pipe.sum_daily_data(starName,data,labSpec)
        #TAG_7_
        #TAG_8_
        #TAG_9_
        plot.plot_daily_data_timeseries(tryData,starName,badD)
        #pipe.plot_daily_data_timeseries(data,starName,bad)
        
        # 
        #we're going to loop through all folders in a stars data and determine if each is 
        #old or new format once thats determined we get the wave/flat/spec data differently
        #then we should not care about the format anymore and will do the pipeline as normal
        #for path, subdirs, files in os.walk(dataPath + starName):
            #walk through each folder this is the main loop over all observations
            #each dir here should only contain one observation
            #print(subdirs)

            
    print(badD)

                