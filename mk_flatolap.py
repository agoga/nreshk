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