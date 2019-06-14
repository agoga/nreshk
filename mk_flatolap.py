def mk_flatolap(lam, flat, idl=''):
    #WORKING MK_FLATOLAP
    ##PORT OF DR. TIM BROWN'S NRES HK CODE 
    ##comments marked brown are Dr. Brown's

    import numpy as np
    import scipy.io as sc
    import sys


    from matplotlib import pyplot as plt

    from astropy.convolution import convolve, Box1DKernel



    ##intializations and hardcoded inputs
    ##TODO take command line input???
    ##mk_flatolap
    lamRan=[380.,420.]
    dLam =0.001
    gOrd=[63,64,65,66]
    nGord= len(gOrd)
    nx=4096
    bounds=[[615,3803],[670,3770],[733,3740],[750,3660]]

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

    #print(dLambx[:,0])
    #print(scale[:,0]) 


    #brown 
    #isolate the desired orders, set contents to zero outside boundaries.



    #TODOmnad;lsad;lsadksadjkl transpose hack blargggggggg
    gFlat = np.transpose(flat[0,gOrd,:])
    sgFlat = np.zeros((nx,nGord),dtype=np.float64)
    #print(np.shape(gFlat))
    #print(np.shape(sgFlat))

    for i in range(len(gOrd)) :
        #bug cant do bounds[i,0] for some reason
        gFlat[0:bounds[i][0],i]=0
        gFlat[bounds[i][1]:,i]=0
        #BOX CAR SMOOTHING INSTEAD OF IDL SMOOTH()
        #https://joseph-long.com/writing/AstroPy-boxcar/
        sgFlat[:,i]=convolve(convolve(convolve(gFlat[:,i]*scale[:,i], Box1DKernel(25)), Box1DKernel(25)), Box1DKernel(25))

    flatOlap = np.zeros(int(nLam),dtype=np.float64)    



    #brown
    #interpolate onto flatolap

    #https://stackoverflow.com/questions/18326714/idl-interpol-equivalent-for-python
    from scipy import interpolate

    for i in range(len(gOrd)) :
        interpfunc = interpolate.interp1d(lam[:,gOrd[i]],sgFlat[:,i], kind='linear', fill_value='extrapolate')
        flatOlap = flatOlap+interpfunc(lamGrid)


    #if any(i > 0 for i in flatOlap) :
     #   print('gotya')

    #brown
    #make output data array -- [lambda,flatolap].  Write it out.
    #output = np.zeros((int(nLam),2),dtype=np.float64)
    #print(np.shape(output))
    #output[:,0] = lamGrid    
    #output[:,1] = flatOlap  

    #plt.figure()
    #plot python data
    #plt.plot(lamGrid, flatOlap, 'k-',color='red')
    #plt.xlabel('wavelength [nm]')


    #PLOT TEST IDL FLAT DATA

    #idlHdu = astropy.io.fits.open(idlfilepath+idlFlatFile)
    if idl != '':
        idlData = idl['flatolap']#idlHdu[0].data[1]#
        plt.plot(lamGrid, idlData, 'k-', color='blue')
        plt.plot(lamGrid, abs(flatOlap-idlData), color='green')
    ##OUTPUTS
    #lamGrid-the x val of all these plots
    #flatOlap-the overlapped flat on these ranges
    #
    #
    #max error BUGGED
    #idlFlat=data[1]
    #maxim = 0
    #for i in range(len(idlFlat)) :
    #    if idlFlat[i] > 0 :
    #        bla = abs(flatOlap[i]-idlFlat[i])
    #        bla = bla/idlFlat[i]
    #        maxim = max(maxim,bla)
    #print('max error')
    #print(maxim)
    return lamGrid, flatOlap