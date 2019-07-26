c=2.99792458e5  #speed of light (km/s)

cahLam=396.847# 396.967 #   #Ca II H line wavelength (nm, vacuum)
cakLam=393.366#393.485 #    #Ca II K line wavelength (nm, vacuum)


#center of blue continuum band (nm, vacuum)
lamB=390.2176-.116#subtraction is the offset from our lab values     
#center of red continuum band (nm, vacuum)
#TODO TAKE FROM VAUGHAN 1978 subtraction is the offset from our lab values 
lamR=400.2204-.116#subtraction is the offset from our lab values  

conWid=2.0         #wavelength width of continuum bands (nm, vacuum)
lineWid=.109       #FWHM of line core window functions (nm, vacuum)
    
sigToFWHM = 2.355

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


def hk_windows(rvcc,lamGrid,cahLam,cakLam,lamB,lamR):
    import numpy as np
    
    #brown
    #make output array
    nLam = len(lamGrid)
    windows = np.zeros((nLam,3),dtype=np.float32)
    z = 1. +rvcc/c
    
    #brown
    #make window functions
    d0 = abs(lamGrid-cahLam*z)/lineWid
    s = (d0<=1.0).nonzero()
    if len(s) > 0:
        windows[s,0]=1.-d0[s] 
    
    d1 = abs(lamGrid-cakLam*z)/lineWid
    s = (d1<=1.0).nonzero()
    if len(s) > 0:
        windows[s,1]=1.-d1[s] 
    
    d2 = abs(lamGrid-lamR*z)*2./conWid
    s = (d2<=1.0).nonzero()
    if len(s) > 0:
        windows[s,2]=1.
    return windows, lamB, lamR

#smarts specific
def smart_hk_windows(rvcc,lamGrid,cahLam,cakLam,lamB,lamR):
    import numpy as np
    
    #brown
    #make output array
    nLam = len(lamGrid)
    windows = np.zeros((nLam,4),dtype=np.float32)
    z = 1. +rvcc/c
    
    #brown
    #make window functions
    d0 = abs(lamGrid-cahLam*z)/lineWid
    s = (d0<=1.0).nonzero()
    if len(s) > 0:
        #print('g')
        #print(1.-d0[s])
        windows[s,0]=1.-d0[s] 
    
    d1 = abs(lamGrid-cakLam*z)/lineWid
    s = (d1<=1.0).nonzero()
    if len(s) > 0:
        windows[s,1]=1.-d1[s] 
    
    d2 = abs(lamGrid-lamR*z)*2./conWid
    s = (d2<=1.0).nonzero()
    if len(s) > 0:
        windows[s,2]=1.
    
    d3 = abs(lamGrid-lamB*z)*2./conWid
    s = (d3<=1.0).nonzero()
    if len(s) > 0:
        windows[s,3]=1.
    return windows, lamB, lamR

#https://stackoverflow.com/questions/11373610/save-matplotlib-file-to-a-directory
def mkdir_p(mypath):
    '''Creates a directory. equivalent to using mkdir -p on the command line'''

    from errno import EEXIST
    from os import makedirs,path

    try:
        makedirs(mypath)
    except OSError as exc: # Python >2.5
        if exc.errno == EEXIST and path.isdir(mypath):
            pass
        else: raise
            

def PlanckFunc (wl , T ) :
    import numpy as np
    import scipy.constants as sc
 #Evaluate the emission intensity for a blackbody of temperature T
 #as a function of wavelength

 #Inputs :
 #wl :: numpy array containing wavelengths [m]
 #T :: temperature [K]
 #Outputs :
 #B :: intensity [W sr ** -1 m** -2]

    wl = np.array(wl) # if the input is a list or a tuple, make itan array
    h = sc.Planck # Planck constant [J s]

    c = sc.c# speed of light [m s**−1]
    kb = sc.k # Boltzmann constant [J K**−1]
    B = ((2* h * c * c ) /( wl **5) ) /( np.exp (( h * c ) /( wl * kb * T ) ) -1)
    return B

def bad_spec_detection_v2(lamGrid, targ):
    import scipy as sc
    
    import numpy as np

    from matplotlib import pyplot as plt
    
    correlation = np.correlate(targ,targ,'full')
    mean = np.mean(correlation[39900:40100])
    maxi = max(correlation[39900:40100])
    print("correlation: " + str(maxi/mean)) 
    
    #plt.figure()
    #plt.plot(range(len(targ)),targ)
    #plt.show()
    #plt.close()
    #plt.figure()
    #plt.plot(range(len(correlation)),correlation)
    #plt.xlim([39900,40100])
    #plt.show()
    #plt.close()
    
    #1.7 might be a better go
    if maxi/mean > 2:
        return True
    else:
        return False
    
def bad_spec_detection(lamGrid, targ):
    
    import scipy as sc
    
    import numpy as np

    from matplotlib import pyplot as plt
    
    #1.5 angstroms allowed diff from the desired spectral features
    allowedError = .24
    
    lowI = np.abs(lamGrid - 392).argmin()
    highI =np.abs(lamGrid - 405).argmin()
    bad = True
    width = 1000
    
    dLam = lamGrid[1] - lamGrid[0]
    
    
    #print('low: ' + str(low) + ' high: ' + str(high))
    gausdTarg = sc.ndimage.filters.gaussian_filter(targ,.05/dLam)
    newLook = gausdTarg[lowI:highI]
    adjLam = lamGrid[lowI:highI]
    #fig = plt.figure()
    #plt.plot(adjLam, newLook, 'k-')
    #plt.show()
    #plt.close()
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