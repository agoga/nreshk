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
    
sigToFWHM = 2.355#used to display a real FWHM value on plots where smoothing occurs

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

#smarts specific hk windows with V-Band included
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
    #v-band inclusion
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
            

##
##This is the massive printing function to create a report of each observation pushed through the
#pipeline.
def pdf_from_data(bGrid, base, oGrid, obs, windows, title, path, descript, flat='', width=1):
    import numpy as np
    from matplotlib import pyplot as plt
    from helpers import mkdir_p
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.gridspec as gridspec
    import scipy as sc
    import helpers as h#for constants
    
    
    calH = h.cahLam#396.847#TODO and make global
    calK = h.cakLam#393.366
    #center of red continuum band (nm, vacuum)
    
    lamR=h.lamR
    
    #if flat was passed as empty, fill it with nans
    if len(flat) == 0:
        flat = np.full(len(bGrid),0)
        
    flatMax = max(flat)    
    
    if flatMax != 0:
        #create a bool array which describes in our grid sections of the flat which higher than .4th of flat max
        #print(flatMax)
        flatSection = flat/flatMax >= .4

        #used to set axis bounds
        mini = min(bGrid[flatSection])
        maxi = max(bGrid[flatSection])
    
    smooth = .01
    
    #colors for hk lines in each plot
    hColor = 'dodgerblue'
    kColor = 'turquoise'
    
    #create a pdf page with 4 rows and 3 columns
    #the bottom 3 rows use all 3 columns but the HK windows output is split inti
    #K feature - H feature - R-Band
    #will need to make 4 wide when V-band is included
    fig, ax = plt.subplots(figsize=(10,10))
    plt.suptitle(title)
    plt.ticklabel_format(useOffset=False)
    with PdfPages(path+descript+"_report.pdf") as curPdf:
        gs = gridspec.GridSpec(4, 3)
        
        #references to each plot
        targPlt = plt.subplot(gs[2,:])
        smoothedPlt = plt.subplot(gs[1,:])
        flatPlt = plt.subplot(gs[3,:])
        kPlt = plt.subplot(gs[0,0])
        hPlt = plt.subplot(gs[0,1])
        rPlt =plt.subplot(gs[0,2])
        
        #please matplotlib don't make my stuff hard to read!
        hPlt.ticklabel_format(useOffset=False)
        kPlt.ticklabel_format(useOffset=False)
        rPlt.ticklabel_format(useOffset=False)
        targPlt.ticklabel_format(useOffset=False)
        flatPlt.ticklabel_format(useOffset=False)
        smoothedPlt.ticklabel_format(useOffset=False)
        
        hPlt.tick_params(axis='both',which= 'major', labelsize=7)
        kPlt.tick_params(axis='both',which= 'major', labelsize=7)
        rPlt.tick_params(axis='both',which= 'major', labelsize=7)
        targPlt.tick_params(axis='both',which= 'major', labelsize=7)
        flatPlt.tick_params(axis='both',which= 'major', labelsize=7)
        smoothedPlt.tick_params(axis='both',which= 'major', labelsize=7)
        
        #titles and lebels
        kPlt.set_xlabel("Wavelength(nm)")
        hPlt.set_ylabel("Irradiance")
        
        hPlt.set_title("Ca-H window")
        kPlt.set_title("Ca-K window")
        rPlt.set_title("Red band window")
        
        targPlt.set_title("Target overlap shifted over reference spectra")
        targPlt.set_xlabel("Wavelength(nm)")
        targPlt.set_ylabel("Irradiance scaled")
        
        smoothedPlt.set_title("Reference and target smoothed by " + str(smooth*h.sigToFWHM) + " nm Kernel")
        smoothedPlt.set_xlabel("Wavelength(nm)")
        smoothedPlt.set_ylabel("Irradiance scaled")
        
        flatPlt.set_title("Flat plot for target")
        flatPlt.set_xlabel("Wavelength(nm)")
        flatPlt.set_ylabel("Irradiance")
        
        #H/K/Red band plots zoomed
        #use exactly the windows that are gotten from hk_windows plus small buffer
        cur=windows[:,0]
        hkWidth=h.lineWid + .005
        rWidth =h.conWid/2 +.05
        
        hPlt.axvline(x=calH,color=hColor)
        hPlt.plot(oGrid[cur!=0],obs[cur!=0],'b-')
        hPlt.set_xlim(calH-hkWidth,calH+hkWidth)
        
        cur=windows[:,1]
        kPlt.axvline(x=calK,color=kColor)
        kPlt.plot(oGrid[cur!=0],obs[cur!=0],'b-')
        kPlt.set_xlim(calK-hkWidth,calK+hkWidth)

        cur=windows[:,2]
        rPlt.plot(oGrid[cur!=0],obs[cur!=0])
        rPlt.set_xlim(lamR-rWidth, lamR+rWidth)
        
        
        #Targolapf plot 
        #get red lines from windows funct too
        rMin = oGrid[cur!=0][0]
        rMax = oGrid[cur!=0][-1]
       #print('rmin: ' + str(rMin) + ' rMax: ' + str(rMax))
        #terrrrible way to get scale fixxxxxxx
        scale = obs[cur!=0]/base[cur!=0]
        avgS = np.mean(scale)
        
        targPlt.axvline(x=calH,color=hColor)
        targPlt.axvline(x=calK,color=kColor)
        targPlt.axvline(x=rMin, color='red')
        targPlt.axvline(x=rMax, color='red')
        
        #if there is a flat to plot then only use the target data when it is greater than .4
        #If we dont do this then on the edges out graph will be divided by a very small number which
        #makes it hard to see the data in the middle that we care about
        if flatMax != 0:
            targPlt.set_xlim([mini,maxi])
            targPlt.plot(bGrid[flatSection],base[flatSection]*avgS,color='lightgray')
            targPlt.plot(oGrid[flatSection],obs[flatSection],'b-')
        else:
            targPlt.set_xlim([391.5,407])
            targPlt.plot(bGrid[base!=0],base[base!=0]*avgS,color='lightgray')
            targPlt.plot(oGrid[obs!=0],obs[obs!=0],'b-')
        
        #smoothed target and lab plot
        dOLam = oGrid[1]-oGrid[0]
        dBLam = bGrid[1]-bGrid[0]
        gdObs = sc.ndimage.filters.gaussian_filter(obs,smooth/dOLam)
        gdBase =  sc.ndimage.filters.gaussian_filter(base,smooth/dBLam)
        
        scale = obs[cur!=0]/base[cur!=0]
        avgS = np.mean(scale)
        
        
        
        smoothedPlt.axvline(x=calH,color=hColor)
        smoothedPlt.axvline(x=calK,color=kColor)
        smoothedPlt.axvline(x=rMin, color='red')
        smoothedPlt.axvline(x=rMax, color='red')
        
        
        #if there is a flat to plot then only use the target data when it is greater than .4
        #If we dont do this then on the edges out graph will be divided by a very small number which
        #makes it hard to see the data in the middle that we care about
        if flatMax != 0:
            smoothedPlt.set_xlim([mini,maxi])
            smoothedPlt.plot(bGrid[flatSection],gdBase[flatSection]*avgS,color='lightgray')
            smoothedPlt.plot(oGrid[flatSection],gdObs[flatSection],'b-')
        else:
            smoothedPlt.set_xlim([391.5,407])
            smoothedPlt.plot(bGrid[gdBase!=0],gdBase[gdBase!=0]*avgS,color='lightgray')
            smoothedPlt.plot(oGrid[gdObs!=0],gdObs[gdObs!=0],'b-')
        
        #flat/other plot
        flatPlt.plot(bGrid[flat!=0], flat[flat!=0], 'k-')
        
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        #plt.show()
        plt.close()
        
        curPdf.savefig(fig)
        
        
#the idea behind this function is to auto-correlate a spectra.
#if the spectra is not complete noise then the correlation function will have a small peek at 0
#if the spectra is noise then the peek should be much more massive since no where besides at no
#offset will the noise overlap.
def bad_spec_detection_v2(lamGrid, targ):
    import scipy as sc
    
    import numpy as np

    from matplotlib import pyplot as plt
    
    correlation = np.correlate(targ,targ,'full')
    
    #TODO maybe take the maximum out of the mean for a different ratio.
    mean = np.mean(correlation[39900:40100])
    maxi = max(correlation[39900:40100])
    print("bad correlation: " + str(maxi/mean)) 
    
    #plt.figure()
    #plt.plot(range(len(targ)),targ)
    #plt.show()
    #plt.close()
    #plt.figure()
    #plt.plot(range(len(correlation)),correlation)
    #plt.xlim([39900,40100])
    #plt.show()
    #plt.close()
    
    
    #TODO for future projects, fiddling with this number will remove or allow different spectra
    #1.7 might be a better go
    if maxi/mean > 2:
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

    #todo clean and place into healpers file
def print_header(header):
    for h in header.keys():
        print(h, header[h])
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
    import os
    return [name for name in os.listdir(a_dir)
            if os.path.isdir(os.path.join(a_dir, name))]