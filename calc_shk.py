# This routine manages the computation of the Ca II H&K activity indes shk.
# Inputs are:
#  mjd = MJD-OBS for the input spectrum.
#  siteid = site from which observation comes.
# rvcc = redshift of target star relative to lab. (km/s)
#  teff = Estimated Teff of target star.
#  extrct(nx,nord) = extracted but not flat-fielded spectra of target. (ADU)
#  lam(nx,nord) = wavelength solution corresp to extrct. (nm)
# Outputs are:
#   shk = (fH + fK)/(fB + fR) where fH, fK are the counts in the cores of
#       the H and K lines, resp., and fB, fR are fluxes in nearby blue and
#       red continuum bands.  Note that shk requires some considerable work
#       to be converted to the more physically useful index R'_HK.
#   eshk = an estimate of the photon+read noise applicable to shk.
# For NRES, fB is not easily accessible, since the relevant echelle order is
# not extracted.  Therefore, fB is approximated as FR*K(Teff), where the 
# function K is the ratio of Planck functions in the red and blue bandpasses.
#END TIM BROWN COMMENTS

#the original calc_shk function has been split into calc_targOlapf and calc_shk to better 
#debug intermediate values 
def calc_targOlapf(lamGrid, lam, extrct, flatOlap, label):
    ##intializations and hardcoded inputs
    import numpy as np
    import scipy.io as sc
    import sys
    import astropy.io.fits
    import os
    import helpers as h
    from matplotlib import pyplot as plt
    from helpers import mkdir_p 
    from astropy.convolution import convolve, Box1DKernel
    from matplotlib.backends.backend_pdf import PdfPages
    
    #good orders?
    gOrd=h.gOrd
    ngord=len(gOrd)

    nx=4096#TODO BAD ADAM
    rdnoi=7.*np.sqrt(5.*5.)           # read noise per resolution element in e-
    resel=.0015                       #resolution element (nm)

    

    #TODO FIX THIS AND ALL SO THAT IT FOLLOWS PYTHON DATA FORMAT INSTEAD OF IDL 
    extrct = np.transpose(extrct)
    lam=np.transpose(lam)



    #brown
    #find the flat field that will be applied to the order-overlapped spectrum.
    #flat field found in above code TODO 
    nLamg=len(lamGrid)

    #brown
    #make the order-overlapped stellar spectrum, divide it by the flat
    nLam = len(lamGrid)
    targOlap=np.zeros(int(nLam),dtype=np.float64)


    from scipy import interpolate
    #https://stackoverflow.com/questions/18326714/idl-interpol-equivalent-for-python
    for i in range(len(gOrd)) :
        lamC = lam[:,gOrd[i]]
        dLamdx = np.gradient(lamC)

        scale = dLamdx/dLamdx[int(nx/2.)]

        #original T.Brown Code
        #flux=interpol(extrct(*,gord(i))*scale,lam(*,gord(i)),lamgrid)
        interpfunc = interpolate.interp1d(lam[:,gOrd[i]],extrct[:,gOrd[i]]*scale, kind='linear', fill_value='extrapolate')
        flux=interpfunc(lamGrid)


        
        #sl=where(lamgrid le min(lam(*,gord(i))),nsl)
        #sh=where(lamgrid ge max(lam(*,gord(i))),nsh)
        #if(nsl gt 0) then flux(sl)=0.d0
        #if(nsh gt 0) then flux(sh)=0.d0
        
        #ADAM CODE OF ABOVE - POTENTIALLY BUGGY
        mini = np.min(lam[:,gOrd[i]])
        maxi = np.max(lam[:,gOrd[i]])
        print('minmax')
        print(mini)
        print(maxi)
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




    fig = plt.figure()
    #print(flatOlap)

    plt.plot(lamGrid, targOlapf, 'k-')
    plt.xlabel('wavelength [nm]')
    plt.ylabel('tragOlapf')
    plt.show()
    plt.close()
    
    #fig = plt.figure()
    #print(flatOlap)

    #plt.plot(lamGrid, flatOlap, 'k-')
    #plt.xlabel('wavelength [nm]')
    #plt.ylabel('flat')
    #plt.show()
    #plt.close()
    
    return targOlapf
    
    

#This is the actual function to calculate the SHK
#Takes the standard lamda grid, and correlated/noise free target data(targOlapf)
#also takes teff for use in the planck function and a radial velocity
#which is not currently used. I have left it in case we decide later to use it
#over cross-correlation
def calc_shk(lamGrid, targOlapf, rvcc, teff=6200.):
    
    from helpers import hk_windows
    import helpers as h#for constants
    
    from astropy.modeling.blackbody import blackbody_lambda
    from astropy import units as u
    from matplotlib import pyplot as plt
    
    
    
    gain=3.4           # e-/ADU
    kk=31.             # factor to make shk into equivalent width (a guess!)
    
    
    
    #this function gives us the regions of our arrays which hold the information we need to sum
    windows = hk_windows(rvcc, lamGrid,h.cahLam,h.cakLam,h.lamB,h.lamR)[0]

    fh=(targOlapf*windows[:,0]).sum()
    fk=(targOlapf*windows[:,1]).sum()
    fr=(targOlapf*windows[:,2]).sum()
    
    
    #plt.figure()
    #cur=windows[:,0]
    #plt.plot(lamGrid[cur!=0],targOlapf[cur!=0])#*cur[cur!=0],'k-')
    #plt.show()
    #plt.close()
    #plt.figure()
    #cur=windows[:,1]
    #plt.plot(lamGrid[cur!=0],targOlapf[cur!=0])#*cur[cur!=0],'k-')
    #plt.show()
    #plt.close()
    #plt.figure()
    #cur=windows[:,2]
    #plt.plot(lamGrid[cur!=0],targOlapf[cur!=0])#*cur[cur!=0],'k-')
    #plt.show()
    #plt.close()
    
    #the SHK calculation with pseudo V-Band
    plFactor = blackbody_lambda(h.lamB*u.nm,teff*u.K).value/blackbody_lambda(h.lamR*u.nm,teff*u.K).value

    fb = fr*plFactor

    num = (fh+fk)*gain
    den = (fr+fb)*gain
    shk = kk*(fh+fk)/(fr+fb)
    #print("shk: "+ str(shk))

    return shk, windows, fb/fr

#smarts specific calc_shk which just calls the special windows function
#to find the v-band window and uses the real v-band data
def smart_calc_shk(lamGrid, targOlapf, rvcc, teff=6200.):
    
    from helpers import smart_hk_windows
    import helpers as h#for constants
    import numpy as np
    
    #from astropy.modeling.blackbody import blackbody_lambda
    from astropy import units as u
    from matplotlib import pyplot as plt
    
    
    gain=3.4           # e-/ADU
    kk=8.*2.40#31.             # factor to make shk into equivalent width (a guess!)
    
    #this function gives us the regions of our arrays which hold the information we need to sum
    windows = smart_hk_windows(rvcc, lamGrid,h.cahLam,h.cakLam,h.lamB,h.lamR)[0]

    print('window shape')
    print(windows.shape)
    print('len ' + str(len(lamGrid)))
    print('len2 ' + str(len(targOlapf)))
    fh=(targOlapf*windows[:,0]).sum()
    fk=(targOlapf*windows[:,1]).sum()
    fr=(targOlapf*windows[:,2]).sum()
    fb=(targOlapf*windows[:,3]).sum()

    print('fr ' +str(fr))
    print('fb ' + str(fb))
    #plt.figure()
    #cur=windows[:,0]
    #plt.plot(lamGrid[cur!=0],targOlapf[cur!=0])#*cur[cur!=0],'k-')
    #plt.show()
    #plt.close()
    #plt.figure()
    #cur=windows[:,1]
    #plt.plot(lamGrid[cur!=0],targOlapf[cur!=0])#*cur[cur!=0],'k-')
    #plt.show()
    #plt.close()
    #plt.figure()
    #cur=windows[:,3]
   # plt.plot(lamGrid[cur!=0],targOlapf[cur!=0])#*cur[cur!=0],'k-')
   # plt.show()
    #plt.close()

    num = (fh+fk)*gain
    den = (fr+fb)*gain
    shk = kk*(fh+fk)/(fr+fb)
    #print("shk: "+ str(shk))

    return shk, windows, fb/fr