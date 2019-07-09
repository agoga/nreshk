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
def calc_targOlapf(lamGrid, lam, extrct, flatOlap, label, idl=''):
    ##intializations and hardcoded inputs
    ##TODO take command line input???
    
    #TODO do we need all this?
    import numpy as np
    import scipy.io as sc
    import sys
    import astropy.io.fits
    import os
    from matplotlib import pyplot as plt
    from helpers import mkdir_p 
    from astropy.convolution import convolve, Box1DKernel
    from matplotlib.backends.backend_pdf import PdfPages
    

    gOrd=[63,64,65,66]
    ngord=len(gOrd)
    #sz=size(lam)#used to find nx TODO remove
    nx=4096#TODO BAD ADAM
    rdnoi=7.*np.sqrt(5.*5.)           # read noise per resolution element in e-
    resel=.0015                       #resolution element (nm)

    

    #TODO FIX THIS AND ALL SO THAT IT FOLLOWS PYTHON DATA FORMAT
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

        #original
        #flux=interpol(extrct(*,gord(i))*scale,lam(*,gord(i)),lamgrid)
        interpfunc = interpolate.interp1d(lam[:,gOrd[i]],extrct[:,gOrd[i]]*scale, kind='linear', fill_value='extrapolate')
        flux=interpfunc(lamGrid)


        
        #sl=where(lamgrid le min(lam(*,gord(i))),nsl)
        #sh=where(lamgrid ge max(lam(*,gord(i))),nsh)
        #if(nsl gt 0) then flux(sl)=0.d0
        #if(nsh gt 0) then flux(sh)=0.d0
        
        #ADAM CODE OF ABOVE - LIKLEY BUGGY
        mini = np.min(lam[:,gOrd[i]])
        maxi = np.max(lam[:,gOrd[i]])

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


    #for i in range(len(gOrd)) :
    #    pyL = extrct[:,gOrd[i]]
    #    idlL = idl['extrct'][gOrd[i],:]
    #    plt.figure()
    #    plt.plot(range(len(pyL)),abs(pyL-idlL))

    #fig = plt.figure()
    #print(flatOlap)

    #plt.plot(lamGrid, targOlapf, 'k-')
    #plt.xlabel('wavelength [nm]')
    #plt.ylabel('tragOlapf')
    #curPdf.savefig(fig)
   # plt.close()
    
    return targOlapf
    
    

    
def calc_shk(lamGrid, targOlapf, rvcc, teff=6200., idl=''):
    
    from helpers import hk_windows
    from helpers import PlanckFunc as planck
    from matplotlib import pyplot as plt
    import helpers as h#for constants
    
    
    gain=3.4           # e-/ADU
    kk=31.             # factor to make shk into equivalent width (a guess!)
    
    
    
    
    windows = hk_windows(rvcc, lamGrid,h.cahLam,h.cakLam,h.lamB,h.lamR)[0]

    fh=(targOlapf*windows[:,0]).sum()
    fk=(targOlapf*windows[:,1]).sum()
    fr=(targOlapf*windows[:,2]).sum()
    
    mini = next((i for i, x in enumerate(windows[:,0]) if x), None)
    maxi = [i for i, e in enumerate(windows[:,0]) if e != 0][-1]
    #print('lamMin H: '+ str(lamGrid[mini]))
    #print('lamMax H: '+ str(lamGrid[maxi]))
    
    mini = next((i for i, x in enumerate(windows[:,1]) if x), None)
    maxi = [i for i, e in enumerate(windows[:,1]) if e != 0][-1]
    #print('lamMin K: '+ str(lamGrid[mini]))
    #print('lamMax K: '+ str(lamGrid[maxi]))
    
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
    ##plt.plot(lamGrid[cur!=0],targOlapf[cur!=0])#*cur[cur!=0],'k-')
    #plt.show()
    #plt.close()
    
    #print('fh')
    #print(abs(fh-idl['fh']))
    #print('fk')
    #print(abs(fk-idl['fk']))
    #print('fr')
    #print(abs(fr-idl['fr']))
    #print('pl')
    #print(abs(plFactor-idl['plfactor']))
    #TODO FIND WORKING PYTHON PLANK FUNCTION
    
    plFactor = planck(h.lamB*10,teff)/planck(h.lamR*10,teff)
    #print('plfactor: '+ str(plFactor))
        #0.977753 #SO HACKED, TODO DOES THIS CHANGE??
    fb = fr*plFactor

    num = (fh+fk)*gain
    den = (fr+fb)*gain
    shk = kk*(fh+fk)/(fr+fb)
    #print("shk: "+ str(shk))
    if idl!='' :
        print("error vs idl; "+str(100*abs(shk-idl['shk'])/idl['shk']))
    return shk, windows

#smarts specific
def smart_calc_shk(lamGrid, targOlapf, rvcc, teff=6200., idl=''):
    
    from helpers import smart_hk_windows
    from helpers import PlanckFunc as planck
    from matplotlib import pyplot as plt
    import helpers as h#for constants
    
    
    gain=3.4           # e-/ADU
    kk=8.*2.40#31.             # factor to make shk into equivalent width (a guess!)
    
    
    
    
    windows = smart_hk_windows(rvcc, lamGrid,h.cahLam,h.cakLam,h.lamB,h.lamR)[0]

    fh=(targOlapf*windows[:,0]).sum()
    fk=(targOlapf*windows[:,1]).sum()
    fr=(targOlapf*windows[:,2]).sum()
    fb=(targOlapf*windows[:,3]).sum()
    
    #mini = next((i for i, x in enumerate(windows[:,0]) if x), None)
    #maxi = [i for i, e in enumerate(windows[:,0]) if e != 0][-1]
    #print('lamMin H: '+ str(lamGrid[mini]))
    #print('lamMax H: '+ str(lamGrid[maxi]))
    
    #mini = next((i for i, x in enumerate(windows[:,1]) if x), None)
    #maxi = [i for i, e in enumerate(windows[:,1]) if e != 0][-1]
    #print('lamMin K: '+ str(lamGrid[mini]))
    #print('lamMax K: '+ str(lamGrid[maxi]))
    
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
    ##plt.plot(lamGrid[cur!=0],targOlapf[cur!=0])#*cur[cur!=0],'k-')
    #plt.show()
    #plt.close()
    
    #print('fh')
    #print(abs(fh-idl['fh']))
    #print('fk')
    #print(abs(fk-idl['fk']))
    #print('fr')
    #print(abs(fr-idl['fr']))
    #print('pl')
    #print(abs(plFactor-idl['plfactor']))
    #TODO FIND WORKING PYTHON PLANK FUNCTION
    
    plFactor = planck(h.lamB*10,teff)/planck(h.lamR*10,teff)
    #print('plfactor: '+ str(plFactor))
        #0.977753 #SO HACKED, TODO DOES THIS CHANGE??
    #fb = fr*plFactor

    num = (fh+fk)*gain
    den = (fr+fb)*gain
    shk = kk*(fh+fk)/(fr+fb)
    #print("shk: "+ str(shk))
    if idl!='' :
        print("error vs idl; "+str(100*abs(shk-idl['shk'])/idl['shk']))
    return shk, windows, fr/fb