    ##STARTING CALC_SHK
    #
    #
def calc_targOlapf(lamGrid, lam, extrct, flatOlap, label, idl=''):
    ##intializations and hardcoded inputs
    ##TODO take command line input???
    ##mk_flatolap

    #oldLam = lamGrid
    #lamGrid = idl['lamgrid']
    import numpy as np
    import scipy.io as sc
    import sys
    import astropy.io.fits
    import os
    from matplotlib import pyplot as plt
    from hk_windows import mkdir_p 
    from astropy.convolution import convolve, Box1DKernel


    gOrd=[63,64,65,66]
    ngord=len(gOrd)
    #sz=size(lam)#used to find nx TODO remove
    nx=4096#TODO BAD ADAM
    rdnoi=7.*np.sqrt(5.*5.)           # read noise per resolution element in e-
    resel=.0015         #resolution element (nm)

    

    #TODO FIX THIS AND ALL SO THAT IT FOLLOWS PYTHON DATA FORMAT
    #extrct= dataHDu1[1].data
    extrct = np.transpose(extrct)
    lam=np.transpose(lam)


    #print(np.shape(lam))
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

       # print('min')
        #print(min(extrct[:,gord[i]]*scale))
        #print('min')
        #print(max(extrct[:,gord[i]]*scale))
        #flux=interpol(extrct(*,gord(i))*scale,lam(*,gord(i)),lamgrid)
        interpfunc = interpolate.interp1d(lam[:,gOrd[i]],extrct[:,gOrd[i]]*scale, kind='linear', fill_value='extrapolate')
        flux=interpfunc(lamGrid)


        #ADAM CODE LIKLEY BUGGY
        mini = np.min(lam[:,gOrd[i]])
        maxi = np.max(lam[:,gOrd[i]])

        sl=(lamGrid<=mini).nonzero()
        sh=(lamGrid>=maxi).nonzero()


        #print(sl)
        if len(sl)>0 :
            for i in sl:
                flux[i]=0.
        if len(sh)>0 :
            for i in sh:
                flux[i]=0.

        #E ADAM CODE
        targOlap=targOlap+flux




    sg=(flatOlap>=.01).nonzero()
    targOlapf = np.zeros(int(nLamg),dtype=np.float64)
    for i in sg :
        targOlapf[i] = targOlap[i]/flatOlap[i]


    #for i in range(len(gOrd)) :
    #    pyL = extrct[:,gOrd[i]]
    #    idlL = idl['extrct'][gOrd[i],:]
    #    plt.figure()
    #    plt.plot(range(len(pyL)),abs(pyL-idlL))

    #plt.figure()
    #print(flatOlap)

    #plt.plot(lamGrid, targOlap, 'k-')
    #plt.xlabel('wavelength [nm]')
    #plt.ylabel('tragOlap')
    
    mkdir_p("images/" + label.split('/')[0])
    fileStr = "images/" + label + "_targOlapf.pdf"
    plt.plot(lamGrid, targOlapf, 'k-')
    plt.xlabel('wavelength [nm]')
    plt.ylabel('tragOlapf')
    plt.tight_layout()
    plt.savefig(fileStr)
    plt.close()

    return targOlapf
    
    
    
def calc_shk(lamGrid, targOlapf, rvcc, teff=6200., idl=''):
    from hk_windows import hk_windows
    from matplotlib import pyplot as plt
    

    gain=3.4           # e-/ADU
    kk=31.             # factor to make shk into equivalent width (a guess!)
    cahLam=396.85# 396.967 #   #Ca II H line wavelength (nm, vacuum)
    cakLam=393.369#393.485 #    #Ca II K line wavelength (nm, vacuum)
    
    #center of blue continuum band (nm, vacuum)
    lamB=390.2176-.116#subtraction is the offset from our lab values     
    #center of red continuum band (nm, vacuum)
    lamR=400.2204-.116#subtraction is the offset from our lab values  
    
    windows = hk_windows(rvcc, lamGrid,cahLam,cakLam,lamB,lamR)[0]

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
    from hk_windows import PlanckFunc as planck
    plFactor = planck(lamB*10,teff)/planck(lamR*10,teff)
    #print('plfactor: '+ str(plFactor))
        #0.977753 #SO HACKED, TODO DOES THIS CHANGE??
    fb = fr*plFactor

    num = (fh+fk)*gain
    den = (fr+fb)*gain
    shk = kk*(fh+fk)/(fr+fb)
    print("shk: "+ str(shk))
    if idl!='' :
        print("error vs idl; "+str(100*abs(shk-idl['shk'])/idl['shk']))
    return shk