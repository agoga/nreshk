def hk_windows(rvcc,lamGrid):
    import numpy as np
    
    
    c=2.99792458e5  #speed of light (km/s)
    cahLam=396.967     #Ca II H line wavelength (nm, vacuum)
    cakLam=393.485     #Ca II K line wavelength (nm, vacuum)
    lamB=390.2176      #center of blue continuum band (nm, vacuum)
    lamR=400.2204      #center of red continuum band (nm, vacuum)
    conWid=2.0         #wavelength width of continuum bands (nm, vacuum)
    lineWid=.109       #FWHM of line core window functions (nm, vacuum)
    
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
    return windows, lamB, lamR