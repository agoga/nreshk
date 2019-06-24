def hk_windows(rvcc,lamGrid,cahLam,cakLam,lamB,lamR):
    import numpy as np
    
    
    c=2.99792458e5  #speed of light (km/s)
        
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