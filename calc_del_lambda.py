

def import_lab_frame_spectra(fluxdir, minwl=None, maxwl=None, resolution=1, residual=False):
    import numpy as np
    import matplotlib as mpl
    from matplotlib import pyplot as plt
    import os
    import os.path
    import re
    import scipy
    import scipy.signal
    import scipy.ndimage.filters
    
    
    wavelength = []
    res_flux = []
    irradiance = []
    
    for f in os.listdir(fluxdir):
        if not re.match('^lm', f): continue
        fpath = os.path.join(fluxdir, f)
        d = np.genfromtxt(fpath, unpack=True)
        wavelength = np.append(wavelength, d[0]) # in nm
        res_flux = np.append(res_flux, d[1]) # from a normalized spectra
        irradiance = np.append(irradiance, d[2]) # from a regular spectra; mu-W / cm^2 / nm
        
    wavelength = np.array(wavelength)
    res_flux = np.array(res_flux)
    irradiance = np.array(irradiance)
    sort = np.argsort(wavelength)
    wavelength = wavelength[sort]
    res_flux = res_flux[sort]
    irradiance = irradiance[sort]

    angstroms = wavelength * 10 # in Angstroms
    
    
    
    if minwl is None:
        minwl = angstroms[0]
    if maxwl is None:
        maxwl = angstroms[-1]
    sel = (angstroms >= minwl) & (angstroms <= maxwl)
    dw = angstroms[1] - angstroms[0]
    if not residual:
        series = irradiance / 10. # nm^-1 => Ang.^-1
    else:
        series = res_flux

    if resolution > 0:
        series = scipy.ndimage.filters.gaussian_filter(series, resolution/dw)

    #print('degrading source: ' + str(resolution/dw))
    
    return angstroms[sel], series[sel]
    
    
#############################################
def tmp_find_del_lam(labGrid, lab, tarGrid, targ, smooth) :

    import scipy as sc
    #import astropy.io.fits
    import numpy as np
    #import os
    #from calc_shk import calc_shk
    #from calc_shk import calc_targOlapf
    #from mk_flatolap import mk_flatolap
    from matplotlib import pyplot as plt
    from astropy.convolution import convolve, Box1DKernel
    from scipy import interpolate
    
    tmpGridScale = 1
    dLam = tarGrid[1] - tarGrid[0]
    
    gausdTarg = sc.ndimage.filters.gaussian_filter(targ,smooth/dLam)
    #print('degrading target: ' + str(smooth/dLam))

    #get the lab spectrum into angs div by 10 on angstrom grid for our purposes
    interpfunc = interpolate.interp1d(labGrid/10, lab, kind='linear')#,fill_value='extrapolate')
    labInterp=interpfunc(tarGrid)



    #ZERO OUT THE EDGES OF OUR LAB SPECTRA
    #from 0 to first nonzero element of targolap
    labInterp[:targ.nonzero()[0][0]]=0
    #from last nonzero element to end
    labInterp[targ.nonzero()[0][-1]:]=0




    #tmp = np.convolve(targOlapf,pikapika,'same')
    #plt.figure(figsize=(12,6))
    #plt.title('convolution function')
    #plt.plot(range(len(tmp)), tmp, 'k-', color='blue')
    #plt.xlabel('delta lamda?')
    #plt.xlabel('f * g')

    #plt.figure(figsize=(12,6))
    #plt.plot(tarGrid, targ, 'k-', color='blue')
    #plt.plot(lamGrid, tmpTarg1, 'k-', color="green")
    #plt.plot(lamGrid, tmpTarg2, 'k-',color='red')
    #plt.plot(tarGrid, labInterp*tmpGridScale, 'k-')
    #plt.axvline(x=393.4, color='red')
    #plt.axvline(x=396.9, color='red')

    #plt.xlabel('wavelength')
    #plt.ylabel('smoothed spectras')

    labInterp[labInterp==0]=np.nan
    targ[targ==0]=np.nan



    rmsx =  np.nanmean(labInterp)
    rmsy =   np.nanmean(targ)


    labInterp[np.isnan(labInterp)]=0
    targ[np.isnan(targ)]=0





    #THIS IS THE CROSS CORRELATION SECTION
    ##
    ##
    #print('valid return: ')
    #print(np.correlate(convoldTarg,pikapika-rmsx,'valid'))
    out = np.correlate(targ[targ!=0]-rmsy,labInterp[targ!=0]-rmsx,'full')

    mval = np.argmax(out)
    minim =  np.argmin(out)
   # if abs(min(out)) > abs(max(out)):
       # mval = minim
        #print('min actually')

    #print((mval-len(out)/2))
    #find the half 
    offset = (mval-(len(out)-1)/2)*(tarGrid[1]-tarGrid[0])
    print('offset: ' + str(offset))
    ##
    ##
    ##



    #NOW PLOTS
   #plt.figure(figsize=(12,6))
   # plt.plot(range(len(out)), out, 'k-', color='blue')
   # plt.xlabel('offset')
   # plt.ylabel('correlation')
   # #tmpMax = np.argmax(out)
    #print('index of maximum: ' + str(tmpMax) + ' and adjusted delLam: ' + str(tmpMax/len(out)))
    #print(out)
    #plt.figure(figsize=(12,6))
    #plt.plot(tarGrid, gausdTarg, 'k-', color='blue')
    #plt.axvline(x=393.369, color='red')
    #plt.axvline(x=396.85, color='red')

    #plt.plot(tarGrid-offset, gausdTarg, 'k-',color='green')
    #SCALE JUST FOR VIEWING
    #plt.plot(tarGrid, labInterp*tmpGridScale, 'k-')
    #plt.show()
    #plt.close()
    return offset,targ,labInterp

    #shk = calc_shk(lamGrid, targOlapf, 17. )
    
##    
##    
##    
##    
##
#grid is the x val grid to plot against
#base is the base data we are comparing against
#adj is the adjusted data whos differences against base we want
#regions is a list of lamda values to find differences in
#
def lamda_zoom(bGrid, base, oGrid, obs, regions, descript, width=1):
    import numpy as np
    from matplotlib import pyplot as plt
    from hk_windows import mkdir_p
    
    numRegions = len(regions)
    #fig index is the pass to subplot which shows shape of figure, this is for a tall skinny figure
    #ex 8 regions we want the index to start at 811 for first figure
    #812 will be second etc.
    figIndex = 100*numRegions + 11
    plt.figure(1,figsize=(12,10*numRegions))
    
    obsLam = []
    baseLam = []
    yBuffer = 50
    for l in regions:
        plt.subplot(figIndex)
        figIndex += 1
        
        #find the index range for the current wavelength with size of given width
        bLowI=int(np.abs(bGrid-l+width/2).argmin())
        bMaxI=int(np.abs(bGrid-l-width/2).argmin())
        #for both since obs should be offset
        oLowI=int(np.abs(oGrid-l+width/2).argmin())
        oMaxI=int(np.abs(oGrid-l-width/2).argmin())
        
        #find index of minimum value in range of base and obs using the index
        baseMin = np.argmin(base[bLowI:bMaxI])
        obsMin = np.argmin(obs[oLowI:oMaxI])

        baseLam.append((bGrid[bLowI+baseMin]))
        obsLam.append(oGrid[oLowI+obsMin])
        

        #print('baseLam: '+ str(baseLam[-1]))
        #print('obsLam: '+ str(obsLam))
        #print('b - o: ' + str(baseLam[-1]-obsLam[-1]))
        plt.xlim(baseLam[-1]-width/2,baseLam[-1]+width/2)

        
        plt.axvline(x=l,color='r')
        plt.axvline(x=obsLam[-1],color='g', linestyle='dashed')
        plt.axvline(x=baseLam[-1],color='k', linestyle='dashed')
        #now do some dirty scaling(not perfect but probably good enough)
        scale = obs[bLowI:bMaxI]/base[bLowI:bMaxI]
        avgS = np.mean(scale)
        
        yMin = min(min(obs[oLowI:oMaxI]),min(avgS*base[bLowI:bMaxI]))
        yMax = max(max(obs[oLowI:oMaxI]),max(avgS*base[bLowI:bMaxI]))
        plt.ylim(yMin-yBuffer,yMax+yBuffer)
        
        #plt.axvline(x=396.847, color='red')
        plt.plot(oGrid, obs,'g-',bGrid, base*avgS, 'k-')
        plt.xlabel('nm')
        plt.ylabel('scaled spectra')
        
    
    #create pdf too 
    mkdir_p("images/" + descript.split('/')[0])
    fileStr = "images/" + descript + "_zoom.pdf"
    plt.tight_layout()
    plt.savefig(fileStr)
    
    plt.close(1)
    
    plt.figure()
    plt.plot(np.asarray(baseLam),np.asarray(baseLam)-np.asarray(obsLam), 'ko')
    plt.xlabel('lab frame lambda')
    plt.ylabel('lab lam - obs lam')
    fileStr = "images/" + descript + "_error_graph.pdf"
    plt.tight_layout()
    plt.savefig(fileStr)
    #plt.show()
    plt.close()
    
    plt.figure(figsize=(12,6))
    for l in baseLam:
        plt.axvline(l, color='red')
        
    #only plot non-zero values    
    plt.plot(bGrid[base!=0],base[base!=0]*avgS,'k-',oGrid[obs!=0],obs[obs!=0],'g-')
    fileStr = "images/" + descript + "_wide_view.pdf"
    plt.tight_layout()
    plt.savefig(fileStr)
    plt.close()