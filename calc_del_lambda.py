

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
    import helpers as h#for constants
    
    
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
#TODO DESCPT
def tmp_find_del_lam(labGrid, lab, tarGrid, targ, smooth) :

    #TODO DO WE NEED ALL THESE
    import scipy as sc
    #import astropy.io.fits
    import numpy as np
    import helpers as h#for constants
    #import os
    #from calc_shk import calc_shk
    #from calc_shk import calc_targOlapf
    #from mk_flatolap import mk_flatolap
    from matplotlib import pyplot as plt
    from astropy.convolution import convolve, Box1DKernel
    from scipy import interpolate
    
    tmpGridScale = 1
    dLam = tarGrid[1] - tarGrid[0]
    #print(np.shape(targ))
    gausdTarg = sc.ndimage.filters.gaussian_filter(targ,smooth/dLam)
    
    dLabLam = labGrid[1] - labGrid[0]
    
    #TODO SHOULD 2.55 be h.sigToFWHM?
    gausedLab = sc.ndimage.filters.gaussian_filter(lab,(dLam/dLabLam)/2.55)

    #get the lab spectrum into angs div by 10 on angstrom grid for our purposes
    interpfunc = interpolate.interp1d(labGrid, gausedLab, kind='linear')#,fill_value='extrapolate')
    labInterp=interpfunc(tarGrid)



    #ZERO OUT THE EDGES OF OUR LAB SPECTRA
    #TODO do this with strict values to remove edge errors affecting correlation
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

    #do not consider 0's while taking mean
    labInterp[labInterp==0]=np.nan
    targ[targ==0]=np.nan
    
    rmsx = np.nanmean(labInterp)
    rmsy =  np.nanmean(targ)

    #place the 0's back
    labInterp[np.isnan(labInterp)]=0
    targ[np.isnan(targ)]=0




    #plt.figure(figsize=(12,6))
    #plt.plot(1000*labInterp[targ!=0]-rmsx,'k-')
    #plt.plot(targ[targ!=0]-rmsy,'g-')
    #plt.show()
    #plt.close()
    
    #THIS IS THE CROSS CORRELATION SECTION
    ##
    ##correlate must have same sized arrays input
    correlation = np.correlate(targ[targ!=0]-rmsy,(labInterp[targ!=0]-rmsx),'full')
    #length of the correlation array is length input array times 2 plus 1
    #if the two arrays are already aligned then the peak of correlation function should be middle
    middle = int((len(correlation)-1)/2)
    
    #width is used for local maximum finding
    width = 50#TODO needs to be a grid based setting
    #needs to be like 50 for smarts and like 250 for nres
    

    
    #max value is the index of the maximum value(local max around middle of array if width used)
    mval = middle-width+np.argmax(correlation[middle-width:middle+width])
    
    
    ##PRE VACAY
#want to make a quadratic to be more precise with 'peak' of correlation
#need to do poly only in certain range around center because wings will take over the fit
    fitWidth = 5
    import numpy.polynomial.polynomial as poly
    xRange = mval + np.arange(2*fitWidth)- fitWidth
    polyFunc = np.polyfit(xRange, correlation[mval-fitWidth:mval+fitWidth],2)
    #print(polyFunc)
    #f=a*x^2 +b*x+c
    #f' = 2*a*x+b = 0 -> x = -b/2a 
    xVal = -polyFunc[1]/(2*polyFunc[0])
    
    #print(xVal)
    #ffit = np.polyval(polyFunc,xRange)
    #plt.figure()
    #plt.xlim(mval-fitWidth*3,mval+fitWidth*3)
    #plt.plot(xRange, ffit,'g-')
    #plt.plot(range(len(correlation)), correlation, 'k-')
    #plt.show()
    #plt.close()
    
    #mval= np.argmax(ffit)
    #the actual lamda offset is how far from middle we are in pixel space times the pixel to grid ratio
    #
    offset = (xVal-middle)*(tarGrid[1]-tarGrid[0])
    #print('offset: ' + str(offset))
    ##


    
    
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
    return offset,targ,labInterp,gausdTarg
   
    
    
##
##TODO descripts
def pdf_from_data(bGrid, base, oGrid, obs, flat, windows, shk, path, descript, width=1):
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
    
    lamR=h.lamR#400.2204-.116#TODO TAKE FROM VAUGHAN 1978 subtraction is the offset from our lab values 
    
    
    smooth = .01
    
    #colors for hk lines in each plot
    hColor = 'dodgerblue'
    kColor = 'turquoise'
    
    
    fig, ax = plt.subplots(figsize=(10,10))
    plt.ticklabel_format(useOffset=False)
    with PdfPages(path+descript+"_report.pdf") as curPdf:
        gs = gridspec.GridSpec(4, 3)
        
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
        
        #titles 
        kPlt.set_xlabel("Wavelength(nm)")
        hPlt.set_ylabel("Irradiance")
        
        hPlt.set_title("Cal-H window")
        kPlt.set_title("Cal-K window")
        rPlt.set_title("Red band window")
        
        targPlt.set_title("Target overlap shifted over lab spectra")
        targPlt.set_xlabel("Wavelength(nm)")
        targPlt.set_ylabel("Irradiance scaled")
        
        smoothedPlt.set_title("Lab and target smoothed by " + str(smooth*h.sigToFWHM) + " nm Kernal")
        smoothedPlt.set_xlabel("Wavelength(nm)")
        smoothedPlt.set_ylabel("Irradiance scaled")
        
        flatPlt.set_title("Flat plot for target with shk: " + str(shk))
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
        targPlt.set_xlim([391.5,407])
        targPlt.axvline(x=calH,color=hColor)
        targPlt.axvline(x=calK,color=kColor)
        targPlt.axvline(x=rMin, color='red')
        targPlt.axvline(x=rMax, color='red')
        targPlt.plot(bGrid[base!=0],base[base!=0]*avgS,color='lightgray')
        targPlt.plot(oGrid[obs!=0],obs[obs!=0],'b-')
        
        #smoothed target and lab plot
        dOLam = oGrid[1]-oGrid[0]
        dBLam = bGrid[1]-bGrid[0]
        gdObs = sc.ndimage.filters.gaussian_filter(obs,smooth/dOLam)
        gdBase =  sc.ndimage.filters.gaussian_filter(base,smooth/dBLam)
        
        scale = obs[cur!=0]/base[cur!=0]
        avgS = np.mean(scale)
        smoothedPlt.set_xlim([391.5,407])
        smoothedPlt.axvline(x=calH,color=hColor)
        smoothedPlt.axvline(x=calK,color=kColor)
        smoothedPlt.axvline(x=rMin, color='red')
        smoothedPlt.axvline(x=rMax, color='red')
        smoothedPlt.plot(bGrid[gdBase!=0],gdBase[gdBase!=0]*avgS,color='lightgray')
        smoothedPlt.plot(oGrid[gdObs!=0],gdObs[gdObs!=0],'b-')
        
        #flat/other plot
        flatPlt.plot(bGrid[flat!=0], flat[flat!=0], 'k-')
        
        plt.tight_layout()
        plt.close()
        
        curPdf.savefig(fig)
    
##    
##    
##    
##    DEPRECATED FUNCTION REMOVE SOON
##
#grid is the x val grid to plot against
#base is the base data we are comparing against
#adj is the adjusted data whos differences against base we want
#regions is a list of lamda values to find differences in
#
def lamda_zoom(bGrid, base, oGrid, obs, regions, path, descript, width=1):
    import numpy as np
    from matplotlib import pyplot as plt
    from hk_windows import mkdir_p
    from matplotlib.backends.backend_pdf import PdfPages

    numRegions = len(regions)
    #fig index is the pass to subplot which shows shape of figure, this is for a tall skinny figure
    #ex 8 regions we want the index to start at 811 for first figure
    #812 will be second etc.
    figIndex = 100*numRegions + 11
    fig = plt.figure(1,figsize=(12,10*numRegions))
    
    obsLam = []
    baseLam = []
    yBuffer = 50
    #with PdfPages('multipage_pdf.pdf') as pdf:
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
    #mkdir_p("images/" + path+ descript.split('/')[0])
    #fileStr = "images/" + path+ descript + "_zoom.pdf"
    plt.tight_layout()
    curPdf.savefig(fig)
    
    plt.close(1)
    
    fig = plt.figure()
    plt.plot(np.asarray(baseLam),np.asarray(baseLam)-np.asarray(obsLam), 'ko')
    plt.xlabel('lab frame lambda')
    plt.ylabel('lab lam - obs lam')
    #fileStr = "images/" +path+ descript + "_error_graph.pdf"
    plt.tight_layout()
    #curPdf.savefig(fig)
    #plt.savefig(fileStr)
    #plt.show()
    plt.close()
    
    fig = plt.figure(figsize=(12,6))
    for l in baseLam:
        plt.axvline(l, color='red')
        
    #only plot non-zero values    
    plt.plot(bGrid[base!=0],base[base!=0]*avgS,'k-',oGrid[obs!=0],obs[obs!=0],'g-')
    fileStr = "images/" + path+descript + "_wide_view.pdf"
    plt.tight_layout()
    #plt.savefig(fileStr)
    #curPdf.savefig(fig)
    plt.close()
    
    fig = plt.figure()
    fileStr = "images/" + path+descript + "_targOlapf.pdf"
    plt.plot(oGrid, obs, 'k-')
    plt.xlabel('wavelength [nm]')
    plt.ylabel('tragOlapf')
    plt.tight_layout()
    #plt.savefig(fileStr)
    #curPdf.savefig(fig)
    plt.close()