# plotting.py 
# Adam Goga
# 
#def plot_timeseries - plots final timeseries using the analyzed_data structure
#def pdf_from_intermediate_data - Produces pdf of analyzable data; zoomed out window 
# fits, raw spectra and flat along with each value found through the pipeline 
import astropy.io.fits 

import numpy as np
import scipy as sc
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.time import Time

import helpers as h#for constants


#this function plots the final time series of SHK values for each star
def plot_timeseries(inData,bad, fig = None, ax = None, multi=None, outputPath=None):
    mjdArray=[o.mjd for o in inData]
    shkValArray=[o.shk for o in inData]

    if fig is None and ax is None:
        fig, ax = plt.subplots(figsize=(12,6))
       
    

    starName = ''
    outputDir = ''
    if inData is not None:
        starName = inData[0].star
        outputDir = inData[0].starDir()
    pltStr = []
    t = []


    for i in range(len(inData)):
        d = inData[i]

        avg=d.average
        s=d.site

        size = h.singleIconSize
        opac = h.singleOpacity
        
        if s in h.siteColors:
            col = h.siteColors[s][0]
        else:
            col = 'y'#some ugly default if we failed to update the site color dict
        if d.mjd in bad and avg == True:
            mark ='x'
        elif avg == True:
            mark = 'o'#circle
        else:
            mark = '^'#triangle

        if mark == 'o':
            size = h.averageIconSize
            opac = h.averageOpacity
        ax.scatter(d.decimalYr.value,d.shk,marker=mark,edgecolors='k', c=col, s=size, alpha=opac)

        if hasattr(d,'bad') and d.bad == True:
            ax.scatter(d.decimalYr.value,d.shk,marker='x',edgecolors='k', c='r', s=size*2/3, alpha=opac)

        #if hasattr(d,'nOrd') and d.nOrd == 67:
        #    ax.scatter(d.decimalYr.value,d.shk,marker='x',edgecolors='k', c='k', s=size*1/3, alpha=opac)
    
    sites =[]
    for key in h.siteColors:
        siteCol = h.siteColors[key][0]
        siteName = h.siteColors[key][1]
        p = mpatches.Patch(color=siteCol, label=siteName)
        sites.append(p)

    plt.legend(handles=sites,prop={'size': 10}, loc="upper left",bbox_to_anchor=(1.04,1))
    

    h.mkdir_p(outputDir)

    plt.style.use('classic')
    ax.ticklabel_format(useOffset=False)
    plt.title('HD '+starName+' magnetic activity time series')
    plt.xlabel('Time(years)')
    plt.ylabel('Adjusted S-index')

    if outputPath is None:
        plt.savefig(outputDir+starName+'_shk_time_series.pdf',bbox_inches='tight')
    else:
        plt.savefig(outputPath+'_shk_time_series.pdf',bbox_inches='tight')
        
    if multi is None:
        plt.show()
        plt.close()



##
##This is the massive printing function to create a report of each observation pushed through the
#pipeline.
#def pdf_from_intermediate_data(bGrid, base, oGrid, obs, windows, title, path, flat='', width=1):
def pdf_from_intermediate_data(bGrid, base, oData, width=1):
    title = oData.pdfTitle()
    windows = oData.window
    oGrid = oData.lamGrid#- oData.offset[0] #Shift the grid by the offset if any
    flat = oData.flat
    obs = oData.targOlapf
    calH = h.cahLam#396.847
    calK = h.cakLam#393.366

    hOffset = oData.offset[0]
    kOffset = oData.offset[1]
    rOffset = oData.offset[2]
    bOffset = oData.offset[3]
    #center of red continuum band (nm, vacuum)

    #create the directories for pdf plotting and save every intermediate data array

    outputDir = oData.outputDir()
    dataPath = oData.data_path()
    reportPath = oData.report_path()

    #probably want to move the data saving to a different function
    h.mkdir_p(outputDir)

    #probably want to move the data saving to a different function
    #IF DEBUG TODO
    np.savez(dataPath, analyzedData=oData)
    #3rd time
    #probably want to move the data saving to a different function

    lamR=h.lamR
    lamB=h.lamB

    #if flat was passed as empty, fill it with nans
    if len(flat) == 0:
        flat = np.full(len(bGrid),0)
        
    flatMax = max(flat)    
    
    if flatMax != 0:
        #create a bool array which describes in our grid sections of the flat which higher than .4th of flat max.
        #edit - basically want all of the stuff TODO check this, was changed from .4 may need to go back
        flatSection = flat/flatMax >= .01

        #used to set axis bounds
        mini = min(bGrid[flatSection])
        maxi = max(bGrid[flatSection])
    
    smooth = .01
    
    #colors for hk lines in each plot
    hColor = 'dodgerblue'
    kColor = 'turquoise'
    

    #comments old, trust no one
    #create a pdf page with 4 rows and 4 columns
    #the bottom 3 rows use all 3 columns but the HK windows output is split inti
    fig, ax = plt.subplots(figsize=(10,10))
    plt.suptitle(title)
    plt.ticklabel_format(useOffset=False)
    with PdfPages(reportPath) as curPdf:
        gs = gridspec.GridSpec(4, 2)
        
        #references to each plot
        targPlt = plt.subplot(gs[1,:])
        smoothedPlt = plt.subplot(gs[2,:])
        flatPlt = plt.subplot(gs[3,:])
        kPlt = plt.subplot(gs[0,0])
        hPlt = plt.subplot(gs[0,1])
        rPlt =plt.subplot(gs[1,1])
        bPlt =plt.subplot(gs[1,0])


        #could turn some of these into loops by making a list of all plots but idk
        #please matplotlib don't make my stuff hard to read!
        hPlt.ticklabel_format(useOffset=False)
        kPlt.ticklabel_format(useOffset=False)
        rPlt.ticklabel_format(useOffset=False)
        bPlt.ticklabel_format(useOffset=False)
        targPlt.ticklabel_format(useOffset=False)
        flatPlt.ticklabel_format(useOffset=False)
        smoothedPlt.ticklabel_format(useOffset=False)
        
        hPlt.tick_params(axis='both',which= 'major', labelsize=7)
        kPlt.tick_params(axis='both',which= 'major', labelsize=7)
        rPlt.tick_params(axis='both',which= 'major', labelsize=7)
        bPlt.tick_params(axis='both',which= 'major', labelsize=7)
        targPlt.tick_params(axis='both',which= 'major', labelsize=7)
        flatPlt.tick_params(axis='both',which= 'major', labelsize=7)
        smoothedPlt.tick_params(axis='both',which= 'major', labelsize=7)
        
        #titles and lebels
        kPlt.set_xlabel("Wavelength(nm)")
        hPlt.set_ylabel("Irradiance")
        
        hPlt.set_title("Ca-H window, offset: " + str(round(hOffset,4))+"nm")
        kPlt.set_title("Ca-K window, offset: " + str(round(kOffset,4))+"nm")
        rPlt.set_title("R band, center: " + str(round(h.lamR,4))+", width: " + str(round(h.conWid,2))+", offset: " + str(round(rOffset,4))+"nm",fontsize=10)
        bPlt.set_title("V band, center: " +  str(round(h.lamB,4))+", width: " + str(round(h.conWid,2))+", offset: " + str(round(bOffset,4))+"nm",fontsize=10)
      
        #@TODO Sept 2021 bug, code would crash here-since we're not making a targplt plot?/??
        #targPlt.set_title("Target overlap over reference spectra")
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

        #cur=windows[:,0] 9/3/2020
        cur=windows[0]
        hkWidth=h.lineWid + .5 #setting how zoomed out the windows will be
        rWidth =h.conWid/2 +.05
        bWidth =h.conWid/2 +.05
        
        yScale = 1.4#how much more than the max we use


        hPlt.axvline(x=calH,color=hColor)
        hPlt.plot(oGrid-hOffset,obs,'b-')
        hPlt.axvline(x=h.cahLam-h.lineWid,color='gray')
        hPlt.axvline(x=h.cahLam+h.lineWid,color='gray')
        hPlt.set_xlim(calH-hkWidth,calH+hkWidth)
        lowY=min(min(obs[cur!=0]),0)
        hPlt.set_ylim(lowY,max(obs[cur!=0])*yScale)
        
        cur=windows[1]
        kPlt.axvline(x=calK,color=kColor)
        kPlt.axvline(x=h.cakLam-h.lineWid,color='gray')
        kPlt.axvline(x=h.cakLam+h.lineWid,color='gray')
        kPlt.plot(oGrid-kOffset,obs,'b-')
        kPlt.set_xlim(calK-hkWidth,calK+hkWidth)
        lowY=min(min(obs[cur!=0]),0)
        kPlt.set_ylim(lowY,max(obs[cur!=0])*yScale)

        cur=windows[3]
        bPlt.plot(oGrid[cur!=0]-bOffset,obs[cur!=0])
        bPlt.set_xlim(lamB-bWidth, lamB+bWidth)

        #get blue lines
        bMin = (oGrid[cur!=0]-bOffset)[0]
        bMax = (oGrid[cur!=0]-bOffset)[-1]
        
        cur=windows[2]
        rPlt.plot(oGrid[cur!=0]-rOffset,obs[cur!=0])
        rPlt.set_xlim(lamR-rWidth, lamR+rWidth)
        
        #get red lines from windows funct too
        rMin = (oGrid[cur!=0]-rOffset)[0]
        rMax = (oGrid[cur!=0]-rOffset)[-1]


        #terrrrible way to get scale TODO fixxxxxxx
        scale = obs[cur!=0]/base[cur!=0]
        avgS = np.mean(scale)
        
        targPlt.axvline(x=calH,color=hColor)
        targPlt.axvline(x=calK,color=kColor)
        targPlt.axvline(x=rMin, color='red')
        targPlt.axvline(x=rMax, color='red')
        targPlt.axvline(x=bMin, color='purple')
        targPlt.axvline(x=bMax, color='purple')

        #if there is a flat to plot then only use the target data when it is greater than .4
        #If we dont do this then on the edges out graph will be divided by a very small number which
        #makes it hard to see the data in the middle that we care about
        if flatMax != 0:
            targPlt.set_xlim([mini,maxi])
            targPlt.plot(bGrid[flatSection],base[flatSection]*avgS,color='lightgray')
            targPlt.plot(oGrid[flatSection]-hOffset,obs[flatSection],'b-')
            #print(mini)
            #print(maxi)
        else:
            targPlt.set_xlim([391.5,407])#BADAM
            targPlt.plot(bGrid[base!=0],base[base!=0]*avgS,color='lightgray')
            targPlt.plot(oGrid[obs!=0]-hOffset,obs[obs!=0],'b-')
        
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

        smoothedPlt.axvline(x=bMin, color='purple')
        smoothedPlt.axvline(x=bMax, color='purple')
        
        #if there is a flat to plot then only use the target data when it is greater than .4
        #If we dont do this then on the edges out graph will be divided by a very small number which
        #makes it hard to see the data in the middle that we care about
        if flatMax != 0:
            smoothedPlt.set_xlim([mini,maxi])
            smoothedPlt.plot(bGrid[flatSection],gdBase[flatSection]*avgS,color='lightgray')
            smoothedPlt.plot(oGrid[flatSection]-hOffset,gdObs[flatSection],'b-')
        else:
            smoothedPlt.set_xlim([391.5,407])
            smoothedPlt.plot(bGrid[gdBase!=0],gdBase[gdBase!=0]*avgS,color='lightgray')
            smoothedPlt.plot(oGrid[gdObs!=0]-hOffset,gdObs[gdObs!=0],'b-')
        
        #flat/other plot
        flatPlt.plot(bGrid[flat!=0], flat[flat!=0], 'k-')
        
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        #plt.show()
        plt.close()
        
        curPdf.savefig(fig)