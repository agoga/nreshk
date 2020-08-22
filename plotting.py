#this function plots the final time series of SHK values for each star
#and is the reason for the lame format of the data array.
#data format
#first array is list of mjds
#data[0].append(mjd)
#second array is list of shks
#data[1].append(shk)
#third array is header data and bool saying whether it is a single observation(false) or average
#data[2].append(tuple((header,False)))
#fourth array is data for the nightly observation function: grid, grid offset, and spectra
#data[3].append(tuple((lamGrid,correlation[0],targOlapf)))
def plot_daily_data_timeseries(inData,starName,bad):
    import numpy as np
    import matplotlib.patches as mpatches
    from matplotlib import pyplot as plt
    from astropy.time import Time
    import astropy.io.fits 


    mjdArray=[o.mjd for o in inData]
    shkValArray=[o.shk for o in inData]

    headerArray = [o.header for o in inData]

    #headerArray=np.asarray(tH)
    boolArray=np.asarray([o.average for o in inData])
    #mjdArray = np.asarray(inData[0])
    #shkValArray = np.asarray(inData[1])
    #headerArray = np.asarray(inData[2])
    
    fig, ax = plt.subplots(figsize=(12,6))
    

    pltStr = []
    for i in range(len(inData)):
        #h[0] is header
        #h[1] is bool saying wheather it's a average observation or single
        h=headerArray[i]
        b=boolArray[i]
        s = h['SITEID']
        tStr =''

        if s == 'lsc':
            tStr +='b'
        elif s == 'cpt':
            tStr +='g'
        elif s== 'elp':
            tStr +='r'
        elif s=='tlv':
            tStr +='k'

        if b == True:
            tStr += 'o'
        else:
            tStr += '^'   
        pltStr.append(tStr)
    
    t= Time(mjdArray, format='mjd')
    t.format = 'decimalyear'
    
    for i in range(len(pltStr)):
        mark = pltStr[i][1]
        col = pltStr[i][0]
        size = 150
        opac = .3
        if mark == 'o':
            size = 100
            opac = 1

        if mjdArray[i] not in bad:    
            plt.scatter(t[i].value,shkValArray[i],marker=mark,c=col, s=size, alpha = opac)
        else:
            plt.scatter(t[i].value,shkValArray[i],marker='x',c=col, s=100, alpha = opac)

    rl = mpatches.Patch(color='red', label='McDonald Obs\'')
    bl = mpatches.Patch(color='blue', label='Cerro Tololo Interamerican Obs\'')
    gl = mpatches.Patch(color='green', label='South African Astro Obs\'')
    kl = mpatches.Patch(color='black', label='Wise Obs\'')
    plt.legend(handles=[rl,bl,gl,kl],prop={'size': 10}, loc=4)
    

    plt.style.use('classic')
    ax.ticklabel_format(useOffset=False)
    plt.title('HD '+starName+' magnetic activity time series')
    plt.xlabel('Time(years)')
    plt.ylabel('Unadjusted S-index')
    plt.savefig('output/'+starName+'/'+starName+'_shk_time_series.pdf')
    plt.show()
    plt.close()



##
##This is the massive printing function to create a report of each observation pushed through the
#pipeline.
def pdf_from_intermediate_data(bGrid, base, oGrid, obs, windows, title, path, descript, flat='', width=1):
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