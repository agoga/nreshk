#This function should not be run every time the pipeline is run.
#It only needs to be run when you have added new flat files to the data
#A two dimensional dictionary is created for quick lookup of flat file 
#names by their site and MJD's and stored in a pickle file to be loaded
#in every time the pipeline is run
def create_flat_dict_file(flatDir,fileName):
    import pickle
    import pprint
    import os
    import astropy.io.fits
    #find every flat file and put into dictionary
    flatDict = dict() 
    flatDict.update({'cpt':dict()})
    flatDict.update({'tlv':dict()})
    flatDict.update({'lsc':dict()})
    flatDict.update({'elp':dict()})


    #this loops through all the flat_....fits files in the flats directory
    #and creates the 2 dimensional dictionary  
    for path, subdirs, files in os.walk(flatDir):
        for name in files:
            if name.startswith("flat_") and name.endswith(".fits"):
                flHDu1 = astropy.io.fits.open(os.path.join(path,name))
                header = flHDu1[0].header
                
                mjd = float(header['MJD-OBS'])
                site = header['SITEID']
                
                siteDict = flatDict[site]
                siteDict.update({mjd:os.path.join(path,name)})
                
    f = open(flatDir+fileName,"wb")
    pickle.dump(flatDict,f)
    f.close()
    
    
#This is the NRES specific function goes through the current folder(which should be one observation)
#and loads, from fits files, the wavelength grid, spectra, header information
#This function is only needed when there is both old and new data from NRES and we don't know which
#is which. When all the data is in the new format then half this function can be removed.
def load_folder_for_pipeline(path,dirs,curFiles):
    import os
    import astropy.io.fits
    
    waveGrid = []
    spec = []
    header = []
    fileName = ''#spec file name
    
    oldFormat = False
    #see if the folder contains a -noflat or -wave file
    #There may be a faster way to determine if this folder is old or new but I don't see a safer way
    for file in curFiles:
        if file.endswith("-noflat.fits") or file.endswith("-wave.fits"):
            oldFormat = True
            break

    #we've now looped through all the files if any are old then we try old format
    if oldFormat:
        #print('old format')
        #load the data in the old format
        for file in curFiles:
            if file.endswith("-noflat.fits"):
                fileName = file
                specHDu = astropy.io.fits.open(os.path.join(path,dirs,file))
            elif file.endswith("-wave.fits"):
                wvHDu = astropy.io.fits.open(os.path.join(path,dirs,file))

        spec = specHDu[0].data[0]
        header = specHDu[0].header
        
        #need the length here because sometimes the wave grid is in the first element
        #of this array but most of the time it's not 
        if len(wvHDu)>1 and type(wvHDu[1]) == astropy.io.fits.hdu.table.BinTableHDU:
        #print('wavehdu')   
        #print(wvHDu.info())
        #print(wvHDu[1].data.dtype)
        #print(type(wvHDu[1].data))
        #print('done')
        #print('stuff ' + str(i))
        #print((wvHDu[1].data)[0])
        
            #there are two different names for this data, it's mostly in the first name but
            #sometimes in the second
            try:
                waveGrid = wvHDu[1].data['WavelenStar'][0]
            except:
                waveGrid=wvHDu[1].data['Wavelength'][0]
                
        elif len(wvHDu)>1:
            waveGrid = (wvHDu[1].data)[0]
        else:
            #print(wvHDu[0].data.dtype)
            waveGrid = (wvHDu[0].data)[0]
            
        specHDu.close()
        wvHDu.close()
    else: #ELSE IT'S NEW FORMAT and super easy
        for file in curFiles:
            if file.endswith(".fits"):
                fileName = file
                oHDu = astropy.io.fits.open(os.path.join(path,dirs,file))
                waveGrid = oHDu[7].data
                spec = oHDu[1].data
                header = oHDu[0].header
                oHDu.close()
   #end big if              
    if len(waveGrid) == 0 or len(spec) == 0 or len(header) == 0:
        print('BIG PROBLEM IN load_folder_for_pipeline')
                    
    return waveGrid, spec, header, fileName, oldFormat



#Takes in the data array with reference spectra and set name, reference spectra and setname used 
#for printing.
#This function goes through the data array for each observation and finds any other obs on that day
#It will combined these observations into one targOLap and run the rest of the pipeline on the 
#averaged data to give a better signal to noise for these observations.
#the output is a longer data array to be used later for printing the final 
#This is how the data structure is created and maintained:
#data[0].append(mjd)
#data[1].append(shk)
#data[2].append(tuple((header,False)))
#data[3].append(tuple((lamGrid,correlation[0],targOlapf)))
#E
#correlation[0] is the delta lamda to place spectra in lab frame and lamGrid is the lab frame
def sum_daily_data(setName,data,labSpec):
    from scipy import interpolate
    import numpy as np
    from calc_shk import calc_shk
    from calc_shk import smart_calc_shk
    from helpers import mkdir_p
    from astropy.time import Time
    from helpers import pdf_from_data
    import helpers as h
    
    rv = 0
    done = []#array to hold which MJD are done

    #make a clean copy
    dailyData = np.asarray(data)
    
    #IF YOU DONT INCLUDE [:] THEN THIS IS A REFERENCE AND LIFE WILL STINK WHEN SORTING
    sortedList = data[0][:]
    datesList = dailyData[0]
    for i in range(len(datesList)):
        curD = dailyData[:,i]
        header = curD[2][0]
        site = header['SITEID']

        #dont double up if we've alrady DONE this day
        curMjd = curD[0]
        if(curMjd in done):
            continue


        #find the closest MJD's to the current
        sortedList.sort(key=lambda x: abs(x - curMjd))

        #uncomment this if you don't want nights with only 1 obs to have circles
        #if there are no observations on this day(using index 1 because 0 is the same obs)
        #if abs(curMjd - sortedList[1]) >= 1:
            #continue
            
        #grab the current day's values
        curTuple = curD[3]

        curLamGrid=curTuple[0]
        curTargOlapf=curTuple[2]
        curDLam=curTuple[1]
        
        
        #need to interpolate because all the observations have slightly different lamda grids
        #put them all onto the referece spectra's grid held in curLamGrid
        interpfunc = interpolate.interp1d(curLamGrid-curDLam, curTargOlapf, kind='linear',fill_value='extrapolate')
        curTarg=interpfunc(curLamGrid)
        
        combinedTarg = curTarg
        
        done.append(curMjd)
        for i in range(len(sortedList)):
            same = sortedList[i]
            
            #don't do the current date
            if same == curMjd:
                continue
                
            #the list is sorted so once we reach a point where the obs is more than 
            #a day apart we're done
            if abs(curMjd - same) >= 1:
                break
            #print(str(curMjd) + " also has " + str(same))
            
            #create an array of all the spectra that have been done, i..e these two are now done
            done.append(same)

            #grab the same day's data(will happen for multiple obs
            sameI = np.argmin(abs(datesList-same))

            sameTuple = dailyData[:,sameI][3]

            sameLamGrid=sameTuple[0]
            sameTargOlapf=sameTuple[2]
            sameDLam=sameTuple[1]


            

            #The issue is the two targetOlaps are on different grids. 
            #Interp them on the same grid then add
            interpfunc = interpolate.interp1d(sameLamGrid-sameDLam, sameTargOlapf, kind='linear',fill_value='extrapolate')
            sameTarg=interpfunc(curLamGrid)

            #this is the big deal, increases signal to noise for each days spectra
            combinedTarg = combinedTarg+sameTarg
        #end combining loop
        
        
        #printing  stuff
        first = str(curMjd).split('.')[0]
        
        #lookup this star's teff
        tempEff = h.tEffLookup[setName.strip('/')]

        #find SHK with new offset to lamda grid
        shkRet = smart_calc_shk(curLamGrid, combinedTarg, rv, teff=tempEff)

        shk = shkRet[0]
        windows = shkRet[1]

        mkdir_p("output/"+setName+"/"+first+'/')
        np.savez("output/"+setName+"/"+first+"/combined_data", targOlapf=combinedTarg,flatOlap=combinedTarg, lamGrid=curLamGrid, adjLamGrid=curLamGrid,windows=windows)


        decimalYr = Time(curMjd,format='mjd')
        decimalYr.format = 'decimalyear'
        title = 'NRES spectra, ' + site +', '+header['DATE-OBS']+' ('+ '{:.6}'.format(decimalYr.value) +'), S='+'{:.4}'.format(shk)
        pdf_from_data(curLamGrid, labSpec,curLamGrid, combinedTarg,windows,title, "output/"+setName+"/"+first+"/","combined",width=.3)

        #print("c " + str(type(curMjd))+ "s " +str(type(shk)))
        #print("adding " + str(curMjd) + " with shk " + str(shk))
        data[0].append(curMjd)
        data[1].append(shk)
        data[2].append(tuple((header,True)))
        data[3].append(tuple((curLamGrid,0,combinedTarg)))

    return data


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
def plot_daily_data_timeseries(data,setName,bad):
    import numpy as np
    import matplotlib.patches as mpatches
    from matplotlib import pyplot as plt
    from astropy.time import Time
    
    
    mjdArray = np.asarray(data[0])
    shkValArray = np.asarray(data[1])
    headerArray = np.asarray(data[2])



    fig, ax = plt.subplots(figsize=(12,6))

    pltStr = []
    for h in headerArray:
        #h[0] is header
        #h[1] is bool saying wheather it's a average observation or single
        s = h[0]['SITEID']
        tStr =''

        if s == 'lsc':
            tStr +='b'
        elif s == 'cpt':
            tStr +='g'
        elif s== 'elp':
            tStr +='r'
        elif s=='tlv':
            tStr +='k'

        if h[1] == True:
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
    plt.title('HD '+setName+' magnetic activity time series')
    plt.xlabel('Time(years)')
    plt.ylabel('Unadjusted S-index')
    plt.savefig('output/'+setName+'/'+setName+'_shk_time_series.pdf')
    plt.show()
    plt.close()