#this is where single run functions such as creating the flat pickle file(flat database) and loading aligning spectra

from pipeline import import_aligning_spectra
from NRES_SHK_Pipeline import create_flat_dict_file
#from matplotlib.backends.backend_pdf import PdfPages
import NRES_SHK_Pipeline as nresshk
import pipeline_testing as t
import helpers as h

#Soon to be settings.xml
#Project folder directory on your computer
homePath = 'D:/My Documents/Coding Projects/nreshk/'

#required paths for pipeline, try not to change ever
dataPath = homePath + 'data/'
outputPath = homePath + 'output/'
flatPath = dataPath + 'flats/'
flatPickle = "flatDict.pkl"

create_flat_dict_file(flatPath,flatPickle)

#need this in memory for the following wrapper
#TODO maybe a better way to do this without lugging the lab spectra around every install
lab = import_aligning_spectra(dataPath+ 'LabSpectra/', resolution=0)#,res*10)

flatDict = nresshk.NRES_SHK_MkFlat(flatPath,flatPickle)

print('done loading lab spectra')

#normal run start up

#MJD dates known to be bad TODO better feature use, maybe not w/ mjd

#skip = [58088.0763592, 58088.0801385, 58162.0245788, 58162.028358, 58272.4211262, 58272.4248934, 58354.2323109, 58354.2381758, 58553.7479017, 58553.7527935, 58821.1965842, 58821.2014756]
badnres=[]#[59195.7275335,59241.8083241]#should be empty list if not skipping any


#a list of stars we MUST run, regardless of if there is already output
#it will not run any other star if this is None
forceRun = ["22049"]#["152391"]#,["100180"]
#forceRun = None  


#manualAdj is for observations which are not being offset properly. combine with 'only' functionality to tweek the adjustment
#this is not really desirable but we'd rather have correct offsets than not and if automatic not working oh well

#window order taken from helpers.analyzedData offsets 9/10/20 Ca H, Ca K, R band, B band
#manAdj = [[58045.1409353],
        #[[-.04,-.04,-.0555,-.044]

        #]]
manAdj = None


#debug mode only has on or off at the moment and will cause lots of printing. Not really implemented and should not be trusted
h.debug = False


#output mode
#0 - all output 
#1 - only final time series pdf's created
#2 - only star N/A
#etc
h.pdfMode = 0


#for testing ONLY specific obs should be None if not testing specific mjd

#None#[58850,58821,58791,58791]
only=None

#TODO impement|
#if using radial velocity to align must have dict of star HD's and their rv"
alignmentDict = []#currently N/A"

#t.test_daily_data_sum()

analyzedData = nresshk.NRES_SHK_Pipeline(dataPath,outputPath,flatDict,lab,badnres,forceRun,manAdj,only)