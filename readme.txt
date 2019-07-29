This pipeline takes raw data from LCO NRES archives and creates a timeseries of the S-Index for each star it is run on. The S-index is a proxy for magnetic activity in stars. Please see the Mount Wilson Observatory HK project for more information.

Authors: Adam Goga - adfgoga@gmail.com
         Ricky Egeland - egeland@ucar.edu
         Travis Metcalfe - travis@spsci.org
         Tim Brown (original pipeline based on Tim Brown's IDL pipeline)

The following is the required folder structure to get the pipeline to run on your data:
nreshk > python files.py
       > data\
            > 22049\ #Star folders must start with a number
            > 17051\ #another example of star folder
            > flats\ #must contain only flat files, the pipeline *should* be able to read flat files in any sub directory
            > LabSpectra\ #the reference spectra folder which will be loaded in for cross correlation
       >output\ #created by the pipeline to show time series information on each star and nightly observation's debug info
      
All spectra and flat files should be in .fits format
