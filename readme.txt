----------------------------------------------------------------------
This pipeline takes raw data from LCO NRES archives and creates a timeseries of the S-Index for each star it is run on. The S-index is a proxy for magnetic activity in stars. Please see the Mount Wilson Observatory HK project for more information.

Authors: Adam Goga - adfgoga@gmail.com
         Ricky Egeland - egeland@ucar.edu
         Travis Metcalfe - travis@spsci.org
         Tim Brown (original pipeline based on Tim Brown's IDL pipeline)


Notes:
Working.ipynb is the main Jupyter notebook which contains the wrapper for the pipeline.

The following is the required folder structure to get the pipeline to run on your data:
nreshk > python files.py
       > data\
            > 22049\ #Star folders must start with a number
            > 17051\ #another example of star folder
            > flats\ #must contain only flat files, the pipeline *should* be able to read flat files in any sub directory
            > LabSpectra\ #the reference spectra folder which will be loaded in for cross correlation
       >output\ #created by the pipeline to show time series information on each star and nightly observation's debug info
      
Star directories should contain sub-directories for each observation, maintaining the directory name contained within the .tar.gz files obtained from the LCO archive. I.e. every .fits file representing an observation has it's own folder. Currently the pipeline needs the reduced 'e91' observations from LCO database not the 'e00' files.

All spectra and flat files should be in .fits format; not the compressed .fits.fz format that comes from the LCO archive.  Use `funpack` to decompress the files.

At this time, stellar effective temperatures are hard-coded into a dictionary named `tEffLookup` in the helpers.py module.  If you are running the pipeline on a star not listed in that dictionary, it must be added first.  (A better configuration method will be implemented in a future version.)



Additional information and references:
Physics today overview - https://physicstoday.scitation.org/doi/10.1063/PT.3.3956
Metcalfe et. al (2016) - https://ui.adsabs.harvard.edu/abs/2016ApJ...826L...2M/abstract
Metcalfe & van Saders - https://ui.adsabs.harvard.edu/abs/2017SoPh..292..126M/abstract


<a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-sa/4.0/88x31.png" /></a><br /><span xmlns:dct="http://purl.org/dc/terms/" href="http://purl.org/dc/dcmitype/Text" property="dct:title" rel="dct:type">NRES SHK Pipeline</span> by <a xmlns:cc="http://creativecommons.org/ns#" href="https://github.com/agoga/nreshk/" property="cc:attributionName" rel="cc:attributionURL">Adam Goga</a> is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/">Creative Commons Attribution-ShareAlike 4.0 International License</a>.
