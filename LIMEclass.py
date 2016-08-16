import fileinput
import numpy as np
import scipy.constants as sc
import headerclass as header 

# Define all the properties for LIME runs.      
     
class model:
    
    def __init__(self,
                 fileout,
                 headerfile,
                 moldatfile,
                 transitions, 
                 thetas=[0],
                 nchan=200,
                 velres=20,
                 pIntensity=1e4,
                 sinkPoints=1e3,
                 antialias=1,
                 sampling=2,
                 lte_only=1,
                 imgres=0.05,
                 distance=54.,
                 pxls=128,
                 unit=1,
                 stellarmass=0.6,
                 phis=None,
                 outputfile=None,
                 binoutputfile=None,
                 gridfile=None,
                 opratio=None,
                 dtemp=None,
                 xmol=None,
                 g2d=None,
                 bvalue=None,
                 btype='absolute',
                 coordsys='cylindrical',
                 dust='jena_thin_e6.tab',
                 directory='../',
                 ):
        
        # Assign all the variables.
        self.fileout = fileout
        self.hdr = header.headerFile(headerfile, coordsys=coordsys)
        self.moldatfile = moldatfile
        self.transitions = transitions
        self.thetas = thetas
        self.phis = phis
        self.nchan = nchan
        self.velres = velres
        self.pIntensity = pIntensity
        self.sinkPoints = sinkPoints
        self.antialias = antialias
        self.sampling = sampling
        self.lte_only = lte_only
        self.imgres = imgres
        self.distance = distance
        self.pxls = pxls
        self.unit = unit
        self.stellarmass = stellarmass
        self.outputfile = outputfile
        self.binoutputfile = binoutputfile
        self.gridfile = gridfile
        self.opratio = opratio
        self.dtemp = dtemp
        self.xmol = xmol
        self.g2d = g2d
        self.bvalue = bvalue
        self.btype = btype
        self.coordsys = coordsys
        self.fileout = fileout
        self.dust = dust
        self.directory = directory
        
        return
                

