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
                 sinkPoints=3e3,
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
                 dtemp=None,
                 xmol=None,
                 g2d=None,
                 bvalue=None,
                 btype='absolute',
                 coordsys='cylindrical',
                 dust='jena_thin_e6.tab',
                 directory='../',
                 nmodels=1,
                 returnNoise=False,
                 ):
        

        # Output configurations.
        
        self.fileout = fileout
        self.outputfile = outputfile
        self.binoutputfile = binoutputfile
        self.gridfile = gridfile
        self.returnNoise = returnNoise
        self.directory = directory
        self.nmodels = nmodels

        # Chemical model properties from header file.
        
        self.coordsys = coordsys
        self.hdr = header.headerFile(headerfile, coordsys=coordsys)
        self.ndim = self.hdr.ndim
        self.ncells = self.hdr.ncells
        self.rin = self.hdr.rin
        self.rout = self.hdr.rout
        self.bvalue = bvalue
        self.btype = btype
       
        # Additional model properties, if None specified, check if provided by
        # the chemical header, otherwise, revert to default value (=None).
        
        if dtemp is None: 
            if 'dtemp' in self.hdr.arrnames:
                self.dtemp = 'dtemp'
            else:
                self.dtemp = None
        else:
            self.dtemp = dtemp
        
        if (type(xmol) is not float and xmol is not None):
            raise TypeError("xmol needs to be a float or 'None'.") 
        else:
            self.xmol = xmol

        if g2d is None:
            if 'gastodust' in self.hdr.arrnames:
                self.g2d = 'gastodust'
            else:
                self.g2d = None
        else:
            if not (type(g2d) is float or type(g2d) is int or type(g2d) is str):
                raise TypeError("g2d needs to be a number (float or int) or an array name (str).")
            self.g2d = g2d

        # TODO: Include ortho-para colliders and molecules.
        # Currently set to deafult.
        self.opratio = None

        # LIME properties.    

        self.pIntensity = pIntensity
        self.sinkPoints = sinkPoints
        self.sampling = sampling
        self.lte_only = lte_only
        self.moldatfile = moldatfile
        self.stellarmass = stellarmass
        self.dust = dust

        # Imaging parameters.
        # Make sure the appropriate variables are lists.        

        if type(transitions) is not list:
            self.transitions = [transitions]
        else:
            self.transitions = transitions

        if type(thetas) is not list:
            self.thetas = [thetas]
        else:
            self.thetas = thetas

        if type(phis) is not list:
            self.phis = [phis]
        else:
            self.phis = phis

        self.nchan = nchan
        self.velres = velres
        self.antialias = antialias

        self.imgres = imgres
        self.distance = distance
        self.pxls = pxls
        self.unit = unit

        return
                

