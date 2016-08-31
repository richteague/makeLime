import os
import warnings
import fileinput
import numpy as np
import scipy.constants as sc
import headerclass as header 

# Define all the properties for LIME runs.
# Check all the variable types here.      
     
class model:
    
    def __init__(self,
                 name='output',
                 headerfile='header.h',
                 moldatfile='molecule.dat',
                 transitions=[1], 
                 inclinations=[0],
                 positionangles=[0],
                 nchan=200,
                 velres=20,
                 pIntensity=1e5,
                 sinkPoints=3e3,
                 antialias=1,
                 sampling=2,
                 lte_only=1,
                 imgres=0.065,
                 distance=54.,
                 pxls=128,
                 unit=0,
                 stellarmass=0.6,
                 outputfile=False,
                 binoutputfile=False,
                 gridfile=False,
                 dens=None,
                 temp=None,
                 dtemp=None,
                 abund=None,
                 g2d=None,
                 doppler=None,
                 dopplertype='absolute',
                 coordsys='cylindrical',
                 dust='jena_thin_e6.tab',
                 directory='../',
                 nmodels=1,
                 returnnoise=False,
                 ):
        
        self.path = os.path.dirname(__file__)
        self.auxfiles = self.path + '/AuxFiles'
        
        # Output configurations.
        
        if name[-5:] == '.fits':
            warnings.warn("Removing '.fits' from model name.")
            self.name = name[:-5]
        else:
            self.name = name
        
        if type(outputfile) is bool: 
            self.outputfile = outputfile
        else:
            raise TypeError("outputfile must be a Boolean.")
            
        if type(binoutputfile) is bool: 
            self.binoutputfile = binoutputfile
        else:
            raise TypeError("binoutputfile must be a Boolean.")            

        if type(gridfile) is bool: 
            self.gridfile = gridfile
        else:
            raise TypeError("gridfile must be a Boolean.") 

        if type(returnnoise) is bool: 
            self.returnnoise = returnnoise
        else:
            raise TypeError("returnnoise must be a Boolean.") 

        if type(directory) is str:
            self.directory = directory
        else:
            raise TypeError('directory must be a path.')
            
        if type(nmodels) is int:
            self.nmodels = nmodels
        elif type(nmodels) is float:
            self.nmodels = int(nmodels)
        else:
            raise TypeError('nmodels must be an integer.')

        if (coordsys == 'cylindrical' or coordsys == 'polar'): 
            self.coordsys = coordsys
        else:
            raise NotImplementedError("Only cylindrical or polar coordinates.")
            
        # Chemical model properties from header file.
        
        if type(headerfile) is str:
            self.hdr = header.headerFile(headerfile, coordsys=self.coordsys)
        else:
            raise TypeError('headerfile must be a path.')
        self.ndim = self.hdr.ndim
        self.ncells = self.hdr.ncells
        self.rin = self.hdr.rin
        self.rout = self.hdr.rout
        
        # Additional model properties, if None specified, check if provided by
        # the chemical header, otherwise, revert to default value (=None).
        
        self.dens = self.checkTypes(dens, 'dens', [str])
        self.abund = self.checkTypes(abund, 'abund', [str, float])
        self.temp = self.checkTypes(temp, 'temp', [str])
        self.dtemp = self.checkTypes(dtemp, 'dtemp', [str, float, None])
        self.g2d = self.checkTypes(g2d, 'gastodust', [str, float, None])
        self.doppler = self.checkTypes(doppler, 'doppler', [str, float, None])
        if (dopplertype == 'absolute' or dopplertype == 'mach'):
            self.dopplertype = dopplertype
        else:
            raise ValueError("dopplertype must be 'absolute' or 'mach'.")
        if type(stellarmass) is float:
            self.stellarmass = stellarmass       
        else:
            raise TypeError("stellarmass must be a float.")
            
        # TODO: Currently set to deafult.
        self.opratio = None

        # LIME properties.    
        
        if type(pIntensity) is float:
            self.pIntensity = pIntensity
        else:
            raise TypeError("pIntensity must be a float.")
        
        if type(sinkPoints) is float:
            self.sinkPoints = sinkPoints
        else:
            raise TypeError("sinkPoints must be a float.")
        
        if self.sinkPoints > self.pIntensity:
            warnings.warn("sinkPoints > pIntensity.")
        
        if (type(sampling) is int and sampling < 3):
            self.sampling = sampling
        else:
            raise ValueError("sampling must be 0, 1 or 2.")
        
        if (type(lte_only) is int and lte_only < 2):
            self.lte_only = lte_only
        elif type(lte_only) is bool:
            self.lte_only = int(lte_only)
        else:
            raise ValueError("lte_only must be True or False.")
        

        if os.path.isfile(self.auxfiles+'/'+moldatfile):
            self.moldatfile = moldatfile
        else:
            raise ValueError("No molecular data file found called %s." % moldatfile)
            
        if os.path.isfile(self.auxfiles+'/'+dust):
            self.dust = dust
        else:
            raise ValueError("No dust opacities file found called %s." % dust)            
    
        # Imaging parameters.
        # Make sure the appropriate variables are lists.        

        if type(transitions) is not list:
            self.transitions = [transitions]
        else:
            self.transitions = transitions

        if type(inclinations) is not list:
            self.thetas = [inclinations]
        else:
            self.thetas = inclinations

        if type(positionangles) is not list:
            self.phis = [positionangles]
        else:
            self.phis = positionangles
        if self.phis != [0]:
            raise NotImplementedError("No non-zero PAs.")

        if (type(nchan) is float or type(nchan) is int):
            self.nchan = float(nchan)
        else:
            raise TypeError("nchan must be a number.")
        
        if (type(velres) is float or type(velres) is int):
            self.velres = float(velres)
        else:
            raise TypeError("velres must be a number.")
        
        if (type(antialias) is float or type(antialias) is int):
            self.antialias = int(antialias)
        else:
            raise TypeError("antialias must be a number")
        if self.antialias > 4:
            warnings.warn("High antialias value of %d, might be slow." % self.antialias)
        
        if type(imgres) is float:
            self.imgres = imgres
        else:
            raise TypeError("imgres must be a float.")
        
        if type(distance) is float:
            self.distance = distance
        else:
            raise TypeError('distance must be a float.')
        
        if (type(pxls) is int or type(pxls) is float):
            self.pxls = int(pxls)
        else:
            raise TypeError('pxls must be a number.')

        if (self.distance * self.imgres * self.pxls * 0.5 < self.rout):
            warnings.warn("Check distance and pixel scaling.")
            warnings.warn("Image has projected distance of %.2f au." % (self.imgres * self.distance * self.pxls))
            warnings.warn("Model has a size of %.2f au." % (2. * self.rout))
        
        if (type(unit) is int and unit in [0, 1, 2, 3]):
            self.unit = unit
        else:
            raise ValueError("unit must be 0, 1, 2 or 3.")

        return
    
    def checkTypes(self, inval, default, types):
        if inval is None:
            if default in self.hdr.arrnames:
                return default
            elif None in types:
                return inval
            else:
                raise TypeError("Cannot find %s in %s." % (default, self.hdr.fn))
        elif type(inval) in types:
            return inval
        else:
            raise TypeError("%s" % self.typestostring(types))
        return
                

    def typestostring(self, types):
        if len(types) == 1:
            string = 'Must be a '
        elif len(types) == 2:
            string = 'Must be either '
        else:
            string = 'Must be one of: '
        for i, t in enumerate(types):
            try:
                string += '%s' % t.__name__
            except:
                string += 'None'
            if i < len(types) - 2:
                string += ', '
            elif i < len(types) - 1:
                string += ' or '
            else:
                string += '.'
        return string
