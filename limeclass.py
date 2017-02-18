"""
Set all the parameters for running LIME models. Part of the makeLIME package.
"""

import os
import headerclass as header
import lamdaclass as rates
import numpy as np


class model:

    def __init__(self, headerfile, moldatfile, **kwargs):

        # This class is called from `tempfolder`.
        # Set up the paths appropriately.
        # TODO: Allow for dust / collisional rates in current directory.

        self.path = os.path.dirname(__file__)
        self.auxfiles = self.path + '/AuxFiles/'
        self.directory = kwargs.get('directory', '../')
        if not os.path.exists(self.directory):
            os.mkdir(self.directory)

        # LIME model parameters, read the manual for more info.
        # https://lime.readthedocs.io/en/latest/usermanual.html#parameters
        # These should mostly be read from the header file.

        if not self.checkexistance(headerfile, '../'):
            raise ValueError('%s not found.' % headerfile)
        self.header = header.headerFile('../' + headerfile)

        self.radius = min(kwargs.get('r_max', 1e10), self.header.rmax)
        self.minScale = max(kwargs.get('r_min', 0.), self.header.rmin)
        if self.minScale >= self.radius:
            raise ValueError('radius < minScale')

        self.pIntensity = float(kwargs.get('pIntensity', 1e5))
        self.sinkPoints = float(kwargs.get('sinkPoints', 3e3))
        if self.sinkPoints > self.pIntensity:
            print("Warning: sinkPoints > pIntensity.")

        self.sampling = int(kwargs.get('sampling', 2))
        if self.sampling not in [0, 1, 2]:
            raise ValueError('sampling must be 0, 1 or 2.')

        self.tcmb = float(kwargs.get('tcmb', 2.725))
        self.moldatfile = self.verifymoldatfile(moldatfile)

        self.dust = kwargs.get('dust', 'jena_thin_e6.tab')
        if not self.checkexistance(self.dust, self.auxfiles):
            raise ValueError('%s not found.' % self.dust)

        # Options to turn on the generation of additional output files. There
        # are four different stages of output. Here we just consider the final
        # one. TODO: How do we deal with all the differeent output files from
        # runs with multiple models?

        self.outputgrid = kwargs.get('outputgrid', False)
        if self.outputgrid and self.outputgrid[-3:] != '.ds':
            self.outputgrid = self.outputgrid + '.ds'

        # self.returnoputs = kwargs.get('outputfile', False)
        # self.returnbputs = kwargs.get('binoutputfile', False)
        # self.returngrids = kwargs.get('gridfile', False)

        self.restart = kwargs.get('restart', None)
        if self.restart is not None:
            raise NotImplementedError()

        self.pregrid = kwargs.get('pregrid', None)
        if self.pregrid is not None:
            raise NotImplementedError()

        self.lte_only = int(kwargs.get('lte_only', 1))
        if self.lte_only > 1:
            raise ValueError('lte_only must be 0 or 1.')
        elif self.lte_only == 0:
            print 'Warning: Running non-LTE model. Will be slow.'

        self.blend = int(kwargs.get('blend', 0))
        if self.blend > 1:
            raise ValueError('blend must be 0 or 1.')

        self.antialias = int(kwargs.get('antialias', 1))
        if self.antialias > 4:
            print "Warning: Large antialias value."

        self.polarization = kwargs.get('polarization', None)
        if self.polarization is not None:
            raise NotImplementedError()

        # If oH2 and pH2 are available collision partners, use them over H2.
        # If no ortho / para ratio is given, assume a value of 3.
        # self.opr is then the n(H2) -> [n(oH2), n(pH2)] rescaling value.

        if (self.oH2 and self.pH2):
            print 'Using oH2 and pH2 as collision partners.'
            self.opr = kwargs.get('opr', 3.)
            self.opr = np.array([self.opr, 1.]) / (1. + self.opr)
            self.collpartIds = [3, 2]
        else:
            print 'Using H2 as single collision partner.'
            self.opr = [1.0]
            self.collpartIds = [1]

        # We then set the molecular and dust weights. Assume that they're
        # equally distributed between them.

        self.nMolWeights = np.array([1.0 for cId in self.collpartIds])
        self.dustWeights = np.array([1.0 for cId in self.collpartIds])

        self.traceRayAlgorithm = int(kwargs.get('traceRayAlgorithm', 0))
        if self.traceRayAlgorithm > 1:
            raise ValueError('traceRayAlgorithm must be 0 or 1.')

        self.nThreads = int(kwargs.get('nThreads', 20))

        # LIME imaging parameters, read the manual for more info.
        # https://lime.readthedocs.io/en/latest/usermanual.html#images

        self.pxls = int(kwargs.get('pxls', 128))
        self.imgres = float(kwargs.get('imgres', 0.065))
        self.distance = float(kwargs.get('distance', 59.))
        self.unit = int(kwargs.get('unit', 1))

        if self.unit not in [0, 1, 2, 3]:
            raise ValueError("unit must be 0, 1, 2 or 3.")

        # We want to remove any directory before the filename,
        # and make sure that there's no extension on the end.

        self.name = kwargs.get('name', self.header.fn[:-2])
        self.name = self.name.split('/')[-1]
        if '.' in self.name:
            self.name = ''.join(self.name.split('.')[:-1])

        self.source_vel = float(kwargs.get('source_vel', 0.0))
        self.source_vel *= 1e3

        self.nchan = int(kwargs.get('nchan', 100))
        self.velres = float(kwargs.get('velres', 50.))

        self.transitions = kwargs.get('transitions', [0])
        if type(self.transitions) is not list:
            self.transitions = [self.transitions]

        # Need to include more than one molecule at a time.

        self.molI = kwargs.get('molI', None)
        if self.molI is not None:
            raise NotImplementedError()

        self.freq = kwargs.get('freq', None)
        if self.freq is not None:
            raise NotImplementedError()

        self.bandwidth = kwargs.get('bandwidth', None)
        if self.bandwidth is not None:
            raise NotImplementedError()

        # LIME image rotation parameters.
        # https://lime.readthedocs.io/en/latest/usermanual.html#image-rotation-parameters

        self.incl = kwargs.get('incl', [0])
        if type(self.incl) is not list:
            self.incl = [self.incl]

        self.posang = kwargs.get('posang', [0])
        if type(self.posang) is not list:
            self.posang = [self.posang]

        self.azimuth = kwargs.get('azimuth', [0])
        if type(self.azimuth) is not list:
            self.azimuth = [self.azimuth]

        # LIME model functions.
        # https://lime.readthedocs.io/en/latest/usermanual.html#model-functions
        # Read https://github.com/richteague/makeLime/blob/master/README.md
        # for details on the format of some of these.
        # TODO: Include non-Keplerian rotation.

        self.dens = kwargs.get('dens', 'dens')
        self.abund = kwargs.get('abund', 'abund')
        self.temp = kwargs.get('temp', 'temp')
        self.dtemp = kwargs.get('dtemp', None)

        self.rescaletemp = kwargs.get('rescaletemp', False)

        self.checkTypes(self.dens, [str])
        self.checkTypes(self.abund, [str, float])
        self.checkTypes(self.temp, [str])
        self.checkTypes(self.dtemp, [str, float, None])

        self.doppler = kwargs.get('doppler', 50.)
        self.dopplertype = kwargs.get('dopplertype', 'absolute')
        self.checkTypes(self.doppler, [str, float])
        if not (self.dopplertype == 'absolute' or self.dopplertype == 'mach'):
            raise ValueError("dopplertype must be 'absolute' or 'mach'.")

        self.stellarmass = float(kwargs.get('stellarmass', 0.7))

        self.magfield = kwargs.get('magfield', None)
        if self.magfield is not None:
            raise NotImplementedError()

        self.g2d = kwargs.get('g2d', None)
        self.checkTypes(self.g2d, [str, float, None])

        self.gridDensity = kwargs.get('gridDensity', None)
        if self.gridDensity is not None:
            raise NotImplementedError()

        # makeLIME specific variables.
        # https://github.com/richteague/makeLime/blob/master/README.md

        self.ndim = self.header.ndim
        self.ncells = self.header.ncells
        self.nmodels = int(kwargs.get('nmodels', 1))
        self.returnnoise = kwargs.get('returnnoise', False)
        self.coordsys = self.header.coordsys

        self.ninc = len(self.incl)
        self.npos = len(self.posang)
        self.nazi = len(self.azimuth)
        self.ntra = len(self.transitions)

        # Check the field of view of the model. Raise a warning if this
        # is not enough to image the whole model.

        if (self.distance * self.imgres * self.pxls * 0.5 < self.radius):
            proj_dist = self.imgres * self.distance * self.pxls
            proj_size = 2. * self.radius
            print "Warning: Check distance and pixel scaling."
            print "\t Image has projected distance of %.2f au," % proj_dist
            print "\t but the model has a size of %.2f au." % proj_size

        # Additional variables. See readme for info.

        self.depletion = float(kwargs.get('depletion', 1.0))
        self.oversample = int(kwargs.get('oversample', 1))

        # If there is no dust, include a dummy density.
        # an opr = 0 (thus we rescale the density to zero).
        # Make this MolWeight = 0 and dustWeight = 1, but remove
        # all other dustweights.

        if not int(kwargs.get('includeContinuum', True)):
            self.opr = np.append(self.opr, 0.)
            self.collpartIds = np.append(self.collpartIds, 4)
            self.nMolWeights = np.append(self.nMolWeights, 0.)
            self.dustWeights = np.append(self.dustWeights * 0., 1.)
            self.dtemp = 0.0

        # niceness. If selected change the niceness of the runs.
        self.niceness = kwargs.get('niceness', False)

        return

    def checkTypes(self, val, types):
        """Check the array names are valid."""
        if not (type(val) in types or val in types):
            raise TypeError(self.typestr(val, types))
        elif (type(val) is str and val not in self.header.anames):
            raise ValueError('{} not found in {}'.format(val, self.header.fn))
        return

    def typestr(self, val, types):
        """Format a string for output."""
        string = '{} should be '.format(val)
        if len(types) == 1:
            string += 'a '
        elif len(types) == 2:
            string += 'either '
        else:
            string += 'one of: '
        for i, t in enumerate(types):
            try:
                string += '{}'.format(t.__name__)
            except:
                string += 'None'
            if i < len(types) - 2:
                string += ', '
            elif i < len(types) - 1:
                string += ' or '
            else:
                string += '.'
        return string

    def checkexistance(self, fn, folder):
        """Check if fn is in the specified folder."""
        return os.path.isfile(folder + fn)

    def verifymoldatfile(self, moldatfile):
        """Verify the collisional rates exist."""
        if not self.checkexistance(moldatfile, self.auxfiles):
            raise ValueError('%s not found.' % moldatfile)
        ratefile = rates.ratefile(self.auxfiles + moldatfile)
        self.H2 = 'H2' in ratefile.partners
        self.oH2 = 'oH2' in ratefile.partners
        self.pH2 = 'pH2' in ratefile.partners
        if (not self.H2 and not self.oH2 and not self.pH2):
            raise ValueError('No appropriate collisonal rates found.')
        return moldatfile
