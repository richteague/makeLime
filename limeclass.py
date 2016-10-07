import os
import headerclass as header
import lamdaclass as rates

# Define all the properties for LIME runs.
# Check all the variable types here.


class model:

    def __init__(self,
                 name='output',
                 headerfile='header.h',
                 moldatfile='molecule.dat',
                 transitions=[1],
                 inclinations=[0.5],
                 positionangles=[1.5708],
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
                 blend=0,
                 opr_cp=None,
                 depletion=1.0,
                 r_inner=None,
                 r_outer=None,
                 oversample=10,
                 ):

        print '\n'
        self.path = os.path.dirname(__file__)
        self.auxfiles = self.path + '/AuxFiles'
        self.oversample = int(oversample)


        # The name will be appended to all outputs. The .fits extension is not
        # wanted so removed.

        if name[-5:] == '.fits':
            print "Removing '.fits' from model name."
            self.name = name[:-5]
        else:
            self.name = name

        # Check that the requested collisional rate and dust opacity files are
        # in the AuxFiles directory.

        if os.path.isfile(self.auxfiles + '/' + moldatfile):
            self.moldatfile = moldatfile
        else:
            raise ValueError("no %s found in AuxFiles." % moldatfile)

        if os.path.isfile(self.auxfiles + '/' + dust):
            self.dust = dust
        else:
            raise ValueError("no %s found in AuxFiles." % dust)

        # This is the output directory. By default it is the starting folder.
        # We check that this is a valid path, otherwise it defaults to '../'.

        if type(directory) is str:
            self.directory = directory
        else:
            raise TypeError('directory must be a path.')
        if not os.path.exists(self.directory):
            print "Warning: output directory not found. Defaulting to '../'."
            self.directory = '../'

        # The output files are booleans. If set to true, they will return
        # the specified output files. The outputfile will combine all the
        # output files run for each model. gridnoise will return one grid
        # per model. TODO: combine binoutputfile.

        if type(outputfile) is bool:
            self.outputfile = outputfile
        else:
            raise TypeError("outputfile must be a Boolean.")

        if type(binoutputfile) is bool:
            self.binoutputfile = binoutputfile
        else:
            raise TypeError("binoutputfile must be a Boolean.")
        if self.binoutputfile:
            raise NotImplementedError()

        if type(gridfile) is bool:
            self.gridfile = gridfile
        else:
            raise TypeError("gridfile must be a Boolean.")

        if type(returnnoise) is bool:
            self.returnnoise = returnnoise
        else:
            raise TypeError("returnnoise must be a Boolean.")

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

        # Chemical model properties from header file. This automatically reads
        # in the number of dimensions (by counting c1arr, c2arr and c3arr),
        # the number of cells and the minimum and maximum radii.

        if type(headerfile) is str:
            self.hdr = header.headerFile('../'+headerfile, coordsys=self.coordsys)
        else:
            raise TypeError('headerfile must be a path.')
        self.ndim = self.hdr.ndim
        self.ncells = self.hdr.ncells
        self.rin = self.hdr.rin
        self.rout = self.hdr.rout

        # If inner and outer radii are specifically supplied, we will
        # overwrite those read from the header file.

        if r_inner is not None:
            self.rin = max(self.rin, r_inner)
        if r_outer is not None:
            self.rout = min(self.rout, r_outer)
        if self.rin >= self.rout:
            raise ValueError("Outer radius less than inner radius.")

        # Additional model properties, if None specified, check if provided by
        # the chemical header, otherwise, revert to default value (=None).

        self.dens = self.checkTypes(dens, 'dens', [str])
        self.abund = self.checkTypes(abund, 'abund', [str, float])
        self.temp = self.checkTypes(temp, 'temp', [str])
        self.dtemp = self.checkTypes(dtemp, 'dtemp', [str, float, None])
        self.g2d = self.checkTypes(g2d, 'gastodust', [str, float, None])
        self.doppler = self.checkTypes(doppler, 'doppler', [str, float, None])

        # The doppler type defines whether the microturbulence is specified in
        # meters per second, 'absolute', or a fraction of the local sound
        # speed, 'mach'.

        if (dopplertype == 'absolute' or dopplertype == 'mach'):
            self.dopplertype = dopplertype
        else:
            raise ValueError("dopplertype must be 'absolute' or 'mach'.")

        # The stellar mass is required for the Keplerian rotation.
        # TODO: Allow the inclusion of an array of velocity profiles.

        if type(stellarmass) is float:
            self.stellarmass = stellarmass
        else:
            raise TypeError("stellarmass must be a float.")

        # Number of grid points for lime. pIntensity should be greater than
        # sinkPoints but not required.

        if type(pIntensity) is float:
            self.pIntensity = pIntensity
        else:
            raise TypeError("pIntensity must be a float.")
        if type(sinkPoints) is float:
            self.sinkPoints = sinkPoints
        else:
            raise TypeError("sinkPoints must be a float.")
        if self.sinkPoints > self.pIntensity:
            print "Warning: sinkPoints > pIntensity."

        # Sampling as in the lim manual. sampling = 2 is spherical.

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
        if not self.lte_only:
            print 'Warning: Running non-LTE model. Will be slow.'

        # Imaging parameters.
        # Make sure the appropriate variables are lists.

        if type(transitions) is not list:
            self.transitions = [transitions]
        else:
            self.transitions = transitions

        # Must we specfiy the (inc, pa, azi) triplet.
        # Due to the old code style, these are theta and phi.
        # TODO: Include azimuthal values, but this is not necessary for
        # azimuthally symmetric models.

        if type(inclinations) is not list:
            self.thetas = [inclinations]
        else:
            self.thetas = inclinations

        if type(positionangles) is not list:
            self.phis = [positionangles]
        else:
            self.phis = positionangles

        # TODO: Include possibiliy of lists of values.
        self.azimuth = 0.0

        # Specify the spectral axis in meters per second.
        # TODO: Allow for a spectral scan.

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
            print "Warning: High antialias value of %d." % self.antialias
            print "\t Ray tracing may be slow..."

        # Specify the field of view. Will complain if the model is not fully
        # imaged with the specified image axes.

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
            proj_dist = self.imgres * self.distance * self.pxls
            proj_size = 2. * self.rout
            print "Warning: Check distance and pixel scaling."
            print "\t Image has projected distance of %.2f au," % proj_dist
            print "\t but the model has a size of %.2f au." % proj_size

        if (type(unit) is int and unit in [0, 1, 2, 3]):
            self.unit = unit
        else:
            raise ValueError("unit must be 0, 1, 2 or 3.")

        if not (blend == 0 or blend == 1):
            self.blend = 0
            print 'Warning: Unknown blend value. Set to 0.'
        else:
            self.blend = blend
        if self.blend:
            print 'Warning: Line blending on, will be slow.'

        # Coillision partners and ortho- / para- ratios.
        # By default, assume H2. If opr_cp is set, rescale 'dens' to
        # ortho-H2 and para-H2.

        ratefile = rates.ratefile(moldatfile)

        if 'H2' in ratefile.partners:
            H2_rates = True
        else:
            H2_rates = False

        if 'oH2' in ratefile.partners:
            oH2_rates = True
        else:
            oH2_rates = False

        if 'pH2' in ratefile.partners:
            pH2_rates = True
        else:
            pH2_rates = False

        if not (opr_cp > 0. or opr_cp is None):
            raise ValueError("opr_cp must be positive or None.")
        if (not H2_rates and not oH2_rates and not pH2_rates):
            raise ValueError("No appropriate collisonal rates found.")

        if (opr_cp is None and not H2_rates):
            if not self.lte_only:
                print 'Warning: No H2 collider density rates found.'
                print '\t Assuming oH2 / pH2 ratio of 3.'
            self.opr_cp = 3.
        elif (type(opr_cp) is float and not oH2_rates and not pH2_rates):
            if not self.lte_only:
                print 'Warning: No ortho-H2 or para-H2 rates found.'
                print '\t Assuming total H2 rates.'
            self.opr_co = None
        else:
            self.opr_cp = opr_cp

        # Include a depletion factor for the molecule.
        # This allows for isotopologues to be included, for example.

        if type(depletion) is float:
            self.rescale_abund = depletion
        else:
            raise TypeError("depletion must be a float.")

        return

    def checkTypes(self, inval, default, types):
        if inval is None:
            if default in self.hdr.arrnames:
                return default
            elif None in types:
                return inval
            else:
                raise TypeError("No %s in %s." % (default, self.hdr.fn))
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
