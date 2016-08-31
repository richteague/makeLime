import os
import numpy as np
from astropy.io import fits
import scipy.constants as sc

''' Functions to create the model.c file. '''

# Generate the model.c file to run with LIME.
# m             : int, model number.
# model         : LIMEclass, parameters.


def makeFile(m, model):

    # Include the headers.
    tempfile = ['#include "lime.h"\n', '#include "math.h"\n',
                '#include "stdio.h"\n', '#include "stdlib.h"\n']
    tempfile.append('#include "%s"\n\n' % model.hdr.fn.split('/')[-1])

    # Include the imaging parameters.
    writeImageParameters(tempfile, m, model)

    # Include the interpolation functions.
    writeFindValue(tempfile, model)

    # Include the model functions.
    writeDensity(tempfile, model)
    writeTemperatures(tempfile, model)
    writeAbundance(tempfile, model)
    writeGastoDust(tempfile, model)
    writeDopplerBroadening(tempfile, model)
    writeVelocityStructure(tempfile, model)

    # Save the output.
    with open('model_%d.c' % m, 'w') as tosave:
        for line in tempfile:
            tosave.write('%s' % line)

    return

# Write an parameter block with the given parameters.
# temp          : str, lines of input file.
# m             : int, model number.
# model         : LIMEclass, parameters.


def writeImageParameters(temp, m, model):

    temp.append('void input(inputPars *par, image *img){\n\n')

    # Include the radiative transfer parameters.

    temp.append('\tpar->radius = %.1f*AU;\n' % model.rout)
    temp.append('\tpar->minScale = %.2f*AU;\n' % model.rin)
    temp.append('\tpar->pIntensity = %.0f;\n' % model.pIntensity)
    temp.append('\tpar->sinkPoints = %.0f;\n' % model.sinkPoints)
    temp.append('\tpar->dust = "%s";\n' % model.dust)
    temp.append('\tpar->moldatfile[0] = "%s";\n' % model.moldatfile)
    temp.append('\tpar->antialias = %d;\n' % model.antialias)
    temp.append('\tpar->sampling = %d;\n' % model.sampling)
    temp.append('\tpar->lte_only = %d;\n' % model.lte_only)
    if model.outputfile is not None:
        temp.append('\tpar->outputfile = "outputfile_%d.out";\n' % m)
    if model.binoutputfile is not None:
        temp.append('\tpar->binoutputfile = "binoutputfile_%d.out";\n' % m)
    if model.gridfile is not None:
        temp.append('\tpar->gridfile = "gridfile_%d.out";\n' % m)
    temp.append('\n')

    # For each permutation, write an image block.

    for i, theta in enumerate(model.thetas):
        for p, phi in enumerate(model.phis):
            for t, trans in enumerate(model.transitions):
                nimg = i * len(model.transitions) * len(model.phis)
                nimg += p * len(model.transitions) + t
                nimg = int(nimg)
                writeImageBlock(temp, nimg, m, theta, phi, trans, model)
    temp.append('}\n\n\n')

    return

# Write an image block with the given parameters.
# temp          : str, lines of input file.
# nimg          : int, image nuumber.
# m             : int, model number.
# theta         : float, inclination.
# phi           : float, position angle. TODO: Change to PA.
# transition    : int, transition number.
# model         : LIMEclass, parameters.


def writeImageBlock(temp, nimg, m, theta, phi, trans, model):

    # Use the filename convention.

    filename = '%s_%.3f_%.3f_%d' % (m, theta, phi, trans)

    temp.append('\timg[%.0f].nchan = %.0f;\n' % (nimg, model.nchan))
    temp.append('\timg[%.0f].velres = %.0f;\n' % (nimg, model.velres))
    temp.append('\timg[%.0f].trans = %.0f;\n' % (nimg, trans))
    temp.append('\timg[%.0f].pxls = %d;\n' % (nimg, model.pxls))
    temp.append('\timg[%.0f].imgres = %.3f;\n' % (nimg, model.imgres))
    temp.append('\timg[%.0f].theta = %.3f;\n' % (nimg, theta))
    temp.append('\timg[%.0f].phi = %.3f;\n' % (nimg, phi))
    temp.append('\timg[%.0f].distance = %.1f*PC;\n' % (nimg, model.distance))
    temp.append('\timg[%.0f].unit = %d;\n' % (nimg, model.unit))
    temp.append('\timg[%.0f].filename = "%s.fits";\n' % (nimg, filename))
    temp.append('\n')
    return

# Write the coordinates for the start of each function.
# temp          : str, lines of input file.
# model         : LIMEclass, parameters.


def writeCoords(temp, model):
    if not (model.coordsys is 'cylindrical' and model.ndim is 2):
        raise NotImplementedError
    if model.coordsys is 'cylindrical':
        if model.ndim is 2:
            temp.append('\tdouble c1 = sqrt(x*x + y*y) / AU;\n')
            temp.append('\tdouble c2 = fabs(z) / AU;\n')
            temp.append('\tdouble c3 = -1.;\n')
    return

# Write the interpolation routine.
# temp          : str, lines of input file.
# model         : LIMEclass, parameters.


def writeFindValue(temp, model):

    # Include the appropriate interpolation routines.

    if not (model.coordsys is 'cylindrical' and model.ndim is 2):
        raise NotImplementedError

    # Find the correct path to the files.
    path = os.path.dirname(__file__)
    path += '/InterpolationRouttines/'
    with open(path+'%dD_%s.c' % (model.ndim, model.coordsys)) as f:
        lines = f.readlines()
    for line in lines:
        line = line.replace('NCELLS', '%d' % model.ncells)
        temp.append(line)
    temp.append('\n\n')
    return

# Write the molecular abundance section.
# temp          : str, lines of input file.
# model         : LIMEclass, parameters.


def writeAbundance(temp, model):
    temp.append('''void abundance(double x, double y, double z, double\
                   *abundance) {\n\n''')
    writeCoords(temp, model)

    if type(model.abund) is str:
        temp.append('\tabundance[0] = findvalue(c1, c2, c3, %s);\n' %
                    model.abund)
    elif type(model.abund) is float:
        temp.append('\tabundance[0] = %.3e;\n' % float(model.abund))
    else:
        raise TypeError()

    temp.append('\tif (abundance[0] < 0.){\n\t\tabundance[0] = 0.;\n\t}\n')
    temp.append('\n}\n\n\n')
    return

# Write the molecular abundance section.
# temp          : str, lines of input file.
# model         : LIMEclass, parameters.


def writeDensity(temp, model):
    temp.append('void density(double x, double y, double z, double *density) {\n\n')
    writeCoords(temp, model)

    if type(model.abund) is str:
        temp.append('\tdensity[0] = findvalue(c1, c2, c3, %s);\n' % model.dens) 
    else:
        raise TypeError()   

    temp.append('\tif (density[0] < 1e-30) {\n\t\tdensity[0] = 1e-30;\n\t}\n\n')
    temp.append('}\n\n\n')
    return 



# Write the molecular abundance section.
# temp          : str, lines of input file.
# model         : LIMEclass, parameters.

def writeDopplerBroadening(temp, model):    
    
    if model.doppler is None:
        model.doppler = 0.0
    
    temp.append('void doppler(double x, double y, double z, double *doppler) {\n\n')
    writeCoords(temp, model)
    
    if model.dopplertype == 'mach':
        temp.append('\tdouble val[2];\n')
        temp.append('\ttemperature(x, y, z, &val[2]);\n')
        
    if type(model.doppler) is str:
        temp.append('\t*doppler = findvalue(c1, c2, c3, %s);\n' % model.doppler) 
    elif type(model.doppler) is float:
        temp.append('\t*doppler = %.3f;\n' % model.doppler)

    if model.dopplertype == 'mach':
        temp.append('\t*doppler *= sqrt(KBOLTZ * val[0] / 2.34 / AMU);\n')
        
    temp.append('\n}\n\n\n')
    return



# Write the gas-to-dust ratio section.
# temp          : str, lines of input file.
# model         : LIMEclass, parameters.
# ming2d        : float, minimum value for gas-to-dust.

def writeGastoDust(temp, model, ming2d=1.):

    if model.g2d is None:
        return
    else: 
        temp.append('void gasIIdust(double x, double y, double z, double *gtd) {\n\n')

    if type(model.g2d) is float:
        temp.append('\t*gtd = %.1f;\n\n' % model.g2d)
    elif type(model.g2d) is str:
        writeCoords(temp, model)
        temp.append('\t*gtd = findvalue(c1, c2, c3, %s);\n\n' % model.g2d)

    temp.append('\tif (*gtd < %.1f) {\n\t\t*gtd = %.1f;\n\t}\n' % (ming2d, ming2d))
    temp.append('}\n\n\n')
    return



# Write the gas and dust temperatures section.
# temp          : str, lines of input file.
# model         : LIMEclass, parameters.

def writeTemperatures(temp, model):

    temp.append('void temperature(double x, double y, double z, double *temperature) {\n\n')
    writeCoords(temp, model)

    if type(model.temp) is str:
        temp.append('\ttemperature[0] = findvalue(c1, c2, c3, %s);\n' % model.temp)
        temp.append('\tif (temperature[0] < 2.73) {\n\t\ttemperature[0] = 2.73;\n\t}\n\n')
    else:
        raise TypeError()
        
    if type(model.dtemp) is float:
        temp.append('\ttemperature[1] = %.3f * temperature[0];\n' % model.dtemp)
        temp.append('\tif (temperature[1] < 2.73){temperature[1] = 2.73;}\n')
    elif type(model.dtemp) is str:
        temp.append('\ttemperature[1] = findvalue(c1, c2, c3, %s);\n' % model.dtemp)
        temp.append('\tif (temperature[1] < 2.73) {\n\t\ttemperature[1] = 2.73;\n\t}\n')
    elif model.dtemp is not None:
        raise TypeError()

    temp.append('}\n\n\n')
    return



# Write the velocity structure section.
# temp          : str, lines of input file.
# model         : LIMEclass, parameters.

def writeVelocityStructure(temp, model):

    temp.append('void velocity(double x, double y, double z, double *velocity) {\n\n')

    if model.stellarmass is not None:
        temp.append('\tvelocity[0] = sqrt(6.67e-11 * %.3f * 2e30 / sqrt(x*x + y*y + z*z));\n' % model.stellarmass)
        temp.append('\tvelocity[0] *= sin(atan2(y,x));\n')
        temp.append('\tvelocity[1] = sqrt(6.67e-11 * %.3f * 2e30 / sqrt(x*x + y*y + z*z));\n' % model.stellarmass)
        temp.append('\tvelocity[1] *= cos(atan2(y,x));\n')
        temp.append('\tvelocity[2] = 0.0;\n')

    else:
        raise NotImplementedError('Must currently use Keplerian rotation.')
        
    temp.append('\n}\n\n\n')

    return


    
# Average over the models run.
# model         : LIMEclass, parameters.
 
def averageModels(model):
    if model.nmodels == 1:
        for t in model.thetas:
            for p in model.phis:
                for j in model.transitions:
                    filename = '%s_%.3f_%.3f_%d.fits' % (model.name, t, p, j)
                    os.system('mv 0_%.3f_%.3f_%d.fits %s' % (t, p, j, filename))
                    writeFitsHeader(filename, model, t, p) 
                    os.system('mv %s %s' % (filename, model.directory))
    else:
        for t in model.thetas:
            for p in model.phis:
                for j in model.transitions:
                    toaverage = np.array([fits.getdata('%d_%.3f_%.3f_%d.fits' % (m, t, p, j), 0) 
                                          for m in range(model.nmodels)])
                    averaged = np.average(toaverage, axis=0)
                    hdulist = fits.open('0_%.3f_%.3f_%d.fits' % (t, p, j))
                    hdulist[0].data = averaged
                    filename = model.name + '_%.3f_%.3f_%d.fits' % (t, p, j)
                    hdulist.writeto(filename)
                    writeFitsHeader(filename, model, t, p)
                    os.system('mv %s %s' % (filename, model.directory))
    return



# Include model data in the final .fits file header.
# filename      : str, filename to edit the properties of.
# model         : LIMEclass, parameters.
# theta         : float, inclination angle.
# phi           : float, position angle TODO: change to real PA.

def writeFitsHeader(filename, model, theta, phi):
       
    data, header = fits.getdata(filename, header=True)
    header['DISTANCE'] = model.distance, 'Distance in parsec.'
    header['CHEMMOD'] = model.hdr.fn, 'Chemical model used.'
    header['INC'] = theta, 'Inclianation in radians.'
    header['PA'] = phi, 'Position angle in radians.'
    header['NMODELS'] = model.nmodels, 'Number of models averaged.' 
    fits.writeto(filename, data, header, clobber=True)    

    print 'Appended header with model properties.'
    return


# Estimate the MCMC noise from the model ensemble.
# model         : LIMEclass, parameters.

def getNoise(model):
    if (model.returnnoise and model.nmodels > 1):
        for t in model.thetas:
            for p in model.phis:
                for j in model.transitions:
                    toaverage = np.array([fits.getdata('%d_%.3f_%.3f_%d.fits' % (m, t, p, j), 0)
                                          for m in range(model.nmodels)])
                    gridnoise = np.std(toaverage, axis=0)
                    hdulist = fits.open('0_%.3f_%.3f_%d.fits' % (t, p, j))
                    hdulist[0].data = gridnoise
                    filename = model.name + '_%.3f_%.3f_%d_noise.fits' % (t, p, j)
                    hdulist.writeto(filename)
                    os.system('mv %s %s' % (filename, model.directory))         
    return 



# Combine the population files from the model ensemble.
# model         : LIMEclass, parameters.

def combinePopfiles(model):
    if model.outputfile:
        popfiles = np.vstack([np.loadtxt('outputfile_%d.out' % m) 
                              for m in range(model.nmodels)]).T
        popfiles[:3] /= sc.au
        np.save('%s%s_popfile' % (model.directory, model.name), popfiles)
    return 


# Combine the binary population files from the model ensemble.
# model         : LIMEclass, parameters.
    
def combineBinPopFiles(model):
    if model.binoutputfile:
        raise NotImplementedError()
    return
    
