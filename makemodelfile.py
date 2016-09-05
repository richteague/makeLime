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
    temp.append('\tpar->blend = %d;\n' % model.blend)

    # Ortho / para ratio of the H2. If an ortho / para ratio is specified, 
    # takes density[0] to be n(oH2) and density[1] to be n(pH2). We set 
    # MolWeights = 1.0 for both densities as wew want to take into account 
    # both collision partners equally:
    #
    #       n(mol) = x(mol) * (w_0 * n(oH2) + w_1 * n(pH2)),
    #
    # thus so long as n(oH2) and n(pH2) are correctly specified in the density
    # function, w_0 = w_1 = 1.0. The same goes for the dust. TODO: Include 
    # other collision partners.

    if model.opr_cp is None:
        temp.append('\tpar->collPartIds[0] = CP_H2;\n')
        temp.append('\tpar->dustWeights[0] = 1.0;\n')
        temp.append('\tpar->nMolWeights[0] = 1.0;\n')
    else:
        temp.append('\tpar->collPartIds[0] = CP_o_H2;\n')
        temp.append('\tpar->collPartIds[1] = CP_p_H2;\n')
        for i in range(2):
            temp.append('\tpar->nMolWeights[%d] = 1.0;\n' % i)
            temp.append('\tpar->dustWeights[%d] = 1.0;\n' % i)

    # Include optional output files.

    if model.outputfile:
        temp.append('\tpar->outputfile = "outputfile_%d.out";\n' % m)
    if model.binoutputfile:
        temp.append('\tpar->binoutputfile = "binoutputfile_%d.out";\n' % m)
    if model.gridfile:
        temp.append('\tpar->gridfile = "gridfile_s%s_%d.out";\n' % (model.name, m))
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

def writeImageBlock(temp, nimg, m, theta, phi, trans, model):

    # temp          : str, lines of input file.
    # nimg          : int, image nuumber.
    # m             : int, model number.
    # theta         : float, inclination.
    # phi           : float, position angle.
    # transition    : int, transition number.
    # model         : LIMEclass, parameters.
    #
    # Write an image block with the given parameters. Use the filename:
    # nmodel_inclination_positionangle_transition.fits

    filename = '%s_%.3f_%.3f_%d' % (m, theta, phi, trans)
    temp.append('\timg[%d].nchan = %d;\n' % (nimg, model.nchan))
    temp.append('\timg[%d].velres = %d;\n' % (nimg, model.velres))
    temp.append('\timg[%d].trans = %d;\n' % (nimg, trans))
    temp.append('\timg[%d].pxls = %d;\n' % (nimg, model.pxls))
    temp.append('\timg[%d].imgres = %.3f;\n' % (nimg, model.imgres))
    temp.append('\timg[%d].incl = %.3f;\n' % (nimg, theta))
    temp.append('\timg[%d].posang = %.3f;\n' % (nimg, phi-np.radians(90.)))
    temp.append('\timg[%d].azimuth = %.3f;\n' % (nimg, model.azimuth))
    temp.append('\timg[%d].distance = %.3f*PC;\n' % (nimg, model.distance))
    temp.append('\timg[%d].unit = %d;\n' % (nimg, model.unit))
    temp.append('\timg[%d].filename = "%s.fits";\n' % (nimg, filename))
    temp.append('\n')
    return

def writeCoords(temp, model):

    # temp          : str, lines of input file.
    # model         : LIMEclass, parameters.
    #
    # Define the model coordinates for each function.

    if not (model.coordsys is 'cylindrical' and model.ndim is 2):
        raise NotImplementedError
    if model.coordsys is 'cylindrical':
        if model.ndim is 2:
            temp.append('\tdouble c1 = sqrt(x*x + y*y) / AU;\n')
            temp.append('\tdouble c2 = fabs(z) / AU;\n')
            temp.append('\tdouble c3 = -1.;\n')
    return

def writeFindValue(temp, model):

    # temp          : str, lines of input file.
    # model         : LIMEclass, parameters.
    #
    # Write the interpolation functions. Current only bilinear interpolation 
    # with cylindrical coords.

    if not (model.coordsys is 'cylindrical' and model.ndim is 2):
        raise NotImplementedError
    path = os.path.dirname(__file__)
    path += '/InterpolationRoutines/'
    with open(path+'%dD_%s.c' % (model.ndim, model.coordsys)) as f:
        lines = f.readlines()
    for line in lines:
        line = line.replace('NCELLS', '%d' % model.ncells)
        temp.append(line)
    temp.append('\n\n')
    return

def writeAbundance(temp, model):

    # temp          : str, lines of input file.
    # model         : LIMEclass, parameters.
    #
    # Include the molecular abundance. Depletion rescales the 'abund' value
    # given in the header file.

    s = 'void abundance(double x, double y, double z,'
    s += ' double *abundance) {\n\n'
    temp.append(s)
    writeCoords(temp, model)

    if type(model.abund) is str:
        temp.append('\tabundance[0] = findvalue(c1, c2, c3, %s);\n' %
                    model.abund)
    elif type(model.abund) is float:
        temp.append('\tabundance[0] = %.3e;\n' % float(model.abund))
    else:
        raise TypeError()

    # Include a rescaling factor for, e.g., ortho / para considerations
    # or rescaling to an isotopologue. Note: this is 'depletion' in limeclass.

    if model.rescale_abund != 1.0:
        temp.append('\tabundance[0] *= %.3e;\n' % model.rescale_abund)

    temp.append('\tif (abundance[0] < 0.){\n\t\tabundance[0] = 0.;\n\t}\n')
    temp.append('\n}\n\n\n')
    return

def writeDensity(temp, model):

    # temp          : str, lines of input file.
    # model         : LIMEclass, parameters.
    #
    # Write the molecular abundances. Currently assumes 'dens' is H2 density.
    # If an opr_cp value is given, we include density[1] and have density[0] 
    # and density[1] as n(oH2) and n(pH2) respectively.

    temp.append('void density(double x, double y, double z, double *density) {\n\n')
    writeCoords(temp, model)

    # Rescale the densities as appropriate: x[0] = x(oH2), x[1] = x(nH2).
    
    x = np.ones(2)
    if model.opr_cp is not None:
        x[0] = model.opr_cp / (1. + model.opr_cp)
        x[1] = 1. / (1. + model.opr_cp)
    
    for i in range(2):
        s = '\tdensity[%d] = %.2f *' % (i, x[i])
        s += ' findvalue(c1, c2, c3, %s);\n' % model.dens
        temp.append(s)
        s = '\tif (density[%d] < 1e-30)' % i
        s += ' {\n\t\tdensity[%d] = 1e-30;\n\t}\n\n' % i
        temp.append(s)
        if model.opr_cp is None:
            break
    temp.append('}\n\n\n')
    return

def writeDopplerBroadening(temp, model):

    # temp          : str, lines of input file.
    # model         : LIMEclass, parameters.
    #
    # Include the mirco-turbulence parameter. This can be either in m/s, 
    # model.dopplertype = 'absolute', or a fraction of the local soundspeed, 
    # model.dopplertyle = 'mach'. This does not include thermal broadening.

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

def writeGastoDust(temp, model, ming2d=1.):

    # temp          : str, lines of input file.
    # model         : LIMEclass, parameters.
    # ming2d        : float, minimum value for gas-to-dust.
    #
    # Write the gas-to-dust ratios. This is relative to n(H2), 
    # which is the sum of all density[i] values.

    if model.g2d is None:
        return
    else:
        s = 'void gasIIdust(double x, double y, double z, double *gtd) {\n\n'
        temp.append(s)

    if type(model.g2d) is float:
        temp.append('\t*gtd = %.1f;\n\n' % model.g2d)
    elif type(model.g2d) is str:
        writeCoords(temp, model)
        temp.append('\t*gtd = findvalue(c1, c2, c3, %s);\n\n' % model.g2d)

    s = '\tif (*gtd < %.1f) {\n\t\t*gtd = %.1f;\n\t}\n' % (ming2d, ming2d)
    temp.append(s)
    temp.append('}\n\n\n')
    return

def writeTemperatures(temp, model):

    # temp          : str, lines of input file.
    # model         : LIMEclass, parameters.
    #
    # Write the gas and dust temperatures. If dtemp is a float, then rescale
    # the gas temperature. If it is a string, use that array name as input.

    s = 'void temperature(double x, double y, double z,'
    s += 'double *temperature) {\n\n'
    temp.append(s)
    writeCoords(temp, model)

    if type(model.temp) is str:
        s = '\ttemperature[0] = findvalue(c1, c2, c3, %s);\n' % model.temp
        temp.append(s)
        s = '\tif (temperature[0] < 2.73) {\n\t\ttemperature[0] = 2.73;\n\t}\n\n'
        temp.append(s)
    else:
        raise TypeError()

    if type(model.dtemp) is float:
        s = '\ttemperature[1] = %.3f * temperature[0];\n' % model.dtemp
        temp.append(s)
        s = '\tif (temperature[1] < 2.73){temperature[1] = 2.73;}\n'
        temp.append(s)
    elif type(model.dtemp) is str:
        s = '\ttemperature[1] = findvalue(c1, c2, c3, %s);\n' % model.dtemp
        temp.append(s)
        s = '\tif (temperature[1] < 2.73) {\n\t\ttemperature[1] = 2.73;\n\t}\n'
        temp.append(s)
    elif model.dtemp is not None:
        raise TypeError()

    temp.append('}\n\n\n')
    return

def writeVelocityStructure(temp, model):

    # temp          : str, lines of input file.
    # model         : LIMEclass, parameters.
    #
    # Write the velocity structure section. Only Keplerian rotation is
    # possible. TODO: Include arrays of velocities.

    s = 'void velocity(double x, double y, double z, double *velocity) {\n\n'
    temp.append(s)

    if model.stellarmass is not None:
        s = '\tvelocity[0] = sqrt(6.67e-11 * %.3f ' % model.stellarmass
        s += '* 2e30 / sqrt(x*x + y*y + z*z));\n' 
        temp.append(s)
        temp.append('\tvelocity[0] *= sin(atan2(y,x));\n')
        s = '\tvelocity[1] = sqrt(6.67e-11 * %.3f ' % model.stellarmass
        s += '* 2e30 / sqrt(x*x + y*y + z*z));\n'
        temp.append(s)
        temp.append('\tvelocity[1] *= cos(atan2(y,x));\n')
        temp.append('\tvelocity[2] = 0.0;\n')
    else:
        raise NotImplementedError('Must currently use Keplerian rotation.')
    temp.append('\n}\n\n\n')
    return

def averageModels(model):

    # model         : LIMEclass, parameters.
    #
    # Average over all the models. Move the combined files to the correct
    # directory.

    if model.nmodels == 1:
        for t in model.thetas:
            for p in model.phis:
                for j in model.transitions:
                    fn = '%s_%.3f_%.3f_%d.fits' % (model.name, t, p, j)
                    os.system('mv 0_%.3f_%.3f_%d.fits %s' % (t, p, j, fn))
                    writeFitsHeader(fn, model, t, p)
                    os.system('mv %s %s' % (fn, model.directory))
    else:
        for t in model.thetas:
            for p in model.phis:
                for j in model.transitions:
                    avg = [fits.getdata('%d_%.3f_%.3f_%d.fits' % (m, t, p, j))
                           for m in range(model.nmodels)]
                    averaged = np.average(avg, axis=0)
                    hdulist = fits.open('0_%.3f_%.3f_%d.fits' % (t, p, j))
                    hdulist[0].data = averaged
                    filename = model.name + '_%.3f_%.3f_%d.fits' % (t, p, j)
                    hdulist.writeto(filename)
                    writeFitsHeader(filename, model, t, p)
                    os.system('mv %s %s' % (filename, model.directory))
    return

def writeFitsHeader(filename, model, theta, phi):

    # filename      : str, filename to edit the properties of.
    # model         : LIMEclass, parameters.
    # theta         : float, inclination angle.
    # phi           : float, position angle.
    #
    # Include model data in the final .fits file header.

    data, header = fits.getdata(filename, header=True)
    header['DISTANCE'] = model.distance, 'Distance in parsec.'
    header['CHEMMOD'] = model.hdr.fn, 'Chemical model used.'
    header['INC'] = theta, 'Inclianation in radians.'
    header['PA'] = phi, 'Position angle in radians.'
    header['NMODELS'] = model.nmodels, 'Number of models averaged.'
    if model.opr_cp is not None:
        header['H2_OPR'] = model.opr_cp, 'H2 ortho to para ratio.'
    fits.writeto(filename, data, header, clobber=True)
    print 'Appended header with model properties.'
    return

def getNoise(model):

    # model         : LIMEclass, parameters.
    #
    # Estimate the MCMC noise from the model ensemble.

    if (model.returnnoise and model.nmodels > 1):
        for t in model.thetas:
            for p in model.phis:
                for j in model.transitions:
                    avg = [fits.getdata('%d_%.3f_%.3f_%d.fits' % (m, t, p, j))
                           for m in range(model.nmodels)]
                    gridnoise = np.std(avg, axis=0)
                    hdulist = fits.open('0_%.3f_%.3f_%d.fits' % (t, p, j))
                    hdulist[0].data = gridnoise
                    filename = model.name + '_%.3f_%.3f_%d_noise.fits' % (t, p, j)
                    hdulist.writeto(filename)
                    os.system('mv %s %s' % (filename, model.directory))
    return

def combinePopfiles(model):

    # model         : LIMEclass, parameters.
    #
    # Combine the population files from the model ensemble.

    if model.outputfile:
        popfiles = np.vstack([np.loadtxt('outputfile_%d.out' % m)
                              for m in range(model.nmodels)]).T
        popfiles[:3] /= sc.au
        np.save('%s%s_popfile' % (model.directory, model.name), popfiles)
    return

def combineBinPopFiles(model):

    # model         : LIMEclass, parameters.
    #
    # Combine the binary population files from the model ensemble.

    if model.binoutputfile:
        raise NotImplementedError()
    return

def combineGridfiles(model):

    # model         : LIMEclass, parameters.
    #
    # Move the grid files. TODO: Can we combine them?

    if model.gridfile:
        os.system('mv gridfile*.out %s' % model.directory)
    return
