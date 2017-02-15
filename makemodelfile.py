"""
Functions to write the model.c file used for LIME.
"""

import os


def makeFile(m, model):
    """Generate a model.c file to compile LIME with."""

    # Include the headers.
    tempfile = ['#include "lime.h"\n', '#include "math.h"\n',
                '#include "stdio.h"\n', '#include "stdlib.h"\n']
    tempfile.append('#include "%s"\n\n' % model.header.fn.split('/')[-1])

    # Include the imaging parameters.
    writeImageParameters(tempfile, m, model)

    # Include the interpolation functions.
    writeFindValue(tempfile, model)

    # Include the model functions.
    writeDensity(tempfile, model)
    writeTemperatures(tempfile, model)
    writeAbundance(tempfile, model)
    writeGastoDust(tempfile, model)
    writeRandomVelocities(tempfile, model)
    writeVelocityStructure(tempfile, model)

    # Save the output.
    with open('model_%d.c' % m, 'w') as tosave:
        for line in tempfile:
            tosave.write('%s' % line)

    return


def writeImageParameters(temp, m, model):
    """
    Include the imaging parameters. These should be the variables in the
    order of the manual. Note that some need to be included.
    """

    temp.append('void input(inputPars *par, image *img){\n\n')
    temp.append('\tpar->radius = %.5f*AU;\n' % model.radius)
    temp.append('\tpar->minScale = %.5f*AU;\n' % (max(model.minScale, 1e-3)))
    temp.append('\tpar->pIntensity = %.d;\n' % model.pIntensity)
    temp.append('\tpar->sinkPoints = %.d;\n' % model.sinkPoints)
    temp.append('\tpar->sampling = %d;\n' % model.sampling)
    temp.append('\tpar->tcmb = %.5f;\n' % model.tcmb)
    temp.append('\tpar->moldatfile[0] = "%s";\n' % model.moldatfile)
    temp.append('\tpar->dust = "%s";\n' % model.dust)
    temp.append('\tpar->lte_only = %d;\n' % model.lte_only)
    temp.append('\tpar->blend = %d;\n' % model.blend)
    temp.append('\tpar->antialias = %d;\n' % model.antialias)
    temp.append('\tpar->traceRayAlgorithm = %d;\n' % model.traceRayAlgorithm)
    temp.append('\tpar->nThreads = %d;\n' % model.nThreads)
    if model.outputgrid:
        temp.append('\tpar->gridOutFiles[3] = "%s";\n' % model.outputgrid)

    # Include the weighting factors included in LIME 1.6.

    abundance_weighting(temp, model)
    dust_weighting(temp, model)
    temp.append('\n')

    # For each permutation, write an image block.

    for i, inc in enumerate(model.incl):
        for p, pa in enumerate(model.posang):
            for a, azi in enumerate(model.azimuth):
                for t, trans in enumerate(model.transitions):
                    nimg = t + (a + p + i) * model.ntra
                    nimg += (p + i) * model.nazi + i * model.npos
                    writeImageBlock(temp, nimg, m, inc, pa, azi, trans, model)
    temp.append('}\n\n\n')
    return


def dust_weighting(temp, model):
    """
    Include the dustWeights. Parameters according to includeDust.
    """
    for i, weight in enumerate(model.dustWeights):
        temp.append('\tpar->dustWeights[{}] = {};\n'.format(i, weight))
    return


def abundance_weighting(temp, model):
    """
    Include the collisional partner IDs and their associated weights. Given
    that we work with an input of n(H2), the weights for all should be 1.0.
    """
    for i, cId in enumerate(model.collpartIds):
        temp.append('\tpar->collPartIds[{}] = {};\n'.format(i, cId))
    for i, weight in enumerate(model.nMolWeights):
        temp.append('\tpar->nMolWeights[{}] = {};\n'.format(i, weight))
    return


def writeImageBlock(temp, nimg, m, inc, pa, azi, trans, model):
    """
    Write an image block with the given parameters. Use the filename:

        nmodel_inc_pa_azi_trans.fits

    to allow for easier parsing of the different models. Parameters are in the
    manual order.
    """
    filename = '%s_%.2f_%.2f_%.2f_%d' % (m, inc, pa, azi, trans)
    imgs = '\timg[%d].' % nimg
    temp.append(imgs+'pxls = %d;\n' % model.pxls)
    temp.append(imgs+'imgres = %.3f;\n' % model.imgres)
    temp.append(imgs+'distance = %.3f*PC;\n' % model.distance)
    temp.append(imgs+'unit = %d;\n' % model.unit)
    temp.append(imgs+'filename = "%s.fits";\n' % filename)
    temp.append(imgs+'source_vel = %.3f;\n' % model.source_vel)
    temp.append(imgs+'nchan = %d;\n' % (model.nchan * model.oversample))
    temp.append(imgs+'velres = %.3f;\n' % (model.velres / model.oversample))
    temp.append(imgs+'trans = %d;\n' % trans)
    temp.append(imgs+'incl = %.3f;\n' % inc)
    temp.append(imgs+'posang = %.3f;\n' % (pa-1.57))
    temp.append(imgs+'azimuth = %.3f;\n\n' % azi)
    return


def writeCoords(temp, model):
    """
    Define the model coordinates for each function.
    """
    if not (model.coordsys is 'cylindrical' and model.ndim is 2):
        raise NotImplementedError
    if model.coordsys is 'cylindrical':
        if model.ndim is 2:
            temp.append('\tdouble c1 = sqrt(x*x + y*y) / AU;\n')
            temp.append('\tdouble c2 = fabs(z) / AU;\n')
            temp.append('\tdouble c3 = -1.;\n')
    return


def writeFindValue(temp, model):
    """
    Include the interpolation functions. TODO: Include 1D models.
    """
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
    """
    Molecular abundances. If 'abund' is a string, assume an array name, but
    if it is a float, assume a homogeneous rescaling. Also include a depletion
    factor through 'depletion'.
    """
    temp.append('void abundance(double x, double y, double z,')
    temp.append(' double *abundance) {\n\n')
    writeCoords(temp, model)
    temp.append('\tabundance[0] = ')
    if type(model.abund) is str:
        temp.append('findvalue(c1, c2, c3, %s);\n' % model.abund)
    elif type(model.abund) is float:
        temp.append('%.3e;\n' % float(model.abund))
    else:
        raise TypeError()
    temp.append('\tabundance[0] *= %.3e;\n' % model.depletion)
    temp.append('\tif (abundance[0] < 0.){\n\t\tabundance[0] = 0.;\n\t}\n')
    temp.append('\n}\n\n\n')
    return


def writeDensity(temp, model):
    """
    Write the molecular densities. Currently assumes 'dens' is H2 density.
    If an opr value is given, we include a second density array and have
    density[0] and density[1] as n(oH2) and n(pH2) respectively.
    """
    temp.append('void density(double x, double y,')
    temp.append('double z, double *density) {\n\n')
    writeCoords(temp, model)
    for i, rescale in enumerate(model.opr):
        temp.append('\tdensity[%d] = %.2f *' % (i, rescale))
        temp.append(' findvalue(c1, c2, c3, %s);\n' % model.dens)
        temp.append('\tif (density[%d] < 1e-30)' % i)
        temp.append(' {\n\t\tdensity[%d] = 1e-30;\n\t}\n\n' % i)
    temp.append('}\n\n\n')
    return


def writeRandomVelocities(temp, model):
    """
    Include the mirco-turbulence parameter. This can be either in m/s,
    model.dopplertype = 'absolute', or a fraction of the local soundspeed,
    model.dopplertyle = 'mach'.
    """
    temp.append('void doppler(double x, double y,')
    temp.append('double z, double *doppler) {\n\n')
    writeCoords(temp, model)
    if type(model.doppler) is str:
        temp.append('\t*doppler = findvalue(c1, c2,')
        temp.append('c3, %s);\n' % model.doppler)
    else:
        temp.append('\t*doppler = %.3f;\n' % model.doppler)
    writeDopplerMach(temp, model)
    temp.append('\n}\n\n\n')
    return


def writeDopplerMach(temp, model):
    """
    Write the doppler function for Mach values.
    """
    if model.dopplertype != 'mach':
        return
    temp.append('\tdouble val[2];\n')
    temp.append('\ttemperature(x, y, z, &val[2]);\n')
    temp.append('\t*doppler *= sqrt(KBOLTZ * val[0] / 2.34 / AMU);\n')
    return


def writeGastoDust(temp, model, ming2d=1.):
    """
    Write the gas-to-dust ratios. This is relative to n(H2), which is the sum
    of all density[i] values.
    """
    if model.g2d is None:
        return
    temp.append('void gasIIdust(double x, double y,')
    temp.append(' double z, double *gtd) {\n\n')
    if type(model.g2d) is float:
        temp.append('\t*gtd = %.1f;\n\n' % model.g2d)
    elif type(model.g2d) is str:
        writeCoords(temp, model)
        temp.append('\t*gtd = findvalue(c1, c2, c3, %s);\n\n' % model.g2d)
    temp.append('\tif (*gtd < %.1f) ' % ming2d)
    temp.append('{\n\t\t*gtd = %.1f;\n\t}\n}\n\n\n' % ming2d)
    return


def writeTemperatures(temp, model):
    """
    Write the gas and dust temperatures.
    """
    temp.append('void temperature(double x, double y, double z,')
    temp.append('double *temperature) {\n\n')
    writeCoords(temp, model)
    writeGasTemperature(temp, model)
    writeDustTemperature(temp, model)
    temp.append('}\n\n\n')
    return


def writeGasTemperature(temp, model):
    """Include the gas temperature."""
    temp.append('\ttemperature[0] = ')
    temp.append('findvalue(c1, c2, c3, %s);\n' % model.temp)
    if model.rescaletemp:
        temp.append('\ttemperature[0] *= %.3f;\n' & model.rescaletemp)
    temp.append('\tif (temperature[0] < 2.73) ')
    temp.append('{\n\t\ttemperature[0] = 2.73;\n\t}\n\n')
    return


def writeDustTemperature(temp, model):
    """Include the dust temperature if appropriate."""
    if type(model.dtemp) is float:
        temp.append('\ttemperature[1] = ')
        temp.append('%.3f * temperature[0];\n' % model.dtemp)
        temp.append('\tif (temperature[1] < 0.)')
        temp.append('{temperature[1] = 0.;}\n')
    elif model.dtemp is not None:
        temp.append('\ttemperature[1] = ')
        temp.append('findvalue(c1, c2, c3, %s);\n' % model.dtemp)
        temp.append('\tif (temperature[1] < 0.0) ')
        temp.append('{\n\t\ttemperature[1] = 0.0;\n\t}\n')
    return


def writeVelocityStructure(temp, model):
    """Include Keplerian rotation."""

    # TODO: Include arrays of velocities.
    # TODO: Is the z*z component required?

    temp.append('void velocity(double x, double y,')
    temp.append('double z, double *velocity) {\n\n')
    temp.append('\tif (sqrt(x*x + y*y + z*z) == 0.0){\n')
    temp.append('\t\tvelocity[0] = 0.0;\n')
    temp.append('\t\tvelocity[1] = 0.0;\n')
    temp.append('\t\tvelocity[2] = 0.0;\n')
    temp.append('\t\t return;\n')
    temp.append('}\n\n')
    temp.append('\tvelocity[0] = sqrt(6.67e-11 * ')
    temp.append('%.3f * 2e30 / ' % model.stellarmass)
    temp.append('sqrt(x*x + y*y + z*z));\n')
    temp.append('\tvelocity[0] *= sin(atan2(y,x));\n')
    temp.append('\tvelocity[1] = sqrt(6.67e-11 * ')
    temp.append('%.3f * 2e30 / ' % model.stellarmass)
    temp.append('sqrt(x*x + y*y + z*z));\n')
    temp.append('\tvelocity[1] *= cos(atan2(y,x));\n')
    temp.append('\tvelocity[2] = 0.0;\n')
    temp.append('\n}\n\n\n')
    return
