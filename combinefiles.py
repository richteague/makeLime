import os
import numpy as np
from astropy.io import fits
import scipy.constants as sc


# Functions to combine the output files from the makeLIME run.

def fmt_fn(prefix, i, p, a, t):
    """Returns the filename"""
    return '{}_{:.2f}_{:.2f}_{:.2f}_{}.fits'.format(prefix, i, p, a, t)


def resampleVelocities(model):
    """Resample the velocities onto the requested velocity axis."""
    if model.oversample == 1:
        return
    files = [fn for fn in os.listdir('./') if fn.endswith('.fits')]
    for fn in files:
        hdu = fits.open(fn)
        data = hdu[0].data
        data = [np.average(data[i:i+model.oversample], axis=0)
                for i in range(data.shape[0])[::model.oversample]]
        hdu[0].data = np.array(data)
        hdu.writeto(fn, clobber=True)
    return


def moveModels(model, prefix='0_', suffix=''):
    """Move the finished models to appropriate end directory."""
    files = [fn for fn in os.listdir('./') if fn.endswith('%s.fits' % suffix)]
    for fn in [fn for fn in files if fn.beginswith(prefix)]:
        os.system('mv %s %s%s' % (fn, model.directory, model.name+fn[1:]))
    return


def averageModels(model):
    """Average over all the models and save to 0_*.fits."""
    if model.nmodels == 1:
        return
    for i in model.incl:
        for p in model.posang:
            for a in model.azimuth:
                for t in model.transitions:
                    fn = fmt_fn(0, i, p, a, t)
                    hdu = fits.open(fn)
                    avg = [fits.getdata(fmt_fn(m, i, p, a, t))
                           for m in range(model.nmodels)]
                    getNoise(avg, i, p, a, t, model)
                    hdu[0].data = np.average(avg, axis=0)
                    hdu.writeto(fn, clobber=True)
    return


def writeFitsHeader(filename, model, inc, pa, azi):
    """Include model data in the final .fits file header."""
    data, header = fits.getdata(filename, header=True)
    header['DISTANCE'] = model.distance, 'Distance in parsec.'
    header['CHEMMOD'] = model.hdr.fn.split('/')[-1], 'Chemical model used.'
    header['INC'] = inc, 'Inclianation in radians.'
    header['PA'] = pa, 'Position angle in radians.'
    header['AZI'] = azi, 'Azimuthal angle in radians.'
    header['NMODELS'] = model.nmodels, 'Number of models averaged.'
    if model.opr_cp is not None:
        header['H2_OPR'] = model.opr_cp, 'H2 ortho to para ratio.'

    # Include the change for over-sampled velocity axis.

    header['NAXIS3'] = data.shape[0]
    header['CDELT3'] = model.velres
    header['CRPIX3'] = (data.shape[0] + 1.) / 2.
    fits.writeto(filename, data, header, clobber=True)
    return


def getNoise(avgmodels, i, p, a, t, model):
    """Estimate the MCMC noise from the model ensemble."""
    if (model.nmodel == 1 or not model.returnnoise):
        return
    noise = np.std(avgmodels, axis=0)
    hdu = fits.open(fmt_fn(1, i, p, a, t))
    hdu.data = noise
    hdu.writeto(fmt_fn(0, i, p, a, t)[:-5]+'_noise.fits')
    return


def combinePopfiles(model):
    """Combine the population files from the model ensemble."""
    if not model.outputfile:
        return
    outputfile = np.vstack([np.loadtxt('outputfile_%d.out' % m)
                            for m in range(model.nmodels)]).T
    outputfile[:3] /= sc.au
    np.save('%s%s_popfile.out' % (model.directory, model.name), outputfile)
    return


def combineBinPopFiles(model):
    """Combine the binary population files from the model ensemble."""
    if model.binoutputfile:
        raise NotImplementedError()
    return


def combineGridfiles(model):
    """Move the grid files."""
    if model.gridfile:
        os.system('mv gridfile*.out %s' % model.directory)
    return
