import os
import numpy as np
import scipy.constants as sc
from astropy.io import fits


def averageModels(model):
    for t in model.thetas:
        for p in model.phis:
            for j in model.transitions:
                toaverage = np.array([fits.getdata('%d_%.3f_%.3f_%d.fits' % (m, t, p, j), 0) 
                                      for m in range(model.nmodels)])
                averaged = np.average(toaverage, axis=0)
                hdulist = fits.open('0_%.3f_%.3f_%d.fits' % (t, p, j))
                hdulist[0].data = averaged
                filename = fileout + '_%.3f_%.3f_%d.fits' % (t, p, j)
                hdulist.writeto(filename)
                fits.setval(filename, 'NMODELS', value='%d' % nmodels, comment='Number of models averaged over.')
                os.system('mv %s %s' % (filename, model.directory))
    return


def getNoise(model):
    for t in model.thetas:
        for p in model.phis:
            for j in model.transitions:
                toaverage = np.array([fits.getdata('%d_%.3f_%.3f_%d.fits' % (m, t, p, j), 0)
                                      for m in range(model.nmodels)])
                gridnoise = np.std(toaverage, axis=0)
                hdulist = fits.open('0_%.3f_%.3f_%d.fits' % (t, p, j))
                hdulist[0].data = gridnoise
                hdulist.writeto(fileout+'_%.3f_%.3f_%d_noise.fits' % (t, p, j))
                os.system('mv %s_%.3f_%.3f_%d_noise.fits %s' % (fileout, t, p, j, model.directory))         
    return 



def combinePopfiles(model):
    popfiles = np.vstack([np.loadtxt('outputfile_%d.out' % m) 
                          for m in range(model.nmodels)]).T
    popfiles[:3] /= sc.au
    i = 0
    while (i < len(fileout) and fileout[i] != '.'):
        i += 1
    np.save('%s%s' % (model.directory, model.fileout), popfiles)
    return 

