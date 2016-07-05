import os
import numpy as np
import scipy.constants as sc
from astropy.io import fits


def averageModels(nmodels, thetas, phis, transitions, fileout,
                  returnnoise=False, directory='../'):

    # Average all the models.

    for t in thetas:
        for p in phis:
            for j in transitions:
                    
                toaverage = np.array([fits.getdata('%d_%.3f_%.3f_%d.fits' % (m, t, p, j), 0) for m in range(nmodels)])
                averaged = np.average(toaverage, axis=0)
                hdulist = fits.open('0_%.3f_%.3f_%d.fits' % (t, p, j))
                hdulist[0].data = averaged
                hdulist.writeto(fileout+'_%.3f_%.3f_%d.fits' % (t, p, j))
                os.system('mv %s_%.3f_%.3f_%d.fits %s' % (fileout, t, p, j, directory))
                
                if returnnoise:
                    getNoise(nmodels, thetas, phis, transitions, fileout, directory)

    return


def getNoise(nmodels, thetas, phis, transitions, fileout, directory='./'):

    # Calculate the dispersion of each set of models.

    for t in thetas:
        for p in phis:
            for j in transitions:
                toaverage = np.array([fits.getdata('%d_%.3f_%.3f_%d.fits' % (m, t, p, j), 0) for m in range(nmodels)])
                gridnoise = np.std(toaverage, axis=0)
                hdulist = fits.open('0_%.3f_%.3f_%d.fits' % (t, p, j))
                hdulist[0].data = gridnoise
                hdulist.writeto(fileout+'_%.3f_%.3f_%d_noise.fits' % (t, p, j))
                os.system('mv %s_%.3f_%.3f_%d_noise.fits %s' % (fileout, t, p, j, directory))         
    
    return 



def combinePopfiles(nmodels, fileout, directory='../'):
        
    # Combine all the popfiles from the averaged models.

    popfiles = np.vstack([np.loadtxt('popfile_%d.out' % m) for m in range(nmodels)]).T
    popfiles[:3] /= sc.au
    np.save('%s%s_popfile' % (directory, fileout), popfiles)

    return 
