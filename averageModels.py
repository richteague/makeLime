import os
import numpy as np
from astropy.io import fits

# Average all the models from LIME. 
# If asked to, create a .fits file of the standard deviation
# of each pixel so that one could understand MCMC noise from
# the gridding. If there's only one model, just rename it.

averageModels(nmodels, thetas, phis, transitions,
              returnnoise=True, directory='../'):

    if nmodels > 1:
        for t in thetas:
            for p in phis:
                for j in transitions:
                
                    # Calculate the averages and move the file to the folder above.
                    
                    toaverage = np.array([fits.getdata('%d_%.3f_%.3f_%d.fits' % (m, t, p, j), 0)])
                    averaged = np.average(toaverage, axis=0)
                    hdulist = fits.open('0_%.3f_%.3f_%d.fits' % (t, p, j))
                    hdulist[0].data = averaged
                    hdulist.writeto(fileout+'_%.3f_%.3f_%d.fits' % (t, p, j))
                    os.system('mv %s_%.3f_%.3f_%d.fits %s' % (fileout, i, t, j, directory))
                    
                    # If appropriate, calculate the grid noise.
                    
                    if returnnoise:
                        gridnoise = np.std(toaverage, axis=0)
                        hdulist = fits.open('0_%.3f_%.3f_%d.fits' % (t, p, j))
                        hdulist[0].data = gridnoise
                        hdulist.writeto(fileout+'_%.3f_%.3f_%d_noise.fits' % (t, p, j))
                        os.system('mv %s_%.3f_%.3f_%d_noise.fits %s' % (fileout, t, p, j, directory))
                        
                        print 'For inclination: %.2f,' % t
                        print 'position angle: %.2f,' % p
                        print 'and transition: %d,' % j
                        print 'we have a grid noise of: %.2f [units].' % (gridnoise.mean())
        
        else:
            for t in thetas:
                for p in phis:
                    for j in transitions:
                        os.system('mv 0_%.3f_%.3f_%d.fits %s' % (t, p, j, directory))
    return

