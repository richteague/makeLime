import os
import numpy as np
import scipy.constants as sc
from astropy.io import fits

# Average all the models from LIME. 
# If asked to, create a .fits file of the standard deviation
# of each pixel so that one could understand MCMC noise from
# the gridding. If there's only one model, just rename it.

def averageModels(nmodels, thetas, phis, transitions, fileout,
                  returnnoise=False, directory='../'):

        for t in thetas:
            for p in phis:
                for j in transitions:
                
                    # Calculate the averages and move the file to the folder above.
                    
                    toaverage = np.array([fits.getdata('%d_%.3f_%.3f_%d.fits' % (m, t, p, j), 0) for m in range(nmodels)])
                    averaged = np.average(toaverage, axis=0)
                    hdulist = fits.open('0_%.3f_%.3f_%d.fits' % (t, p, j))
                    hdulist[0].data = averaged
                    hdulist.writeto(fileout+'_%.3f_%.3f_%d.fits' % (t, p, j))
                    os.system('mv %s_%.3f_%.3f_%d.fits %s' % (fileout, t, p, j, directory))
                    
                    # If appropriate, calculate the grid noise.
                    
                    if returnnoise:
                        gridnoise = np.std(toaverage, axis=0)
                        hdulist = fits.open('0_%.3f_%.3f_%d.fits' % (t, p, j))
                        hdulist[0].data = gridnoise
                        hdulist.writeto(fileout+'_%.3f_%.3f_%d_noise.fits' % (t, p, j))
                        os.system('mv %s_%.3f_%.3f_%d_noise.fits %s' % (fileout, t, p, j, directory))         

                        # Update this text!
 
                        print 'For inclination: %.2f,' % t
                        print 'position angle: %.2f,' % p
                        print 'and transition: %d,' % j
                        print 'we have a mean fractional grid noise of: %.5e [units].' % (np.nanmean(gridnoise/averaged))
                        print 'Through averaging we reduce this to: %.5e [units].' % (np.nanmean(gridnoise/averaged)/np.sqrt(nmodels))
    
    return


def combinePopfiles(nmodels, fileout, directory='../'):
        
    # Combine all the popfiles from the averaged models.

    popfiles = np.vstack([np.loadtxt('popfile_%d.out' % m) for m in range(nmodels)]).T
    popfiles[:3] /= sc.au
    np.save('%s%s_popfile' % (directory, fileout), popfiles)

    return 
