import os
import sys
import time
import numpy as np
import scipy.constants as sc
from astropy.io import fits
import writeLIME as wl

# Runs LIME for a given header file and specified parameters.
# If nmodels > 1, will run LIME nmodels times and average the results.
# This is good to:
#   a) remove 'popcorning' in lower density regions,
#   b) reduce the Monte-Carlo noise.
# We can include an ortho / para ratio for the H2 (main collisional
# partner). This is NOT YET passed to wrieLIME which takes this into
# account. It is specified such that: orthoratio = n(oH2) / n(H2).


def runModels(chemheader, fileout, transitions, stellarmass, mach,
              pIntensity, sinkPoints, antialias, lte_only, blend, 
              nchan, velres, pxls, imgres, thetas, phi, distance, unit,
              datfile, nmodels=1, modelfile='model_template.c',
              dustfile=None, orthoratio=None, equaltemp=True):

    # Time the duration of run to print out at the end.
    
    t0 = time.time()
    print 'Running: ', chemheader
    
    # Make sure there are no negative models.

    nmodels = max(1,nmodels)

    # Create a temporary folder to work in.
    # We clean this up once finished.

    foldername = 'tempfolder'
    if not os.path.isdir(foldername):
        os.makedirs(foldername)
    else:
        suffix = 1
        while os.path.isdir(foldername+'%d' % suffix):
            suffix += 1
        foldername = foldername+'%d' % suffix
        os.makedirs(foldername)
    os.chdir(foldername)

    # Import the chemistry header file, template file,
    # dust opacities and collisional rates.

    os.system('cp ../%s .' % chemheader)
    os.system('cp ~/LimeFiles/%s .' % modelfile)
    if dustfile is None:    
        dustfile = 'jena_thin_e6.tab'
    os.system('cp ~/LimeFiles/%s .' % dustfile)
    os.system('cp ~/LimeFiles/%s .' % datfile.lower())

    # For each iteration of the model, create a new model_X.c file.

    for m in range(nmodels):
        print 'Running model %d...' % m
        
        wl.generateModelFile(chemheader, '%d' %  m, transitions, stellarmass,
                             mach, pIntensity, sinkPoints,
                             dustfile, antialias, lte_only, blend,  nchan,
                             velres, pxls, imgres, thetas, phi, distance,
                             unit, datfile, modelfile=modelfile, 
                             equaltemp=equaltemp)

        # Run the model.
        
        os.system('nohup lime -n -f -p 20 %d.c >nohup_%d.out 2>&1 &'
                  % (m, m))
        time.sleep(120)

    # Make sure all the models have run through.
    # This is achieved by checking that there are no more .x files in 
    # the current directory.
    # Currently, if the run fails, the .x file is left and not removed.
    # Need to implement some way to check that this is not the case.
    
    remaining = -1
    while len([fn for fn in os.listdir(os.curdir) if fn.endswith('.x')]) > 0:
        newremaining = len([fn for fn in os.listdir(os.curdir) if fn.endswith('.x')])
        if newremaining != remaining:
            print 'Waiting on %d models to run.' % newremaining
            remaining = newremaining
        time.sleep(60*remaining)
    print 'All instances complete.'
    
    # Average all the files with the same transition.
    # Implement something which calculates the standard devitation of
    # the models so we know what the averaging has done.

    for t in transitions:
        for i in thetas:
            readdata = np.array([fits.getdata('%d_%.3f_%d.fits' % (m, i, t), 0)
                                 for m in range(nmodels)])
            averaged = np.average(readdata, axis=0)          
            hdulist = fits.open('0_%.3f_%d.fits' % (i,t))
            hdulist[0].data = averaged
            hdulist.writeto(fileout+'_%.3f_%d.fits' % (i, t))
            os.system('mv %s_%.3f_%d.fits ../' % (fileout, i, t))
            np.save('../%s_%.3f_%d_stddev' % (fileout, i, t), averaged)

    # Exit the folder and delete it and everything in it.
    
    os.chdir('../')
    os.system('rm -rf %s' % foldername)

    # Print the total time.
    totaltime = time.time() - t0
    if totaltime < 60.:
        print 'Finished in %.2f seconds.' % totaltime
    elif (totaltime/60.) < 60.:
        print 'Finished in %.2f minutes.' % (totaltime/60.)
    else:
        print 'Finished in %.2f hours.' % (totaltime/3600.)

    print '\n\n'    
    return
