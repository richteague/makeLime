import os
import sys
import time
import writeLIME as wl
import averageModels as aM

# Given some input, run LIME models. Includes the option to run serveral
# and average to increase the resolution in the outer disk regions.

# chemheader (string) - 
#   File containing the header model. See writeHEADER for more details.
# fileout (string) - 
#   Prefix for the output .fits files.
# transitions (list of ints) -
#   Transitions to image the molecule at. Note that the transition number
#   is T-1 where T is the transition label in the LAMDA data file.
# stellarmass (float) - 
#   Mass of the central star used to calculate the Keplerian rotation.
# mach (float) -
#   Mach number of the micro-turbulence. This will be included in addition
#   to thermal broadening.
#   TODO: Include a way to have position specific values.
# pIntensity, sinkPoints, antialias, lte_only, blend, nchan,
# velres, pxls, imgres, distance, unit (ints and floats) - 
#   All values for the raytracing of LIME. See the LIME manual for info on these.
#   http://lime.readthedocs.io/en/v1.5/usermanual.html#parameters
# thetas (list of floats) - 
#   List of inclinations in radians to image the disk at. 0 is face on.
# phis (list of floats) - 
#   List of position angles in radians to image the disk at.
#   TODO: Check that these give sane results.
# datfile (string) - 
#   Collisional rate file to use. Make sure this is in LAMDA format.
# nmodels (int, optional) -
#   Number of models to run and average over.
# modelfile (string, optional) - 
#   Template file used to build the model.c file.
# dustfile (string, optional) - 
#   Dust opacities to use. Default are the Ossenlopf & Henning (1994).
# orthoratio (None, float or string, optional) -  
# orthoratio == None (default);
#   Assume {density} is the total H2 density and let density[0]
#   be the total H2 density.
# orthoratio == float;
#   Assume {density} is the total H2 density and scale
#   density[0] and density[1] to ortho and para H2 respectively.
# orthoratio == string;
#   Assume {density} is the oH2 density and {string} is the 
#   pH2 density, setting density[0] and density[1] respectively.
# popfile (boolean, optional) - 
#   Return the level population file. Will return for model 0.
#   TODO: Include a way to average over several models.
# returnnoise (boolean, optional) - 
#   If true and nmodels > 1, will return a .fits file containing the standard
#   deviation of each pixel from the array of models to give an idea of the MCMC
#   noise in the grid.
# directory (string, optional) - 
#   Directory to place all final .fits files. Default is the current
#   working directory.


def runModels(chemheader, fileout, transitions, stellarmass, mach,
              pIntensity, sinkPoints, antialias, lte_only, blend, 
              nchan, velres, pxls, imgres, distance, unit, thetas, phis,
              datfile, nmodels=1, modelfile=None,
              dustfile=None, orthoratio=None, popfile=False,
              returnnoise=False, directory='../'):


    # Record the time.
    t0 = time.time()
    print 'Running: ', chemheader
    
    
    # Check all the input values and change if necessary.
    nmodels = max(1,nmodels)
    if modelfile is None:
        modelfile = 'model_template.c'
    if dustfile is None: 
        dustfile = 'jena_thin_e6.tab'
    if unit > 3:
        unit = 0


    # Create a temporary folder to work in. We clean this up once 
    # finished to remove intermediate files. Into this folder we 
    # copy all the appropriate files.
    # TODO: Have a changable 'LimeFiles' folder.
    
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
    os.system('cp %s%s .' % (directory, chemheader))
    os.system('cp ~/LimeFiles/%s .' % modelfile)
    os.system('cp ~/LimeFiles/%s .' % dustfile)
    os.system('cp ~/LimeFiles/%s .' % datfile.lower())


    # For each iteration of the model, create a new model_X.c file
    # and run it. All output is saved into the nohup_X.out file which
    # is then deleted. Wait for 2 minutes before starting a new one.

    for m in range(nmodels):
        wl.generateModelFile(chemheader, m, transitions, stellarmass,
                             mach, pIntensity, sinkPoints,
                             dustfile, antialias, lte_only, blend,  nchan,
                             velres, pxls, imgres, thetas, phi, distance,
                             unit, datfile, modelfile=modelfile, 
                             equaltemp=equaltemp, popfile=popfile)
        os.system('nohup lime -n -f -p 20 %d.c >nohup_%d.out 2>&1 &'
                  % (m, m))
        time.sleep(120)


    # Make sure all the models have run through. This is achieved by checking 
    # that there are no more .x files in the directory.
    # TODO: Work out how to quit if a model segment faults.
    
    remaining = -1
    while len([fn for fn in os.listdir(os.curdir) if fn.endswith('.x')]) > 0:
        newremaining = len([fn for fn in os.listdir(os.curdir) if fn.endswith('.x')])
        if newremaining != remaining:
            print 'Waiting on %d models to run.' % newremaining
            remaining = newremaining
        time.sleep(60*remaining)
    print 'All instances complete.'
    
    
    # Average all the models and move to the appropriate folder.
    aM.averageModels(nmodels, 
                     thetas, 
                     phis, 
                     transitions, 
                     returnnoise=returnnoise,
                     directory=directory)

    
    # Exit the folder and delete it.
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
