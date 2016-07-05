''' Script to run a batch of Lime jobs. '''

def makeUniqueFolder(fname='tempfolder', path='./'):

    # Make a unique folder.

    import os 
    if not os.path.isdir(path+fname):
        os.makedirs(path+fname)
    else:
        suffix = 1
        while os.path.isdir(path+fname+'%d' % suffix):
            suffix += 1
        fname = fname+'%d' % suffix
        os.makedirs(path+fname)
    return fname


def seconds2hms(seconds):

    # Convert seconds to hours, minutes, seconds.

    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    return '%d:%02d:%02d' % (h, m, s)

def runLime(chemheader, moldatfile, fileout, thetas, phis, transitions, nchan, velres, nmodels=1, pIntensity=1e4, sinkPoints=1e3, dust='jena_thin_e6.tab',
            antialias=1, sampling=2, outputfile=None, binoutputfile=None, gridfile=None, lte_only=1, imgres=0.05, distance=54., pxls=128,
            unit=0, coordsys='cylindrical', ndim=2, opratio=None, dtemp=None, xmol=None, g2d=None, bvalue=50., btype='absolute', stellarmass=0.6,
            cleanup=True, waittime=120, directory='../'):

    # Build, run and average multiple LIME models for a given chemical model.

    import os
    import sys
    import time
    import makemodel as make

    t0 = time.time()

    # Create the temporary folder to work in.
    fname = makeUniqueFolder()

    # Move all appropriate files there.
    path = os.path.dirname(__file__) + '/AuxFiles'
    os.chdir(fname)
    os.system('cp ../%s .' % chemheader)
    os.system('cp %s/%s .' % (path, moldatfile))
    os.system('cp %s/%s .' % (path, dust))


    # For each iteration, run a model with a pause of waittime seconds.
    for m in range(nmodels):
        print 'Running model %d of %d.' % (m+1, nmodels)

        # Make sure the file outputs are not over written.
        if outputfile is not None:
            toutputfile = outputfile + '_%d' % m 
        else:
            toutputfile = None

        if binoutputfile is not None:
            tbinoutputfile = binoutputfile + '_%d' % m 
        else:
            tbinoutputfile = None

        if gridfile is not None:
            tgridfile = gridfile + '_%d' % m 
        else:
            tgridfile = None

        # Make the model.c file.
        make.makeModelFile(chemheader=chemheader, moldatfile=moldatfile, thetas=thetas, phis=phis, transitions=transitions, nchan=nchan, velres=velres,
                           pIntensity=pIntensity, sinkPoints=sinkPoints, dust=dust, antialias=antialias, sampling=sampling, outputfile=toutputfile, 
                           binoutputfile=binoutputfile, gridfile=gridfile, lte_only=lte_only, imgres=imgres, distance=distance, pxls=pxls, unit=unit,
                           coordsys=coordsys, ndim=ndim, opratio=opratio, dtemp=dtemp, xmol=xmol, g2d=g2d, bvalue=bvalue, btype=btype,
                           stellarmass=stellarmass, modelnumber=m)

        # Run the file.
        os.system('nohup lime -n -f -p 20 model_%d.c >nohup_%d.out 2>&1 &' % (m, m))
        waittime = max(10., waittime)
        time.sleep(waittime)


    # Make sure all the models have run.
    remaining = -1
    while len([fn for fn in os.listdir(os.curdir) if fn.endswith('.x')]) > 0:
        newremaining = len([fn for fn in os.listdir(os.curdir) if fn.endswith('.x')])
        if newremaining != remaining:
            print 'Waiting on %d models to run.' % newremaining
            remaining = newremaining
        time.sleep(60*remaining)
    print 'All instances complete.'


    # If more than one model is run, average them.
    # Move the resulting file to the output directory.
    if fileout[-5:] == '.fits':
        fileout = fileout[:-5]
    if nmodels > 1:
        import averagemodels as avg
        avg.averageModels(nmodels, thetas, phis, transitions, fileout,
                          returnnoise=False, directory=directory)
    else:
        for t in thetas:
            for p in phis:
                for j in transitions:
                    os.system('mv *.fits %s/%s_%.3f_%.3f_%d.fits' % (directory, fileout, t, p, j))
    

    # Move the model out and then, if required, clear the folder.
    if cleanup:
        print 'Cleaning up temporary folders.'
        os.system('rm -rf %f' % fname)


    # Print the total time.
    print 'Finished in %s.\n\n' % seconds2hms(time.time() - t0)
    
    return

