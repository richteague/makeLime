''' Script to run a batch of Lime jobs. '''

import os
import sys
import time
import makemodelfile as make
import LIMEclass as lime
import averagemodels as avg


# Make a unique folder.
def makeUniqueFolder(fname='tempfolder', path='./'):
    if not os.path.isdir(path+fname):
        os.makedirs(path+fname)
    else:
        suffix = 1
        while os.path.isdir(path+fname+'%d' % suffix):
            suffix += 1
        fname = fname+'%d' % suffix
        os.makedirs(path+fname)
    return fname


# Convert seconds to hours, minutes, seconds.
def seconds2hms(seconds):
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    return '%d:%02d:%02d' % (h, m, s)


# Run LIME models.
def runLime(headerfile, moldatfile, thetas, transitions, nchan, velres, 
            fileout='temporary.fits', phis=[0], nmodels=1, pIntensity=1e4, 
            sinkPoints=1e3, dust='jena_thin_e6.tab',
            antialias=1, sampling=2, outputfile=None,
            binoutputfile=None, gridfile=None, lte_only=1, imgres=0.05, 
            distance=54., pxls=128,
            unit=0, coordsys='cylindrical', opratio=None, dtemp=None, 
            xmol=None, g2d=None, bvalue=50., btype='absolute', 
            stellarmass=0.6, cleanup=True, waittime=20, directory='../'):


    if opratio is not None:
        raise NotImplementedError('Still need work on this.')

    # Make sure the lists of imaging parameters are lists.

    if type(thetas) is not list:
        thetas = [thetas]
    if type(transitions) is not list:
        transitions = [transitions]
    if type(phis) is not list:
        phis = [phis]
    if coordsys not in ['cylindrical', 'polar']:
        raise ValueError("""Wrong coordsys value.
                            Must be either 'cylindrical' or 'polar'.""")
    if btype not in ['absolute', 'mach']:
        raise ValueError("""Wrong btype value.
                            Must be either 'absolute' or 'mach'.""")
    if opratio is not None:
        raise NotImplementedError


    t0 = time.time()

    # Create the temporary folder to work in and move there.
    fname = makeUniqueFolder()
    os.chdir(fname)
    path = os.path.dirname(__file__) + '/AuxFiles'
    os.system('cp ../%s .' % headerfile)
    os.system('cp %s/%s .' % (path, moldatfile))
    os.system('cp %s/%s .' % (path, dust))


    # Generate a LIME model class instance.
    model = lime.model(fileout=fileout,
                       headerfile=headerfile,
                       moldatfile=moldatfile,
                       transitions=transitions, 
                       thetas=thetas,
                       nchan=nchan,
                       velres=velres,
                       pIntensity=pIntensity,
                       sinkPoints=sinkPoints,
                       antialias=antialias,
                       sampling=sampling,
                       lte_only=lte_only,
                       imgres=imgres,
                       distance=distance,
                       pxls=pxls,
                       unit=unit,
                       stellarmass=stellarmass,
                       outputfile=outputfile,
                       binoutputfile=binoutputfile,
                       gridfile=gridfile,
                       opratio=opratio,
                       dtemp=dtemp,
                       phis=phis,
                       xmol=xmol,
                       g2d=g2d,
                       bvalue=bvalue,
                       btype=btype,
                       coordsys=coordsys,
                       dust=dust,
                       directory=directory,
                       )
                       
    print '\nAssuming input is %dD-%s coordinates.' % (model.hdr.ndim, model.coordsys)
    print 'Found minScale = %.2f au and radius = %.2f au.\n' % (model.hdr.rin, model.hdr.rout) 

    # For each iteration, run a model with a pause of waittime seconds.
    
    print '\n'
    for m in range(nmodels):
        print 'Running model %d of %d.' % (m+1, nmodels)
        make.makeFile(m, model)
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
    if len([fn for fn in os.listdir('./') if fn.endswith('.fits')]) < nmodels:
        print 'Not all models were successfully run. Aborting without clean-up.\n'
        return
    else:
        print 'All instances complete.\n'

    # If more than one model is run, average them and move them.
    
    if fileout[-5:] == '.fits':
        fileout = fileout[:-5]
    if nmodels > 1:
        avg.averageModels(model)
        if outputfile is not None:
            avg.combinePopfiles(model)
    else:
        for t in thetas:
            for p in phis:
                for j in transitions:
                    os.system('mv *.fits %s/%s_%.3f_%.3f_%d.fits' % (directory, fileout, t, p, j))



    # Move the model out and then, if required, clear the folder.
    os.chdir('../')
    if cleanup:
        print 'Cleaning up temporary folders.'
        os.system('rm -rf %s' % fname)


    # Print the total time.
    print 'Finished in %s.\n\n' % seconds2hms(time.time() - t0)

    return
