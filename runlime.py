''' Script to run a batch of Lime jobs. '''

import os
import time
import makemodelfile as make
import limeclass as lime


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
def run(name='tempmodelname',
        headerfile='header.h',
        moldatfile='molecularrates.dat',
        transitions=[0],
        inclinations=[0],
        positionangles=[0],
        nchan=200,
        velres=20,
        pIntensity=1e5,
        sinkPoints=3e3,
        antialias=1,
        sampling=2,
        lte_only=1,
        imgres=0.065,
        distance=54.,
        pxls=128,
        unit=1,
        stellarmass=0.6,
        outputfile=False,
        binoutputfile=False,
        gridfile=False,
        dens=None,
        temp=None,
        dtemp=None,
        abund=None,
        g2d=None,
        doppler=None,
        dopplertype='absolute',
        coordsys='cylindrical',
        dust='jena_thin_e6.tab',
        directory='../',
        nmodels=1,
        returnnoise=False,
        cleanup=True,
        blend=1,
        opr_cp=None,
        r_inner=None,
        r_outer=None,
        depletion=1.0,
        waittime=10.):

    # Start the clock to time the running of models.
    t0 = time.time()

    # Create the temporary folder to work in and move there.
    fname = makeUniqueFolder()
    os.chdir(fname)

    # Generate a LIME model class instance.
    model = lime.model(name=name,
                       headerfile=headerfile,
                       moldatfile=moldatfile,
                       transitions=transitions,
                       inclinations=inclinations,
                       positionangles=positionangles,
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
                       dens=dens,
                       temp=temp,
                       dtemp=dtemp,
                       abund=abund,
                       g2d=g2d,
                       doppler=doppler,
                       dopplertype=dopplertype,
                       coordsys=coordsys,
                       dust=dust,
                       directory=directory,
                       nmodels=nmodels,
                       returnnoise=returnnoise,
                       blend=blend,
                       opr_cp=opr_cp,
                       depletion=depletion,
                       r_inner=r_inner,
                       r_outer=r_outer,
                       )

    # Copy the appropriate files into the working directory.
    path = os.path.dirname(__file__) + '/AuxFiles'
    os.system('cp ../%s .' % headerfile)
    os.system('cp %s/%s .' % (path, moldatfile))
    os.system('cp %s/%s .' % (path, dust))

    print '\n'
    print 'Input is %dD-%s coordinates.' % (model.ndim, model.coordsys)
    print 'Model bounds are %.2f and %.2f au' % (model.rin, model.rout)

    # For each iteration, run a model with a pause of waittime seconds.
    print '\n'
    for m in range(nmodels):
        print 'Running model %d of %d.' % (m+1, nmodels)
        make.makeFile(m, model)
        cmd = 'nohup lime -n -f -p 20 model_%d.c >nohup_%d.out 2>&1 &' % (m, m)
        os.system(cmd)
        waittime = max(10., waittime)
        time.sleep(waittime)

    # Make sure all the models have run.
    # Check the number of *.x files to guess how many are running.
    # If nohup_X.out contains segmentation fault, then quit.

    remaining = -1
    print '\n'
    while len([fn for fn in os.listdir('./') if fn.endswith('.x')]) > 0:
        nremaining = len([fn for fn in os.listdir('./') if fn.endswith('.x')])
        if nremaining != remaining:
            print 'Waiting on %d model(s) to run.' % nremaining
            remaining = nremaining
        time.sleep(10*remaining)

        # Check for segmentation faults.
        # If the number of nohup_X.out files with 'core dumped' in them
        # equals the number of lime_$$.x files left, assumped all have quit.

        nohups = [fn for fn in os.listdir('./') if fn.startswith('nohup')]
        coresdumped = 0
        for nh in nohups:
            if 'core dumped' in open(nh).read():
                coresdumped += 1
        if coresdumped == nremaining:
            print 'Found %d segmentation faults.' % coresdumped
            break

    if len([fn for fn in os.listdir('./') if fn.endswith('.fits')]) < nmodels:
        print 'Not all models were successfully run.'
        print 'Aborting without clean-up.\n'
        return
    else:
        print 'All instances complete.\n'

    # Combine the model ensemble.
    make.averageModels(model)
    make.combineGridfiles(model)
    make.combinePopfiles(model)
    make.combineBinPopFiles(model)
    make.getNoise(model)

    # Clean up.
    os.chdir('../')
    if cleanup:
        print 'Cleaning up temporary folders.'
        os.system('rm -rf %s' % fname)

    # Print the total time.
    print 'Finished in %s.\n\n' % seconds2hms(time.time() - t0)

    return
