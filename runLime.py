''' Script to run a batch of Lime jobs. '''


def runLime(chemheader, moldatfile, thetas, phis, transitions, nchan, velres, nmodels=1, pIntensity=1e4, sinkPoints=1e3, dust='jena_thin_e6.tab',
            antialias=1, sampling=2, outputfile=None, binoutputfile=None, gridfile=None, non_lte=1, imgres=0.035, distance=54., pxls=128,
            unit=0, coordsys='cylindrical', ndim=2, opratio=None, dtemp=None, xmol=None, g2d=None, bval=50., btype='absolute', stellarmass=0.6):

    # Create the folder structure to hold it all in.
    # Move appropriate files there.
    # For each model, create the model.c file and run it.
    # Once all are run, average them.
    # Cleanup.

    return
