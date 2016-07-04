''' Make model.c file. '''

#import inputparameters as ip

# Generate the model.c file to run with LIME.
def makeModelFile(chemheader, moldatfile, thetas, phis, transitions, nchan, velres, 
                  pIntensity=1e4, sinkPoints=1e3, dust='jena_thin_e6.tab', antialias=1, sampling=2,
                  outputfile=None, binoutputfile=None, gridfile=None, non_lte=1, imgres=0.035, distance=54., pxls=128, unit=0,
                  coordsys='cyclindrical', ndim=2, opratio=None, dtemp=None, xmol=None, g2d=None, bvalue=50., btype='absolute',
                  stellarmass=0.6, modelnumber=0):

    tempfile = ['#include "lime.h"\n', '#include "math.h"\n', '#include "stdio.h"\n', '#include "stdlib.h"\n\n']

    # Calculate the radius, minScale and ncells from the chemheader.
    radius = 200.
    minScale = 2.
    ncells = 50000


    # Include the imaging parameters.
    import imaging as im
    im.writeImageParameters(tempfile, radius=radius, minScale=minScale, moldatfile=moldatfile, modelnumber=modelnumber, 
                            thetas=thetas, phis=phis, transitions=transitions, nchan=nchan, velres=velres,
                            pIntensity=pIntensity, sinkPoints=sinkPoints, dust=dust, antialias=antialias, 
                            sampling=sampling, outputfile=outputfile, binoutputfile=binoutputfile,
                            gridfile=gridfile, non_lte=non_lte, imgres=imgres, distance=distance, 
                            pxls=pxls, unit=unit)


    # Include the model functions.
    import physicalstructure as ps
    ps.writeDensity(tempfile, opratio, coordsys, ndim)
    ps.writeTemperatures(tempfile, dtemp, coordsys, ndim)
    ps.writeAbundance(tempfile, xmol, opratio, coordsys, ndim)
    ps.writeGastoDust(tempfile, g2d, coordsys, ndim)
    ps.writeDopplerBroadening(tempfile, bvalue, btype, coordsys, ndim)
    ps.writeVelocityStructure(tempfile, stellarmass, coordsys, ndim)
    

    with open('temp.c', 'w') as tosave:
        for line in tempfile:
            tosave.write('%s' % line)

    return 
