''' Make model.c file. '''

#import inputparameters as ip

# Generate the model.c file to run with LIME.
def makeModelFile(chemheader, moldatfile, thetas, phis, transitions, nchan, velres, 
                  pIntensity=1e4, sinkPoints=1e3, dust='jena_thin_e6.tab', antialias=1, sampling=2,
                  outputfile=None, binoutputfile=None, gridfile=None, lte_only=1, imgres=0.05, distance=54., pxls=128, unit=0,
                  coordsys='cylindrical', opratio=None, dtemp=None, xmol=None, g2d=None, bvalue=50., btype='absolute',
                  stellarmass=0.6, modelnumber=0):

    tempfile = ['#include "lime.h"\n', '#include "math.h"\n', '#include "stdio.h"\n', '#include "stdlib.h"\n']
    
    # Include the header.
    if chemheader[-2:] == '.h':
        chemheader = chemheader[:-2]
    tempfile.append('#include "%s.h"\n\n' % chemheader)

    # Calculate the radius, minScale and ncells from the chemheader.
    import chemicalheader as ch
    minScale, radius, ncells, ndim = ch.valsfromheader(chemheader)
    if modelnumber == 0:
        print '\nAssuming input is %dD-%s coordinates.' % (ndim, coordsys)
        print 'Found minScale = %.2f au and radius = %.2f au.\n' % (minScale, radius) 
    
    # Include the imaging parameters.
    import imaging as im
    im.writeImageParameters(tempfile, radius=radius, minScale=minScale, moldatfile=moldatfile, modelnumber=modelnumber, 
                            thetas=thetas, phis=phis, transitions=transitions, nchan=nchan, velres=velres,
                            pIntensity=pIntensity, sinkPoints=sinkPoints, dust=dust, antialias=antialias, 
                            sampling=sampling, outputfile=outputfile, binoutputfile=binoutputfile,
                            gridfile=gridfile, lte_only=lte_only, imgres=imgres, distance=distance, 
                            pxls=pxls, unit=unit)


    # Include the interpolation functions.
    import interpolation as interp
    interp.writeFindValue(tempfile, ncells=ncells, coordsys=coordsys, ndim=ndim)


    # Include the model functions.
    import physicalstructure as ps
    ps.writeDensity(tempfile, coordsys, ndim, opratio)
    ps.writeTemperatures(tempfile, coordsys, ndim, dtemp)
    ps.writeAbundance(tempfile, coordsys, ndim, xmol, opratio)
    ps.writeGastoDust(tempfile, coordsys, ndim, g2d)
    ps.writeDopplerBroadening(tempfile, coordsys, ndim, bvalue, btype)
    ps.writeVelocityStructure(tempfile, coordsys, ndim, stellarmass)
    

    # Save the output.
    with open('model_%d.c' % modelnumber, 'w') as tosave:
        for line in tempfile:
            tosave.write('%s' % line)

    return 
