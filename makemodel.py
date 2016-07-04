''' Make model.c file. '''

#import inputparameters as ip

# Generate the model.c file to run with LIME.
def makeModelFile(chemheader, thetas, phis, transitions, nchan, velres, imgres=0.035, distance=54., pxls=128, unit=0,
                  coordsys='cyclindrical', ndim=2, opratio=None, dtemp=None, xmol=None, g2d=None, bvalue=50., btype='absolute',
                  stellarmass=0.6, modelnumber=0):

    tempfile = ['#include "lime.h"\n', '#include "math.h"\n', '#include "stdio.h"\n', '#include "stdlib.h"\n\n']

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
