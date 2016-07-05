''' Interpolation routines for LIME '''

import fileinput

def writeCoords(temp, coordsys='cylindrical', ndim=2):

    if not (coordsys is 'cylindrical' and ndim is 2):
        raise NotImplementedError

    if coordsys is 'cylindrical':
        if ndim is 2:
            temp.append('\tdouble c1 = sqrt(x*x + y*y) / AU;\n')
            temp.append('\tdouble c2 = fabs(z) / AU;\n')
            temp.append('\tdouble c3 = -1.;\n')
    return

def writeFindValue(temp, ncells, coordsys='cyclindrical', ndim=2):

    # Include the appropriate interpolation routines. 
    if not (coordsys is 'cylindrical' and ndim is 2):
        raise NotImplementedError 

    # Find the correct path to the files.
    import os
    with open(os.path.dirname(__file__)+'/InterpolationRoutines/%dD_%s.c' % (ndim, coordsys)) as f:
        lines = f.readlines()
    for line in lines:
        line = line.replace('NCELLS', '%d' % ncells)
        temp.append(line)
    temp.append('\n\n')


    return
