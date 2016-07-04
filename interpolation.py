''' Interpolation routines for LIME '''

import fileinput

def writeCoords(temp, coordsys='cyclindrical', ndim=2):

    if not (coordsys is 'cyclindrical' and ndim is 2):
        raise NotImplementedError

    if coordsys is 'cyclindrical':
        if ndim is 2:
            temp.append('\tdouble cone = sqrt(x*x + y*y) / AU;\n')
            temp.append('\tdouble ctwo = atan2(y, x);\n')
            temp.append('\tdouble cthree = 0.;\n\n')

    return

def writeFindValue(temp, ncells, coordsys='cyclindrical', ndim=2):

    # Include the appropriate interpolation routines. 
    if not (coordsys is 'cyclindrical' and ndim is 2):
        raise NotImplementedError 


    # Grab the interpolation codes.
    with open('InterpolationRoutines/interp_%dD_%s.txt' % (ndim, coordsys)) as f:
        lines = f.readlines()
    for line in lines:
        temp.append(line+'\n')
    temp.append('\n\n')


    # Grab the cell finding codes.
    with open('InterpolationRoutines/findcell_%dD_%s.txt' % (ndim, coordsys)) as f:
        lines = f.readlines()
    for line in lines:
        temp.append(line+'\n')
    temp.append('\n\n')


    # Include the wrapper.
    temp.append('double findvalue(double cone, double ctwo, double cthree, cont double arr[%d]){\n' % ncells) 
    temp.append('double value;\n')
    temp.append('int aidx, int bidx, int cidx, didx;\n')
    temp.append('int eidx, int fidx, int gidx, didx;\n')    

    if ndim is 2:
        temp.append('findcell(cone, ctwo, &aidx, &bidx, &cidx, &didx);\n')
        temp.append('if (aidx >= 0){\n')
        temp.append('{value = linterpolate(cone, ctwo, aidx, bidx, cidx, didx, arr);}')
        temp.append('else {value = -1.;}\n')
        temp.append('return value;')

    elif ndim is 3:
        temp.append('findcell(cone, ctwo, cthree, &aidx, &bidx, &cidx, &didx, &eidx, &fidx, &gidx, &hidx);\n')
        temp.append('if (aidx >= 0){\n')
        temp.append('{value = linterpolate(cone, ctwo, cthree, aidx, bidx, cidx, didx, eidx, fidx, gidx, hidx, arr);}')
        temp.append('else {value = -1.;}\n')
        temp.append('return value;')
        
    else:
        raise ValueError

    temp.append('}\n\n')

    return
