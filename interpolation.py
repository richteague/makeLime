''' Interpolation routines for LIME '''

import fileinput

def writeCoords(temp, coordsys='cylindrical', ndim=2):

    if not (coordsys is 'cylindrical' and ndim is 2):
        raise NotImplementedError

    if coordsys is 'cylindrical':
        if ndim is 2:
            temp.append('\tdouble cone = sqrt(x*x + y*y) / AU;\n')
            temp.append('\tdouble ctwo = atan2(y, x);\n')
            temp.append('\tdouble cthree = 0.;\n\n')

    return

def writeFindValue(temp, ncells, coordsys='cyclindrical', ndim=2):

    # Include the appropriate interpolation routines. 
    if not (coordsys is 'cylindrical' and ndim is 2):
        raise NotImplementedError 


    # Find the correct path to the files.
    import os
    path =  os.path.dirname(__file__)

    # Grab the interpolation codes.
    with open(path+'/InterpolationRoutines/interp_%dD_%s.txt' % (ndim, coordsys)) as f:
        lines = f.readlines()
    for line in lines:
        line = line.replace('NCELLS', '%d' % ncells)
        temp.append(line)
    temp.append('\n\n')


    # Grab the cell finding codes.
    with open(path+'/InterpolationRoutines/findcell_%dD_%s.txt' % (ndim, coordsys)) as f:
        lines = f.readlines()
    for line in lines:
        line = line.replace('NCELLS', '%d' % ncells)
        temp.append(line)
    temp.append('\n\n')


    # Include the wrapper.
    temp.append('double findvalue(double cone, double ctwo, double cthree, const double arr[%d]){\n\n' % ncells) 
    temp.append('\tdouble value;\n')
    temp.append('\tint aidx, bidx, cidx, didx;\n')

    if ndim is 2:
        temp.append('\n\tfindcell(cone, ctwo, &aidx, &bidx, &cidx, &didx);\n')
        temp.append('\tif (aidx >= 0) {\n')
        temp.append('\t\tvalue = linterpolate(cone, ctwo, aidx, bidx, cidx, didx, arr);\n')
        temp.append('\t} else {\n\t\tvalue = -1.;\n\t}\n\n')
        temp.append('\treturn value;\n')

    elif ndim is 3:
        temp.append('\tint eidx, fidx, gidx, hidx;\n')    
        temp.append('\n\tfindcell(cone, ctwo, cthree, &aidx, &bidx, &cidx, &didx, &eidx, &fidx, &gidx, &hidx);\n')
        temp.append('\tif (aidx >= 0) {\n')
        temp.append('\t\tvalue = linterpolate(cone, ctwo, cthree, aidx, bidx, cidx, didx, eidx, fidx, gidx, hidx, arr);\n')
        temp.append('\t} else {\n\t\tvalue = -1.;\n\t}\n\n')
        temp.append('\treturn value;\n')
        
    else:
        raise ValueError

    temp.append('\n}\n\n\n')

    return
