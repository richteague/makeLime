import fileinput
import numpy as np
import scipy.constants as sc

# Create a header file for a LIME model.c file.
def makeHeader(chemfile, name=None, opr=1., headervals=None):

    # Default file formatting.
    if headervals is None:
        headervals = np.array(['rvals', 'zvals', 'dens', 'temp', 
                               'dtemp', 'size', 'd2g', 'abund'])

    # Read in the data.
    with np.errstate(divide='ignore'):
        data = np.loadtxt(chemfile, skiprows=3).T

    # Manipulate the data into LIME format.        
    arrname = ['rvals', 'zvals', 'dens', 'temp', 'dtemp', 'abund', 'g2d']
    towrite = np.zeros((len(arrname),data.shape[1]))    
    towrite[0] = data[headervals == 'rvals'][0]
    towrite[1] = data[headervals == 'zvals'][0]
    towrite[2] = data[headervals == 'dens'][0] * 1e3 / sc.m_p / 2.35
    towrite[3] = data[headervals == 'temp'][0]
    towrite[4] = data[headervals == 'dtemp'][0]
    towrite[5] = data[headervals == 'abund'][0] * opr * 1e6 / towrite[2]
    towrite[6] = data[headervals == 'd2g'][0]**-1.
    towrite[6] = np.where(towrite[6] < 50., 50., towrite[6])

    # Mask all the infinite values and all values where there is no density.
    towrite = np.where(np.isfinite(towrite), towrite, 0.)
    for i in range(3,7):
        towrite[i] = np.where(towrite[2] > 0., towrite[i], 0.)

    # Function to write the headers.
    def getHeaderString(array, name):
        tosave = 'const static double %s[%d] = {' % (name, array.size)
        for val in array:
            tosave += '%.3e, ' % val
        tosave = tosave[:-2] + '};\n'
        return tosave

    # Write all the arrays as a C array.    
    hstring = ''
    for i, array in enumerate(towrite):
        hstring += getHeaderString(array, arrname[i])
            
    # Write the output file.
    if name == None:
        name = chemfile[:-4]
    elif name[-2:] == '.h':
        name = name[:-2]
    with open('%s.h' % name, 'w') as hfile:
        hfile.write('%s' % hstring)
    
    print 'Written file to %s.h, have fun!' % name

    return 

