import fileinput
import numpy as np
import scipy.constants as sc


# Reads the inner and outer radial points to define LIME's 
# computational domain. Reads in the number of cells for the
# interpolate routine. Assumes that rvals is the first line 
# and zvals is the second.

def valsfromheader(headername):
    with open('../'+headername) as f:
        header = f.readlines()
    i = 18
    while header[0][i] != ']':
        if header[0][i] == '[':
            j = i
        i += 1
    ncells= int(header[0][j+1:i])
    while header[0][i] != ',':
        if header[0][i] == '{':
            j = i
        i += 1
    rin = float(header[0][j+1:i])
    i = -2
    while header[0][i] != ',':
        i -= 1
    rout = float(header[0][i+2:-3])
    i = -2
    while header[1][i] != ',':
        i -= 1
    rout = np.hypot(float(header[1][i+2:-3]), rout)
    return rin, rout, ncells


# Reads the array names from the header file.
# Calls the parsename() function for each array in the header.
    
def arrsfromheader(headername):    
    with open('../'+headername) as f:
        header = f.readlines()
    return np.array([parsename(line) for line in header])


# Parses the array name from a C array declaration.

def parsename(line):
    i = 0
    while line[i] != '[':
        i += 1
    j = i
    while line[i] != ' ':
        i -= 1
    return line[i+1:j]

# Create a LIME compatible header from Dima's chemical model.
# TODO: This is specific. Make a more general version.

def makeHeader(chemfile, headername=None, writevals=True, returnvals=False):

    # Format the data into LIME appropriate units.
    # If returnvals == True then return NCELLS, RIN and ROUT.

    with np.errstate(divide='ignore'):
        data = np.loadtxt(chemfile, skiprows=3).T
        data[2] /= 2.33 * sc.m_p * 1e3 / 2.
        data[7] = np.where(data[2] != 0.0, data[7]/data[2], 0.0)
        data[2] *= 1e6
        data[6] = np.where(data[6] != 0.0, 1./data[6], 100.)
        data    = np.where(~np.isfinite(data), 0.0, data)

    if writevals:
    
        # Only do this if we want to write.
    
        def getHeaderString(array, name):
            tosave = 'const static double %s[%d] = {' % (name, array.size)
            for val in array:
                tosave += '%.3e, ' % val
            tosave = tosave[:-2] + '};\n'
            return tosave

        # Make sure these are the same as in the template file.
        # Note that we skip the average grainsize array.
        
        arrnames = ['rvals', 'zvals', 'dens', 'temp', 'dtemp', 'size', 'gastodust', 'abund']
        hstring = ''
        for i, array in enumerate(data):
            if i == 5:
                continue
            else:
                hstring += getHeaderString(array, arrnames[i])
                
        # Write the output file.
        
        if headername == None:
            headername = chemfile[:-4]
        elif headername[-2:] == '.h':
            headername = headername[:-2]
        with open('%s.h' % headername, 'w') as hfile:
            hfile.write('%s' % hstring)

    if returnvals == False:
        return
    else:
        return data[0].size, np.hypot(data[0], data[1]).min(), np.hypot(data[0], data[1]).max()

