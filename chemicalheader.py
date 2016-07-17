import fileinput
import numpy as np
import scipy.constants as sc
import re

# Reads the inner and outer radial points to define LIME's 
# computational domain. Reads in the number of cells for the
# interpolate routine. 

def valsfromheader(headername, coordsys='cylindrical'):

    tempname = headername
    if tempname[-2:] != '.h':
        tempname += '.h'
    with open('../'+tempname) as f:
        header = f.readlines()

    # Find the numer of cells.
    i = 0
    while header[0][i] != ']':
        if header[0][i] == '[':
            j = i
        i += 1
    ncells= int(header[0][j+1:i])

    # Calculate the minimum and maximum values.
    # Find the indices of the c1arr and c2arr.
    
    arrnames = arrsfromheader(headername)
    c1 = arrnames.index('c1arr')
    c2 = arrnames.index('c2arr')
    if 'c3arr' in arrnames:
        ndim = 3
    else:
        ndim = 2

    if coordsys is 'cylindrical':
        while header[c1][i-1] != '{':
            i += 1
        rvals = header[c1][i:-3] 
        rvals = np.array([float(v) for v in re.split(', ', rvals)])
        i = 0
        while header[c2][i-1] != '{':
            i += 1
        zvals = header[c2][i:-3]
        zvals = np.array([float(v) for v in re.split(', ', zvals)])        
        rin = np.nanmin(np.hypot(rvals, zvals))
        rout = np.nanmax(np.hypot(rvals, zvals))

    elif coordsys is 'polar':
        while header[c1][i-1] != '{':
            i += 1
        rvals = header[c1][i:-3] 
        rvals = np.array([float(v) for v in re.split(', ', rvals)])
        rin = np.nanmin(rvals)
        rout = np.nanmax(rvals)

    else:
        raise ValueError

    return rin, rout, ncells, ndim


# Reads the array names from the header file.
# Calls the parsename() function for each array in the header.
    
def arrsfromheader(headername):
    tempname = headername
    if tempname[-2:] != '.h':
        tempname += '.h'
    with open('../'+tempname) as f:
        header = f.readlines()
    return [parsename(line) for line in header]


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

def makeHeader(chemfile, headername=None, writevals=True, returnvals=False):

    # Format the data into LIME appropriate units.
    # If returnvals == True then return NCELLS, RIN and ROUT.

    # For a gridded model, make sure that we clip areas.
    # We choose regions where the dust temperature, fourth column, is 0.

    data = np.loadtxt(chemfile, skiprows=3).T
    data = np.array([np.where(data[3] != 0, param, 0.0) for param in data])

    # Make the conversions to LIME appropriate units.
    # Main collider density is in [m^-3]. 
    # Relative abundance is with respect to the main collider density.
    # Temperatures are all in [K].

    with np.errstate(divide='ignore'):
        data[2] /= 2.33 * sc.m_p * 1e3
        data[7] = np.where(data[2] != 0.0, data[7]/data[2], 0.0)
        data[2] *= 1e6
        data[6] = np.where(data[6] != 0.0, 1./data[6], 100.)
        data    = np.where(~np.isfinite(data), 0.0, data)

    if writevals:
    
        # Function to write the arrays.   
 
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

