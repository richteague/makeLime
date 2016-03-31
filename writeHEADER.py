import numpy as np
import scipy.constants as sc

# Create a header file for a LIME model.c file. Assumes Dima's output.


def makeHeader(chemfile, headername=None, writevals=True, returnvals=False):

    # Format the data into LIME appropriate units.
    # If returnvals == True then return NCELLS, RIN and ROUT.

    with np.errstate(divide='ignore'):
        data = np.loadtxt(chemfile, skiprows=3).T
        data[2] /= 2.33 * sc.m_p * 1e3 / 2.
        data[7] = np.where(data[2] != 0.0, data[7]/data[2], 0.0)
        data[2] *= 1e6
        data[6] = np.where(data[6] != 0.0, 1./data[6], 100.)
        data = np.where(~np.isfinite(data), 0.0, data)

    if writevals:

        # Only do this if we want to write to file.

        def getHeaderString(array, name):
            tosave = 'const static double %s[%d] = {' % (name, array.size)
            for val in array:
                tosave += '%.3e, ' % val
            tosave = tosave[:-2] + '};\n'
            return tosave

        # Make sure these are the same as in the template file.
        # Note that we skip the average grainsize array.

        arrnames = ['rvals', 'zvals', 'dens', 'temp', 'dtemp',
                    'size', 'gastodust', 'abund']
        hstring = ''
        for i, array in enumerate(data):
            if i == 5:
                continue
            else:
                hstring += getHeaderString(array, arrnames[i])

        # Write the output file.

        if headername is None:
            headername = chemfile[:-4]
        elif headername[-2:] == '.h':
            headername = headername[:-2]
        with open('%s.h' % headername, 'w') as hfile:
            hfile.write('%s' % hstring)

    if returnvals is False:
        return
    else:
        rmin = np.hypot(data[0], data[1]).min()
        rmax = np.hypot(data[0], data[1]).max()
        ncells = data[0].size
        return ncells, rmin, rmax
