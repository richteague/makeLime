import numpy as np

# Class to read in the header files used.


class headerFile:

    def __init__(self, filename, coordsys='cylindrical'):

        if coordsys == ('cylindrical' or coordsys == 'polar'):
            self.coordsys = coordsys
        else:
            raise NotImplementedError("Only 'cyclindrical' or 'polar'.")
        self.fn = filename
        if self.fn[-2:] != '.h':
            raise ValueError('headerfile must have a *.h extention.')
        # Read in the header.

        with open(self.fn) as f:
            self.hdr = f.readlines()
        self.ncells = self.parsencells()
        self.arrnames = [self.parsename(ln) for ln in self.hdr]

        c1 = self.arrnames.index('c1arr')
        c2 = self.arrnames.index('c2arr')
        if 'c3arr' in self.arrnames:
            self.ndim = 3
        else:
            self.ndim = 2

        i = 0
        if self.coordsys == 'cylindrical':
            while self.hdr[c1][i-1] != '{':
                i += 1
            self.rvals = self.hdr[c1][i:-3]
            self.rvals = np.array([float(v) for v in self.rvals.split(', ')])
            i = 0
            while self.hdr[c2][i-1] != '{':
                i += 1
            self.zvals = self.hdr[c2][i:-3]
            self.zvals = np.array([float(v) for v in self.zvals.split(', ')])
            self.rin = np.nanmin(np.hypot(self.rvals, self.zvals))
            self.rout = np.nanmax(np.hypot(self.rvals, self.zvals))

        elif self.coordsys == 'polar':
            while self.hdr[c1][i-1] != '{':
                i += 1
            self.rvals = self.hdr[c1][i:-3]
            self.rvals = np.array([float(v) for v in self.rvals.split(', ')])
            self.rin = np.nanmin(self.rvals)
            self.rout = np.nanmax(self.rvals)

        else:
            raise NotImplementedError("Only 'cylindrical' or 'polar'.")

        return

    # Parse the number of cells.
    def parsencells(self):
        i = 0
        while self.hdr[0][i] != ']':
            if self.hdr[0][i] == '[':
                j = i
            i += 1
        return int(self.hdr[0][j+1:i])

    # Parse the array name from a line.
    def parsename(self, ln):
        i = 0
        while ln[i] != '[':
            i += 1
        j = i
        while ln[i] != ' ':
            i -= 1
        return ln[i+1:j]
