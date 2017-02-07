"""
Read in the LAMDA collisional rate file and parse information.
"""

import numpy as np


class ratefile:

    # Read in a LAMDA rate file. Very rough for now.

    def __init__(self, fn):
        with open(fn) as f:
            self.filein = f.readlines()

        self.molecule = self.filein[1].strip()
        self.mu = float(self.filein[3].strip())
        self.nlevels = int(self.filein[5].strip())
        self.ntransitions = int(self.filein[8+self.nlevels].strip())
        self.npartners = self.filein[11+self.nlevels+self.ntransitions].strip()
        self.npartners = int(self.npartners)

        # Read in the partner names and the bounding line values.
        self.partners = []
        self.linestarts = []
        self.lineends = []
        linestart = 12 + self.nlevels+self.ntransitions
        for n in range(self.npartners):
            self.linestarts.append(linestart)
            names = np.array(['H2', 'pH2', 'oH2', 'e', 'H', 'He', 'H+'])
            name = names[int(self.filein[linestart+1][0])-1]
            lineend = linestart + 9 + int(self.filein[linestart+3].strip())
            self.partners.append(name)
            self.lineends.append(lineend)
            linestart = lineend

        # Read in the energy level structure.
        self.levels = self.filein[7:7+self.nlevels]
        self.levels = np.array([[n for n in levelsrow.strip().split()]
                                for levelsrow in self.levels]).T

        # Read in the radiative transitions.
        trans_start = 10 + self.nlevels
        trans_end = trans_start + self.ntransitions
        self.transitions = self.filein[trans_start:trans_end]
        self.transitions = np.array([[float(n) for n in t.strip().split()]
                                     for t in self.transitions]).T

        # Split into appropriate arrays.
        self.deltaE = self.levels[1]
        self.weights = self.levels[2]
        self.J = self.levels[3]
        self.EinsteinA = self.transitions[3]
        self.frequencies = self.transitions[4] * 1e9
        self.E_upper = self.transitions[5]

        return


def RADEXnames(key):
    """Switch between RADEX names and numbers."""
    if type(key) is list:
        return [RADEXnames(k) for k in key]
    names = np.array(['H2', 'pH2', 'oH2', 'electrons', 'H', 'He', 'H+'])
    if type(key) is not str:
        return names[int(key)-1]
    else:
        for n in range(names.size):
            if names[n] == key:
                return n+1
