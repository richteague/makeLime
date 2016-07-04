# Write the physical structure.

'''
    Write the collider density. Three options for opratio:
    opratio = None  : assumes 'dens' is n(H2) array.
    opratio = float : assumes an opratio = n(oH2)/n(pH2).
    opratio = str   : assumes 'dens' is n(oH2) array and opratio is n(pH2) array.
'''
def writeDensity(temp, opratio=None): 
    temp.append('void density(double x, double y, double z, *double density){\n'
    temp.append('density[0] = findvalue(cone, ctwo, cthree, dens);\n')
    if type(opratio) is float:
        temp.append('density[0] *= (opratio / (1. + opratio));\n') 
        temp.append('density[1] = findvalue(cone, ctwo, cthree, dens);\n')
        temp.append('density[1] /= (1. + opratio);\n') 
        temp.append('if (density[1] < 1e-30){density[1] = 1e-30;}\n')
    elif type(opratio) is str:
        temp.append('density[1] = findvalue(cone, ctwo, cthree, %s);\n' % opratio)
        temp.append('if (density[1] < 1e-30){density[1] = 1e-30;}\n')
    else:
        raise TypeError("opratio must be either None, a float or a string.")
    temp.append('if (density[0] < 1e-30){density[0] = 1e-30;}\n')
    temp.append('}\n')
    return 


'''
    Write the temperature profiles. Three options for dtemp:
    dtemp = None    : assumes equal gas and dust temperatures.
    dtemp = float   : assumes dust temp is a homogeneous rescaling, 'dtemp', of the gas temp.
    dtemp = str     : assumes dust temp given in 'dtemp' array.
'''
def writeTemperatures(temp, dtemp=None):
    temp.append('void temperature(double x, double y, double z, double *temperature){\n')
    temp.append('temperature[0] = findvalue(cone, ctwo, cthree, temp);\n')
    temp.append('if (temperature[0] < 2.73){temperature[0] = 2.73;}\n')
    if type(dtemp) is float:
        temp.append('temperature[1] = %.3f * temperature[0];\n' % dtemp)
        temp.append('if (temperature[1] < 2.73){temperature[1] = 2.73;}\n')
    elif type(dtemp) is str:
        temp.append('temperature[1] = findvalue(cone, ctwo, cthree, %s);\n' % dtemp)
        temp.append('if (temperature[1] < 2.73){temperature[1] = 2.73;}\n')
    else:
        raise TypeError("dtemp must be either None, a float or a string.")
    temp.append('}\n')
    return


def writeAbundance(temp, xmol=None, opratio=None):
    temp.append('void abundance(double x, double y, double z, double *density){\n')
    if xmol is None:
        temp.append('abundance[0] = findvalue(cone, ctwo, cthree, abund);\n')
    else:
        temp.append('abundance[0] = %.3f;\n' % float(xmol))
    if type(opratio) is float:
        temp.append('abundance[0] *= (1. + opratio) / opratio;\n')
    elif type(opratio) is str:
        temp.append('double val[2];\n')
        temp.append('density(x, y, z, &val[2]);\n')
        temp.append('abundance[0] *= val[1]/val[0] + 1.;\n')  
    else:
        raise TypeError("opratio must be either None, a float or a string.")
    temp.append('if (abundance[0] < 0.){abundance[0] = 0.;}\n')  
    temp.append('}\n')
    return
