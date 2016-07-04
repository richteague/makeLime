''' Write the physical structure functions. '''

def writeDensity(temp, opratio=None): 

    # Collider density.
    temp.append('void density(double x, double y, double z, *double density){\n')
    temp.append('density[0] = findvalue(cone, ctwo, cthree, dens);\n')

    # Ortho/para ratio correction.
    if type(opratio) is float:
        temp.append('density[0] *= (opratio / (1. + opratio));\n') 
        temp.append('density[1] = findvalue(cone, ctwo, cthree, dens);\n')
        temp.append('density[1] /= (1. + opratio);\n') 
        temp.append('if (density[1] < 1e-30){density[1] = 1e-30;}\n')
    elif type(opratio) is str:
        temp.append('density[1] = findvalue(cone, ctwo, cthree, %s);\n' % opratio)
        temp.append('if (density[1] < 1e-30){density[1] = 1e-30;}\n')
    elif opratio is not None:
        raise TypeError("opratio must be either None, a float or a string.")
    temp.append('if (density[0] < 1e-30){density[0] = 1e-30;}\n')
    temp.append('}\n')
    return 


def writeTemperatures(temp, dtemp=None):

    # Gas temperature.
    temp.append('void temperature(double x, double y, double z, double *temperature){\n')
    temp.append('temperature[0] = findvalue(cone, ctwo, cthree, temp);\n')
    temp.append('if (temperature[0] < 2.73){temperature[0] = 2.73;}\n')

    # Dust temperature.
    if type(dtemp) is float:
        temp.append('temperature[1] = %.3e * temperature[0];\n' % dtemp)
        temp.append('if (temperature[1] < 2.73){temperature[1] = 2.73;}\n')
    elif type(dtemp) is str:
        temp.append('temperature[1] = findvalue(cone, ctwo, cthree, %s);\n' % dtemp)
        temp.append('if (temperature[1] < 2.73){temperature[1] = 2.73;}\n')
    elif dtemp is not None:
        raise TypeError("dtemp must be either None, a float or a string.")
    temp.append('}\n')

    return


def writeAbundance(temp, xmol=None, opratio=None):

    # Molecular abundance.
    temp.append('void abundance(double x, double y, double z, double *abundance){\n')
    if xmol is None:
        temp.append('abundance[0] = findvalue(cone, ctwo, cthree, abund);\n')
    else:
        temp.append('abundance[0] = %.3e;\n' % float(xmol))

    # Ortho/para correction.
    if type(opratio) is float:
        temp.append('abundance[0] *= (1. + opratio) / opratio;\n')
    elif type(opratio) is str:
        temp.append('double val[2];\n')
        temp.append('density(x, y, z, &val[2]);\n')
        temp.append('abundance[0] *= val[1]/val[0] + 1.;\n')  
    elif opratio is not None:
        raise TypeError("opratio must be either None, a float or a string.")
    temp.append('if (abundance[0] < 0.){abundance[0] = 0.;}\n')  
    temp.append('}\n')

    return


def writeDusttoGas(temp, g2d=None):

    # Gas to dust ratio.
    if g2d is None:
        return
    else 
        temp.append('void gasIIdust(double x, double y, double z, double *gtd){\n')
    if type(g2d) is float:
        temp.append('*gtd = %.3e;\n' % g2d)
    elif type(g2d) is str:
        temp.append('*gtd = findvalue(cone, ctwo, cthree, %s);\n' % g2d)
    temp.append('if *gtd < 1){*gtd = 0.;}\n')
    temp.append('}\n')

    return
