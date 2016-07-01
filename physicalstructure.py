# Write the physical structure.

'''
    Write the collider density. Three options for opratio:
    opratio = None  : assumes 'dens' is n(H2) array.
    opratio = float : assumes an opratio = n(oH2)/n(pH2).
    opratio = str   : assumes 'dens' is n(oH2) array and opratio is n(pH2) array.
'''
def writeDensity(template, opratio=None):

    # Find the appropriate line.
    for l, line in enumerate(template):
        parsedline = ''.join(line.split())
            if parsedline.startswith('//Density'):
                break

    # Add in appropriate lines.
    if opratio is None:
        template.insert(l+1, 'density[0] = findvalue(cone, ctwo, cthree, dens);\n')
        template.insert(l+2, 'if (density[0] < 1e-30){density[0] = 1e-30;}\n')
    if type(opratio) is float:
        template.insert(l+1, 'density[0] = findvalue(cone, ctwo, cthree, dens);\n')
        template.insert(l+2, 'density[1] = findvalue(cone, ctwo, cthree, dens);\n')
        template.insert(l+3, 'density[0] *= (opratio / (1. + opratio));\n') 
        template.insert(l+4, 'density10] /= (1. + opratio);\n') 
        template.insert(l+5, 'if (density[0] < 1e-30){density[0] = 1e-30;}\n')
        template.insert(l+6, 'if (density[1] < 1e-30){density[1] = 1e-30;}\n')
    elif type(opratio) is str:
        template.insert(l+1, 'density[0] = findvalue(cone, ctwo, cthree, dens);\n')
        template.insert(l+2, 'density[1] = findvalue(cone, ctwo, cthree, %s);\n' % opratio)
        template.insert(l+5, 'if (density[0] < 1e-30){density[0] = 1e-30;}\n')
        template.insert(l+6, 'if (density[1] < 1e-30){density[1] = 1e-30;}\n')
    else:
        raise TypeError("opratio must be either None, a float or a string.")

    return template



'''
    Write the temperature profiles. Three options for dtemp:
    dtemp = None    : assumes equal gas and dust temperatures.
    dtemp = float   : assumes dust temp is a homogeneous rescaling, 'dtemp', of the gas temp.
    dtemp = str     : assumes dust temp given in 'dtemp' array.
'''
def writeTemperatures(template, dtemp=None):

    # Find the appropriate line.
    for l, line in enumerate(template):
        parsedline = ''.join(line.split())
            if parsedline.startswith('//Temperature'):
                break

    # Add in the appropriate lines
    template.insert(l+1, 'temperature[0] = findvalue(cone, ctwo, cthree, temp);\n')
    template.insert(l+2, 'if (temperature[0] < 2.73){temperature[0] = 2.73;}\n')
    if type(dtemp) is float:
        template.insert(l+3, 'temperature[1] = %.3f * temperature[0];\n' % dtemp)
        template.insert(l+4, 'if (temperature[1] < 2.73){temperature[1] = 2.73;}\n')
    elif type(dtemp) is str:
        template.insert(l+3, 'temperature[1] = findvalue(cone, ctwo, cthree, %s);\n' % dtemp)
        template.insert(l+4, 'if (temperature[1] < 2.73){temperature[1] = 2.73;}\n')
    else:
        raise TypeError("dtemp must be either None, a float or a string.")

    return template



def writeAbundance(template, xmol=None, opratio=None):

    # Find the appropriate line.
    for l, line in enumerate(template):
        parsedline = ''.join(line.split())
            if parsedline.startswith('//Abundance'):
                break

    # Add in the appropriate lines
    if xmol is None:
        template.insert(l+1, 'abundance[0] = findvalue(cone, ctwo, cthree, abund);\n')
    else:
        template.insertl+1, 'abundance[0] = %.3f;\n' % float(xmol))
    template.insert(l+2, 'if (abundance[0] < 0.){abundance[0] = 0.;}\n')
    
    # Add in ortho / para correction. 
    if type(opratio) is float:
        template.insert(l+2, 'abundance[0] *= (1. + opratio) / opratio;\n')
    elif type(opratio) is str:
        template.insert(l+2, 'double val[2];\n')
        template.insert(l+3, 'density(x, y, z, &val[2]);\n')
        template.insert(l+4, 'abundance[0] *= val[1]/val[0] + 1.;\n')  
    else:
        raise TypeError("opratio must be either None, a float or a string.")

    return
