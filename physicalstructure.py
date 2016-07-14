''' Write the physical structure functions. '''

import interpolation as interp


# Write the molecular abundances.
def writeAbundance(temp, xmol=None, opratio=None, coordsys='cylindrical', ndim=2):

    temp.append('void abundance(double x, double y, double z, double *abundance) {\n\n')    
    interp.writeCoords(temp, coordsys, ndim)

    if xmol is None:
        temp.append('\tabundance[0] = findvalue(c1, c2, c3, abund);\n') 
    else:
        temp.append('\tabundance[0] = %.3e;\n' % float(xmol))

    if type(opratio) is int:
        opratio = float(opratio)

    if type(opratio) is float:
        temp.append('\tabundance[0] *= (1. + opratio) / opratio;\n\n')

    elif type(opratio) is str:
        temp.append('\tdouble val[2];\n')
        temp.append('\tdensity(x, y, z, &val[2]);\n')
        temp.append('\tabundance[0] *= val[1]/val[0] + 1.;\n\n')  

    elif opratio is not None:
        raise TypeError("opratio must be either None, a float or a string.")

    temp.append('\tif (abundance[0] < 0.){\n\t\tabundance[0] = 0.;\n\t}\n')  
    temp.append('\n}\n\n\n')

    return

# Write the main collider density.
def writeDensity(temp, opratio=None, coordsys='cylindrical', ndim=2): 

    temp.append('void density(double x, double y, double z, double *density) {\n\n')
    interp.writeCoords(temp, coordsys, ndim)
    
    temp.append('\tdensity[0] = findvalue(c1, c2, c3, dens);\n')

    if type(opratio) is int:
        opratio = float(opratio)

    if type(opratio) is float:
        temp.append('\tdensity[0] *= (opratio / (1. + opratio));\n') 
        temp.append('\tdensity[1] = findvalue(c1, c2, c3, dens);\n')
        temp.append('\tdensity[1] /= (1. + opratio);\n') 
        temp.append('\tif (density[1] < 1e-30) {\n\tdensity[1] = 1e-30;\n\t}\n\n')

    elif type(opratio) is str:
        temp.append('\tdensity[1] = findvalue(c1, c2, c3, %s);\n' % opratio)
        temp.append('\tif (density[1] < 1e-30) {\n\t\tdensity[1] = 1e-30;\n\t}\n\n')

    elif opratio is not None:
        raise TypeError("opratio must be either None, a float or a string.")

    temp.append('\tif (density[0] < 1e-30) {\n\t\tdensity[0] = 1e-30;\n\t}\n\n')
    temp.append('}\n\n\n')

    return 


# Write the doppler broadening component.
def writeDopplerBroadening(temp, bvalue=0., btype='absolute', coordsys='cylindrical', ndim=2):

    #Check that the broadening type is correct.
    if not (btype is 'absolute' or btype is 'mach'):
        raise ValueError("btype should be either 'absolute' or 'mach'.")
    
    temp.append('void doppler(double x, double y, double z, double *doppler) {\n\n')
    interp.writeCoords(temp, coordsys, ndim)
    
    if btype.lower() == 'absolute':
        # Constant value.
        temp.append('\t*doppler = %.1f;\n' % bvalue)
    else:
        # Constant mach number.
        temp.append('\tdouble val[2];\n')
        temp.append('\ttemperature(x, y, z, &val[2]);\n')
        temp.append('\t%.2f * sqrt(KBOLTZ * val[0] / 2.34 / AMU);\n')
        
    temp.append('\n}\n\n\n')

    return


# Write the gas-to-dust ratio.
def writeGastoDust(temp, g2d=None, ming2d=1., coordsys='cylindrical', ndim=2):

    if g2d is None:
        return
    else: 
        temp.append('void gasIIdust(double x, double y, double z, double *gtd) {\n\n')

    if (type(g2d) is float or type(g2d) is int):
        temp.append('\t*gtd = %.1f;\n\n' % g2d)

    elif type(g2d) is str:
        interp.writeCoords(temp, coordsys, ndim)
        temp.append('\t*gtd = findvalue(c1, c2, c3, %s);\n\n' % g2d)

    temp.append('\tif (*gtd < %.1f) {\n\t\t*gtd = %.1f;\n\t}\n' % (ming2d, ming2d))
    temp.append('}\n\n\n')

    return


# Write the gas and dust temperatures.
def writeTemperatures(temp, dtemp=None, coordsys='cylindrical', ndim=2):

    temp.append('void temperature(double x, double y, double z, double *temperature) {\n\n')
    interp.writeCoords(temp, coordsys, ndim)

    temp.append('\ttemperature[0] = findvalue(c1, c2, c3, temp);\n')
    temp.append('\tif (temperature[0] < 2.73) {\n\t\ttemperature[0] = 2.73;\n\t}\n\n')

    if type(dtemp) is float:
        temp.append('\ttemperature[1] = %.3f * temperature[0];\n' % dtemp)
        temp.append('\tif (temperature[1] < 2.73){temperature[1] = 2.73;}\n')

    elif type(dtemp) is str:
        temp.append('\ttemperature[1] = findvalue(c1, c2, c3, %s);\n' % dtemp)
        temp.append('\tif (temperature[1] < 2.73) {\n\t\ttemperature[1] = 2.73;\n\t}\n')

    elif dtemp is not None:
        raise TypeError("dtemp must be either None, a float or a string.")

    temp.append('}\n\n\n')

    return


# Write the bulk velocity structure.
def writeVelocityStructure(temp, stellarmass=None, coordsys='cylindrical', ndim=2):

    temp.append('void velocity(double x, double y, double z, double *velocity) {\n\n')

    if stellarmass is not None:
        temp.append('\tvelocity[0] = sqrt(6.67e-11 * %.3f * 2e30 / sqrt(x*x + y*y + z*z));\n' % stellarmass)
        temp.append('\tvelocity[0] *= sin(atan2(y,x));\n')
        temp.append('\tvelocity[1] = sqrt(6.67e-11 * %.3f * 2e30 / sqrt(x*x + y*y + z*z));\n' % stellarmass)
        temp.append('\tvelocity[1] *= cos(atan2(y,x));\n')
        temp.append('\tvelocity[2] = 0.0;\n')

    else:
        temp.append('\tvelocity[0] = findvalue(c1, c2, c3, velx);\n')
        temp.append('\tvelocity[1] = findvalue(c1, c2, c3, vely);\n')
        temp.append('\tvelocity[2] = findvalue(c1, c2, c3, velz);\n')
        
    temp.append('\n}\n\n\n')

    return
