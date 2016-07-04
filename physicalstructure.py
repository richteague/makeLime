''' Write the physical structure functions. '''

import interpolation as interp

# Write the doppler broadening component.
def writeDopplerBroadening(temp, bvalue=0., btype='absolute', coordsys='cyclindrical', ndim=2):
    
    #Check that the broadening type is correct.
    if not (btype.lower() is 'absolute' or btype.lower() is 'mach'):
    raise ValueError("btype should be either 'absolute' or 'mach'.")
    
    temp.append('void doppler(double x, double y, double z, double *doppler){\n')
    interp.writeCoords(temp, coordsys, ndim)
    
    if btype.lower() == 'absolute':
        # Constant value.
        temp.append('*doppler = %.3f;\n' % bvalue)
    else:
        # Constant mach number.
        temp.append('double val[2];\n')
        temp.append('temperature(x, y, z, &val[2]);\n')
        temp.append('%.5f * sqrt(KBOLTZ * val[0] / 2.34 / AMU);\n')
        
    temp.append('}\n')

    return


# Write the bulk velocity structure.
def writeVelocityStructure(temp, stellarmass=None, coordsys='cyclindrical', ndim=2):

    temp.append('void velocity(double x, double y, double z, double *velocity){\n')
    interp.writeCoords(temp, coordsys, ndim)

    if stellarmass is not None:
        temp.append('velocity[0] = sqrt(6.67e-11 * %.3f * 2e30 / sqrt(x*x + y*y + z*z));\n' % stellarmass)
        temp.append('velocity[0] *= sin(atan2(y,x));\n')
        temp.append('velocity[1] = sqrt(6.67e-11 * %.3f * 2e30 / sqrt(x*x + y*y + z*z));\n' % stellarmass)
        temp.append('velocity[1] *= cos(atan2(y,x));\n')
        temp.append('velocity[2] = 0.0;\n')

    else:
        temp.append('velocity[0] = findvalue(cone, ctwo, cthree, velx);\n')
        temp.append('velocity[1] = findvalue(cone, ctwo, cthree, vely);\n')
        temp.append('velocity[2] = findvalue(cone, ctwo, cthree, velz);\n')
        
    temp.append('}\n')

    return


# Write the main collider density.
def writeDensity(temp, opratio=None, coordsys='cyclindrical', ndim=2): 

    temp.append('void density(double x, double y, double z, *double density){\n')
    interp.writeCoords(temp, coordsys, ndim)
    
    temp.append('density[0] = findvalue(cone, ctwo, cthree, dens);\n')

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


# Write the gas and dust temperatures.
def writeTemperatures(temp, dtemp=None, coordsys='cyclindrical', ndim=2):

    temp.append('void temperature(double x, double y, double z, double *temperature){\n')
    interp.writeCoords(temp, coordsys, ndim)

    temp.append('temperature[0] = findvalue(cone, ctwo, cthree, temp);\n')
    temp.append('if (temperature[0] < 2.73){temperature[0] = 2.73;}\n')

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


# Write the molecular abundances.
def writeAbundance(temp, xmol=None, opratio=None):

    temp.append('void abundance(double x, double y, double z, double *abundance){\n')    
    interp.writeCoords(temp, coordsys, ndim)

    if xmol is None:
        temp.append('abundance[0] = findvalue(cone, ctwo, cthree, abund);\n')
    
    else:
        temp.append('abundance[0] = %.3e;\n' % float(xmol))

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


# Write the gas-to-dust ratio.
def writeGastoDust(temp, coordsys='cyclindrical', ndim=2, g2d=None):

    if g2d is None:
        return
    else 
        temp.append('void gasIIdust(double x, double y, double z, double *gtd){\n')

    interp.writeCoords(temp, coordsys, ndim)

    if type(g2d) is float:
        temp.append('*gtd = %.3e;\n' % g2d)

    elif type(g2d) is str:
        temp.append('*gtd = findvalue(cone, ctwo, cthree, %s);\n' % g2d)

    temp.append('if *gtd < 1){*gtd = 0.;}\n')
    temp.append('}\n')

    return

