
# Write the velocity structure of the LIME model.c file.

# Write the doppler broadening component.
def writeDopplerBroadening(temp, bvalue=0., btype='absolute'):

    # Check that the broadening type is correct.
    if not (btype.lower() is 'absolute' or btype.lower() is 'mach'):
        raise ValueError("btype should be either 'absolute' or 'mach'.")

    temp.append('void doppler(double x, double y, double z, double *doppler){\n')

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


def writeVelocityStructure(temp, stellarmass=None, velocityheader=None):

    # Check that one of the two options are specified.
    if (stellarmass is None and velocityheader is None):
        raise ValueError("One of stellarmass or velocityheader must be specified.")

    temp.append('void velocity(double x, double y, double z, double *velocity){\n')   
 
    if stellarmass is not None:
        # Parametric Keplerian rotation.
        temp.append('velocity[0] = sqrt(6.67e-11 * %.3f * 2e30 / sqrt(x*x + y*y + z*z));\n' % stellarmass)
        temp.append('velocity[0] *= sin(atan2(y,x));\n')
        temp.append('velocity[1] = sqrt(6.67e-11 * %.3f * 2e30 / sqrt(x*x + y*y + z*z));\n' % stellarmass)
        temp.append('velocity[1] *= cos(atan2(y,x));\n')
        temp.append('velocity[2] = 0.0;\n')
    else:
        # Provided velocity structure.
        temp.append('velocity[0] = findvalue(cone, ctwo, cthree, velx);\n')
        temp.append('velocity[1] = findvalue(cone, ctwo, cthree, vely);\n')
        temp.append('velocity[2] = findvalue(cone, ctwo, cthree, velz);\n')

    temp.append('}\n')

    return
