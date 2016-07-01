
# Write the velocity structure of the LIME model.c file.

# Write the doppler broadening component.
def writeDopplerBroadening(template, bvalue=0., btype='absolute'):

    # Check that the broadening type is correct.
    if not (btype.lower() is 'absolute' or btype.lower() is 'mach'):
        raise ValueError("btype should be either 'absolute' or 'mach'.")

    # Find the loaction to insert broadening.
    for l, line in enumerate(template):
        parsedline = ''.join(line.split())
        if parsedline.startswith('//Turbulence'):
            break

    # Write appropriate lines.
    if btype.lower() == 'absolute':
        template.insert(l+1, '*doppler = %.3f;\n' % bvalue)
    else:
        template.insert(l+1, 'double val[2];\n')
        template.insert(l+2, 'temperature(x, y, z, &val[2]);\n')
        template.insert(l+3, '%.3f * sqrt(KBOLTZ * val[0] / 2.34 / AMU);\n')

    return


# Write the bulk velocity structure.
# If stellar mass is specified, assume Keplerian rotation.
# If velocityheader is specified, assume the velocity structure is that.
def writeVelocityStructure(template, stellarmass=None, velocityheader=None):

    # Check that one of the two options are specified.
    if (stellarmass is None and velocityheader is None):
        raise ValueError("One of stellarmass or velocityheader must be specified.")

    # Find the location to insert the velocity structure.
    for l, line in enumerate(template):
        parsedline = ''.join(line.split())
        if parsedline.startswith('//Velocity'):
            break
    
    # Write appropriate lines.
    # Note that this is only for (r,z) coordinates at the moment.
    if stellarmass is not None:
        template.insert(l+1, 'velocity[0] = sqrt(6.67e-11 * %.3f * 2e30 / sqrt(x*x + y*y + z*z));\n' % stellarmass)
        template.insert(l+2, 'velocity[0] *= sin(atan2(y,x));\n')
        template.insert(l+3, 'velocity[1] = sqrt(6.67e-11 * %.3f * 2e30 / sqrt(x*x + y*y + z*z));\n' % stellarmass)
        template.insert(l+4, 'velocity[1] *= cos(atan2(y,x));\n')
        template.insert(l+5, 'velocity[2] = 0.0;\n')
    else:
        for ll, lline in enumerate(template):
            parsedline = ''.join(line.split())
            if parsedline.startswith('//Include'):
                template.insert(ll+1, '#include "%s"\n' % velocityheader)
                break
        template.insert(l+1, 'velocity[0] = findvalue(sqrt(x*x + y*y)/AU, fabs(z)/AU, velx);\n')
        template.insert(l+2, 'velocity[1] = findvalue(sqrt(x*x + y*y)/AU, fabs(z)/AU, vely);\n')
        template.insert(l+3, 'velocity[2] = findvalue(sqrt(x*x + y*y)/AU, fabs(z)/AU, velz);\n')

    return
