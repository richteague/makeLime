# All the functions to write a LIME file.


def imageblock(nimg, nchan, velres, trans, pxls, imgres, theta,
               phi, distance, unit, filename=None):

    # Write an image block with the given parameters.
    # Make sure that 'nimg' is unique per model file!

    if filename is None:
        filename = 'model_%d' % nimg
    elif filename[-5:] == '.fits':
        filename = filename[:-5]

    lines = ['' for i in range(10)]
    lines[0] = 'img[%.0f].nchan = %.0f;' % (nimg, nchan)
    lines[1] = 'img[%.0f].velres = %.0f;' % (nimg, velres)
    lines[2] = 'img[%.0f].trans = %.0f;' % (nimg, trans)
    lines[3] = 'img[%.0f].pxls = %.0f;' % (nimg, pxls)
    lines[4] = 'img[%.0f].imgres = %.3f;' % (nimg, imgres)
    lines[5] = 'img[%.0f].theta = %.3f;' % (nimg, theta)
    lines[6] = 'img[%.0f].phi = %.3f;' % (nimg, phi)
    lines[7] = 'img[%.0f].distance = %.1f*PC;' % (nimg, distance)
    lines[8] = 'img[%.0f].unit = %.0f;' % (nimg, unit)
    lines[9] = 'img[%.0f].filename = "%s.fits";' % (nimg, filename)
    toinsert = ''
    for line in lines:
        toinsert += line + '\n'
    return toinsert


def inputparameters(pIntensity, sinkPoints, dustfile, molfile, antialias,
                    lte_only, blend, rin, rout):

    # Write the input parameters block for the model file.
    # rin and rout are the minimum and maximum (polar) radial points
    # from the chemical model.

    lines = ['' for i in range(10)]
    lines[0] = 'par->radius = %.2f*AU;' % (rout)
    lines[1] = 'par->minScale = %.2f*AU;' % (rin)
    lines[2] = 'par->pIntensity = %.0f;' % (pIntensity)
    lines[3] = 'par->sinkPoints = %.0f;' % (sinkPoints)
    lines[4] = 'par->dust = "%s";' % (dustfile)
    lines[5] = 'par->moldatfile[0] = "%s";' % (molfile)
    lines[6] = 'par->antialias = %.0f;' % (antialias)
    lines[7] = 'par->sampling = 2;'
    lines[8] = 'par->lte_only = %.0f;' % (lte_only)
    lines[9] = 'par->blend = %.0f;' % (blend)

    toinsert = ''
    for line in lines:
        toinsert += line + '\n'

    return toinsert


def generateModelFile(chemheader, fileout, transitions, stellarmass, rin, rout,
                      ncells, pIntensity, sinkPoints, dustfile, antialias,
                      lte_only, blend, nchan, velres, pxls, imgres, thetas,
                      phi, distance, unit, collisionfile,
                      modelfile='model_template.c'):

    # Read in the template with which to generate the model.c file.
    # Make sure this template is in your directory.
    with open(modelfile) as f:
        template = f.readlines()

    # Write the stellar mass for the Keplerian rotation.
    if chemheader[-2:] == '.h':
        chemheader = chemheader[:-2]
    for l, line in enumerate(template):
        parsedline = ''.join(line.split())
        if parsedline.startswith('//ChemicalModel'):
            template.insert(l+1, '#include "%s.h"' % chemheader)
            break

    # Write the stellar mass for the Keplerian rotation.
    for l, line in enumerate(template):
        parsedline = ''.join(line.split())
        if parsedline.startswith('//StellarMass'):
            template.insert(l+1, '#define MSTAR %.2f' % stellarmass)
            break

    # Include the number of grid cells for the findcell() function.
    for l, line in enumerate(template):
        parsedline = ''.join(line.split())
        if parsedline.startswith('//GridCells'):
            template.insert(l+1, '#define NCELLS %.d\n' % ncells)
            break

    # Radiative transfer properties.
    for l, line in enumerate(template):
        parsedline = ''.join(line.split())
        if parsedline.startswith('//InputParameters'):
            template.insert(l+1, inputparameters(pIntensity,
                                                 sinkPoints,
                                                 dustfile,
                                                 collisionfile,
                                                 antialias,
                                                 lte_only,
                                                 blend,
                                                 rin,
                                                 rout))
            break

    # For each transition and inclination specified, add in an image block.
    for i, theta in enumerate(thetas):
        for t, trans in enumerate(transitions):
            nimg = i*len(transitions) + t
            for l, line in enumerate(template):
                parsedline = ''.join(line.split())
                if parsedline.startswith('//%0.f' % nimg):
                    template.insert(l+1, imageblock(nimg,
                                                    nchan,
                                                    velres,
                                                    trans,
                                                    pxls,
                                                    imgres,
                                                    theta,
                                                    phi,
                                                    distance,
                                                    unit))
                    template.insert(l+2, '// %.0f\n' % (nimg+1))
                    break

    # Save as a model.c file to run.
    # Remove the '.c' if it has been added.
    if fileout[-2:] == '.c':
        fileout = fileout[:-2]
    with open('%s.c' % fileout, "w") as tempfile:
        for line in template:
            tempfile.write("%s" % line)

    return
