# Write an image block with the given parameters.

def writeImageParameters(temp, radius, minScale, moldatfile, thetas, phis, transitions, nchan,
                         velres, pIntensity=1e4, sinkPoints=1e3, dust='jena_thin_e6.tab', antialias=1,
                         sampling=2, outputfile=None, binoutputfile=None, gridfile=None, lte_only=1,
                         imgres=0.05, distance=54., pxls=128, unit=0, modelnumber=0):


    temp.append('void input(inputPars *par, image *img){\n\n')

    # Include the radiative transfer parameters.
    temp.append('\tpar->radius = %.1f*AU;\n' % radius)
    temp.append('\tpar->minScale = %.2f*AU;\n' % minScale)
    temp.append('\tpar->pIntensity = %.0f;\n' % pIntensity)
    temp.append('\tpar->sinkPoints = %.0f;\n' % sinkPoints)
    temp.append('\tpar->dust = "%s";\n' % dust)
    temp.append('\tpar->moldatfile[0] = "%s";\n' % moldatfile)
    temp.append('\tpar->antialias = %d;\n' % antialias)
    temp.append('\tpar->sampling = %d;\n' % sampling)
    temp.append('\tpar->lte_only = %d;\n' % lte_only)
    if outputfile is not None:
        temp.append('\tpar->outputfile = "%s";\n' % outputfile)
    if binoutputfile is not None:
        temp.append('\tpar->binoutputfile = "%s";\n' % binoutputfile)
    if gridfile is not None:
        temp.append('\tpar->gridfile = "%s";\n' % gridfile)
    temp.append('\n')


    # For each permutation, write an image block.
    for i, theta in enumerate(thetas):
        for p, phi in enumerate(phis):
            for t, trans in enumerate(transitions):
                nimg = int(i*len(thetas)*len(phis) + p*len(phis) + t)
                writeImageBlock(temp, nimg, modelnumber, theta, phi, trans, nchan,
                                velres, imgres, distance, pxls, unit)


    temp.append('}\n\n\n')

    return


def writeImageBlock(temp, nimg, modelnumber, theta, phi, trans, nchan, velres,
                    imgres=0.035, distance=54., pxls=128, unit=0):
   
    # Filename following convention. 
    filename = '%s_%.3f_%.3f_%d' % (modelnumber, theta, phi, trans)

    temp.append('\timg[%.0f].nchan = %.0f;\n' % (nimg, nchan))
    temp.append('\timg[%.0f].velres = %.0f;\n' % (nimg, velres))
    temp.append('\timg[%.0f].trans = %.0f;\n' % (nimg, trans))
    temp.append('\timg[%.0f].pxls = %d;\n' % (nimg, pxls))
    temp.append('\timg[%.0f].imgres = %.3f;\n' % (nimg, imgres))
    temp.append('\timg[%.0f].theta = %.3f;\n' % (nimg, theta))
    temp.append('\timg[%.0f].phi = %.3f;\n' %  (nimg, phi))
    temp.append('\timg[%.0f].distance = %.1f*PC;\n' % (nimg, distance))
    temp.append('\timg[%.0f].unit = %d;\n' % (nimg, unit))
    temp.append('\timg[%.0f].filename = "%s.fits";\n' % (nimg, filename))
    temp.append('\n') 
    
    return 
    
