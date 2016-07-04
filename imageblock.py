# Write an image block with the given parameters.

def writeImageParameters(temp, modelnumber, thetas, phis, transitions,
                         nchan, velres, imgres=0.035, distance=54., pxls=128, unit=0):


    temp.append('void input(inputPars *par, image *img){\n')

    # For each permutation, write an image block.
    for i, theta in enumerate(thetas):
        for p, phi in enumerate(phis):
            for t, trans in enumerate(transitions):
                nimg = int(i*len(transitions) + p*len(phis) + t)
                    writeImageBlock(temp, modelnumber, theta, phi, trans, nchan,
                                    velres, imgres, distance, pxls, unit)


    temp.append('}\n\n')

    return


def writeImageBlock(temp, nimg, modelnumber, theta, phi, trans, nchan, velres.
                    imgres=0.035, distance=54., pxls=128, unit=0):
   
    # Filename following convention. 
    filename = '%s_%.3f_%.3f_%d' % (model, theta, phi, trans)

    temp.append('img[%.0f].nchan = %.0f;\n' % (nimg, nchan))
    temp.append('img[%.0f].velres = %.0f;\n' % (nimg, velres))
    temp.append('img[%.0f].trans = %.0f;\n' % (nimg, trans))
    temp.append('img[%.0f].pxls = %.0f;\n' % (nimg, pxls))
    temp.append('img[%.0f].imgres = %.3f;\n' % (nimg, imgres))
    temp.append('img[%.0f].theta = %.3f;\n' % (nimg, theta))
    temp.append('img[%.0f].phi = %.3f;\n' %  (nimg, phi))
    temp.append('img[%.0f].distance = %.1f*PC;\n' % (nimg, distance))
    temp.append('img[%.0f].unit = %.0f;\n' % (nimg, unit))
    temp.append('img[%.0f].filename = "%s.fits";\n' % (nimg, filename))
    temp.append('\n\n') 
    
    return 
    
