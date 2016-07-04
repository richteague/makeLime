# Write an image block with the given parameters.

def imageblock(temp, model, nimg, nchan, velres, trans, pxls, imgres,
               theta, phi, distance, unit):
    
    filename = '%s_%.3f_%.3f_%d' % (model, theta, phi, trans)
    
    lines = ['' for i in range(10)]
    lines[0] = 'img[%.0f].nchan = %.0f;\n' % (nimg, nchan)
    lines[1] = 'img[%.0f].velres = %.0f;\n' % (nimg, velres)
    lines[2] = 'img[%.0f].trans = %.0f;\n' % (nimg, trans)
    lines[3] = 'img[%.0f].pxls = %.0f;\n' % (nimg, pxls)
    lines[4] = 'img[%.0f].imgres = %.3f;\n' % (nimg, imgres)
    lines[5] = 'img[%.0f].theta = %.3f;\n' % (nimg, theta) 
    lines[6] = 'img[%.0f].phi = %.3f;\n' %  (nimg, phi)
    lines[7] = 'img[%.0f].distance = %.1f*PC;\n' % (nimg, distance)
    lines[8] = 'img[%.0f].unit = %.0f;\n' % (nimg, unit)
    lines[9] = 'img[%.0f].filename = "%s.fits";\n' % (nimg, filename)   
    toinsert = ''
    for line in lines:
        temp.insert(line)
    
    return 
    
