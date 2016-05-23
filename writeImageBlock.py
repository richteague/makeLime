    
# Write an image block with the given parameters.
# For the averaging component, we use the naming convention:
#   filename = modelnumber_theta_phi_trans.fits
# Note that model is a string.

def imageblock(template, model, nimg, nchan, velres, trans, pxls, imgres,
               theta, phi, distance, unit):
    
    # Generate the image block to insert. 
    filename = '%s_%.3f_%.3f_%d' % (model, theta, phi, trans)
    
    lines = ['' for i in range(10)]
    lines[0] = 'img[%.0f].nchan = %.0f;' % (nimg, nchan)
    lines[1] = 'img[%.0f].velres = %.0f;' % (nimg, velres)
    lines[2] = 'img[%.0f].trans = %.0f;' % (nimg, trans)
    lines[3] = 'img[%.0f].pxls = %.0f;' % (nimg, pxls)
    lines[4] = 'img[%.0f].imgres = %.3f;' % (nimg, imgres)
    lines[5] = 'img[%.0f].theta = %.3f;' % (nimg, theta) 
    lines[6] = 'img[%.0f].phi = %.3f;' %  (nimg, phi)
    lines[7] = 'img[%.0f].distance = %.1f*PC;' % (nimg, distance)
    lines[8] = 'img[%.0f].unit = %.0f;' % (nimg, unit)
    lines[9] = 'img[%.0f].filename = "%s.fits";' % (nimg, filename)   
    toinsert = ''
    for line in lines:
        toinsert += line + '\n'
        
    # Find the correct place to insert and insert.
    for l, line in enumerate(template):
        parsedline = ''.join(line.split())
        if parsedline.startswith('//%d' % nimg):
            template.insert(l+1, toinsert)
            template.insert(l+2, '// %d\n' % (nimg+1))
    
    return 
    
