import fileinput
import numpy as np

# All the functions to write a LIME file.

def valsfromheader(headername):

    # Assumes that rvals is the first line and zvals is the second.

    with open('../'+headername) as f:
        header = f.readlines()

    i = 18
    while header[0][i] != ']':
        if header[0][i] == '[':
            j = i
        i += 1
    ncells= int(header[0][j+1:i])
    while header[0][i] != ',':
        if header[0][i] == '{':
            j = i
        i += 1
    rin = float(header[0][j+1:i])
    i = -2
    while header[0][i] != ',':
        i -= 1
    rout = float(header[0][i+2:-3])
    i = -2
    while header[1][i] != ',':
        i -= 1
    rout = np.hypot(float(header[1][i+2:-3]), rout)

    return rin, rout, ncells

def arrsfromheader(headername):
    with open('../'+headername) as f:
        header = f.readlines()
    return np.array([parsename(line) for line in header])

def parsename(line):
    i = 0
    while line[i] != '[':
        i += 1
    j = i
    while line[i] != ' ':
        i -= 1
    return line[i+1:j]

def imageblock(model, nimg, nchan, velres, trans, pxls, imgres, theta,
               phi, distance, unit):
    
    # Write an image block with the given parameters.
    # Make sure that 'nimg' is unique per model file!
    filename = '%s_%.3f_%d' % (model, theta, trans)
    
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

# Write a model.c file for LIME. 
def generateModelFile(chemheader, model, transitions, stellarmass,
                      mach, pIntensity, sinkPoints,
                      dustfile, antialias, lte_only, blend, nchan,
                      velres, pxls, imgres, thetas, phi, distance, unit,
                      collisionfile, orthoratio=None,
                      modelfile='model_template.c', equaltemp=True):
                      
    # Read in the template with which to generate the model.c file.
    # Make sure this template is in your directory.
    with open(modelfile) as f:
        template = f.readlines()
    
    # Determin the inner and outer radii and the number of cells.
    rin, rout, ncells = valsfromheader(chemheader)

    # If appropriate, include the difference in dust and gas temperatures.
    if 'dtemp' in arrsfromheader(chemheader) and not equaltemp:
        print 'Assuming different dust and gas temperatures.'
        
        dstring = ('''temperature[1]=findvalue(hypot(x,y)/AU,fabs(z)/AU,dtemp);\n
                   if(temperature[1]<2.7){\n
                   temperature[1] = 2.7;\n
                   }\n''') 
     
        for l, line in enumerate(template):
            parsedline = ''.join(line.split())
            if parsedline.startswith('//DustTemperature'):
                template.insert(l+1, dstring)
                break

    else:
        print 'Assuming equal dust and gas temperatures.'

    # Include the gas to dust mass ratio.
    for l, line in enumerate(template):
        parsedline = ''.join(line.split())
        if parsedline.startswith('//GasToDust'):
            if 'g2d' in arrsfromheader(chemheader):
                gtd = ('''*gtd=findvalue(hypot(x,y)/AU,fabs(z)/AU,g2d);\n
                         if(*gtd<50.){\n
                         *gtd = 50.;\n
                         }\n''')
            else:
                gtd = ('''*gtd = 100.;\n''')
            template.insert(l+1, gtd)    
            break
    
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
            
    # Define the mach number for the turbulent broadening.
    for l, line in enumerate(template):
        parsedline = ''.join(line.split())
        if parsedline.startswith('//Mach'):
            template.insert(l+1, '#define MACH %.3f\n' % mach)
            break
    
    # Assume no o:p ratio and a single H2 density.
    if orthoratio is None:
        print 'Assuming {dens} is total H2 density.'
            
    # Assume a homogenous o:p ratio and a single total H2 density.                  
    elif type(orthoratio) is float:
        print ('''Assuming {dens} is total H2 density
                  with a %.2f:1.00 o/p ratio.''' 
                  % (orthoratio/(1.-orthoratio)))
         
        pH2string = ('''density[1] = (1. - ORTHORATIO) * density[0];\n
                       density[0] *= ORTHORATIO;\n''')

        for l, line in enumerate(template):
            parsedline = ''.join(line.split())
            if parsedline.startswith('//OrthoRatio'):
                template.insert(l+1, '#define ORTHORATIO %.3f\n' % orthoratio)
                break

        for l, line in enumerate(template):
            parsedline = ''.join(line.split())
            if parsedline.startswith('//ParaH2Density'):
                template.insert(l+1, pH2string)
                break
                
        for l, line in enumerate(template):
            parsedline = ''.join(line.split())
            if parsedline.startswith('//OrthoCorrection'):
                template.insert(l+1, 'abundance[0] /= ORTHORATIO;\n')
                break
                 
    # Assume two separate densities, oH2 and pH2.
    elif type(orthoratio) is str:
        print ('''Assuming {dens} is oH2 density and {%s} is pH2 density.''' 
               % orthoratio)

        pH2string = ('''density[1]=findvalue(hypot(x,y)/AU,fabs(z)/AU,%s);\n
                       if(density[1]<1e-30){\n
                       density[1] = 1e-30;\n
                       }\n''' % orthoratio)
                       
        for l, line in enumerate(template):
            parsedline = ''.join(line.split())
            if parsedline.startswith('//ParaH2Density'):
                template.insert(l+1, pH2string)
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
                    template.insert(l+1, imageblock(model,
                                                    nimg, 
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
    if model[-2:] == '.c':
        model = model[:-2]
    with open('%s.c' % model, "w") as tempfile:
        for line in template:
            tempfile.write("%s" % line)       
                                                
    return
