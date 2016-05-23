
# Write the input parameters block for the model file.
# rin and rout are the minimum and maximum (polar) radial points
# from the chemical model.  
    
def inputparameters(pIntensity, sinkPoints, dustfile, molfile, antialias,
                    lte_only, blend, rin, rout, popfile=False):

    # Generate the parameter block to insert.
    if popfile:
        lines = ['' for i in range(11)]
    else:
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
    
    if popfile:
        lines[10] = 'par->outputfile = "popfile";'
    
    toinsert = ''
    for line in lines:
        toinsert += line + '\n'
    
    # Find the correct place and insert.
    for l, line in enumerate(template):
        parsedline = ''.join(line.split())
        if parsedline.startswith('//InputParameters'):
            template.insert(l+1, toinsert)
            break
    
    return 
