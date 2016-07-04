## imaging.py 

Functions to write the `input` section of the `model.c` file. This requies `thetas`, `phis` and `transitions` to be a list, even if there is only one entry. This is because it will loop over all permutations of these variables.

---

```python
writeImageParameters(temp, radius, minScale, moldatfile, thetas, phis, transitions, nchan,
                     velres, pIntensity=1e4, sinkPoints=1e3, dust='jena_thin_e6.tab', antialias=1,
                     sampling=2, outputfile=None, binoutputfile=None, gridfile=None, non_lte=1,
                     imgres=0.035, distance=54., pxls=128, unit=0, modelnumber=0)
```

The `radius` and `minScale` can be read from the `chemheader` file, see [writeSomething](google.com) All the information about the `LIME` specific values can be found in the [manual](https://lime.readthedocs.io/en/v1.5/usermanual.html).

---

```python
writeImageBlock(tempfile, nimg, modelnumber, theta, phi, trans, nchan, velres, 
                imgres=None, distance=140., pxls=128, unit=0)
```

Writes the specifics for each image block. This should be called with a unique `nimg` value. We adopt the naming convention: `modelnumber_theta_phi_trans.fits` in order to help averaging the values at the end. `theta`, `%.3f`, is the inclination, `phi` is postition angle (not yet!), `%.3f`, and `trans`, `%d`, is the transition specified in the LAMDA rate files.
