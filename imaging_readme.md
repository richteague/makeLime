### Filename Convention

We adopt the naming convention: `modelnumber_theta_phi_trans.fits` in order to help averaging the values at the end. 
`theta` is the inclination,
`phi` is postition angle (not yet!) and `trans` is the transition specified in the LAMDA rate files.

---

### imageblock.py

```python
imageblock(tempfile, nmig, modelnumber, theta, phi, trans, nchan, velres, imgres=None, distance=140., pxls=128, unit=0)
```

`nimg` - Integer of the image block.
`modelnumber` - Number of the model where `modelnumber` <= the number of models to be averaged over.
`theta` - Inclination angle in radians.
`phi` - Rotation angle in radians. Do not use!
`trans` - Transition to trace, note LIME uses zero-indexing while the LAMDA file does not!
`nchan` - Number of channels for the output image.
`velres` - Velocity resolution in [m/s] for the output image.
`imgres = None` - Pixel scaling of the output image in [arcseconds/pixel].
`distance = 140.` - Distance of the observation in [pc].
