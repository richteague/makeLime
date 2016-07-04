### Filename Convention

We adopt the naming convention: `modelnumber_theta_phi_trans.fits` in order to help averaging the values at the end. 
`theta` is the inclination,
`phi` is postition angle (not yet!) and `trans` is the transition specified in the LAMDA rate files.

---

### imageblock.py

```python
imageblock(tempfile, nimg, modelnumber, theta, phi, trans, nchan, velres, imgres=None, distance=140., pxls=128, unit=0)
```

Writes the specifics for the image block. For more information about the `LIME` specific values, read the [manual](https://lime.readthedocs.io/en/v1.5/usermanual.html). This should be called with a unique `nimg` value. Unless both `imgres` and `distance` are set, they are adjusted so that the image fills the given pixel size.


`nimg` - Integer of the image block.

`modelnumber` - Number of the model where `modelnumber` <= the number of models to be averaged over.

`theta` - Inclination angle in radians.

`phi` - Rotation angle in radians. Do not use!

`trans` - Transition to trace, note LIME uses zero-indexing while the LAMDA file does not!

`nchan` - Number of channels for the output image.

`velres` - Velocity resolution in [m/s] for the output image.

`imgres = None` - Pixel scaling of the output image in [arcseconds/pixel].

`distance = 140.` - Distance of the observation in [pc].

`pxls = 128` - Number of pixels per side of the output image.

`unit = 0` - Unit of the output, 0 = Kelving, 1 - Jy/pix.
