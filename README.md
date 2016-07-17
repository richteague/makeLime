# makeLime

Scripts to make running LIME models more efficient. To run:

```python
from makeLime import runLime
runLime.runLime('chemheader.h', 'collisionalrates.dat', 'fileoutname', [0.34], [0], stellarmass=0.6)
```

---

### Chemical Model

As a minimum, a pre-calculated chemical model this must be included as a C header file, `chemheader.h`. The header file must contain the coordinate positions of each point in `c1arr` and `c2arr` where the former is the radial point (in either cylindrical or polar coordinates), and the latter is either the altitude in the cylindrical case, or theta in the polar case. With each point there must be an associated H2 number density, `dens`, gas kinetic temperature `temp` and molecular relative abundance with respect to H2, `abund`. 

If the header file is set up in such a format the scrips should be able to extract all the necessary information automatically.

---

### Coordinate Systems

Currently this can deal with both cylindrical and polar coordinates in both 2D (assuming azimuthal symmetry) or 3D. In include an azimuthal component, this must be done with the `c3arr` array. For the interpolation to work the points must be ordered such that we cycle first through `c2`, then `c1` before finally `c3` if specified. 

---

### Additional Parameters

All inputs which LIME accepts can be included into `chemheader.h` to make a more advanced model. 

__Dust Temperature__: This is controlled by the `dtemp` variable. If no dust temperature is given, equal gas and dust temperatures are assumed. If a float is specified, this is assumed to be a disk-wide rescaling of the gas temperature. Finally, if a string is given, the script assumes this is the array in `chemheader.h` which holds the dust temperatures.

__Doppler Broadening__: Controlled by both `bvalue` and `btype`. By default a 50 m/s, disk-wide value is assumed. `btype` can be either `absolute`, showing `bvalue` gives the linewidth in [m/s], while `mach` means `bvalue` is describing a constant fraction of the local sound speed. For example, `bvalue = 0.3` and `bvalue = 'mach'` assumes turbulent broadening equivalent to 0.3 the local soundspeed. On the other hand, `bvalue = 50.` and `bvalue = 'absolute'` assumes a constant 50 m/s broadening everywhere.


---
---

## averageModels.py

Averages the ensemble of models run if `nmodels > 1`. Has basic options to return noise maps to understand the MC noise.

---

```python
averageModels(nmodels, thetas, phis, transitions, fileout, returnnoise=False, directory='../')
```

All parameters as described in `makemodel.py`.

---

```python
combinePopfiles(nmodels, fileout, directory='../')
```

If multiple models are run, all producing outputfiles, this combines them into a single file.

---

```python
getNoise(nmodels, thetas, phis, transitions, fileout, directory='./')
```

Returns the standard deviation map for the ensemble of models. This is useful to understand how the MC noise affects the data.

---
---

## makemodel.py

Creates `model.c` using the functions specified below. _*NOT YET FINISHED*_

---

```python
makeModelFile(chemheader, moldatfile, thetas, phis, transitions, 
              nchan, velres, pIntensity=1e4, sinkPoints=1e3, dust='jena_thin_e6.tab', 
              antialias=1, sampling=2, outputfile=None, binoutputfile=None, 
              gridfile=None, non_lte=1, imgres=0.035, distance=54., pxls=128, unit=0,
              coordsys='cylindrical', ndim=2, opratio=None, dtemp=None, xmol=None, 
              g2d=None, bvalue=50., btype='absolute', stellarmass=0.6, modelnumber=0)
```
---
---

## physicalstructure.py

Write the functions for the physical model. All functions have to be passed `coordsys` and `ndim` which is used to call the correct interpolation functions. 

---

```python
writeAbundance(tempfile, xmol=None, opratio=None, coordsys='cylindrical', ndim=2)
```

writes the relative abundance of the molecule with respect to the main collider `density[0]` as specified by `writeDensity`. `xmol` can be either a float for homogeneous rescaling of the main collider density, or a string of the name of the array in `chemheader`.

If ortho and para states of the main collider are specified, a correction must be made to the abundance. The same options for `opratio` are as in `writeDensity`.

---

```python
writeDensity(tempfile, opratio=None, coordsys='cylindrical', ndim=2)
```

Writes the main collider density, by deafult assumes this is the H2 density. If `opratio` is specified, then can consider the case for separate oH2 and pH2 densities. If `opratio` is a float we assume that opratio = n(oH2) / n(pH2) and set `density[0]` as n(oH2) and `density[1]` as n(pH2). Alternatively, if `opratio` is a string, we assume that `chemheader` contains the array `dens` for n(oH2) and the array `opratio` for n(pH2).

---

```python
writeDopplerBroadening(tempfile, bval=0, btype='absolute', coordsys='cylindrical', ndim=2)
```

Where `bval` is the value, `btype` is either `'absolute'` for [m/s] or `'mach'` for a factor of the local sound speed.

---

```python
writeGastoDust(tempfile, g2d=None, coordsys='cylindrical', ndim=2)
```

Include the gas-to-dust ratio. If `g2d` is left as None, a standard 100:1 ratio is assumed. 
If `g2d` is a float, use a homoegenous value throughout the disk or specify the array name as a string. `ming2d` is the minimum value it can take.

---

```python
writeTemperatures(tempfile, dtemp=None, coordsys='cylindrical', ndim=2)
```

Includes the gas and dust temperatures. If `dtemp` is `None` then assumes equal gas and dust temperatures. The dust temperature can be assumed to be a homogeneous rescaling of the gas temperature through a float, while a string says that `chemheader` will contain a separate array of the dust temperatures.

---

```python
writeVelocityStructure(tempfile, stellarmass=None, coordsys='cylindrical', ndim=2)
```

Two options are available. If `stellarmass` is specified then Keplerian rotation is used (this is not cyclindrical rotation). If `stellarmass` is not specified, it assumes the `v_x`, `v_y` and `v_z` components are stored in the arrays `velx`, `vely` and `velz` respectively in `chemheader`.

---

---

## imaging.py 

Functions to write the `input` section of the `model.c` file. This requies `thetas`, `phis` and `transitions` to be a list, even if there is only one entry. This is because it will loop over all permutations of these variables.

---

```python
writeImageParameters(temp, radius, minScale, moldatfile, thetas, phis, transitions, nchan,
                     velres, pIntensity=1e4, sinkPoints=1e3, dust='jena_thin_e6.tab', antialias=1,
                     sampling=2, outputfile=None, binoutputfile=None, gridfile=None, non_lte=1,
                     imgres=0.035, distance=54., pxls=128, unit=0, modelnumber=0)
```

The `radius` and `minScale` can be read from the `chemheader` file, see `tobeMade`. All the information about the `LIME` specific values can be found in the [manual](https://lime.readthedocs.io/en/v1.5/usermanual.html).

---

```python
writeImageBlock(tempfile, nimg, modelnumber, theta, phi, trans, nchan, velres, 
                imgres=None, distance=140., pxls=128, unit=0)
```

Writes the specifics for each image block. This should be called with a unique `nimg` value. We adopt the naming convention: `modelnumber_theta_phi_trans.fits` in order to help averaging the values at the end. `theta`, `%.3f`, is the inclination, `phi` is postition angle (not yet!), `%.3f`, and `trans`, `%d`, is the transition specified in the LAMDA rate files.

---

---

## interpolation.py

Writes the interpolation routines for various data structures. Currently we only have 2D cyclindrical coordinates.

---

```python
writeCoords(tempfile, coordsys='cylindrical', ndim=2)
```

For each function in `physicalstructure.py`, writes the appropriate coordinate transform.

---

```python
writeFindValues(tempfile, ncells, coordsys='cylindrical', ndim=2)
```

Copies in the approprirate cell finding and interpolation routines. These are stored in `makeLine/InterpolationRoutines/`. `ncells` can be read in from the `chemheader.h` file.
