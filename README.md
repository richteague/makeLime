# makeLime

Script to make running LIME models more efficient. Assuming the chemical models were converted with the `makeheader.py` routine then to run a simple model:

```python
from makeLime import runlime as lime
lime.run(name='modelname', 
         headerfile='chemheader.h',
         ratefile='collisionalrates.dat', 
         inclinations=[0.34],
         transitions=[3,2,1])
```

This will run a model based on the physical structure stored in `chemheader.h`, assuming the collisional rates `collisionalrates.dat` (which are stored in the `AuxFiles` folder), at an inclination of 0.34 radians, and for the transitions 3, 2 and 1 as listed in the `collisionalrates.dat` file. 

---

### Chemical Model

As a minimum, a pre-calculated chemical model this must be included as a C header file, `chemheader.h`. The header file must contain the (currently) 2D coordinate positions of each point in `c1arr` and `c2arr` where the former is the radial point (in either cylindrical or polar coordinates), and the latter is either the altitude in the cylindrical case, or theta in the polar case. With each point there must be an associated H2 number density, `dens`, gas kinetic temperature `temp` and molecular relative abundance with respect to H2, `abund`. 

If the header file is set up in such a format the scrips should be able to extract all the necessary information automatically. Using the `makeheader.py` routine, this is done automatically:

```
>>> python makeheader.py chemicalmodel.out
```

---

### Additional LIME Parameters

The example above is minimal. All parameters which LIME accepts are also accepted as arguments in `runLime.run`. Read the [LIME manual](https://lime.readthedocs.io/en/v1.5/) for an idea of what is possible. 

When called, `runLIME` will search the chemical header file for arrays containing the density, gas temperature, dust temperature, molecular relative abundance, gas-to-dust ratio and turbulent velocity structures. Unless directed, it will assume the default names: `dens`, `temp`, `dtemp`, `abund`, `gastodust` and `doppler` respectively. 

If the specified array name cannot be found for `dtemp`, `gastodust` or `doppler`, the value will revert to their default values. These are equal gas and dust temperatures, a model wide gas-to-dust ratio of 100 and Doppler broadening of 40 m/s respectively.

In addition to array names and default settings, `dtemp`, `abund`, `gastodust` and `doppler` can also be a float. This specifies a disk-wide value for that property. For example, this can be used to compare a CO structure from a chemical model and one where x(CO) = 1e-4 throughout the disk by setting `abund=1e-4`.

For the turbulent motions set by `doppler`, the default `dopplertype` is assumed to be `absolute`, meaning the Doppler broadening is given in meters per second. Alternatively this can be `mach`, where it is a fraction of the local soundspeed. The bulk motion of the gas is assumed to be Keplerian rotation around a central star of mass specified with `stellarmass` in solar masses.

### Averaging Models

If `nmodels > 1`, then several models are averaged together. This creates a smoother outer disk which usually suffers from gridding noise. If `returnnoise` is set, the standard devation of each voxel through the ensemble of models is returned as an additional file. If `outputfile` or `binoutputfile` are selected, these will be combined if more than one model is run.



