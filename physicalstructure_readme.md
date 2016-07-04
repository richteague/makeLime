### velocitystructure.py

##### Doppler Broadening
```python
writeDopplerBroadening(tempfile, bval=0, btype='absolute')
```

Where `bval` is the value, `btype` is either `'absolute'` for [m/s] or `'mach'` for a factor of the local sound speed.


##### Velocity
```python
writeVelocityStructure(tempfile, stellarmass=None)
```

Two options are available. If `stellarmass` is specified then Keplerian rotation is used (this is not cyclindrical rotation).
If `stellarmass` is not specified, it assumes the `v_x`, `v_y` and `v_z` components are stored in the arrays `velx`, `vely` and `velz`
respectively in `chemheader`.

---

### physicalstructure.py

##### Density
```python
writeDensity(tempfile, opratio=None)
```

Writes the main collider density, by deafult assumes this is the H2 density.
If `opratio` is specified, then can consider the case for separate oH2 and pH2 densities. 
If `opratio` is a float we assume that opratio = n(oH2) / n(pH2) and set `density[0]` as n(oH2) and `density[1]` as n(pH2).
Alternatively, if `opratio` is a string, we assume that `chemheader` contains the array `dens` for n(oH2) and the array `opratio` for n(pH2).

##### Temperature
```python
writeTemperatures(tempfile, dtemp=None)
```

Includes the gas and dust temperatures. If `dtemp` is `None` then assumes equal gas and dust temperatures.
The dust temperature can be assumed to be a homogeneous rescaling of the gas temperature through a float, while a string says
that `chemheader` will contain a separate array of the dust temperatures.

##### Abundance
```python
writeAbundance(tempfile, xmol=None, opratio=None)
```

writes the relative abundance of the molecule with respect to the main collider `density[0]` as specified by `writeDensity`.
`xmol` can be either a float for homogeneous rescaling of the main collider density, or a string of the name of the array in
`chemheader`.

If ortho and para states of the main collider are specified, a correction must be made to the abundance. The same options for
`opratio` are as in `writeDensity`.


##### Gas to Dust Ratio
```python
writeGastoDust(tempfile, g2d=None)
```

Include the gas-to-dust ratio. If `g2d` is left as None, a standard 100:1 ratio is assumed. 
If `g2d` is a float, use a homoegenous value throughout the disk or specify the array name as a string.
