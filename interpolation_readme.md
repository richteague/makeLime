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
