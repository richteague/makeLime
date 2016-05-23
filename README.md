# makeLime

Scripts to make running LIME models more efficient.

## Current Scripts.

`scriptLIME.py` -
Script that brings it all together. For a given number of models and parameters, will generate the model files with `writeLIME`, run them, if appropriate, average them with `averageModels`, and finally move them to a given directory. 

`averageModels.py` -
If more than one models are run, it will average over them to increase the resolution in the outer regions of the disk. If requested, it will also produce .fits files containing the standard deviation of each pixel for the array of models run. Once averaged, it will move the files to the given directory.

`writeLIME.py` -
Given a C header file containing the chemical model and the required parameters, will create a model.c file to be run. It can create several image blocks for lists of inclinations, postition angles and transitions.
__TO DO:__ Need to account for different ortho- and para-H2 colliders. 

`writeImageBlock.py` -
Called by `writeLIME` to create an image block for each permutation of inclination, position angle ane transition.

`writeParameters.py` -
Called by `writeLIME` to write the radiative transfer parameters.

`writeHeader.py` -
Tools to work with the C headers. Has a script to convert `ALCHEMIC` output files into C headers and functions to read values and names which is needed for the other scripts.
__TO DO:__ Make a more general header writing function.

`runLIME.py` -
Need to write this to run multiple models simultaneously and then save some statistics about the averaging.
