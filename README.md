# makeLime

## Progress.

`writeLIME.py` - Can write and insert an arbitrary number of image blocks and an input parameters block into a given model_template.c file.

`writeHEADER.py` - Can read in an ALCHEMIC model file and write it as a C header file. Need to expand this to write an arbitrary set of arrays. Potentially make this so that it writes everything and just assumes default values if they are not specificed. This would save on having several template files.

`runLIME.py` - Need to write this to run multiple models simultaneously and then save some statistics about the averaging.
