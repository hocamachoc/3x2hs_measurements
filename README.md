# 3x2hs_measurements

## Galaxy-galaxy lensing

### Some useful links

#### Y1 Data

- [1xY6xKLUSUlVTjLYSJJO8WnEkV41XS7QK](https://drive.google.com/drive/folders/1xY6xKLUSUlVTjLYSJJO8WnEkV41XS7QK)

#### ggl slides 

- https://docs.google.com/presentation/d/1yJbyqhzFdVgLQmAzzcR6UD5SBFyA39lGLrWM5cFj3IY/edit#slide=id.gb6330c5c3e_0_23

### File description

- ggl_cls-flask.py is the script used to measure the $C_\ell$ from Flask mocks caught from NERSC (path = global/homes/l/ljfaga/nersc_namaster/src/ggl_cls.py)

- ggl_cls-data.py script for measuring data $C_\ell$ based on ggl_cls-flask.py script.

- maps/ contains the maps from data (same maps as those on google drive)

## Current issues

- The data $C_\ell$ are very weird; the E and B modes have the same order of magnitude.
This problem occurs both in ggl cls and shear cls. 
Plenty of tests were performed and the measuring pipeline works in other situations.
Maybe the maps on google drive aren't the final version? 

- I can't reproduce the maps on the drive running the catalog processing pipeline in Nersc
Maybe the code isn't in the final version

### See ggl slides for more details
