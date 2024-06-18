## What are these files for?
* triso_build.py --- constructing the geometry, settings, tallies, and materials.xml files. This file will also generate a plot of the geometry from 'xy' viewpoint and at z=0, saved at geometry.png
* zernike_animation.py --- generating a GIF of the Zernike tally. NOTE: This file is rudimentary and partially hardcoded. This file requires there to be only one tally in the simulation. Further, the dimensions of the compact must be respecified. Also, this file will take a decently long time, and will generate a large number of warnings from cell id conflicts, which I am not entirely sure as to why. 

## How to run
1. python triso_build.py 
2. openmc .   OR   nohup mpirun -np 16 openmc . >&logfile&
2.1 the second is for running in the background so you can exit the terminal
3. python zernike_animation.py
4. Done!
