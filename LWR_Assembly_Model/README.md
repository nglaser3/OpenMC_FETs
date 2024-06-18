## Contents of Folder
This folder contains .py files to construct a rudimentary LWR fuel assembly, construct a zernike animation for a single fuel rod, and a script to run everything. 

## How to Run
- First, edit the input parameters (if desired) at the very top of assembly_build.py
- If you want everything run all at once, and to open the .mp4 file:
	- In the command line, run: source run.sh
- Otherwise:
	- rm *.xml *.h5 *.out --- this removes all lingering files used / outputted by previous simulation
	- python *build.py --- this constructs the xml files for the openmc sim
	- openmc . --- this runs openmc using the .xml files
	- python *plot.py --- this constructs the .mp4 video of the zernike expansion
	- Done! To open the .mp4 file you can use either open XXX.mp4, where XXX is the name of the file, or simply open with your systems finder

## Work to do in future
- Add a seperate input.py file that only containts variable definitions so no editting is done to script files themselves  
