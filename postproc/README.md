-------------------------------------------------------------------

### How to for the post processing of the NEWA production runs

Date: 25.10.2019

Authors:  The NEWA Consortium

-----------------------------------------------------------------

Required python3 modules that need to be installed a priori:
numpy, netCDF4, pyproj, wrf

Needed files: 	`newa2timeseries.py` (main script)
		`newa2nc.py` (module loaded inside main script)

usage:

`python3 newa2timeseries.py 	$PATHTOWRFOUTFILES`

**Other comments:**

* The setup assumes daily output files from WRF with using three nested domains.
Just the innermost domain is post-processed (D03). I.e. the script takes 
(wrfout_d03*) as input that are located in $PATHTOWRFOUTFILES. 
In a regular WRF setup that is the directory where wrf.exe is executed. 

* The script loops over all files and ignores the last one (in our setup case the one that just contained the
single time step of 00:00:00 of the last day). 

* The file prefix can be defined by a system variable "NEWARUNNAME" that the python3
script reads for the system or defined by the parameter prefix in newa2timeseries.py
