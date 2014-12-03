ISc
===

Python script to simplify computation of 3D-RISM initial state corrected (ISc) hydration free energy.

Authors
-------
Maksim Mišin <mishin1991@gmail.com>

Usage examples
--------------

As an input script takes a pdb file with single solute. For example (methane.pdb):
```text
ATOM      1  C1  MOL     1       3.537   1.423   0.000  1.00  0.00
ATOM      2  H1  MOL     1       4.089   2.224   0.496  1.00  0.00
ATOM      3  H2  MOL     1       4.222   0.611  -0.254  1.00  0.00
ATOM      4  H3  MOL     1       2.759   1.049   0.669  1.00  0.00
ATOM      5  H4  MOL     1       3.077   1.810  -0.912  1.00  0.00
TER
END
```

1) 298.15 K methane hydration free energy calculation:

	$ python rism3d_isc.py methane.pdb

2) 350 K calculation with tip3p water
	
	$ python rism3d_isc.py methane.pdb -t 350 --wmodel TP3

Prerequisites
-------------

The script requires:

* Python 2.7 or later: http://www.python.org/
* AmberTools13 or later: http://ambermd.org/


Get some help
-------------

    $ python rism3d_isc.py -h
	usage: rism3d_isc.py [-h] [-t TEMPERATURE] [--clean_up CLEAN_UP]
		             [--timeout TIMEOUT] [--write_h] [--write_c]
		             [--buffer_distance BUFFER_DISTANCE]
		             [--tollerance TOLLERANCE] [--grdsp GRDSP GRDSP GRDSP]
		             [--polar_decomp] [--wmodel WMODEL] [--closure CLOSURE]
		             [--rism1d RISM1D]
		             molec.pdb

	Run 3D-RISM single point calculation at given temperature. Uses GAFF force
	field and AM1-BCC charges.

	positional arguments:
	  molec.pdb             Input file. Must be in pdb format acceptable by
		                Antechamber. Must have pdb extension.

	optional arguments:
	  -h, --help            show this help message and exit
	  -t TEMPERATURE, --temperature TEMPERATURE
		                Temperature in K at which simulation will be run
		                [298.15]
	  --clean_up CLEAN_UP   How should auxiliary files be treated: 0 - delete
		                nothing; 1 - delete some [default]; 2 - delete all but
		                input, results, and log.
	  --timeout TIMEOUT     Minutes after which 3D-RISM calculation will be
		                killed. Use 0 for no timeout. [30]. Only works on
		                Unix-like system.
	  --write_h             Write total correlation function produced in 3D-RISM
		                calculation
	  --write_c             Write direct correlation function produced in 3D-RISM
		                calculation
	  --buffer_distance BUFFER_DISTANCE
		                Minimum distance between the solute and the edge of
		                the solvent box in A for 3D-RISM calculation [30]
	  --tollerance TOLLERANCE
		                Maximum residual values for 3D-RISM solution
		                convergenc. [1E-10]
	  --grdsp GRDSP GRDSP GRDSP
		                Linear grid spacings for x, y and z dimensions. Should
		                be separated with spaces [0.3 0.3 0.3]
	  --polar_decomp        Decomposes solvation free energy into polar and non-
		                polar components
	  --wmodel WMODEL       Water model available in
		                $AMBERHOME/dat/rism1d/model/{wmodel}.mdl [SPC]
	  --closure CLOSURE     Brdige closure which will be used in both 1D-RISM and
		                3D-RISM simmulations. Note that ISc and PMVc
		                corrections that are printed out in results file are
		                meaningless with non-HNC closures [HNC]
	  --rism1d RISM1D       Type of 1D-RISM theory. Only DRISM has been
		                extensively tested [DRISM]

Notes
-----
The script has been tested only on Ubuntu, but it should work on most Linux distributions and on Mac OS X. To make it Windows-friendly one would probably need to change the names of executable programs and add them to PATH.



