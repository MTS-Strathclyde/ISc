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
    usage: rism3d_isc.py [-h] [-c MOLCHARGE] [--multiplicity MULTIPLICITY]
                         [-t TEMPERATURE] [--clean_up CLEAN_UP]
                         [--timeout TIMEOUT] [--write_g] [--write_c] [--write_u]
                         [--write_asymp] [--noasympcorr] [--buffer BUFFER]
                         [--solvbox SOLVBOX SOLVBOX SOLVBOX]
                         [--tolerance [TOLERANCE [TOLERANCE ...]]]
                         [--grdsp GRDSP GRDSP GRDSP] [--polar_decomp]
                         [--charge_f CHARGE_F] [--smodel SMODEL]
                         [--closure [CLOSURE [CLOSURE ...]]]
                         [--verbose3d VERBOSE3D] [--maxstep3d MAXSTEP3D]
                         [--rism1d RISM1D] [--xvv XVV] [--therm THERM] [-p PRMTOP]
                         [--dir_name DIR_NAME] [--scale_chg SCALE_CHG]
                         molec.pdb

    Run 3D-RISM single point calculation at given temperature. Uses GAFF force
    field and AM1-BCC charges.

    positional arguments:
      molec.pdb             Input file. Must be in pdb format acceptable by
                            Antechamber. Must have pdb extension.

    optional arguments:
      -h, --help            show this help message and exit
      -c MOLCHARGE, --molcharge MOLCHARGE
                            Charge of the solute [0]
      --multiplicity MULTIPLICITY
                            Multiplicity of the solute [1]
      -t TEMPERATURE, --temperature TEMPERATURE
                            Temperature in K at which simulation will be run
                            [298.15]
      --clean_up CLEAN_UP   How should auxiliary files be treated: 0 - delete
                            nothing; 1 - delete some [default]; 2 - delete all but
                            input, results, and log.
      --timeout TIMEOUT     Minutes after which 3D-RISM calculation will be
                            killed. Use 0 for no timeout. [30]. Only works on
                            Unix-like system.
      --write_g             Write radial distribution functions produced in 3D-
                            RISM calculation
      --write_c             Write direct correlation function produced in 3D-RISM
                            calculation
      --write_u             Write solute solvent potential energy grid
      --write_asymp         Write asymptotics of total and direct correlation
                            fuctions in real space.
      --noasympcorr         Thermodynamics of 3D-RISM is calculated without long-
                            range asymptotics.
      --buffer BUFFER       Minimum distance between the solute and the edge of
                            the solvent box in A for 3D-RISM calculation [25]
      --solvbox SOLVBOX SOLVBOX SOLVBOX
                            Size of the x, y, and z dimensions of the box in
                            Angstroms. Specifying this parameter overrides buffer.
      --tolerance [TOLERANCE [TOLERANCE ...]]
                            Maximum residual values for 3D-RISM solution
                            convergence. If many closures a list of closures can
                            be supplied. [1E-5]
      --grdsp GRDSP GRDSP GRDSP
                            Linear grid spacings for x, y and z dimensions. Should
                            be separated with spaces. Units: A. [0.5 0.5 0.5]
      --polar_decomp        Decomposes solvation free energy into polar and non-
                            polar components
      --charge_f CHARGE_F   Supply a charge file to the antechamber. A file should
                            contain a list of atomic partial charges appearing in
                            the same ordear as are atoms in tge pdb file. One row
                            shouldn't contain more than 8 charges.
      --smodel SMODEL       Solvent model available in
                            $AMBERHOME/dat/rism1d/model/{smodel}.mdl [SPC]
      --closure [CLOSURE [CLOSURE ...]]
                            Brdige closure which will be used in both 1D-RISM and
                            3D-RISM simmulations. Either HNC, PSEn or KH (n in
                            PSEn should be an integer) or a list of them for
                            sequential convergence. If solvent rism1d is
                            necessary, only last closure will be used for it.
                            [HNC]
      --verbose3d VERBOSE3D
                            Verbosity of 3D-RISM calculation. 0 - print nothing; 1
                            - print iterations; 2 - print all. [2]
      --maxstep3d MAXSTEP3D
                            Maximum number of iterations in 3D-RISM calculation
                            [500] .
      --rism1d RISM1D       Type of 1D-RISM theory. Only DRISM has been
                            extensively tested [DRISM]
      --xvv XVV             Submit existing xvv file. The 1D-RISM calculation will
                            not be run and all related parameters will be ignored.
                            The output therm file of 1D-RISM calculations must be
                            submited along the xvv file. Solvent density will be
                            taken from it. Otherwise, water solvent density will
                            be used for Isc calculation.
      --therm THERM         Output of a 1D-RISM calculation. Note that the
                            temperature option will be ignored.
      -p PRMTOP, --prmtop PRMTOP
                            Submit existing prmtop file.
      --dir_name DIR_NAME   Custom name for produced calculation directory. The
                            default one is: {mol_name}_{temperature}.
      --scale_chg SCALE_CHG
                            Scale all charges predicted by the antechamber by a
                            certain value. Only applicable if you don't provide
                            your own prmtop file through --prmtop. [1.0]
                        

Notes
-----
The script has been tested only on Ubuntu, but it should work on most Linux distributions and on Mac OS X. To make it Windows-friendly one would probably need to change the names of executable programs and add them to PATH.



