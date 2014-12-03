#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 16:04:49 2014

@author: Maksim Misin (mishin1991@gmail.com)

Compute hydration free energie with rism3d.singlpnt using GAFF force field
and AM1-BCC charges. The script prepares topology, runs 1D-RISM and 
3D-RISM calculation and computes free energies using PMVc and ISc corrections.
The output is written to separate resutls.txt file.

As an input takes pdb file compatible with antechamber. For example:
ATOM      1  C1  MOL     1       3.537   1.423   0.000  1.00  0.00
ATOM      2  H1  MOL     1       4.089   2.224   0.496  1.00  0.00
ATOM      3  H2  MOL     1       4.222   0.611  -0.254  1.00  0.00
ATOM      4  H3  MOL     1       2.759   1.049   0.669  1.00  0.00
ATOM      5  H4  MOL     1       3.077   1.810  -0.912  1.00  0.00
TER
END

To run the simmulation simply type:
python run_3drism_ambertools.py molecule.pdb

The script requires working installations of python2.7 and AmberTools 12+

For more information run:
python rism3d_isc.py -h


    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division

import sys
import argparse
import subprocess
import distutils.spawn
import shutil
import os
import glob
import datetime
import time
import threading
import signal


## Non RISM globals ##
__version__ = '2014.1'

REQUIRED_EXECUTABLES = ['antechamber', 'parmchk', 'tleap', 'rism3d.snglpnt',
                        'rism1d']

IS_UNIX = os.name == 'posix'

## RISM-related globals ##

SOLV_SUCEPT_SCRPT = """#!/bin/csh -f

cat > water_{temp}.inp <<EOF
&PARAMETERS
	THEORY='{rism1d}', CLOSUR='{closure}',           !Theory
	NR=16384, DR=0.025,                    !Grid
	OUTLST='xCGT', routup=384, toutup=0,   !Output
	NIS=20, DELVV=0.3, TOLVV=1.e-12,       !MDIIS
	KSAVE=-1, KSHOW=1, maxstep=10000,      !Check pointing and iterations
	SMEAR=1, ADBCOR=0.5,                   !Electrostatics
	TEMPER={temp}, DIEps={diel},           !bulk solvent properties
	NSP=1
/
	&SPECIES                               !SPC water
	DENSITY={conc}d0,
	MODEL="$AMBERHOME/dat/rism1d/model/{wmodel}.mdl"
/
EOF

rism1d water_{temp} > water_{temp}.out || goto error

"""

RUNLEAP = """source leaprc.gaff
loadamberprep {name}.prepin
check MOL
loadamberparams {name}.frcmod
SaveAmberParm MOL {name}.prmtop {name}.incrd
SavePdb MOL {name}.pdb
quit
"""

RESULTS_NAME = 'results.txt'

RESULTS = """dGhyd(HNC)= {hnc} kcal/mol
dGhyd(GF)= {gf} kcal/mol
PMV= {pmv} AA^3

dGhyd(ISc)= {isc} kcal/mol
dGhyd(PMVc)= {pmvc} kcal/mol

"""


def process_command_line(argv):
    """Processes arguments

    Parameters
    ----------
    argv : list
        Command line arguments.

    Returns
    -------
    out : argparse.Namespace
        Namespace of command line arguments.
    """
    parser = argparse.ArgumentParser(description="""Run 3D-RISM single point
            calculation at given temperature. Uses GAFF force field and
            AM1-BCC charges.""")
    #Positional args
    parser.add_argument('file', metavar='molec.pdb',
                        help="""Input file. Must be in pdb format
                        acceptable by Antechamber. Must have pdb
                        extension.""")
    #Optional args
    parser.add_argument('-t', '--temperature',
                        help="""Temperature in K at which simulation will be
                        run [298.15]""", default=298.15, type=float)
    parser.add_argument('--clean_up',
                        help=""" How should auxiliary files be treated:
                        0 - delete nothing;
                        1 - delete some [default];
                        2 - delete all but input, results, and log.
                        """, default=1, type=int)
    parser.add_argument('--timeout',
                        help=""" Minutes after which 3D-RISM calculation
                        will be killed. Use 0 for no timeout. [30]. Only works
                        on Unix-like system.
                        """, default=30, type=float)
    parser.add_argument('--write_h',
                        help="""Write total correlation function produced
                        in 3D-RISM calculation""",
                        action='store_true')
    parser.add_argument('--write_c',
                        help="""Write direct correlation function produced
                        in 3D-RISM calculation""",
                        action='store_true')
    parser.add_argument('--buffer_distance',
                        help="""Minimum distance between the solute and the
                        edge of the solvent box in A for 3D-RISM
                        calculation [30]""",
                        default=30, type=float)
    parser.add_argument('--tollerance',
                        help=""" Maximum residual values for 3D-RISM solution
                        convergenc. [1E-10]""",
                        default=1e-10, type=float)
    parser.add_argument('--grdsp',
                        help="""Linear grid spacings for x, y and z
                        dimensions. Should be separated with spaces
                        [0.3 0.3 0.3]""",
                        default=(.3, .3, .3), nargs=3)
    parser.add_argument('--polar_decomp',
                        help="""Decomposes solvation free energy into polar
                        and non-polar components""",
                        action='store_true')
    parser.add_argument('--wmodel',
                        help="""Water model available in
                        $AMBERHOME/dat/rism1d/model/{wmodel}.mdl [SPC]""",
                        default="SPC")
    parser.add_argument('--closure',
                        help="""Brdige closure which will be used in both
                        1D-RISM and 3D-RISM simmulations. Note that ISc and
                        PMVc corrections that are printed out in results
                        file are meaningless with non-HNC closures [HNC]""",
                        default="HNC")
    parser.add_argument('--rism1d',
                        help="""Type of 1D-RISM theory. Only DRISM has been
                        extensively tested [DRISM]""",
                        default="DRISM")
    return parser.parse_args(argv)


def water_dielectric_const(T):
    """Return water dielectric constant for temperature 253.15K < T < 383.15K.
    Uses interpolation equation (eq. 9) for static dielectri constant found in
    the doucment by The International Association for the Properties of
    Water and Steam from 2011
    <http://www.iapws.org/relguide/LiquidWater.pdf>`__.
    Pressure = 0.1 MPa

    Parameters
    ----------
    T : float
        Temperature in K

    Returns
    -------
    e : float
        Water dielectric constant at T

    Examples
    --------
    >>> round(water_dielectric_const(273.15), 3)
    87.927
    >>> round(water_dielectric_const(298.15), 3)
    78.375
    >>> round(water_dielectric_const(375), 3)
    55.266
    """
    T_star = T/300.0
    coefs = [-43.7527, 299.504, -399.364, 221.327]
    exp_f = [-0.05, -1.47, -2.11, -2.31]
    e = 0
    for i in range(4):
        e += coefs[i]*T_star**(exp_f[i])
    return e


def water_concentration(T):
    """Return water concentration for temperature range 253.15K < T < 383.15K.
    Uses interpolation equation (eq. 2) for specific volume found in
    the doucment by The International Association for the Properties of
    Water and Steam from 2011
    <http://www.iapws.org/relguide/LiquidWater.pdf>`__.
    Pressure = 0.1 MPa

    Parameters
    ----------
    T : float
        Temperature in K

    Returns
    -------
    conc : float
        Water conentration at T in mol/l

    Examples
    --------
    >>> round(water_concentration(273.15), 3)
    55.498
    >>> round(water_concentration(298.15), 3)
    55.343
    """
    p0 = 10.0**5    # Pa
    R = 8.31464     # J/mol/K
    Tr = 10.0
    Ta = 593.0
    Tb = 232.0
    a = [1.93763157E-2,
         6.74458446E+3,
        -2.22521604E+5,
         1.00231247E+8,
        -1.63552118E+9,
         8.32299658E+9]
    b = [5.78545292E-3,
        -1.53195665E-2,
         3.11337859E-2,
        -4.23546241E-2,
         3.38713507E-2,
        -1.19946761E-2]
    n = [None, 4., 5., 7., 8., 9.]
    m = [1., 2., 3., 4., 5., 6.]
    def alpha(T):
        return Tr/(Ta - T)
    def beta(T):
        return Tr/(T - Tb)
    coef = a[0] + b[0]*beta(T)**m[0]
    for i in range(1, 6):
        coef += a[i]*alpha(T)**n[i] + b[i]*beta(T)**m[i]
    v0 = R*Tr/p0*coef  # m3/mol
    return 1/(v0*1000)    # mol/L


class RunCmd(threading.Thread):
    """ Will only work on Unix-like systems. And sometimes it will not work
    even on there. """
    def __init__(self, cmd, timeout, cwd='.'):
        """ Run subprocess for a fixed ammount of time and then kill it.

        Parameters
        ----------
        cmd : list
            Command to execute.

        timeout : float
            Time in minutes after which process will be killed. Pass 0 to
            let process run indefinitely.

        cwd : string, default .
            Directory in which process will be run
        """
        threading.Thread.__init__(self)
        self.cmd = cmd
        if timeout == 0:
            self.timeout = None
        else:
            self.timeout = timeout*60 #convert to seconds
        self.cwd = cwd

    def run(self):
        self.p = subprocess.Popen(self.cmd, preexec_fn=os.setsid,
                           stdout=subprocess.PIPE,
                           cwd=self.cwd)
        self.p.wait()

    def run_and_timeout(self):
        """Run subprocess

        Returns
        -------
        out : tuple
            stout and stderr outputs
        """
        self.start()
        self.join(self.timeout)

        if self.is_alive():
            os.killpg(self.p.pid, signal.SIGTERM)
            self.join()
            print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            print(self.cmd)
            print('Has run out of time')
            print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            raise RuntimeError("3D-RISM calc didn't finish in time.")
        return self.p.communicate()


class RISM3D_Singlpnt(object):
    """ A class used to assist setting up 3D-RISM calculation.

    Init is used to specify temperature as well as non-standard names for
    topology (prmtop) or water susceptibility (xvv) files.

    The calculation details like closure or tolerance are defined in
    setup_calculation method.
    """
    def __init__(self, name, T, logfile, prmtop_name=None, xvv_name=None):
        """ Create a class for running rism3d.snglpnt.

        Parameters
        ----------
        name : string
            Full path to pdb file without extension

        T : float
            A calculation temperature

        logfile : A writable file object
            A file to which calculation std. output will be written

        prmtop_name : string, default None
            A name of topology file for calculation without extension. If
            it is not specified defaults to name.prmtop

        xvv_name : string, default None
            A name of susceptibility file for this calculation. Defaults
            to water_T.xvv, where T is calculation temperature rounded to
            two digits.

        """
        self.name = name
        self.T = T
        self.p, self.no_p_name = os.path.split(name)
        self.logfile = logfile
        if prmtop_name:
            self.prmtop_name = prmtop_name
        else:
            self.prmtop_name = '{}.prmtop'.format(self.no_p_name)
        if xvv_name:
            self.xvv_name = xvv_name
        else:
            self.xvv_name = 'water_{temp}.xvv'.format(temp=self.T)
        self.run_flags_list = None

    def setup_calclation(self, closure='hnc', write_h=True,
                         write_c=True,
                         write_u=False, buffer_distance=30.0,
                         grdspc=(0.3, 0.3, 0.3),
                         tollerance=1e-10, polar_decomp=False):
        """ Setup calculation rism3d.snglpnt. calculation.

        More details on each of the parameter can be found in AmberTools
        manual RISM section.

        Parameters
        ----------
        closure : string, default hnc
            Allowed closure values are kh, hnc, pseN. Here N is an
            integer.

        write_h : boolean, default True
            Specifies whether program will write total correlation
            functions.

        write_c : boolean, default True
            Specifies wheter program will write direct correlation
            functions.

        write_u : boolean, default False
            Specifies wheter program will write potential energy
            grid.

        buffer_distance : float, default 30.0
            Minimum distance between the solute and the edge of
            the solvent box in A.

        grdsp: array-like (should contain 3 floats), default (0.3, 0.3, 0.3)
            Comma separated linear grid spacings for x, y and z dimensions.

        tollerance: float, default 1e-10
            Maximum residual values for solution convergence.

        polar_decomp: boolean, default False
            Decomposes solvation free energy into polar and non-polar
            components
        """
        if len(grdspc) == 3:
            grdspc = ','.join(map(str, grdspc))
        tollerance = str(tollerance)
        self.run_flags_list = ['rism3d.snglpnt',
             '--pdb', '{}.pdb'.format(self.no_p_name),
             '--prmtop', self.prmtop_name,
             '--closure', closure,
             '--xvv', self.xvv_name,
             '--buffer', '{:.6f}'.format(buffer_distance), #distance between solute
                                                      #and the edge of solvent box
             '--grdspc', grdspc,
             '--tolerance', tollerance]
        if write_h:
            self.run_flags_list.extend(['--huv',
                                         'h_{}'.format(self.no_p_name)])
        if write_c:
            self.run_flags_list.extend(['--cuv',
                                         'c_{}'.format(self.no_p_name)])
        if write_u:
            self.run_flags_list.extend(['--uuv',
                                         'u_{}'.format(self.no_p_name)])
        if polar_decomp:
            self.run_flags_list.extend(['--polarDecomp'])

    def run_calculation_and_log(self, timeout=30):
        """Run 3D-RISM single point calculation and log.

        Parameters
        ----------
        timeout : float, defult 30
            When HNC calculations get stuck they tend
            to run for a long amount of time. If calculation is
            running more than 30min most likely it is stuck.
            This option records 3D-RISM caclulation PID and
            kills it after supplied number of minutes. Doesn't work
            on windows. And in some other cases as well.
        """
        start_time = time.time()
        if IS_UNIX:  # use timeout
            run3drism = RunCmd(self.run_flags_list, timeout, cwd=self.p)
            rism_out = run3drism.run_and_timeout()[0]
        else:    # windows
            rism_out = subprocess.check_output(self.run_flags_list, cwd=self.p)
        self.logfile.write(rism_out)
        self.logfile.flush()
        #write timestamp and close
        end_time = time.time()
        self.logfile.write(str(datetime.datetime.now()) + '\n')
        runtime = end_time - start_time
        self.logfile.write(str(round(runtime)))
        self.logfile.flush()
        self.logfile.close()


def prepare_3drism_calc(name, T=298.15, wmodel="SPC", rism1d="DRISM",
                        closure="HNC"):
    """Generate topology and xvv file at given temperature. Parse and get
    water compressibility from water thermodynamic output. Return data needed
    for subsequent calculations.

    Parameters
    ----------
    name : string
        Full path to pdb file without extension

    T : float, default 298.15
        A calculation temperature

    wmodel : string, default SPC
        Water model available in $AMBERHOME/dat/rism1d/model/{wmodel}.mdl

    rism1d : string, default DRISM
        Type of 1D-RISM theory. Only DRISM has been extensively tested

    closure : string, default HNC
        Brdige closure which will be used in both 1D-RISM simmulation

    Returns
    -------
    out: tuple
        Returns (logfile, prmtop_name, compres).
        logfile is a writable file object containing std. out. prmtop_name
        is the name of prepared topology file. compres is the float value of
        water compressibility [10e-4/MPa]
    """
    p, no_p_name = os.path.split(name)
    log_name = '{}.log'.format(name)
    logfile = open(log_name, 'wb')
    logfile.write(str(datetime.datetime.now()))     # timestamp
    #Firstly we use antechamber to recognize atom and bonding types, and
    #generate topology
    ante_out = subprocess.check_output(['antechamber',
                     '-i', '{}.pdb'.format(no_p_name),
                     '-fi', 'pdb',
                     '-o', '{}.prepin'.format(no_p_name), #output file
                     '-fo', 'prepi',   #output format describing each residue
                     '-c', 'bcc',      #charge method  (AM1-BCC)
                     '-s', '2',    #status info ; 2 means verbose
                     '-nc', '0',   #Net molecule charge
                     '-m', '1'],   #Multiplicity
                     cwd=p)
    logfile.write(ante_out)
    #Run parmchk to generate missing gaff force field parameters
    parm_out = subprocess.check_output(['parmchk',
                     '-i', '{}.prepin'.format(no_p_name),
                     '-f', 'prepi',
                     '-o', '{}.frcmod'.format(no_p_name)], #file with missing FF params
                     cwd=p)
    logfile.write(parm_out)
    logfile.flush()
    #Run tleap to generate topology and coordinates for the molecule
    leap_input_name = os.path.join(p, 'runleap.in')
    with open(leap_input_name, 'wb') as f:
        f.write(RUNLEAP.format(name=no_p_name))
    leap_out = subprocess.check_output(['tleap', '-f', 'runleap.in'], cwd=p)
    logfile.write(leap_out)
    logfile.flush()
    #Generate water susceptibility file
    xvv_script_name_no_p = 'water_{}_script.sh'.format(T)
    xvv_script_name = os.path.join(p, xvv_script_name_no_p)
    diel = round(water_dielectric_const(T), 3)
    conc = round(water_concentration(T), 3)
    succ_srcirpt = SOLV_SUCEPT_SCRPT.format(temp=T, diel=diel, conc=conc,
                                            wmodel=wmodel, rism1d=rism1d,
                                            closure=closure)
    with open(xvv_script_name, 'wb') as f:
        f.write(succ_srcirpt)
    xvv_out = subprocess.check_output(['bash', xvv_script_name_no_p], cwd=p)
    logfile.write(xvv_out)
    logfile.flush()
    #Get compressibility from 1D-RISM therm. output
    water_therm_name = 'water_{}.therm'.format(T)
    water_therm_p = os.path.join(p, water_therm_name)
    with open(water_therm_p, 'rb') as f:
        w_lines = f.readlines()
    compres = float(w_lines[2].split()[-1])
    prmtop_name = '{}.prmtop'.format(no_p_name)
    return logfile, prmtop_name, compres


def prepare_calc_directory(mol_path, T):
    """Copy pdb file into the directory with the same name. If such directory
    doesn't exist it will try to create it.

    Parameters
    ----------
    mol_path: string
        Path to solute
    T: float
        A calculation temperature

    Returns
    -------
    name: string
        Full path to pdb file without extension
    """
    pdb_path, name_without_path = os.path.split(mol_path)
    dir_name = os.path.join(pdb_path, name_without_path[:-4] + '_' + str(T))
    try:
        os.mkdir(dir_name)
    except OSError, e:
        if e.errno == 17:
            pass
        else:
            raise e
    name = os.path.join(dir_name, name_without_path)
    shutil.copy(mol_path, name)
    return name[:-4]


def clean_up(name, T, level):
    """Delete junk.

    Parameters
    ----------
    name : string
        Full path to pdb file without extension

    T: float
        A calculation temperature

    level : {0, 1, 2}
        0 - delete nothing; 
        1 - delete ANTECHAMBER*, all water but .sh and .therm,
        .frcmod, .prepin, NEWPDB.PDB, PREP.INF, ATOMTYPE.INF, runleap.in
        sqm*, leap.log;
        2 - delete ALL but RESULTS_NAME and logfile.
    """
    p, no_p_name = os.path.split(name)
    water_name = 'water_{}'.format(T)
    to_del1_glob = ['ANTECHAMBER*', 'sqm*', 'water*vv*']
    to_del1_files = [no_p_name + '.prepin', no_p_name + '.frcmod',
                    water_name + '.inp', water_name + '.out',
                    water_name + '.sav', 'ATOMTYPE.INF',
                    'leap.log', 'NEWPDB.PDB', 'PREP.INF', 'runleap.in']
    will_be_deleted_list = []
    if level == 1:
        for wildcard in to_del1_glob:
            will_be_deleted_list.extend(glob.glob(os.path.join(p, wildcard)))
        will_be_deleted_list.extend([os.path.join(p, f) for f in \
                                                            to_del1_files])
    if level == 2:
        all_files = os.listdir(p)
        all_files.remove(RESULTS_NAME)
        log_name = '{}.log'.format(no_p_name)
        all_files.remove(log_name)
        will_be_deleted_list.extend([os.path.join(p, f) for f in all_files])
    for f in will_be_deleted_list:
        os.unlink(f)


def write_results(name, T, compres):
    """ Parses log file and writes free energies and corrections to
    results.txt.

    Parameters
    ----------
    name: string
        Full path to pdb file without extension

    T: float
        A calculation temperature

    compres: float
        Water compressibility obtained from 1D-RISM output [10e-4/MPa]
        (note that experimental compressibility is meaningless; here we are
        using HNC compressibility to correct gigantic pressure which water has
        according to 1D-RISM calculation).
    """
    p, _ = os.path.split(name)
    log_name = '{}.log'.format(name)
    with open(log_name, 'rb') as f:
        for line in f:
            if line[0:11] == "rism_exchem":
                hnc = float(line.split()[1])
            if line[0:11] == "rism_exchGF":
                gf = float(line.split()[1])
            if line[0:11] == "rism_volume":
                pmv = float(line.split()[1])
    density = water_concentration(T)*6.0221413E-4 #number density 1/A3
    rho_c_k0 = 1 - 1e-20/(density*1.3806488E-23*T*compres) #density*c(k=0)
    cor = 1.9872041E-3*T/2*rho_c_k0*density*pmv
    fix_cor = -density*1.9872041E-3*T*pmv
    isc = hnc + cor + fix_cor # initial state correction [kcal/mol]
    pmvc = hnc + cor # partial molar volume correction [kcal/mol]

    results = RESULTS.format(hnc=hnc, gf=gf, pmv=pmv, isc=isc, pmvc=pmvc)
    #Write results
    with open(os.path.join(p, RESULTS_NAME), 'wb') as f:
        f.write(results)
    return isc


def main(argv):
    args = process_command_line(argv)
    for executable in REQUIRED_EXECUTABLES:
        if not distutils.spawn.find_executable(executable):
            raise NameError("{} is not found!".format(executable))
    if not 253.15 <= args.temperature <= 383.15:
        raise ValueError("Temperature is outside of allowed range.")
    print('Starting HFE calculation for {} at T={} K'.format(args.file,
                                                            args.temperature))
    name = prepare_calc_directory(args.file, args.temperature)
    print('Running AM1-BCC and 1D-RISM calculations...')
    logfile, prmtop_name, compres = prepare_3drism_calc(name,
                            args.temperature, args.wmodel, args.rism1d,
                            args.closure)
    rism_calc = RISM3D_Singlpnt(name, args.temperature,
                                logfile, prmtop_name=prmtop_name)
    rism_calc.setup_calclation(args.closure, args.write_h, args.write_c,
                               buffer_distance=args.buffer_distance,
                               grdspc=args.grdsp,
                               tollerance=args.tollerance,
                               polar_decomp=args.polar_decomp)
    print('Running 3D-RISM calculation...')
    rism_calc.run_calculation_and_log(args.timeout)
    isc = write_results(name, args.temperature, compres)
    print('Calculation has finished')
    print('ISc dG*(hyd)={} kcal/mol'.format(isc))
    print('Detailed output can be found in {} and in {}.log'\
                                        .format(RESULTS_NAME, name))
    clean_up(name, args.temperature, args.clean_up)



if __name__ == '__main__':
    main(sys.argv[1:])



