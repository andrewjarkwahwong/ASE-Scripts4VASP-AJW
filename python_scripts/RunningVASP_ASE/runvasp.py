# Input File for ASE
# Authors: Andrew Jark-Wah Wong and Jin Lin
# Date June 2024

# Inputs:
# 1. K-Points
# 2. INCAR inputs
# 3. POTCAR and Exchange Correlation Functional

import numpy as np
from collections import OrderedDict
from ase.io.vasp import read_vasp
from ase.calculators.vasp import Vasp

#%%
def con2Atoms(poscar):
    atoms = read_vasp(poscar)
    
    return atoms
atoms =con2Atoms("POSCAR")
#%%

## -- Define K-Point Grid from Monkhorst-Pack Method Here -- ##
# Note always test convergence between # of irreducible K-Points and energy for each new surface slab.
'''
Write KPOINT here
'''
KPOINT=(7,7,1)
#%%
"""
All vasp parameter ASE can read:
    
float_keys = [
    'aexx',       # Fraction of exact/DFT exchange
    'aggac',      # Fraction of gradient correction to correlation
    'aggax',      # Fraction of gradient correction to exchange
    'aldac',      # Fraction of LDA correlation energy
    'amin',       #
    'amix',       #
    'amix_mag',   #
    'bmix',       # tags for mixing
    'bmix_mag',   #
    'deper',      # relative stopping criterion for optimization of eigenvalue
    'ebreak',     # absolute stopping criterion for optimization of eigenvalues (EDIFF/N-BANDS/4)
    'emax',       # energy-range for DOSCAR file
    'emin',       #
    'enaug',      # Density cutoff
    'encut',      # Planewave cutoff
    'encutfock',  # FFT grid in the HF related routines
    'hfscreen',   # attribute to change from PBE0 to HSE
    'potim',      # time-step for ion-motion (fs)
    'nelect',     # total number of electrons
    'param1',     # Exchange parameter 
    'param2',     # Exchange parameter 
    'pomass',     # mass of ions in am
    'sigma',      # broadening in eV
    'time',       # special control tag
    'weimin',     # maximum weight for a band to be considered empty
    'zab_vdw',    # vdW-DF parameter
    'zval',       # ionic valence
#The next keywords pertain to the VTST add-ons from Graeme Henkelman's group at UT Austin
    'jacobian',   # Weight of lattice to atomic motion
    'ddr',        # (DdR) dimer separation
    'drotmax',    # (DRotMax) number of rotation steps per translation step
    'dfnmin',     # (DFNMin) rotational force below which dimer is not rotated
    'dfnmax',     # (DFNMax) rotational force below which dimer rotation stops
    'stol',       # convergence ratio for minimum eigenvalue
    'sdr',        # finite difference for setting up Lanczos matrix and step size when translating
    'maxmove',    # Max step for translation for IOPT > 0
    'invcurve',   # Initial curvature for LBFGS (IOPT = 1)
    'timestep',   # Dynamical timestep for IOPT = 3 and IOPT = 7
    'sdalpha',    # Ratio between force and step size for IOPT = 4
#The next keywords pertain to IOPT = 7 (i.e. FIRE)
    'ftimemax',   # Max time step
    'ftimedec',   # Factor to dec. dt
    'ftimeinc',   # Factor to inc. dt
    'falpha',     # Parameter for velocity damping
    'falphadec',  # Factor to dec. alpha
]

exp_keys = [
    'ediff',      # stopping-criterion for electronic upd.
    'ediffg',     # stopping-criterion for ionic upd.
    'symprec',    # precession in symmetry routines
#The next keywords pertain to the VTST add-ons from Graeme Henkelman's group at UT Austin
    'fdstep',     # Finite diference step for IOPT = 1 or 2
]

string_keys = [
    'algo',       # algorithm: Normal (Davidson) | Fast | Very_Fast (RMM-DIIS)
    'gga',        # xc-type: PW PB LM or 91
    'prec',       # Precission of calculation (Low, Normal, Accurate)
    'system',     # name of System
    'tebeg',      #
    'teend',      # temperature during run
]

int_keys = [
    'ialgo',      # algorithm: use only 8 (CG) or 48 (RMM-DIIS)
    'ibrion',     # ionic relaxation: 0-MD 1-quasi-New 2-CG
    'icharg',     # charge: 0-WAVECAR 1-CHGCAR 2-atom 10-const
    'idipol',     # monopol/dipol and quadropole corrections
    'iniwav',     # initial electr wf. : 0-lowe 1-rand
    'isif',       # calculate stress and what to relax
    'ismear',     # part. occupancies: -5 Blochl -4-tet -1-fermi 0-gaus >0 MP
    'ispin',      # spin-polarized calculation
    'istart',     # startjob: 0-new 1-cont 2-samecut
    'isym',       # symmetry: 0-nonsym 1-usesym 2-usePAWsym
    'iwavpr',     # prediction of wf.: 0-non 1-charg 2-wave 3-comb
    'ldauprint',  # 0-silent, 1-occ. matrix written to OUTCAR, 2-1+pot. matrix written
    'ldautype',   # L(S)DA+U: 1-Liechtenstein 2-Dudarev 4-Liechtenstein(LDAU)
    'lmaxmix',    #
    'lorbit',     # create PROOUT
    'maxmix',     #
    'ngx',        # FFT mesh for wavefunctions, x
    'ngxf',       # FFT mesh for charges x
    'ngy',        # FFT mesh for wavefunctions, y
    'ngyf',       # FFT mesh for charges y
    'ngz',        # FFT mesh for wavefunctions, z
    'ngzf',       # FFT mesh for charges z
    'nbands',     # Number of bands
    'nblk',       # blocking for some BLAS calls (Sec. 6.5)
    'nbmod',      # specifies mode for partial charge calculation
    'nelm',       # nr. of electronic steps (default 60)
    'nelmdl',     # nr. of initial electronic steps
    'nelmin',
    'nfree',      # number of steps per DOF when calculting Hessian using finitite differences
    'nkred',      # define sub grid of q-points for HF with nkredx=nkredy=nkredz
    'nkredx',      # define sub grid of q-points in x direction for HF
    'nkredy',      # define sub grid of q-points in y direction for HF
    'nkredz',      # define sub grid of q-points in z direction for HF
    'npar',       # parallelization over bands
    'nsim',       # evaluate NSIM bands simultaneously if using RMM-DIIS
    'nsw',        # number of steps for ionic upd.
    'nupdown',    # fix spin moment to specified value
    'nwrite',     # verbosity write-flag (how much is written)
    'smass',      # Nose mass-parameter (am)
    'vdwgr',      # extra keyword for Andris program
    'vdwrn',      # extra keyword for Andris program
    'voskown',    # use Vosko, Wilk, Nusair interpolation
#The next keywords pertain to the VTST add-ons from Graeme Henkelman's group at UT Austin
    'ichain',     # Flag for controlling which method is being used (0=NEB, 1=DynMat, 2=Dimer, 3=Lanczos)
                  # if ichain > 3, then both IBRION and POTIM are automatically set in the INCAR file
    'iopt',       # Controls which optimizer to use.  for iopt > 0, ibrion = 3 and potim = 0.0
    'snl',        # Maximum dimentionality of the Lanczos matrix
    'lbfgsmem',   # Steps saved for inverse Hessian for IOPT = 1 (LBFGS)
    'fnmin',      # Max iter. before adjusting dt and alpha for IOPT = 7 (FIRE) 
]

bool_keys = [
    'addgrid',    # finer grid for augmentation charge density
    'laechg',     # write AECCAR0/AECCAR1/AECCAR2
    'lasph',      # non-spherical contributions to XC energy (and pot for VASP.5.X)
    'lasync',     # overlap communcation with calculations
    'lcharg',     #
    'lcorr',      # Harris-correction to forces
    'ldau',       # L(S)DA+U
    'ldiag',      # algorithm: perform sub space rotation
    'ldipol',     # potential correction mode
    'lelf',       # create ELFCAR
    'lhfcalc',    # switch to turn on Hartree Fock calculations
    'loptics',    # calculate the frequency dependent dielectric matrix
    'lpard',      # evaluate partial (band and/or k-point) decomposed charge density
    'lplane',     # parallelisation over the FFT grid
    'lscalapack', # switch off scaLAPACK
    'lscalu',     # switch of LU decomposition
    'lsepb',      # write out partial charge of each band seperately?
    'lsepk',      # write out partial charge of each k-point seperately?
    'lthomas',    #
    'luse_vdw',   # Invoke vdW-DF implementation by Klimes et. al 
    'lvhar',      # write Hartree potential to LOCPOT (vasp 5.x)
    'lvtot',      # create WAVECAR/CHGCAR/LOCPOT
    'lwave',      #
#The next keywords pertain to the VTST add-ons from Graeme Henkelman's group at UT Austin
    'lclimb',     # Turn on CI-NEB
    'ltangentold', # Old central difference tangent
    'ldneb',      # Turn on modified double nudging
    'lnebcell',   # Turn on SS-NEB
    'lglobal',    # Optmizize NEB globally for LBFGS (IOPT = 1)
    'llineopt',   # Use force based line minimizer for translation (IOPT = 1)
]

list_keys = [
    'dipol',      # center of cell for dipol
    'eint',       # energy range to calculate partial charge for
    'ferwe',      # Fixed band occupation (spin-paired)
    'ferdo',      # Fixed band occupation (spin-plarized)
    'iband',      # bands to calculate partial charge for
    'magmom',     # initial magnetic moments
    'kpuse',      # k-point to calculate partial charge for
    'ropt',       # number of grid points for non-local proj in real space
    'rwigs',      # Wigner-Seitz radii
    'ldauu',      # ldau parameters, has potential to redundant w.r.t. dict
    'ldaul',      # key 'ldau_luj', but 'ldau_luj' can't be read direct from
    'ldauj',      # the INCAR (since it needs to know information about atomic
                  # species. In case of conflict 'ldau_luj' gets written out
                  # when a calculation is set up
]

special_keys = [
    'lreal',      # non-local projectors in real space
]

dict_keys = [
    'ldau_luj',   # dictionary with L(S)DA+U parameters, e.g. {'Fe':{'L':2, 'U':4.0, 'J':0.9}, ...}
]

keys = [
    # 'NBLOCK' and KBLOCK       inner block; outer block
    # 'NPACO' and APACO         distance and nr. of slots for P.C.
    # 'WEIMIN, EBREAK, DEPER    special control tags
]
"""
#%%

## -- INCAR INPUTS -- ##

'''
Write INCAR tags here
'''
vasp=OrderedDict(ibrion=2,
                 nsw=500,
                 isif=2,
                 isym=0,
                 kpts=KPOINT,
                 lreal='AUTO',
                 ediff=1e-5,
                 ediffg=-0.5e-1,
                 encut=450.,
                 algo='Very_Fast',
                 symprec=1e-10,
                 lwave=False,
                 lcharg=False,
                 lvtot=False,
                 ldipol=True,
                 idipol=3,
                 dipol = (0, 0,0.5), #Defines the center of the slab. Note only direct coordinates in the z direction specified since dipole correction is in the z-direct
                 potim=0.5,
                 nelm=400,
                 ismear=2,
                 ispin=2,
                 prec = 'NORMAL',
                 voskown = 1,
                 sigma = 0.2
                 )

## -- Exchange Correlation Functional and Pseudopotenitals (POTCAR) -- ##

'''
Exchange correlation Functional (xc) and "non-standard" pseudopotential are specified here
calc = Vasp(xc='PBE', setups={'Be': '_sv','Mg': '_sv_GW','Ca': '_sv_GW','Ba': '_sv','Al': '_sv_GW','Nd': '_3'})
Other xc value and parameters set can be found in https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html

If pseduopotential is not specified above, the standard pseudopotential will be used. This may NOT be appropriate for your atom of interest.
Please consult the VASP wiki here regarding this: https://www.vasp.at/wiki/index.php/Available_pseudopotentials#potpaw.54
'''
calc = Vasp(**vasp,xc='PBE',setups={'Li': '_sv','Na': '_sv','K': '_sv','Rb':'_sv','Cs':'_sv'})
atoms.calc=calc
'''
runvasp
'''
atoms.get_potential_energy()


