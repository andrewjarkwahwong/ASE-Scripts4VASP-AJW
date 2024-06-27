# VASP_ASE_AJW
ASE and bash scripts to assist in VASP5.4.4 and VASP6 calculations through my Ph.D.
Author: Andrew Jark-Wah Wong and Jin Li

Email: ajwongphd@gmail.com

The purpose of this repository is to share Python, specifically ASE, and bash scripts that made my life easier as a DFT researcher.
The intended audience is any DFT research that either started running a few calculations or are a veteran that would like scripts that improve their quality of life. Mainly, I wanted to introduce scripts from the Atomic Simulation Environment (ASE) and show its potential through a few examples. Disclaimer is that I am not a great coder but i've learned enough coding to know how beneficial it is so that is why I started this little project.

# What is Provided in the Repo:

## A Detailed Guide & Tutorial:
Please first refer to the LaTex PDF before running any of the scripts. It's a brief pdf is provided in "guide_PDF" to provide a comprehensive guide of my tips and tricks. It's more of a messy and brief guide on how to install ASE and utilize my scripts. More detailed documentation should be referred to https://wiki.fysik.dtu.dk/ase/ as the ASE environment wiki as it is well documented. 

However, my guide will give you a walk through of my experiences installing ASE, how to use ASE to run a VASP calculation, and how to utilize my scripts in the context of catalysis. 

## Python Scripts:
Python scripts are provided in python_scripts folder. The following scripts are currently provided.

### Running VASP through ASE (Folder: RunningVASP_ASE)
runvasp.py: The ASE input file to run VASP

asevasp: a SLURM submission script example.


### Generating POSCARs of surface slabs and simple adsorbate
Bare_FCC_Facet.py: Creates bare fcc surface facets. (Folder: BareSurface)

addadsorbate_monoatomic_fcc.py: Places a monoatomic adsorbate on a surface slab. (Folder: AddAdsorbate)

loop_addadsorbate_monoatomic_fcc.py: Utilizes a for loop to place a different single adsorbate on different surfaces. (Folder: AddAdsorbate)

## Bash Scripts:
Soon to provide these scripts in the bash_scripts folder that improved my quality-of-life as a DFT research.

## A Great Figure:
Refer to the Figures folder for a pleasant surprise.

If you have any recommendations or questions, feel free to reach out at ajwongphd@gmail.com!
