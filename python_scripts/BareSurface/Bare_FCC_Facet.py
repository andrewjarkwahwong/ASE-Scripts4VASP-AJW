# -*- coding: utf-8 -*-
"""
Facet_creator.py
Author: Andrew Jark-Wah Wong

Purpose: Creates a POSCAR of a bare fcc metal surface facet through ASE

"""

from ase import Atoms
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.build import fcc111, add_adsorbate, fcc100, fcc110
from ase.visualize import view
from ase.io import write 
import os 
import shutil

## --- USER INPUTS --- ##
# Define the late-transition-metal of interest
metal = 'Cu'

# Replace fcc111 with fcc110 or fcc100
facet_int = fcc100 #fcc(111), fcc(110), and fcc(100)

#Specify Slab size, layers, and Vacuum 
surface_size = (3,3,5) 
vacuum_size = 7.5

slab = facet_int(metal, surface_size, vacuum_size) #Make surface that's 3x3 with 7.5 A each side (15 A total)

## --- Code --- ### 
# Define a function that determines which Cu atoms to freeze based on Z
# CHANGE THIS VALUE FOR DIFFERENT FACETS MUST BE DEFINED TO FREEZE BOTTOM 3 LAYERS of a 5 layer slab
if facet_int == fcc111:
    z = 13
elif facet_int == fcc110:
    z = 11
elif facet_int == fcc100:
    z = 12
else: 
    print('Facet not defined/avaiable')
    

    
def should_freeze(atom):
    if atom.symbol == metal and atom.z < z: #13 for 111, 11 for 110, 12 for 100
        return True
    else:
        return False

'''
Converts 'Ase.Atoms' objects to VASP POSCAR/CONTCAR filesfor compatibility

Args:
    atoms              Converted 'Ase.Atoms' object of the surface given
    
write POSCAR file with vasp format
'''
def Atoms2con(filename,atoms,format1):
    write(filename=filename, images=atoms, format=format1)
    
# ... (the imports and function definitions remain unchanged)
'''
Removes the adsorbate from the given surface (in our case C2H4)

Args:
    surface_atoms           'Ase.Atoms' object of the surface which has adsorbate to be removed from the surface (possibly)
    
Returns:
    slab_constrained        'Ase.Atoms' object of the bare surface with atoms 2 Angstroms below the surface to be constrained (frozen)
'''


# Create a Boolean array that specifies which atoms to freeze
mask = [should_freeze(a) for a in slab]

# Use the modified mask in the FixAtoms command 
constraint = FixAtoms(mask=mask)

# View the Slab as a POSCAR output (here you can save the file from the popup widntow to whatever and drop it into the cluster)
slab.set_constraint(constraint)
Atoms2con('POSCAR', slab, 'vasp')
view(slab)
write(filename='POSCAR',images=slab, format='vasp')
#choose your own file name
filename = f"{metal}_{facet_int.__name__}_bare.vasp"
Atoms2con(filename, slab, 'vasp')
directory_name = os.path.splitext(filename)[0]
if not os.path.exists(directory_name):
    os.makedirs(directory_name)
shutil.move(filename, os.path.join(directory_name, filename))
os.rename(os.path.join(directory_name, filename), os.path.join(directory_name, "POSCAR"))
