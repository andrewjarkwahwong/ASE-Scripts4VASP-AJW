# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from ase import Atoms
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.build import fcc111, add_adsorbate, fcc100, fcc110
from ase.visualize import view
from ase.io import write 

# Create Au Facet
Metal = 'Cu'

# Replace fcc111 with fcc110 or fcc100
slab = fcc100(Metal, size=(4, 4, 5), vacuum=7.5) #Make surface that's 3x3 with 7.5 A each side (15 A total)

#Define adsorbate (CO2- which is bent)
#d_CO = 1.1551 #These are values I got from optimized CO2- in gas as an initial guess
#h_CO = 0.447
#molecule = Atoms('CO2', positions=[(0., 0., 0.), (0, d_CO, h_CO),(0,-d_CO,h_CO)]) #Create CO2-


# Add adsorbate onto surface
#h = 2 #height above the site
#add_adsorbate(slab, molecule, h, 'hollow') ls


# Sites for FCC111: ontop, bridge, fcc, or hcp
# Sites for FCC110: ontop, longbridge, shortbridge, hollow
# Sites for FCC100: ontop, bridge, hollow


# Define a function that determines which Cu atoms to freeze based on Z
# CHANGE THIS VALUE FOR DIFFERENT FACETS MUST BE DEFINED TO FREEZE BOTTOM 3 LAYERS
z= 12 #13 for 111, 11 for 110, 12 for 100
def should_freeze(atom):
    if atom.symbol == 'Cu' and atom.z < z: #13 for 111, 11 for 110, 12 for 100
        return True
    else:
        return False

# Create a Boolean array that specifies which atoms to freeze
mask = [should_freeze(a) for a in slab]


# Use the modified mask in the FixAtoms command 
constraint = FixAtoms(mask=mask)

# Change Atom type 
# Sub surface vs Surface atoms 


# View the Slab as a POSCAR output (here you can save the file from the popup widntow to whatever and drop it into the cluster)
slab.set_constraint(constraint)
view(slab)
write(filename='POSCAR',images=slab, format='vasp')

#print('Adsorption energy:', e_slab + e_N2 - slab.get_potential_energy())