# Single Adsorbate Placer on fcc metal
# Authors: Jin Li and Andrew Wong

# For a monoatomic adsorbates on a fcc metal facet of your choice
# Inputs: 
# -Facet
# -Adsorbate
# -Adsorption Sites
# -Metal
# -slab size, layer, and vacuum
# -constrained layers
# Outputs:
# -a .vasp file of your slab with adsorbate (POSCAR) placed in a directory

# Code layout:
"""
1. Packages
2. Definitions
3. Inputs
4. The actual code
"""

## -- Packages -- ##

from ase.build import fcc111,fcc110,fcc100, add_adsorbate
from ase.io import read, write
from ase.visualize import view
import numpy as np
from ase.constraints import FixAtoms
import os 
import shutil

## -- Definitions -- ##

# Creates bare surface
def baresurface(surfaceelement,surface_size,vacuum):
    """
    Create a bare FCC (111),(110), or (100) surface using the specified element, surface size, and vacuum.

    Parameters:
    surfaceelement (str): The chemical symbol of the surface element (e.g., 'Cu', 'Au', 'Pt').
    surface_size (tuple): A tuple of three integers specifying the size of the surface in the form (a, b, c),
                          where 'a' and 'b' are the number of surface unit cells in the x and y directions,
                          and 'c' is the number of atomic layers in the z direction.
    vacuum (float): The vacuum thickness in angstroms to be added on both sides of the slab along the z direction.

    Returns:
    slab (ase.Atoms): An Atoms object representing the generated FCC (111) surface.
    """
    # Define the Cu (111) surface
    if facet_int == '111':
        slab = fcc111(surfaceelement, size=surface_size, vacuum=vacuum,periodic=True)
    elif facet_int == '110':
        slab = fcc110(surfaceelement, size=surface_size, vacuum=vacuum,periodic=True)
    elif facet_int == '100':
        slab = fcc100(surfaceelement, size=surface_size, vacuum=vacuum,periodic=True)
    else:
        print('incorrect syntax when selecting slab')
    
    return slab

# %% Constrain slab layers
def constrain_slab(atoms, z_cutoff=2.):
    mask = []
    scaled_positions = atoms.get_scaled_positions()
    unit_cell_height = np.linalg.norm(atoms.cell[2])

    if atoms.cell[2, 2] > 0:
        max_height = max(position[2] for position in scaled_positions)
        threshold = max_height - z_cutoff / unit_cell_height
        for position in scaled_positions:
            if position[2] < threshold:
                mask.append(True)
            else:
                mask.append(False)
    else:
        raise RuntimeError('Tried to constrain a slab that points in neither '
                           'the positive nor negative z directions, so we do '
                           'not know which side to fix')

    atoms.constraints += [FixAtoms(mask=mask)]
    return atoms


# %% add adsorbate on the constrained slab
def addadsorbate(slab, adsorbate, height, position):
    """
    Add an adsorbate atom to a given slab.

    Parameters:
    slab (ase.Atoms): The slab (surface) to which the adsorbate will be added.
    adsorbate (str): The chemical symbol of the adsorbate element (e.g., 'Pt', 'O', 'C').
    height (float): The initial distance (in angstroms) between the adsorbate and the surface along the z direction.
    position (str): The position of the adsorbate on the surface ('ontop', 'bridge', 'hollow', or 'fcc').

    Returns:
    adslab (ase.Atoms): The slab with the adsorbate added.
    """
    
    # Add the adsorbate to the slab
    add_adsorbate(slab=slab, adsorbate=adsorbate, height=height, position=position)
    return slab

# %%
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

## -- User Inputs -- ## 
#Specify metal
metal = 'Au'

#specify fcc(111),(110),or (100)
facet_int='100' #'111','110', or '100'

#Specify Slab parameters
slab = baresurface(metal, (3, 3, 5), 7.5) # This specifies a y x y slab with z layers and v vacuum each side of the slab.

#Specify Adsorbate and the adsorption site
adsorbate = 'Mg'
position = 'ontop' #fcc(111):ontop, bridge, fcc, hcp
                 #fcc(110):ontop, longbridge, shortbridge, hollow
                 #fcc(100):ontop, bridge, hollow

## Constrain the sub-surface atoms of the slab (CHECK Z_CUTOFF VALUE)
constrained_slab = constrain_slab(slab, z_cutoff=3.) # 3 for fcc(111) 5 layer
                                                     # 2 for fcc(110) 5 layer
                                                     # 3 for fcc(100) 5 layer 


## -- Code Runs -- ##
#The Code will output open ase gui to visulaize the structure and save the output as a .vasp file (it's the same as a POSCAR but different name).

# Add the adsorbate
adslab = addadsorbate(constrained_slab, adsorbate, 3, position)
Atoms2con('POSCAR', adslab, 'vasp')
view(adslab)        
#choose your own file name
filename = f"{metal}_{facet_int}_{adsorbate}_{position}.vasp"
Atoms2con(filename, adslab, 'vasp')
# view(adslab)
directory_name = os.path.splitext(filename)[0]
if not os.path.exists(directory_name):
    os.makedirs(directory_name)
shutil.move(filename, os.path.join(directory_name, filename))
os.rename(os.path.join(directory_name, filename), os.path.join(directory_name, "POSCAR"))
                                                                                                              
