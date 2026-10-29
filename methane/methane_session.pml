# Builds the methane session. Run from this folder:
#   pymol methane_session.pml     (or File > Run Script in PyMOL)
reinitialize
bg_color white
load methane.xyz, mol
show sticks, mol
set stick_radius, 0.12
color gray40, elem H
color black, elem C

# one localized C-H bond orbital (the C-H1 bond, along (1,1,1)); blue +, red -
load CH_bond_orbital.cube, bond
isosurface bond_pos, bond, 0.06
isosurface bond_neg, bond, -0.06
color blue, bond_pos
color red, bond_neg
set transparency, 0.3, bond_pos
set transparency, 0.3, bond_neg

# the delocalized canonical orbitals the bond orbital is built from
load canonical_a1.cube, a1
isosurface a1_pos, a1, 0.06
isosurface a1_neg, a1, -0.06
color blue, a1_pos
color red, a1_neg
set transparency, 0.3, a1_pos
set transparency, 0.3, a1_neg

load canonical_t2_a.cube, t2
isosurface t2_pos, t2, 0.06
isosurface t2_neg, t2, -0.06
color blue, t2_pos
color red, t2_neg
set transparency, 0.3, t2_pos
set transparency, 0.3, t2_neg

# start with only the bond orbital; toggle the canonical ones in the object panel
disable a1
disable a1_pos
disable a1_neg
disable t2
disable t2_pos
disable t2_neg
disable bond
# put the C-H1 bond (the one the orbital sits on) across the screen
orient (id 1 or id 2)
zoom mol, 2
save methane_session.pse
