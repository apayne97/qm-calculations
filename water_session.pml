# Builds the water density/HOMO/LUMO session. Run from this folder:
#   pymol water_session.pml     (or File > Run Script in PyMOL)
reinitialize
bg_color white
load water.xyz, mol
show sticks, mol
set stick_radius, 0.12
color gray40, elem H
color red, elem O

load water_density.cube, dens
isosurface dsurf, dens, 0.05
color slate, dsurf
set transparency, 0.5, dsurf

load water_HOMO.cube, homo
isosurface homo_pos, homo, 0.05
isosurface homo_neg, homo, -0.05
color blue, homo_pos
color red, homo_neg

load water_LUMO.cube, lumo
isosurface lumo_pos, lumo, 0.05
isosurface lumo_neg, lumo, -0.05
color forest, lumo_pos
color orange, lumo_neg

# start with only the density shown; toggle the others in the object panel
disable homo_pos
disable homo_neg
disable lumo_pos
disable lumo_neg
disable homo
disable lumo
disable dens
orient mol
# orient looks down the axis perpendicular to the molecular plane, which hides
# the HOMO's two lobes behind each other; turn to an oblique side view
turn y, 90
turn y, 40
zoom mol, 4
set transparency, 0.3, homo_pos
set transparency, 0.3, homo_neg
set transparency, 0.3, lumo_pos
set transparency, 0.3, lumo_neg
save water_session.pse
