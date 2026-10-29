# Builds the water dimer session. Run from this folder:
#   pymol water_dimer_session.pml     (or File > Run Script in PyMOL)
reinitialize
bg_color white
load water_dimer.xyz, dimer
show sticks, dimer
set stick_radius, 0.12
color gray40, elem H
color red, elem O
# the hydrogen bond (donor H to acceptor O)
distance hbond, id 5, id 1
color black, hbond
set dash_gap, 0.25
set dash_radius, 0.03
hide labels, hbond

load water_dimer_density.cube, dens
isosurface dsurf, dens, 0.05
color slate, dsurf
set transparency, 0.5, dsurf

# difference density: rho(AB) - rho(A) - rho(B); blue = gained, red = lost on binding
load water_dimer_diff.cube, diff
isosurface diff_pos, diff, 0.004
isosurface diff_neg, diff, -0.004
color blue, diff_pos
color red, diff_neg
set transparency, 0.3, diff_pos
set transparency, 0.3, diff_neg

# start with only the difference density; toggle the total density in the object panel
disable dens
disable dsurf
orient dimer
zoom dimer, 1.5
save water_dimer_session.pse
