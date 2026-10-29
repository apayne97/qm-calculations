# Builds the carbon p-orbital session. Run from this folder:
#   pymol C_p_session.pml     (or File > Run Script in PyMOL)
reinitialize
bg_color white
load C.xyz, atom
show spheres, atom
set sphere_scale, 0.25, atom
color gray40, atom

# the p orbital: blue = positive lobe, red = negative lobe (sign is just phase)
load C_p_orbital.cube, porb
isosurface p_pos, porb, 0.05
isosurface p_neg, porb, -0.05
color blue, p_pos
color red, p_neg
set transparency, 0.3, p_pos
set transparency, 0.3, p_neg

# the lobe axis n = (0.7871, 0.5417, 0.2950) from carbon.py, drawn in Angstrom
pseudoatom axis_a, pos=[-1.968, -1.354, -0.738]
pseudoatom axis_b, pos=[1.968, 1.354, 0.738]
distance axis, axis_a, axis_b
hide labels, axis
set dash_gap, 0.15
set dash_radius, 0.02
color black, axis
hide everything, axis_a
hide everything, axis_b

# total density, off at first
load C_density.cube, dens
isosurface dsurf, dens, 0.05
color slate, dsurf
set transparency, 0.6, dsurf
disable dens
disable dsurf
# put the lobe axis across the screen
orient axis_a or axis_b
zoom axis_a or axis_b, 3
save C_p_session.pse
