# H def2-tzvpd basis functions, one scene per function. Run from this folder.
reinitialize
bg_color white
load H_atom.xyz, nucleus
show spheres, nucleus
set sphere_scale, 0.05, nucleus
color black, nucleus
load H_1s_0.cube, m_H_1s_0
isosurface H_1s_0_pos, m_H_1s_0, 0.33835
color blue, H_1s_0_pos
set transparency, 0.3, H_1s_0_pos
load H_2s_1.cube, m_H_2s_1
isosurface H_2s_1_pos, m_H_2s_1, 0.076661
color blue, H_2s_1_pos
set transparency, 0.3, H_2s_1_pos
load H_3s_2.cube, m_H_3s_2
isosurface H_3s_2_pos, m_H_3s_2, 0.032232
color blue, H_3s_2_pos
set transparency, 0.3, H_3s_2_pos
load H_2pz_5.cube, m_H_2pz_5
isosurface H_2pz_5_pos, m_H_2pz_5, 0.12828
color blue, H_2pz_5_pos
set transparency, 0.3, H_2pz_5_pos
isosurface H_2pz_5_neg, m_H_2pz_5, -0.12828
color red, H_2pz_5_neg
set transparency, 0.3, H_2pz_5_neg
load H_3pz_8.cube, m_H_3pz_8
isosurface H_3pz_8_pos, m_H_3pz_8, 0.026105
color blue, H_3pz_8_pos
set transparency, 0.3, H_3pz_8_pos
isosurface H_3pz_8_neg, m_H_3pz_8, -0.026105
color red, H_3pz_8_neg
set transparency, 0.3, H_3pz_8_neg
disable *
orient nucleus
turn x, -75
turn y, 25
disable *
enable nucleus
enable H_1s_0_pos
zoom nucleus, buffer=0.954
scene H_1s_0, store, H 1s
disable *
enable nucleus
enable H_2s_1_pos
zoom nucleus, buffer=1.69
scene H_2s_1, store, H 2s
disable *
enable nucleus
enable H_3s_2_pos
zoom nucleus, buffer=2.63
scene H_3s_2, store, H 3s
disable *
enable nucleus
enable H_2pz_5_pos
enable H_2pz_5_neg
zoom nucleus, buffer=1.57
scene H_2pz_5, store, H 2pz
disable *
enable nucleus
enable H_3pz_8_pos
enable H_3pz_8_neg
zoom nucleus, buffer=3.61
scene H_3pz_8, store, H 3pz
scene H_1s_0, recall
save H_basis.pse
