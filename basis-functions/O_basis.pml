# O def2-tzvpd basis functions, one scene per function. Run from this folder.
reinitialize
bg_color white
load O_atom.xyz, nucleus
show spheres, nucleus
set sphere_scale, 0.05, nucleus
color black, nucleus
load O_1s_0.cube, m_O_1s_0
isosurface O_1s_0_pos, m_O_1s_0, 2.8428
color blue, O_1s_0_pos
set transparency, 0.3, O_1s_0_pos
load O_2s_1.cube, m_O_2s_1
isosurface O_2s_1_pos, m_O_2s_1, 0.94079
color blue, O_2s_1_pos
set transparency, 0.3, O_2s_1_pos
load O_3s_2.cube, m_O_3s_2
isosurface O_3s_2_pos, m_O_3s_2, 0.19937
color blue, O_3s_2_pos
set transparency, 0.3, O_3s_2_pos
load O_4s_3.cube, m_O_4s_3
isosurface O_4s_3_pos, m_O_4s_3, 0.099733
color blue, O_4s_3_pos
set transparency, 0.3, O_4s_3_pos
load O_5s_4.cube, m_O_5s_4
isosurface O_5s_4_pos, m_O_5s_4, 0.049991
color blue, O_5s_4_pos
set transparency, 0.3, O_5s_4_pos
load O_6s_5.cube, m_O_6s_5
isosurface O_6s_5_pos, m_O_6s_5, 0.024187
color blue, O_6s_5_pos
set transparency, 0.3, O_6s_5_pos
load O_2pz_8.cube, m_O_2pz_8
isosurface O_2pz_8_pos, m_O_2pz_8, 0.32176
color blue, O_2pz_8_pos
set transparency, 0.3, O_2pz_8_pos
isosurface O_2pz_8_neg, m_O_2pz_8, -0.32176
color red, O_2pz_8_neg
set transparency, 0.3, O_2pz_8_neg
load O_3pz_11.cube, m_O_3pz_11
isosurface O_3pz_11_pos, m_O_3pz_11, 0.094092
color blue, O_3pz_11_pos
set transparency, 0.3, O_3pz_11_pos
isosurface O_3pz_11_neg, m_O_3pz_11, -0.094092
color red, O_3pz_11_neg
set transparency, 0.3, O_3pz_11_neg
load O_4pz_14.cube, m_O_4pz_14
isosurface O_4pz_14_pos, m_O_4pz_14, 0.040993
color blue, O_4pz_14_pos
set transparency, 0.3, O_4pz_14_pos
isosurface O_4pz_14_neg, m_O_4pz_14, -0.040993
color red, O_4pz_14_neg
set transparency, 0.3, O_4pz_14_neg
load O_5pz_17.cube, m_O_5pz_17
isosurface O_5pz_17_pos, m_O_5pz_17, 0.016301
color blue, O_5pz_17_pos
set transparency, 0.3, O_5pz_17_pos
isosurface O_5pz_17_neg, m_O_5pz_17, -0.016301
color red, O_5pz_17_neg
set transparency, 0.3, O_5pz_17_neg
load O_3dz2_20.cube, m_O_3dz2_20
isosurface O_3dz2_20_pos, m_O_3dz2_20, 0.28051
color blue, O_3dz2_20_pos
set transparency, 0.3, O_3dz2_20_pos
isosurface O_3dz2_20_neg, m_O_3dz2_20, -0.28051
color red, O_3dz2_20_neg
set transparency, 0.3, O_3dz2_20_neg
load O_4dz2_25.cube, m_O_4dz2_25
isosurface O_4dz2_25_pos, m_O_4dz2_25, 0.10765
color blue, O_4dz2_25_pos
set transparency, 0.3, O_4dz2_25_pos
isosurface O_4dz2_25_neg, m_O_4dz2_25, -0.10765
color red, O_4dz2_25_neg
set transparency, 0.3, O_4dz2_25_neg
load O_5dz2_30.cube, m_O_5dz2_30
isosurface O_5dz2_30_pos, m_O_5dz2_30, 0.035504
color blue, O_5dz2_30_pos
set transparency, 0.3, O_5dz2_30_pos
isosurface O_5dz2_30_neg, m_O_5dz2_30, -0.035504
color red, O_5dz2_30_neg
set transparency, 0.3, O_5dz2_30_neg
load O_4f0_36.cube, m_O_4f0_36
isosurface O_4f0_36_pos, m_O_4f0_36, 0.1928
color blue, O_4f0_36_pos
set transparency, 0.3, O_4f0_36_pos
isosurface O_4f0_36_neg, m_O_4f0_36, -0.1928
color red, O_4f0_36_neg
set transparency, 0.3, O_4f0_36_neg
disable *
orient nucleus
turn x, -75
turn y, 25
disable *
enable nucleus
enable O_1s_0_pos
zoom nucleus, buffer=0.57
scene O_1s_0, store, O 1s
disable *
enable nucleus
enable O_2s_1_pos
zoom nucleus, buffer=0.715
scene O_2s_1, store, O 2s
disable *
enable nucleus
enable O_3s_2_pos
zoom nucleus, buffer=1.13
scene O_3s_2, store, O 3s
disable *
enable nucleus
enable O_4s_3_pos
zoom nucleus, buffer=1.5
scene O_4s_3, store, O 4s
disable *
enable nucleus
enable O_5s_4_pos
zoom nucleus, buffer=2.09
scene O_5s_4, store, O 5s
disable *
enable nucleus
enable O_6s_5_pos
zoom nucleus, buffer=3.08
scene O_6s_5, store, O 6s
disable *
enable nucleus
enable O_2pz_8_pos
enable O_2pz_8_neg
zoom nucleus, buffer=1.12
scene O_2pz_8, store, O 2pz
disable *
enable nucleus
enable O_3pz_11_pos
enable O_3pz_11_neg
zoom nucleus, buffer=1.82
scene O_3pz_11, store, O 3pz
disable *
enable nucleus
enable O_4pz_14_pos
enable O_4pz_14_neg
zoom nucleus, buffer=2.8
scene O_4pz_14, store, O 4pz
disable *
enable nucleus
enable O_5pz_17_pos
enable O_5pz_17_neg
zoom nucleus, buffer=4.75
scene O_5pz_17, store, O 5pz
disable *
enable nucleus
enable O_3dz2_20_pos
enable O_3dz2_20_neg
zoom nucleus, buffer=1.23
scene O_3dz2_20, store, O 3dz^2
disable *
enable nucleus
enable O_4dz2_25_pos
enable O_4dz2_25_neg
zoom nucleus, buffer=1.89
scene O_4dz2_25, store, O 4dz^2
disable *
enable nucleus
enable O_5dz2_30_pos
enable O_5dz2_30_neg
zoom nucleus, buffer=3.42
scene O_5dz2_30, store, O 5dz^2
disable *
enable nucleus
enable O_4f0_36_pos
enable O_4f0_36_neg
zoom nucleus, buffer=1.54
scene O_4f0_36, store, O 4f+0
scene O_1s_0, recall
save O_basis.pse
