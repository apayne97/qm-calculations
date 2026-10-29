# qm-calculations

Small DFT test calculation on water, plus visualization files. Made while
working out what the OMol25 level of theory means (see slides 56–57 of the
2026.10.06 Schrödinger job talk).

## Level of theory

- Functional: ωB97M-V (range-separated hybrid meta-GGA, VV10 nonlocal dispersion)
- Basis set: def2-TZVPD (58 basis functions for water; 40 on O, 9 on each H)
- Integration grid: (99, 590)
- Code: PySCF 2.14.0 (not ORCA, so energies may differ slightly from ORCA's)
- Geometry: water, O at (0, 0, 0.117790), H at (0, ±0.755453, -0.471161) Å
- Result: E = -76.4334 Ha, SCF converged; HOMO -0.4280 Ha, LUMO 0.0400 Ha

For comparison, SPICE uses ωB97M-D3(BJ)/def2-TZVPPD (74 basis functions for
water). The two basis sets differ mainly on hydrogen (an extra p shell and a
d shell in TZVPPD) and on one d shell for S and Cl.

## Files

| File | What it is |
|---|---|
| `water.py` | Runs the calculation and writes the cube files |
| `plots.py` | Plots oxygen's def2-TZVPD basis functions along z |
| `water.xyz` | Geometry, for showing atoms in PyMOL |
| `water_density.cube` | Total electron density |
| `water_HOMO.cube`, `water_LUMO.cube` | Frontier orbitals |
| `oxygen_basis.png` | Output of `plots.py` |
| `water_session.pse` | PyMOL session: atoms, density, HOMO, LUMO (only the density is shown at first; toggle the others in the object panel) |
| `water_session.pml` | Script that rebuilds the session; run it from this folder |

## Reproducing

Needs `pyscf`, `numpy`, `matplotlib`. The scripts write to the current
directory, so run them from this folder:

```
pip install pyscf numpy matplotlib
python water.py
python plots.py
```

## Viewing in PyMOL

Open `water_session.pse`, or run `water_session.pml` from this folder to rebuild
it. The `.pse` embeds the cube data (~14 MB). To do it by hand:

PyMOL reads cube files as maps only, so load the xyz for the atoms. The 0.05
isovalue is a starting point; adjust it until the surface looks right. Orbital
sign (red/blue) is just phase.

```
load water.xyz, mol
show sticks, mol
load water_density.cube, dens
isosurface dsurf, dens, 0.05
set transparency, 0.5, dsurf
load water_HOMO.cube, homo
isosurface homo_pos, homo, 0.05
isosurface homo_neg, homo, -0.05
color blue, homo_pos
color red, homo_neg
```

Avogadro and IboView also open the cube files directly and let you slide the
isovalue.

## Notes

- `plots.py` only plots s, pz, dz² and f along one axis. The 1s spike dominates
  the s panel, so the other s functions look flat; use a log scale or plot
  s2–s6 separately to see them.
- To get cube files from ORCA instead, run `orca_plot` on the `.gbw` file.

## 3D basis functions (`basis-functions/`)

What the def2-TZVPD basis functions look like in space, for O and H
separately. One PyMOL session per atom, with one scene per function
(`O_basis.pse`, `H_basis.pse`); step through the scenes in PyMOL's scene
panel or with the scene keys. `oxygen_basis_3d.png` and `hydrogen_basis_3d.png`
are overview grids of every scene.

- One component per shell is shown (s, pz, dz², f0). The other components of a
  shell are rotations of the same shape.
- Blue is positive, red is negative. Each function is contoured at 25% of its
  own maximum |value|, so the tight 1s and the diffuse functions can be compared
  by shape but not by absolute size at a common isovalue. The black cross is the
  nucleus.
- Labels like `3pz` are PySCF's running count of shells with that angular
  momentum, not principal quantum numbers. The later (larger) ones in each
  family are the more diffuse functions.
- Oxygen: 6 s, 4 p, 3 d, 1 f shells (40 functions in all). Hydrogen: 3 s and
  2 p shells (9 functions in all).
- `basis3d.py` regenerates the cube files and the `.pml` scripts (needs
  `pyscf`); `O_basis.pml` / `H_basis.pml` rebuild the sessions and write the
  `.pse` files. Run them from the `basis-functions/` folder.

## Water dimer (`dimer/`)

A counterpoise-corrected interaction energy for a water dimer, at the same
level of theory as the monomer above (ωB97M-V/def2-TZVPD, PySCF).

- `water_dimer.py` builds the dimer (O···O = 2.91 Å, donor H pointing at the
  acceptor O, O–H···O = 180°, rigid monomers, **not optimized**) and runs five
  calculations: the dimer, each monomer in its own basis, and each monomer in
  the dimer basis (partner replaced by ghost atoms, `ghost-O` / `ghost-H`).
- E_int = E(AB) − E(A) − E(B). The counterpoise (CP) version uses the monomer
  energies from the dimer basis.
- Result: E_int(raw) = −4.30, E_int(CP) = −4.22 kcal/mol, BSSE = 0.08 kcal/mol.
  The BSSE is small because def2-TZVPD already has diffuse functions. An
  optimized dimer would bind by roughly −5 kcal/mol.
- `scan.py` repeats this at several O···O scale factors (0.9 to 2.0) and writes
  `scan.csv`. Moving the donor along the O···O axis is a simplified version of
  the "distance scaling" evaluation in OMol25.
- Scan result (CP E_int, kcal/mol): 2.62 Å −2.01, 2.91 Å −4.22, 3.20 Å −3.90,
  3.64 Å −2.64, 4.37 Å −1.27, 5.82 Å −0.40. Repulsion dominates on the short
  side, binding peaks near 2.9 Å, and it decays slowly (dipole–dipole) rather
  than vanishing, so zero is only reached far out. BSSE shrinks with distance
  (0.11 to 0.004 kcal/mol).
- Cube files: `water_dimer_density.cube` is the dimer density;
  `water_dimer_diff.cube` is ρ(AB) − ρ(A) − ρ(B), with the monomer densities
  computed in the dimer basis so the subtraction is on the same footing.
  `water_dimer_session.pse` (rebuilt by `water_dimer_session.pml`) shows the
  difference density at ±0.004: blue is density gained on binding, red is
  density lost.
- Run from the `dimer/` folder; needs `pyscf` and `numpy`. The full set takes
  about 3 minutes on a laptop.

## Carbon atom p orbital (`carbon/`)

One occupied p orbital of the carbon atom (triplet ground state, 1s² 2s² 2p²;
unrestricted ωB97M-V/def2-TZVPD, 37 basis functions), written out as an
explicit formula.

- `carbon.py` runs the calculation, takes alpha orbital 3, prints its
  coefficients and the p-shell primitives, and checks the closed formula below
  against PySCF at 50 random points (they agree to about 1e-8). The output is
  in `carbon_output.txt`.
- The atom is not spherical in a DFT calculation with integer occupations: the
  two occupied alpha p orbitals are degenerate (−0.3602 Ha), so the program
  picks an arbitrary rotation of them. This orbital's lobe axis is
  n = (0.7871, 0.5417, 0.2950), not x, y or z. The x : y : z coefficient ratio
  is the same on all three p shells, which is how you can tell.
- Orbital: ψ(r) = sqrt(3/4π) (n·r) [1.357 e^(−34.70 r²) + 1.448 e^(−7.958 r²)
  + 1.268 e^(−2.378 r²) + 0.8007 e^(−0.8143 r²) + 0.2614 e^(−0.2889 r²)
  + 0.04705 e^(−0.1006 r²)], with r in bohr. Each bracket coefficient is
  a_s · d_k · N_k (shell weight × contraction coefficient × primitive norm).
  All six are positive, so the radial part has no node.
- Left out of the formula: a non-p admixture with coefficients up to 9e-4
  (mostly f functions), numerical noise from the finite integration grid.
- `C_p_session.pse` (rebuilt by `C_p_session.pml`) shows the two lobes
  (blue +, red −) and the lobe axis as a dashed line; the total density is
  there but switched off.

## Methane: an sp3 bond orbital (`methane/`)

Tetrahedral methane (C–H 1.087 Å), restricted ωB97M-V/def2-TZVPD, 73 basis
functions, E = −40.4910 Ha. Hybridization only exists inside a molecule, so this
is the sp3 case that the isolated carbon atom cannot show.

- `methane.py` runs the calculation, localizes the four valence orbitals (Boys;
  the C 1s core is left alone), and analyzes one C–H bond orbital. Output is in
  `methane_output.txt`.
- Canonical orbitals: one a1 (−0.8128 Ha) and three degenerate t2 (−0.4932 Ha),
  delocalized over the whole molecule. Neither looks like a bond.
- Each localized bond orbital is exactly **25% a1 + 75% t2** (squared mixing
  coefficients). That is where "sp3" comes from: a1 has the symmetry of carbon's
  s function and t2 of its p functions, so one s plus three p make four
  equivalent orbitals pointing at the corners of a tetrahedron.
- Population analysis of the chosen bond orbital depends on the method
  (`analysis_methods.py`). Löwdin: 51% on carbon, and on carbon s 20.1%, p 75.8%,
  d 4.1%. Mulliken: 61% on carbon, and on carbon s 27.7%, p 71.4%, d 0.9%.
  Ideal sp3 is 25% s, 75% p. The d part is the usual polarization role (a small
  shape correction); I did not test that by removing the d shells.
- The carbon p functions all point along (0.577, 0.577, 0.577), the C→H bond
  direction, to better than 0.01 degrees.
- Carbon s shell coefficients in the bond orbital: −0.047, −0.101 (tight shells),
  +0.160, +0.160, +0.074 (diffuse). The sign change with distance is a radial node,
  built from basis functions that have none. Carbon p shell weights along
  n: 0.282, 0.202, 0.166. Bonded-H s shells: 0.202, 0.275, −0.013.
- Closed form: ψ ≈ Σ over s shells of c_s Y₀₀ R_s(r) + √(3/4π) (n·r) Σ over p
  shells of a_p R_p(r), with R(r) = Σ_k d_k N_k e^(−α_k r²), summed over carbon and
  the bonded H (the other three H are small). It reproduces PySCF's value to a
  maximum error of 8.6e-3 against a peak amplitude of 0.388. The leftover is the
  omitted d and f functions.
- `methane_session.pse` (rebuilt by `methane_session.pml`) shows the bond
  orbital: a big blue lobe along the C–H bond and a small red back lobe behind
  carbon. The canonical a1 and t2 orbitals are in the same file, switched off.

## Methane dimer: does dispersion matter? (`methane-dimer/`)

Counterpoise-corrected interaction energies for two methanes at C···C = 3.80 Å
(rigid, staggered with three H facing each other, not optimized), def2-TZVPD.
The dimer is bound almost entirely by dispersion.

| Method | E_int, CP (kcal/mol) |
|---|---|
| HF | +0.28 |
| B3LYP | +0.32 |
| ωB97M-V with VV10 off | +0.07 |
| ωB97M-V | −0.42 |

- HF and B3LYP have no dispersion, so they find the pair slightly repulsive.
- Turning VV10 off in ωB97M-V takes the interaction energy from −0.42 to +0.07, so
  VV10 supplies about 0.5 kcal/mol here. "VV10 off" is only a diagnostic: the
  other functional parameters were fit with VV10 on.
- For scale, CCSD(T)/CBS for the optimized methane dimer is about −0.5 kcal/mol
  (S22 value, from memory); this geometry isn't optimized, so don't expect an
  exact match.
- In PySCF, `mf.nlc = False` is what disables VV10. `mf.nlc = ""` defers to the
  functional name and leaves it on. `methane_dimer.py` runs everything
  (about 8 minutes, mostly the two ωB97M-V cases).
