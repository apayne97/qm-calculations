"""Methane (tetrahedral) at wB97M-V/def2-TZVPD: localize the four valence
orbitals into C-H bond orbitals and look at one in terms of C s/p functions.
Run from this folder:  python methane.py   (needs pyscf, numpy)
"""
import numpy as np
from pyscf import gto, dft, lo
from pyscf.tools import cubegen

R_CH = 1.087  # Angstrom
d = R_CH / np.sqrt(3)
atoms = [("C", (0, 0, 0)), ("H", (d, d, d)), ("H", (d, -d, -d)),
         ("H", (-d, d, -d)), ("H", (-d, -d, d))]
mol = gto.M(atom=atoms, basis="def2-tzvpd", verbose=0)
mf = dft.RKS(mol)
mf.xc = "wb97m_v"
mf.nlc = "VV10" if not getattr(mf, "nlc", None) else mf.nlc
mf.grids.atom_grid = (99, 590)
e = mf.kernel()
assert mf.converged
print(f"N basis functions: {mol.nao}   E = {e:.6f} Ha")
print("occupied canonical orbital energies (Ha):", np.round(mf.mo_energy[:5], 4))
print("  -> 1 core (C 1s), then a1 (1 orbital), then t2 (3 degenerate orbitals)")

# --- localize the 4 valence orbitals (leave the C 1s core alone) ---------
S = mf.get_ovlp()
mo_val = mf.mo_coeff[:, 1:5]
loc = lo.Boys(mol, mo_val)
loc.init_guess = "atomic"
c_loc = loc.kernel()
U = mo_val.T @ S @ c_loc            # rows: canonical (a1,t2,t2,t2); cols: localized
print("\nlocalized orbitals as mixtures of the canonical a1 and t2 orbitals")
print("(squared coefficients; columns = localized orbital 1..4):")
print(np.round(U ** 2, 3))
print("a1 share of each localized orbital:", np.round(U[0] ** 2, 3),
      " t2 share:", np.round((U[1:] ** 2).sum(0), 3))

# --- Lowdin analysis: how much of each localized orbital is on which atom ---
w, v = np.linalg.eigh(S)
S_half = v @ np.diag(np.sqrt(w)) @ v.T
labels = mol.ao_labels(fmt=False)       # (atom index, element, shell label, comp)
ang = lambda lab: "spdfg".index(lab[2][-1]) if lab[2][-1] in "spdfg" else None
def ang_of(lab):
    # label like '2px', '3dxy', '4f-2', '1s'
    s = lab[2]
    for ch in s:
        if ch in "spdfg":
            return "spdfg".index(ch)
ls = np.array([ang_of(l) for l in labels])
atom_of = np.array([l[0] for l in labels])

best = None
for j in range(4):
    cj = S_half @ c_loc[:, j]
    wH = [(cj[atom_of == a] ** 2).sum() for a in range(1, 5)]
    if best is None or max(wH) > best[0]:
        best = (max(wH), j, int(np.argmax(wH)) + 1)
_, J, H = best
cj_low = S_half @ c_loc[:, J]
print(f"\nlocalized orbital {J} is the C-H{H} bond orbital")
print(f"  Lowdin weight on C: {(cj_low[atom_of == 0] ** 2).sum():.3f}   on H{H}: "
      f"{(cj_low[atom_of == H] ** 2).sum():.3f}   on the other three H: "
      f"{sum((cj_low[atom_of == a] ** 2).sum() for a in range(1, 5) if a != H):.3f}")
wl = [(cj_low[(atom_of == 0) & (ls == l)] ** 2).sum() for l in range(3)]
tot = sum(wl)
print(f"  on carbon: s {wl[0]/tot:.1%}  p {wl[1]/tot:.1%}  d {wl[2]/tot:.1%}"
      f"   (sp3 would be 25% s, 75% p)")
print(f"  p/s weight ratio on carbon: {wl[1]/wl[0]:.2f}   (sp3: 3.00)")

# --- the carbon p part: direction compared with the C->H bond -------------
c = c_loc[:, J]
ao_start = mol.ao_loc_nr()
c_shells = [ib for ib in range(mol.nbas) if mol.bas_atom(ib) == 0 and mol.bas_angular(ib) == 1]
vec = np.array([c[ao_start[ib]:ao_start[ib + 1]] for ib in c_shells])
nrm = np.linalg.norm(vec, axis=1)
direction = vec[0] / nrm[0]
for v_, n_ in zip(vec, nrm):
    assert np.allclose(v_ / n_, direction, atol=1e-3) or np.allclose(v_ / n_, -direction, atol=1e-3)
bond = np.array(atoms[H][1]) / np.linalg.norm(atoms[H][1])
print("\ncarbon p shells all point along n =", np.round(direction, 4))
print("  C->H bond direction              =", np.round(bond, 4))
print(f"  angle between them: {np.degrees(np.arccos(abs(direction @ bond))):.2f} degrees")

# --- exact reconstruction as shells: s: Y00 c Σ d N e^(-a r²), p: sqrt(3/4π)(v·r) ... ----
rng = np.random.default_rng(0)
pts = rng.normal(scale=1.2, size=(200, 3))   # Bohr
ang_ = 1.8897261246
def part(ib):
    A = np.array(atoms[mol.bas_atom(ib)][1]) * ang_
    rel = pts - A
    r2 = (rel ** 2).sum(1)
    l = mol.bas_angular(ib)
    ex = mol.bas_exp(ib); dd = mol.bas_ctr_coeff(ib)[:, 0]
    N = np.array([gto.gto_norm(l, a) for a in ex])
    radial = (dd[None] * N[None] * np.exp(-ex[None] * r2[:, None])).sum(1)
    cc = c[ao_start[ib]:ao_start[ib + 1]]
    if l == 0:
        return cc[0] * np.sqrt(1 / (4 * np.pi)) * radial
    if l == 1:
        return np.sqrt(3 / (4 * np.pi)) * (rel @ cc) * radial
    return 0.0
psi_sp = sum(part(ib) for ib in range(mol.nbas) if mol.bas_angular(ib) <= 1)
psi_full = mol.eval_gto("GTOval_sph", pts) @ c
print(f"\nfull orbital: max |psi| = {abs(psi_full).max():.3f}")
print(f"s+p formula vs full orbital: max error = {abs(psi_sp - psi_full).max():.2e}")
print("  (the missing piece is the d and f functions)")
for lname, sel in (("C-centered", 0), ("H%d-centered" % H, H)):
    s_part = sum(part(ib) for ib in range(mol.nbas) if mol.bas_atom(ib) == sel and mol.bas_angular(ib) <= 1)
    print(f"  {lname} s+p part: rms {np.sqrt((s_part ** 2).mean()):.3f}")

# --- tables for the write-up ----------------------------------------------
print("\nCARBON s shells (c = coefficient on the shell; then primitives alpha, d*N):")
for ib in range(mol.nbas):
    if mol.bas_atom(ib) == 0 and mol.bas_angular(ib) == 0:
        cc = c[ao_start[ib]]
        ex = mol.bas_exp(ib); dd = mol.bas_ctr_coeff(ib)[:, 0]
        print(f"  shell {ib}: c = {cc: .5f}")
print("CARBON p shells (weight a = |v|, signed along n):")
sgn = np.sign(vec[0] @ direction)
for ib, v_ in zip(c_shells, vec):
    print(f"  shell {ib}: a = {np.sign(v_ @ direction) * np.linalg.norm(v_): .5f}")
print(f"H{H} shells:")
for ib in range(mol.nbas):
    if mol.bas_atom(ib) == H and mol.bas_angular(ib) <= 1:
        cc = c[ao_start[ib]:ao_start[ib + 1]]
        print(f"  shell {ib} l={mol.bas_angular(ib)}: c = {np.round(cc, 5)}")

cubegen.orbital(mol, "CH_bond_orbital.cube", c, nx=80, ny=80, nz=80)
cubegen.orbital(mol, "canonical_a1.cube", mf.mo_coeff[:, 1], nx=80, ny=80, nz=80)
cubegen.orbital(mol, "canonical_t2_a.cube", mf.mo_coeff[:, 2], nx=80, ny=80, nz=80)
with open("methane.xyz", "w") as f:
    f.write("5\nmethane\n")
    for el, (x, y, z) in atoms:
        f.write(f"{el} {x:.6f} {y:.6f} {z:.6f}\n")
np.save("localized_J_H.npy", np.array([J, H]))
