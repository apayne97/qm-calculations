"""Carbon atom (triplet ground state) at wB97M-V/def2-TZVPD, then one occupied
alpha p orbital written out as an explicit formula and checked numerically.
Run from this folder:  python carbon.py   (needs pyscf, numpy)
"""
import numpy as np
from pyscf import gto, dft
from pyscf.tools import cubegen

mol = gto.M(atom="C 0 0 0", basis="def2-tzvpd", spin=2, charge=0, verbose=0)
mf = dft.UKS(mol)
mf.xc = "wb97m_v"
mf.nlc = "VV10" if not getattr(mf, "nlc", None) else mf.nlc
mf.grids.atom_grid = (99, 590)
e = mf.kernel()
assert mf.converged
print(f"N basis functions: {mol.nao}   E = {e:.6f} Ha")
print("alpha occupied:", int(mf.mo_occ[0].sum()), "  beta occupied:", int(mf.mo_occ[1].sum()))
print("alpha orbital energies (Ha):", np.round(mf.mo_energy[0][:5], 4))
print("beta  orbital energies (Ha):", np.round(mf.mo_energy[1][:5], 4))

# alpha orbitals 0,1 = 1s,2s ; 2,3 = occupied p (degenerate) ; 4 = empty p
IDX = 3
c = mf.mo_coeff[0][:, IDX]
print(f"\nalpha orbital {IDX}: energy {mf.mo_energy[0][IDX]:.4f} Ha, occupation {mf.mo_occ[0][IDX]:.0f}")

# --- basis shells -------------------------------------------------------
print("\nshells (angular momentum l, number of primitives):")
ao_start = mol.ao_loc_nr()
for ib in range(mol.nbas):
    print(f"  shell {ib:2d}: l={mol.bas_angular(ib)}  nprim={mol.bas_nprim(ib)}"
          f"  AOs {ao_start[ib]}-{ao_start[ib+1]-1}")

# --- the p part of the orbital: ψ = Σ_shell (cx x + cy y + cz z) R_shell(r) ----
p_shells = [ib for ib in range(mol.nbas) if mol.bas_angular(ib) == 1]
vec = np.array([c[ao_start[ib]:ao_start[ib + 1]] for ib in p_shells])  # (shell, xyz)
# PySCF orders real p functions as px, py, pz
norms = np.linalg.norm(vec, axis=1)
direction = vec[0] / norms[0]
for v, n in zip(vec, norms):
    assert np.allclose(v / n, direction, atol=1e-6), "p shells not co-linear"
print("\ndirection n of the lobe axis (unit vector):", np.round(direction, 4))
print("radial weights a_s per p shell (coefficient vector length):", np.round(norms, 5))
p_aos = np.concatenate([np.arange(ao_start[ib], ao_start[ib + 1]) for ib in p_shells])
other = np.delete(np.arange(mol.nao), p_aos)
print("largest coefficient outside the p shells:", f"{abs(c[other]).max():.1e}")

print("\np-shell primitives (exponent alpha, contraction coefficient d, primitive norm N):")
for s, ib in enumerate(p_shells):
    ex = mol.bas_exp(ib)
    d = mol.bas_ctr_coeff(ib)[:, 0]
    print(f"  p shell {s + 1} (weight a = {norms[s]:.5f}):")
    for a_, d_ in zip(ex, d):
        print(f"      alpha = {a_:12.6f}   d = {d_: .6f}   N = {gto.gto_norm(1, a_): .6f}")

# --- check ψ(r) = sqrt(3/4π) (n·r) Σ_s a_s Σ_k d_sk N_k exp(-alpha_sk r²) ----
rng = np.random.default_rng(0)
pts = rng.normal(scale=1.5, size=(50, 3))                      # Bohr
r2 = (pts ** 2).sum(1)
radial = np.zeros(len(pts))
for s, ib in enumerate(p_shells):
    ex = mol.bas_exp(ib); d = mol.bas_ctr_coeff(ib)[:, 0]
    nrm = np.array([gto.gto_norm(1, a_) for a_ in ex])
    radial += norms[s] * (d[None, :] * nrm[None, :] * np.exp(-ex[None, :] * r2[:, None])).sum(1)
psi_formula = np.sqrt(3 / (4 * np.pi)) * (pts @ direction) * radial
ao = mol.eval_gto("GTOval_sph", pts)
c_p = np.zeros_like(c); c_p[p_aos] = c[p_aos]
psi_p = ao @ c_p        # PySCF value of the p part only
psi_full = ao @ c       # whole orbital, including the tiny non-p admixture
print("\nmax |formula - PySCF p part| over 50 random points:",
      f"{abs(psi_formula - psi_p).max():.2e}")
print("max |non-p part| (tiny f admixture from the finite grid):",
      f"{abs(psi_full - psi_p).max():.2e}")
print("typical |psi| at those points:", f"{abs(psi_full).mean():.2e}")
assert np.allclose(psi_formula, psi_p, atol=1e-10)

# --- cubes and xyz for PyMOL ---------------------------------------------
cubegen.orbital(mol, "C_p_orbital.cube", c, nx=80, ny=80, nz=80)
cubegen.density(mol, "C_density.cube", mf.make_rdm1()[0] + mf.make_rdm1()[1],
                nx=80, ny=80, nz=80)
with open("C.xyz", "w") as f:
    f.write("1\ncarbon atom\nC 0.0 0.0 0.0\n")
