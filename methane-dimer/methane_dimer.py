"""Counterpoise-corrected methane dimer interaction energy with and without
dispersion. Run from this folder:  python methane_dimer.py [R_CC ...]
Geometry: two rigid methanes (C-H 1.087 A) with three H facing each other
(staggered, D3d), related by inversion, so the two monomers are equivalent and
only one ghost-atom calculation is needed. Not optimized.
"""
import sys, time
import numpy as np
from pyscf import gto, dft, scf

HA2KCAL = 627.5095
BASIS = "def2-tzvpd"
R_CH = 1.087


def monomer():
    z_apex = -R_CH
    z_tri = R_CH / 3
    rad = R_CH * np.sqrt(8) / 3
    atoms = [("C", np.zeros(3)), ("H", np.array([0, 0, z_apex]))]
    for k in range(3):
        a = np.radians(120 * k)
        atoms.append(("H", np.array([rad * np.cos(a), rad * np.sin(a), z_tri])))
    return atoms


def dimer(r_cc):
    a = monomer()
    # B = inversion of A through the midpoint (0, 0, r_cc/2)
    b = [(el, np.array([0, 0, r_cc]) - xyz) for el, xyz in a]
    return a, b


def fmt(atoms, ghost=()):
    return "\n".join(
        f"{'ghost-' if i in ghost else ''}{el} {x[0]:.6f} {x[1]:.6f} {x[2]:.6f}"
        for i, (el, x) in enumerate(atoms))


def energy(atoms, method, ghost=()):
    mol = gto.M(atom=fmt(atoms, ghost), basis=BASIS, verbose=0)
    if method == "HF":
        mf = scf.RHF(mol)
    else:
        mf = dft.RKS(mol)
        mf.grids.atom_grid = (99, 590)
        if method == "B3LYP":
            mf.xc = "b3lyp"
        else:
            mf.xc = "wb97m_v"
            if method == "wB97M-V":
                mf.nlc = "VV10" if not getattr(mf, "nlc", None) else mf.nlc
            else:   # "wB97M-V, VV10 off": nlc=False disables it (nlc="" would defer to the xc name)
                mf.nlc = False
    mf.kernel()
    assert mf.converged
    return mf.e_tot


METHODS = ["HF", "B3LYP", "wB97M-V, VV10 off", "wB97M-V"]
if __name__ == "__main__":
    radii = [float(x) for x in sys.argv[1:]] or [3.8]
    print("R_CC(A)  method               E_int raw  E_int CP  (kcal/mol)")
    for r in radii:
        a, b = dimer(r)
        ab = a + b
        for m in METHODS:
            t = time.time()
            e_ab = energy(ab, m)
            e_a = energy(a, m)                       # monomers are equivalent
            e_a_g = energy(ab, m, ghost=range(5, 10))
            raw = (e_ab - 2 * e_a) * HA2KCAL
            cp = (e_ab - 2 * e_a_g) * HA2KCAL
            print(f"{r:6.2f}   {m:20s} {raw:8.3f} {cp:9.3f}   ({time.time()-t:.0f}s)", flush=True)
