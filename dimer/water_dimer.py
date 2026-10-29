"""Counterpoise-corrected interaction energy of a water dimer, plus the
difference density rho(AB) - rho(A) - rho(B).

Same level of theory as ../water.py: wB97M-V/def2-TZVPD (OMol25), PySCF.
Run from this folder:  python water_dimer.py
The geometry is built, not optimized, and the monomers stay rigid, so E_int
is only roughly the textbook water-dimer value (about -5 kcal/mol).
"""
import numpy as np
from pyscf import gto, dft
from pyscf.tools import cubegen

HA2KCAL = 627.5095
BASIS = "def2-tzvpd"
R_OO = 2.91  # O...O distance, Angstrom


def water(o, bisector, normal, r=0.9572, ang=104.52):
    """O at `o`; both H placed in the plane spanned by bisector and normal."""
    half = np.radians(ang / 2)
    b = bisector / np.linalg.norm(bisector)
    n = normal / np.linalg.norm(normal)
    h1 = o + r * (np.cos(half) * b + np.sin(half) * n)
    h2 = o + r * (np.cos(half) * b - np.sin(half) * n)
    return [("O", o), ("H", h1), ("H", h2)]


def dimer_atoms(r_oo=R_OO):
    x = np.array([1.0, 0, 0]); y = np.array([0, 1.0, 0]); z = np.array([0, 0, 1.0])
    # Donor D sits at +x with one H (H1) pointing straight at the acceptor O.
    # Its bisector is rotated by half the H-O-H angle off -x, so H1 lies on -x.
    half = np.radians(104.52 / 2)
    bis_d = np.cos(half) * x + np.sin(half) * y
    donor = water(r_oo * x, bis_d, y)
    # Rebuild so H1 is exactly along -x: place H1 explicitly, H2 from the angle.
    o_d = r_oo * x
    h1 = o_d - 0.9572 * x
    h2 = o_d + 0.9572 * (np.cos(np.radians(104.52)) * (-x) + np.sin(np.radians(104.52)) * z)
    donor = [("O", o_d), ("H", h1), ("H", h2)]
    # Acceptor A at the origin; its bisector points away from the donor,
    # tilted ~55 deg off the O...O axis, with its plane perpendicular to
    # the donor's plane (xz), i.e. the acceptor plane is xy.
    tilt = np.radians(125)
    bis_a = np.cos(tilt) * x + np.sin(tilt) * y
    acceptor = water(np.zeros(3), bis_a, z)
    return acceptor, donor


def fmt(atoms, ghost_idx=()):
    lines = []
    for i, (el, xyz) in enumerate(atoms):
        sym = f"ghost-{el}" if i in ghost_idx else el
        lines.append(f"{sym} {xyz[0]:.6f} {xyz[1]:.6f} {xyz[2]:.6f}")
    return "\n".join(lines)


def run(atoms, ghost_idx=(), label=""):
    mol = gto.M(atom=fmt(atoms, ghost_idx), basis=BASIS, verbose=0)
    mf = dft.RKS(mol)
    mf.xc = "wb97m_v"
    mf.nlc = "VV10" if not getattr(mf, "nlc", None) else mf.nlc
    mf.grids.atom_grid = (99, 590)
    e = mf.kernel()
    assert mf.converged, f"{label} did not converge"
    print(f"{label:28s} E = {e:.8f} Ha  ({mol.nao} basis functions)", flush=True)
    return mol, mf, e


def interaction(r_oo=R_OO, save=False):
    a, d = dimer_atoms(r_oo)
    ab = a + d
    A, B = (0, 1, 2), (3, 4, 5)
    mol_ab, mf_ab, e_ab = run(ab, label="dimer AB")
    # Monomers in their own basis (uncorrected).
    _, _, e_a = run(a, label="A, own basis")
    _, _, e_b = run(d, label="B, own basis")
    # Monomers in the dimer basis: partner replaced by ghost atoms.
    _, mf_ag, e_a_g = run(ab, ghost_idx=B, label="A, dimer basis (CP)")
    _, mf_bg, e_b_g = run(ab, ghost_idx=A, label="B, dimer basis (CP)")
    out = dict(
        r=r_oo,
        e_raw=(e_ab - e_a - e_b) * HA2KCAL,
        e_cp=(e_ab - e_a_g - e_b_g) * HA2KCAL,
    )
    out["bsse"] = out["e_cp"] - out["e_raw"]
    if save:
        with open("water_dimer.xyz", "w") as f:
            f.write(f"6\nwater dimer, O...O = {r_oo} A\n")
            for el, xyz in ab:
                f.write(f"{el} {xyz[0]:.6f} {xyz[1]:.6f} {xyz[2]:.6f}\n")
        # Same AO ordering in all three calculations, so the density
        # matrices can be subtracted directly.
        dm = mf_ab.make_rdm1() - mf_ag.make_rdm1() - mf_bg.make_rdm1()
        cubegen.density(mol_ab, "water_dimer_density.cube", mf_ab.make_rdm1(),
                        nx=100, ny=80, nz=80)
        cubegen.density(mol_ab, "water_dimer_diff.cube", dm,
                        nx=100, ny=80, nz=80)
    return out


if __name__ == "__main__":
    r = interaction(save=True)
    print(f"\nE_int uncorrected   = {r['e_raw']:7.3f} kcal/mol")
    print(f"E_int counterpoise  = {r['e_cp']:7.3f} kcal/mol")
    print(f"BSSE                = {r['bsse']:7.3f} kcal/mol")
