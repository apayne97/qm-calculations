"""Same C-H bond orbital as methane.py, with s/p/d makeup on carbon from two
population analyses (Lowdin and Mulliken).
Run from this folder:  python analysis_methods.py
"""
import numpy as np
from pyscf import gto, dft, lo

R = 1.087 / np.sqrt(3)
atoms = [("C", (0, 0, 0)), ("H", (R, R, R)), ("H", (R, -R, -R)),
         ("H", (-R, R, -R)), ("H", (-R, -R, R))]
mol = gto.M(atom=atoms, basis="def2-tzvpd", verbose=0)
mf = dft.RKS(mol); mf.xc = "wb97m_v"
mf.nlc = "VV10" if not getattr(mf, "nlc", None) else mf.nlc
mf.grids.atom_grid = (99, 590)
mf.kernel()
S = mf.get_ovlp()
loc = lo.Boys(mol, mf.mo_coeff[:, 1:5]); loc.init_guess = "atomic"
c_all = loc.kernel()
labels = mol.ao_labels(fmt=False)
atom_of = np.array([l[0] for l in labels])
ls = np.array(["spdfg".index([ch for ch in l[2] if ch in "spdfg"][0]) for l in labels])

w, v = np.linalg.eigh(S)
S_half = v @ np.diag(np.sqrt(w)) @ v.T
# pick the orbital most on H1
J = int(np.argmax([((S_half @ c_all[:, j])[atom_of == 1] ** 2).sum() for j in range(4)]))
c = c_all[:, J]

def report(name, pop):
    onC = pop[atom_of == 0].sum()
    wl = [pop[(atom_of == 0) & (ls == l)].sum() for l in range(3)]
    t = sum(wl)
    print(f"{name:10s} carbon holds {onC:5.1%}; on carbon  s {wl[0]/t:5.1%}  "
          f"p {wl[1]/t:5.1%}  d {wl[2]/t:5.1%}")

report("Lowdin", (S_half @ c) ** 2)
report("Mulliken", c * (S @ c))
print("(Mulliken splits the overlap between two atoms 50/50, Lowdin first\n"
      " orthogonalizes the basis; neither is 'right', and Mulliken can even\n"
      " give negative populations.)")
