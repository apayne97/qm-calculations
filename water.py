import numpy as np
from pyscf import gto, dft
from pyscf.tools import cubegen
mol = gto.M(atom="""O 0 0 0.117790
H 0 0.755453 -0.471161
H 0 -0.755453 -0.471161""", basis="def2-tzvpd", spin=0, charge=0, verbose=3)
print("N basis functions:", mol.nao)
mf = dft.RKS(mol); mf.xc = "wb97m_v"
mf.nlc = "VV10" if not getattr(mf, "nlc", None) else mf.nlc
mf.grids.atom_grid = (99, 590)
e = mf.kernel()
print("E =", e, "converged", mf.converged)
dm = mf.make_rdm1()
cubegen.density(mol, "water_density.cube", dm, nx=80, ny=80, nz=80)
homo = mol.nelectron//2 - 1
cubegen.orbital(mol, "water_HOMO.cube", mf.mo_coeff[:, homo], nx=80, ny=80, nz=80)
cubegen.orbital(mol, "water_LUMO.cube", mf.mo_coeff[:, homo+1], nx=80, ny=80, nz=80)
print("orbital energies HOMO/LUMO (Ha):", mf.mo_energy[homo], mf.mo_energy[homo+1])
