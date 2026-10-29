import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pyscf import gto
mol = gto.M(atom="O 0 0 0", basis="def2-tzvpd", spin=2, verbose=0)
labels = mol.ao_labels()
z = np.linspace(0, 5, 600)
pts = np.column_stack([np.zeros_like(z), np.zeros_like(z), z*1.8897259886])  # Angstrom->bohr
ao = mol.eval_gto("GTOval_sph", pts)
groups = {"s": lambda l: l.split()[2].endswith("s"), "p (pz)": lambda l: "pz" in l, "d (dz^2)": lambda l: "dz^2" in l, "f (f0)": lambda l: "f+0" in l or "f0" in l}
fig, axs = plt.subplots(1, 4, figsize=(16, 3.8), sharex=True)
for ax, (name, f) in zip(axs, groups.items()):
    idx = [i for i, l in enumerate(labels) if f(l)]
    for i in idx: ax.plot(z, ao[:, i], label=labels[i].strip())
    ax.set_title(f"Oxygen {name} functions"); ax.set_xlabel("distance from O along z (Å)"); ax.axhline(0, c="k", lw=.5)
axs[0].set_ylabel("basis function value"); axs[0].set_ylim(-2, 8)
fig.suptitle("def2-TZVPD basis functions on oxygen"); fig.tight_layout(); fig.savefig("oxygen_basis.png", dpi=140)
print(len(labels), labels[:12])
