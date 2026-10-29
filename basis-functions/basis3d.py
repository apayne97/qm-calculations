"""Write cube files for the def2-TZVPD basis functions of O and H (one
component per shell) and PyMOL scripts that build one session per atom, with
one scene per basis function. Run from this folder."""
import numpy as np
from pyscf import gto
from pyscf.tools import cubegen

BASIS, N, FRAC = "def2-tzvpd", 56, 0.25   # grid points per side, contour as fraction of max|phi|
ATOMS = {"O": 2, "H": 1}                  # element -> spin (only needed to build the mol)
PICK = ("s", "pz", "dz^2", "f+0")         # one component per shell

def shell_name(label):                     # '0 O 3pz    ' -> '3pz'
    return label.split()[2]

for el, spin in ATOMS.items():
    mol = gto.M(atom=f"{el} 0 0 0", basis=BASIS, spin=spin, verbose=0)
    open(f"{el}_atom.xyz", "w").write(f"1\n{el}\n{el} 0.0 0.0 0.0\n")
    labels = mol.ao_labels()
    # radial profile along z to size each box
    zs = np.linspace(0, 40, 4000)
    prof = mol.eval_gto("GTOval_sph", np.column_stack([np.zeros_like(zs), np.zeros_like(zs), zs]))
    funcs = []
    for i, lab in enumerate(labels):
        nm = shell_name(lab)
        if nm.endswith(PICK):
            funcs.append((i, nm))
    pml = [f"# {el} {BASIS} basis functions, one scene per function. Run from this folder.",
           "reinitialize", "bg_color white", f"load {el}_atom.xyz, nucleus",
           "show spheres, nucleus", "set sphere_scale, 0.05, nucleus", "color black, nucleus"]
    scenes = []
    for i, nm in funcs:
        safe = nm.replace("^", "").replace("+", "")
        v = np.abs(prof[:, i]); r_ext = zs[np.where(v > 0.02 * v.max())[0].max()]
        r_lvl = zs[np.where(v > FRAC * v.max())[0].max()]      # extent of the drawn surface, bohr
        zoom_buf = 0.529 * r_lvl * 1.1 + 0.5                   # Angstrom
        margin = float(min(max(1.2 * r_ext, 2.0), 30.0))
        cube = cubegen.Cube(mol, N, N, N, margin=margin)
        field = (mol.eval_gto("GTOval_sph", cube.get_coords())[:, i]).reshape(N, N, N)
        fn = f"{el}_{safe}_{i}.cube"; cube.write(field, fn, comment=f"{el} {nm} {BASIS}")
        lvl = FRAC * np.abs(field).max()
        obj = f"{el}_{safe}_{i}"
        pml += [f"load {fn}, m_{obj}", f"isosurface {obj}_pos, m_{obj}, {lvl:.5g}",
                f"color blue, {obj}_pos", f"set transparency, 0.3, {obj}_pos"]
        has_neg = field.min() < -lvl
        if has_neg:
            pml += [f"isosurface {obj}_neg, m_{obj}, {-lvl:.5g}", f"color red, {obj}_neg",
                    f"set transparency, 0.3, {obj}_neg"]
        scenes.append((obj, has_neg, nm, zoom_buf))
    pml += ["disable *", "orient nucleus", "turn x, -75", "turn y, 25"]
    for obj, has_neg, nm, buf in scenes:
        pml += ["disable *", "enable nucleus", f"enable {obj}_pos"] + ([f"enable {obj}_neg"] if has_neg else [])
        pml += [f"zoom nucleus, buffer={buf:.3g}", f"scene {obj}, store, {el} {nm}"]
    pml += [f"scene {scenes[0][0]}, recall", f"save {el}_basis.pse"]
    open(f"{el}_basis.pml", "w").write("\n".join(pml) + "\n")
    print(el, len(funcs), "functions:", [s[2] for s in scenes])
