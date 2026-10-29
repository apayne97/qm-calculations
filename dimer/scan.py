"""Distance scan of the water dimer (rigid monomers, donor translated along
the O...O axis). Scale factor 1.0 is the geometry in water_dimer.py. Writes
scan.csv (one row per point, so a crash keeps earlier points).
Run from this folder:  python scan.py
"""
import water_dimer as w

SCALES = [0.9, 1.0, 1.1, 1.25, 1.5, 2.0]

with open("scan.csv", "w") as f:
    f.write("scale,r_oo,e_int_raw,e_int_cp,bsse\n")
    for s in SCALES:
        r = w.interaction(w.R_OO * s)
        f.write(f"{s},{r['r']:.4f},{r['e_raw']:.4f},{r['e_cp']:.4f},{r['bsse']:.4f}\n")
        f.flush()
