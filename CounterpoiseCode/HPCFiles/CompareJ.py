"""Compare the J_diagnosis.npz results of two or more dimers side by side.

Usage:  python CompareJ.py <root folder> <dimer 1> <dimer 2> [...]
  e.g.  python CompareJ.py $EPHEMERAL/outputs_cp/Y6_20000ps_5A Mol116_Lattice0_Mol162_Lattice0_20000ps Mol192_Lattice0_Mol203_Lattice0_20000ps
"""
import os
import sys
import numpy as np

np.set_printoptions(precision=2, suppress=True, linewidth=120)

root = sys.argv[1]
dimers = sys.argv[2:]
data = {d: np.load(os.path.join(root, d, "J_diagnosis.npz")) for d in dimers}
H, L = 2, 3                                   # positions of HOMO and LUMO in the arrays


def short(name):
    """Mol116_Lattice0_Mol162_Lattice0_20000ps -> 116-162 (L0)"""
    p = name.split("_")
    return f"{p[0][3:]}-{p[2][3:]} (L{p[3][7:]})"


# ---------- summary table: one column per dimer ----------
rows = {
    "nbasis A, B, AB":                lambda r: str(r["nbasis"]),
    "nmo    A, B, AB":                lambda r: str(r["nmo"]),
    "J HOMO-HOMO (meV)":              lambda r: f"{r['J_eff'][H, H]:.1f}",
    "J LUMO-LUMO (meV)":              lambda r: f"{r['J_eff'][L, L]:.1f}",
    "A: HOMO - HOMO-1 (meV)":         lambda r: f"{1000*(r['E_A'][H] - r['E_A'][H-1]):.0f}",
    "B: HOMO - HOMO-1 (meV)":         lambda r: f"{1000*(r['E_B'][H] - r['E_B'][H-1]):.0f}",
    "A HOMO - B HOMO (meV)":          lambda r: f"{1000*(r['E_A'][H] - r['E_B'][H]):.0f}",
    "Pair HOMO - HOMO-1 (meV)":       lambda r: f"{1000*(r['E_AB'][H] - r['E_AB'][H-1]):.0f}",
    "2|J HOMO-HOMO| (meV)":           lambda r: f"{2*abs(r['J_eff'][H, H]):.0f}",
    "Pair LUMO+1 - LUMO (meV)":       lambda r: f"{1000*(r['E_AB'][L+1] - r['E_AB'][L]):.0f}",
    "2|J LUMO-LUMO| (meV)":           lambda r: f"{2*abs(r['J_eff'][L, L]):.0f}",
    "S HOMO-HOMO":                    lambda r: f"{r['S_MO'][H, H]:.4f}",
    "S LUMO-LUMO":                    lambda r: f"{r['S_MO'][L, L]:.4f}",
}

w = 22
print(f"{'':28s}" + "".join(f"{short(d):>{w}s}" for d in dimers))
for name, f in rows.items():
    print(f"{name:28s}" + "".join(f"{f(data[d]):>{w}s}" for d in dimers))

# ---------- full coupling matrices ----------
for d in dimers:
    r = data[d]
    labels = list(r["labels"])
    print(f"\nJ_eff (meV) for {short(d)}: rows = orbital of A, columns = orbital of B")
    print(f"{'':6s}" + "".join(f"{l:>9s}" for l in labels))
    for lab, row in zip(labels, r["J_eff"]):
        print(f"{lab:6s}" + "".join(f"{x:9.1f}" for x in row))