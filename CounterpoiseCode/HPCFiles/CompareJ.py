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
        "A HOMO: total weight":           lambda r: f"{r['wA_HOMO_total']:.4f}",
    "B HOMO: total weight":           lambda r: f"{r['wB_HOMO_total']:.4f}",
    "A HOMO: in pair HOMO + HOMO-1":  lambda r: f"{r['wA_HOMO'][3] + r['wA_HOMO'][4]:.3f}",
    "B HOMO: in pair HOMO + HOMO-1":  lambda r: f"{r['wB_HOMO'][3] + r['wB_HOMO'][4]:.3f}",
    "A LUMO: in pair LUMO + LUMO+1":  lambda r: f"{r['wA_LUMO'][5] + r['wA_LUMO'][6]:.3f}",
    "B LUMO: in pair LUMO + LUMO+1":  lambda r: f"{r['wB_LUMO'][5] + r['wB_LUMO'][6]:.3f}",
    "Site energy diff, HOMO (meV)":   lambda r: f"{1000*(r['site_HOMO'][0] - r['site_HOMO'][1]):.0f}",
    "Predicted HOMO splitting (meV)": lambda r: f"{1000*np.sqrt((r['site_HOMO'][0]-r['site_HOMO'][1])**2 + 4*(r['J_eff'][H,H]/1000)**2):.0f}",
    "A occupied found in pair empty": lambda r: f"{r['occA_in_empty']:.3f}",
    "B occupied found in pair empty": lambda r: f"{r['occB_in_empty']:.3f}",
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

#This produces matrix M:
#M_ij = <phi_i^A|H|phi_j^B> , in MeV
#ie it contains all the possible couplings between the six frontier orbitals
#HOMO-HOMO is Th, LUMO-LUMO is Te

#in two level picture, HOMO and HOMO-1 are split by:
#sqrt(deltaE^2 +4 T^2) (T=Th in this case)
#so the splitting must be greater than 2J
# ---------- where A's and B's frontier orbitals end up among the pair's orbitals ----------
for d in dimers:
    r = data[d]
    names = []
    for p in r["pair_window"]:
        if p < 0:
            names.append(f"H{p}")
        elif p == 0:
            names.append("HOMO")
        elif p == 1:
            names.append("LUMO")
        else:
            names.append(f"L+{p-1}")
    print(f"\nWeight of each orbital in the PAIR's orbitals, for {short(d)}")
    print(f"{'':8s}" + "".join(f"{n:>7s}" for n in names))
    for key in ["wA_HOMO", "wB_HOMO", "wA_LUMO", "wB_LUMO"]:
        print(f"{key:8s}" + "".join(f"{x:7.3f}" for x in r[key]))

for d in dimers:
    r = data[d]
    labels = list(r["labels"])
    names = []
    for p in r["pair_window"]:
        if p < 0:
            names.append(f"H{p}")
        elif p == 0:
            names.append("HOMO")
        elif p == 1:
            names.append("LUMO")
        else:
            names.append(f"L+{p-1}")
    for mol, key, site in [("A", "compA", "siteA_window"), ("B", "compB", "siteB_window")]:
        print(f"\n{short(d)}: fraction of molecule {mol}'s orbitals (rows) in each PAIR orbital (columns)")
        print(f"{'':6s}" + "".join(f"{n:>7s}" for n in names) + f"{'site (eV)':>11s}")
        for lab, row, e in zip(labels, r[key], r[site]):
            print(f"{lab:6s}" + "".join(f"{x:7.3f}" for x in row) + f"{e:11.3f}")

# ---------- what the pair's HOMO is made of ----------
def orb_name(p):
    if p < 0:
        return f"H{p}"
    elif p == 0:
        return "HOMO"
    elif p == 1:
        return "LUMO"
    else:
        return f"L+{p-1}"

for d in dimers:
    r = data[d]
    print(f"\nPair HOMO of {short(d)} (energy {r['E_pairHOMO']:.3f} eV)")
    print(f"  fraction in A's occupied orbitals: {r['pairHOMO_Aocc']:.3f}   in B's occupied orbitals: {r['pairHOMO_Bocc']:.3f}")
    print("  largest from A: " + ", ".join(f"{orb_name(p)} {w:.3f}" for p, w in zip(r["pairHOMO_topA"], r["pairHOMO_wA"])))
    print("  largest from B: " + ", ".join(f"{orb_name(p)} {w:.3f}" for p, w in zip(r["pairHOMO_topB"], r["pairHOMO_wB"])))