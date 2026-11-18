"""
gen_features_csv.py  -  geometric features for every dimer in one snapshot.

Usage (on CX3 or Mac, in the msci env):
    python ~/MSci/common/gen_features_csv.py Y6_20000ps_5A

Reads (all written by prepare_snapshot when the run was set up):
    init.file, fin.file, finlattice.file, dist.txt   the dimer list and closest C-C distance
    MoleculesPBC/Mol*_Lattice*.xyz                   monomer geometries

Output (in the run folder):
    Ts_features.csv   one row per dimer: init, fin, fin_lattice + features
                      (join to Ts_compare.csv on init, fin, fin_lattice)
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

# ---------- geometry: read a monomer and get its plane / long axis ----------

def read_xyz(path, heavy_only=True):
    xyz = pd.read_csv(path, sep=r'\s+', header=None, names=['El', 'X', 'Y', 'Z'])
    if heavy_only:
        xyz = xyz[xyz.El != 'H']        # capping H's would distort the plane fit
    return xyz[['X', 'Y', 'Z']].values


def frame(coords):
    """Centroid, long axis and plane normal from a principal-axis (SVD) fit.
    Rows of vt are ordered by how spread out the atoms are along them:
    vt[0] = longest direction (long axis), vt[2] = least spread (plane normal)."""
    centroid = coords.mean(axis=0)
    _, _, vt = np.linalg.svd(coords - centroid)
    return centroid, vt[0], vt[2]


def angle(u, v):
    """Angle between two axes in degrees, 0-90. The sign of an axis from the
    SVD is arbitrary, so we use |u.v|."""
    return np.degrees(np.arccos(np.clip(abs(np.dot(u, v)), 0, 1)))


def geometry(run_dir, init, fin, fin_lattice):
    folder = os.path.join(run_dir, 'MoleculesPBC')
    cA, longA, nA = frame(read_xyz(os.path.join(folder, f'Mol{init}_Lattice0.xyz')))
    cB, longB, nB = frame(read_xyz(os.path.join(folder, f'Mol{fin}_Lattice{fin_lattice}.xyz')))
    d = cB - cA
    perp = abs(np.dot(d, nA))                     # separation along A's plane normal
    return {
        'centroid_dist':   np.linalg.norm(d),     # Å
        'plane_angle':     angle(nA, nB),         # 0 = parallel planes (cofacial)
        'long_axis_angle': angle(longA, longB),   # 0 = backbones lie the same way
        'stack_dist':      perp,                  # Å, ~3.5 for a π-stack
        'slip':            np.sqrt(max(np.dot(d, d) - perp**2, 0)),  # Å, in-plane offset
    }


def load_pairs(run_dir):
    path = os.path.join(run_dir, 'Ts_compare.csv')
    if not os.path.exists(path):
        sys.exit(f"Can't find {path} - run gen_comp_csv.py first")
    return pd.read_csv(path, usecols=['init', 'fin', 'fin_lattice', 'dist'])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run', help='run name, e.g. Y6_20000ps_5A')
    a = ap.parse_args()

    run_dir = os.path.expanduser(f'~/MSci/runs/{a.run}')
    pairs = load_pairs(run_dir)

    rows = []
    for p in pairs.itertuples(index=False):
        try:
            g = geometry(run_dir, p.init, p.fin, p.fin_lattice)
        except FileNotFoundError as e:
            print(f"  skipping Mol{p.init} Mol{p.fin} Lattice{p.fin_lattice}: {e}")
            continue
        rows.append({'init': p.init, 'fin': p.fin, 'fin_lattice': p.fin_lattice,
                     'dist': p.dist, **g})

    feats = pd.DataFrame(rows).sort_values('dist')
    out_csv = os.path.join(run_dir, 'Ts_features.csv')
    feats.to_csv(out_csv, index=False, float_format='%.6f')
    print(f"Wrote {len(feats)} of {len(pairs)} dimers to {out_csv}")


if __name__ == '__main__':
    main()