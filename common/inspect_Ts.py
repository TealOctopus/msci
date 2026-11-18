"""
inspect_Ts.py - find the dimers where MOO and CP agree least (or most).

Usage (on CX3, in the msci env):
    python inspect_Ts.py Y6_20000ps_5A                  # 10 worst, holes
    python inspect_Ts.py Y6_20000ps_5A -n 5 -c e       # 5 worst, electrons
    python inspect_Ts.py Y6_20000ps_5A -n 5 --best     # 5 most accurate
    python inspect_Ts.py Y6_20000ps_5A --min_cp 5      # only pairs with |J_CP| > 5 meV
"""
import argparse
import os

import numpy as np
import pandas as pd


# ---------- 1. read the csv and add the fractional difference ----------

def load_compare(run):
    run_dir = os.path.expanduser(f'~/MSci/runs/{run}')
    df = pd.read_csv(os.path.join(run_dir, 'Ts_compare.csv'))
    for c in ['h', 'e']:
        cp = df[f'T{c}_CP'].abs()
        moo = df[f'T{c}_MOO'].abs()
        # (|MOO| - |CP|) / |CP| : +0.5 means MOO is 50% too big, -0.5 means 50% too small
        df[f'frac_diff_{c}'] = (moo - cp) / cp
    return df, run_dir


# ---------- 2. geometry: read a monomer and get its plane / long axis ----------

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
        'centroid_dist': np.linalg.norm(d),       # Å
        'plane_angle':   angle(nA, nB),           # 0 = parallel planes (cofacial)
        'long_axis_angle': angle(longA, longB),   # 0 = backbones lie the same way
        'stack_dist':    perp,                    # Å, ~3.5 for a pi-stack
        'slip':          np.sqrt(max(np.dot(d, d) - perp**2, 0)),  # Å, in-plane offset
    }


# ---------- 3. pick the N worst / best and describe them ----------

def show(df, run_dir, n=10, carrier='h', best=False, min_cp=1.0):
    col = f'frac_diff_{carrier}'
    # ignore pairs whose CP coupling is too small for the fractional difference to mean much
    keep = 1000 * df[f'T{carrier}_CP'].abs() > min_cp
    print(f'\nKeeping {keep.sum()} of {len(df)} pairs with |T{carrier}_CP| > {min_cp:g} meV')
    ranked = (df[keep].assign(abs_frac=df[col].abs())
                      .sort_values('abs_frac', ascending=best).head(n))

    rows = []
    for r in ranked.itertuples():
        g = geometry(run_dir, r.init, r.fin, r.fin_lattice)
        rows.append({
            'init': r.init, 'fin': r.fin, 'fin_lattice': r.fin_lattice,
            'dist': r.dist,
            f'T{carrier}_CP_meV':  1000 * getattr(r, f'T{carrier}_CP'),
            f'T{carrier}_MOO_meV': 1000 * getattr(r, f'T{carrier}_MOO'),
            col: getattr(r, col),
            **g,
        })
    out = pd.DataFrame(rows)

    which = 'BEST' if best else 'WORST'
    label = 'holes' if carrier == 'h' else 'electrons'
    print(f'\n{which} {n} pairs for {label} (dist, centroid_dist, stack_dist, slip in Å; angles in degrees)')
    print(out.round(3).to_string(index=False))
    return out


def worst(df, run_dir, n=10, carrier='h', min_cp=1.0):
    return show(df, run_dir, n, carrier, best=False, min_cp=min_cp)


def best(df, run_dir, n=10, carrier='h', min_cp=1.0):
    return show(df, run_dir, n, carrier, best=True, min_cp=min_cp)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('run', help='e.g. Y6_20000ps_5A')
    ap.add_argument('-n', type=int, default=10, help='number of pairs to show')
    ap.add_argument('-c', '--carrier', choices=['h', 'e'], default='h')
    ap.add_argument('--best', action='store_true', help='show the most accurate pairs instead')
    ap.add_argument('--min_cp', type=float, default=1.0,
                    help='only rank pairs with |J_CP| above this, in meV (default 1; 0 = no filter)')
    a = ap.parse_args()

    df, run_dir = load_compare(a.run)
    show(df, run_dir, a.n, a.carrier, a.best, a.min_cp)