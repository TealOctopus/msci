"""
compare_Ts.py  -  compare MOO and counterpoise (CP) transfer integrals for one snapshot.

Usage (on the HPC, in the msci conda env):
    python ~/MSci/common/compare_Ts.py Y6_20000ps_5A
    python ~/MSci/common/compare_Ts.py Y6_20000ps_5A --moo path/to/moo.csv --cp path/to/cp.csv

Defaults:
    MOO csv : ~/MSci/runs/<run>/Ts_MOO.csv   (written by RunAllPairs_OneSnapshot.sh)
    CP csv  : ~/MSci/runs/<run>/Ts_CP.csv    (written by CalcJs.py)

Outputs (in the run folder):
    Ts_compare.csv   one row per dimer with both methods side by side
    Ts_compare.png   |J_MOO| vs |J_CP| for holes and electrons (log-log)

Notes:
    - Both codes give J in eV; this script converts to meV.
    - The sign of J depends on the arbitrary phase of each monomer orbital,
      so the comparison uses |J|. Sign agreement is printed for information only.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

KEYS = ['init', 'fin', 'fin_lattice']


def load(path, label):
    if not os.path.exists(path):
        sys.exit(f"Can't find the {label} results: {path}")
    df = pd.read_csv(path)
    df[KEYS] = df[KEYS].astype(int)
    return df


def stats(x, y):
    """x = CP (reference), y = MOO. Both |J| in meV."""
    ok = (x > 0) & (y > 0)
    return {
        'N':              len(x),
        'Pearson r':      np.corrcoef(x, y)[0, 1],
        'Pearson r (log10)': np.corrcoef(np.log10(x[ok]), np.log10(y[ok]))[0, 1],
        'MAE (meV)':      np.mean(np.abs(y - x)),
        'RMSE (meV)':     np.sqrt(np.mean((y - x) ** 2)),
        'median MOO/CP':  np.median(y[ok] / x[ok]),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run', help='run name, e.g. Y6_20000ps_5A')
    ap.add_argument('--moo', help='MOO csv (default: runs/<run>/Ts_MOO.csv)')
    ap.add_argument('--cp',  help='CP csv  (default: runs/<run>/Ts_CP.csv)')
    ap.add_argument('--noplot', action='store_true')
    a = ap.parse_args()

    run_dir = os.path.expanduser(f'~/MSci/runs/{a.run}')
    moo = load(a.moo or os.path.join(run_dir, 'Ts_MOO.csv'), 'MOO')
    cp  = load(a.cp  or os.path.join(run_dir, 'Ts_CP.csv'),  'CP')

    # match dimers on (init, fin, fin_lattice); outer join so we can see what's missing
    m = moo.merge(cp, on=KEYS, how='outer', indicator=True)
    only_moo = m[m._merge == 'left_only']
    only_cp  = m[m._merge == 'right_only']
    both     = m[m._merge == 'both'].drop(columns='_merge').copy()

    print(f"\nRun {a.run}: {len(moo)} MOO pairs, {len(cp)} CP pairs, {len(both)} in both")
    for name, d in [('MOO', only_moo), ('CP', only_cp)]:
        if len(d):
            print(f"  only in {name} (probably failed in the other):")
            for r in d.itertuples():
                print(f"    Mol{r.init}  Mol{r.fin}  Lattice{r.fin_lattice}")
    if len(both) < 2:
        sys.exit("Not enough matched pairs to compare.")

    for c in ['Th_MOO', 'Te_MOO', 'Th_CP', 'Te_CP']:
        both[c + '_meV'] = 1000 * both[c]

    print()
    summary = {}
    for carrier, lab in [('h', 'holes (HOMO-HOMO)'), ('e', 'electrons (LUMO-LUMO)')]:
        x = both[f'T{carrier}_CP_meV'].abs().values
        y = both[f'T{carrier}_MOO_meV'].abs().values
        s = stats(x, y)
        s['same sign (info)'] = np.mean(np.sign(both[f'T{carrier}_CP']) == np.sign(both[f'T{carrier}_MOO']))
        summary[lab] = s
    print(pd.DataFrame(summary).round(3).to_string())

    out_csv = os.path.join(run_dir, 'Ts_compare.csv')
    both.sort_values('dist' if 'dist' in both else KEYS).to_csv(out_csv, index=False, float_format='%.6f')
    print(f"\nWrote {out_csv}")

    if a.noplot:
        return
    import matplotlib
    matplotlib.use('Agg')            # no display on the HPC
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.6))
    for ax, carrier, title, lab in zip(axes, ['h', 'e'], ['Holes', 'Electrons'], summary):
        x = both[f'T{carrier}_CP_meV'].abs()
        y = both[f'T{carrier}_MOO_meV'].abs()
        c = both['dist'] if 'dist' in both else None
        sc = ax.scatter(x, y, c=c, s=18, cmap='viridis')
        pos = pd.concat([x, y])
        pos = pos[pos > 0]
        lo, hi = pos.min() / 2, pos.max() * 2
        ax.plot([lo, hi], [lo, hi], 'k--', lw=0.8, label='y = x')
        ax.set(xscale='log', yscale='log', xlim=(lo, hi), ylim=(lo, hi),
               xlabel='|J| counterpoise (meV)', ylabel='|J| MOO (meV)',
               title=f"{title}  (r_log = {summary[lab]['Pearson r (log10)']:.2f})")
        ax.legend(loc='upper left', frameon=False)
        if c is not None:
            fig.colorbar(sc, ax=ax, label='closest C–C distance (Å)')
    fig.suptitle(a.run)
    fig.tight_layout()
    out_png = os.path.join(run_dir, 'Ts_compare.png')
    fig.savefig(out_png, dpi=150)
    print(f"Wrote {out_png}")


if __name__ == '__main__':
    main()
