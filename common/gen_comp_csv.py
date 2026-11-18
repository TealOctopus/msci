"""
gen_comp_csv.py  -  match MOO and counterpoise (CP) transfer integrals for one snapshot
and write them side by side to a csv. No statistics or plots (see analyse_Ts.py).

Usage (on the HPC, in the msci conda env):
    python ~/MSci/common/gen_comp_csv.py Y6_20000ps_5A
    python ~/MSci/common/gen_comp_csv.py Y6_20000ps_5A --moo path/to/moo.csv --cp path/to/cp.csv

On Mac:
ln -s "$HOME/Documents/Uni/Year 4/Project/MSci" ~/MSci
Run this

Output (in the run folder):
    Ts_compare.csv   one row per dimer with both methods side by side (J in eV and meV)
"""
import argparse
import os
import sys

import pandas as pd

KEYS = ['init', 'fin', 'fin_lattice']
#dimer identified by init_mol, fin_mol and fin_lattice, so uses (init, fin, fin_lattice)
#as keys to merge the two file


def load(path, label):
    if not os.path.exists(path):
        sys.exit(f"Can't find the {label} results: {path}")
    df = pd.read_csv(path)
    #reads cv file into a df
    df[KEYS] = df[KEYS].astype(int)
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run', help='run name, e.g. Y6_20000ps_5A')
    #first arg given saved as run (this is name as run)
    ap.add_argument('--moo', help='MOO csv (default: runs/<run>/Ts_MOO.csv)')
    ap.add_argument('--cp',  help='CP csv  (default: runs/<run>/Ts_CP.csv)')
    #these are optional, in case runs are stored somewhere else
    a = ap.parse_args()

    run_dir = os.path.expanduser(f'~/MSci/runs/{a.run}')
    moo = load(a.moo or os.path.join(run_dir, 'Ts_MOO.csv'), 'MOO')
    cp  = load(a.cp  or os.path.join(run_dir, 'Ts_CP.csv'),  'CP')
    #locate the files storing the Ts

    # both files may have 'dist'; check they agree, then keep only one copy
    # 'dist' may be in MOO only (older CP runs) or in both (after the CalcJs.py change)
    if 'dist' in moo and 'dist' in cp:
        # both have it: check they agree, then keep one copy so the merge doesn't
        # rename them dist_x / dist_y
        chk = moo[KEYS + ['dist']].merge(cp[KEYS + ['dist']], on=KEYS, suffixes=('_moo', '_cp'))
        bad = (chk.dist_moo - chk.dist_cp).abs() > 1e-3          # Å
        if bad.any():
            print(f"  WARNING: dist differs between MOO and CP for {bad.sum()} pairs:")
            print(chk[bad].to_string(index=False))
        cp = cp.drop(columns='dist')
    # if only one file has 'dist' there's no name clash, so nothing to do

    # match dimers on (init, fin, fin_lattice); outer join so we can see what's missing
    m = moo.merge(cp, on=KEYS, how='outer', indicator=True)
    #outer keeps dimer that has been in either file

    only_moo = m[m._merge == 'left_only']
    only_cp  = m[m._merge == 'right_only']
    both     = m[m._merge == 'both'].drop(columns='_merge').copy()
    #these lines arent v. important, they just keep track of where
    #the values are coming from, and if a comparison is missing either
    #MOO or CSV, we can quickly identify what job failed


    print(f"\nRun {a.run}: {len(moo)} MOO pairs, {len(cp)} CP pairs, {len(both)} in both")
    for name, d in [('MOO', only_moo), ('CP', only_cp)]:
        if len(d):
            print(f"  only in {name} (probably failed in the other):")
            for r in d.itertuples():
                print(f"    Mol{r.init}  Mol{r.fin}  Lattice{r.fin_lattice}")
    #lists likely failed jobs

    for c in ['Th_MOO', 'Te_MOO', 'Th_CP', 'Te_CP']:
        both[c + '_meV'] = 1000 * both[c]
    #unit conversion

    out_csv = os.path.join(run_dir, 'Ts_compare.csv')
    both.sort_values('dist' if 'dist' in both else KEYS).to_csv(out_csv, index=False, float_format='%.6f')
    #this just sorts the df, either by dist or if that fails, by keys
    print(f"Wrote {out_csv}")


if __name__ == '__main__':
    main()