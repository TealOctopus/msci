from CounterpoiseJarvVersion import CalcJ
import os
import re
import sys
import numpy as np
from multiprocessing import Pool

run_name = sys.argv[1]  # e.g. Y6_20000ps_5A
data_file = os.path.join(os.environ['EPHEMERAL'], 'outputs_cp', run_name)
subfiles = sorted(d for d in os.listdir(data_file) if os.path.isdir(os.path.join(data_file, d)))


def one_dimer(file):
    """Calculate Th, Te for one dimer folder. Returns None if it fails."""
    mols = [int(s) for s in re.findall(r'\d+', file)]   # [init, 0, fin, fin_lattice, t]
    try:
        Th, Te = CalcJ(data_file, file)
    except Exception as e:
        print(f"FAILED {file}: {e}", flush=True)
        return None
    return [mols[0], mols[2], mols[3], mols[4], Th, Te]

def load_dists(run_name):
    """{(init, fin, fin_lattice): closest C-C distance in Å}, from the files prepare_snapshot wrote."""
    run_dir = os.path.expanduser(f'~/MSci/runs/{run_name}')
    init = np.loadtxt(os.path.join(run_dir, 'init.file'),       dtype=int, ndmin=1)
    fin  = np.loadtxt(os.path.join(run_dir, 'fin.file'),        dtype=int, ndmin=1)
    fl   = np.loadtxt(os.path.join(run_dir, 'finlattice.file'), dtype=int, ndmin=1)
    dist = np.loadtxt(os.path.join(run_dir, 'dist.txt'),        ndmin=1)
    return dict(zip(zip(init, fin, fl), dist))


if __name__ == '__main__':
    nproc = int(os.environ.get('NCPUS', 1))   # PBS sets NCPUS from the select line
    print(f"{len(subfiles)} dimers, {nproc} workers", flush=True)
    with Pool(nproc) as pool:
        results = pool.map(one_dimer, subfiles)   # keeps the original order
    rows = [r for r in results if r is not None]

    # add the closest C-C distance for each dimer (NaN if it can't be found)
    dists = load_dists(run_name)
    rows = [r[:3] + [dists.get(tuple(r[:3]), np.nan)] + r[3:] for r in rows]

    out = os.path.expanduser(f'~/MSci/runs/{run_name}/Ts_CP.csv')
    np.savetxt(out, np.array(rows), delimiter=',', comments='',
               fmt=['%d', '%d', '%d', '%.8f', '%d', '%.8f', '%.8f'],
               header='init,fin,fin_lattice,dist,t_ps,Th_CP,Te_CP')
    print(f"Wrote {len(rows)} dimers to {out}")