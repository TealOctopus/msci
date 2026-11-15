import sys, os
import numpy as np
import pandas as pd
import generate_CP_GJF as cp_gjf

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, '..', 'common'))
import gen_dimer_and_mols as gd

METHOD = ('LC-wHPBE/6-31G(d,p) IOp(3/107=0140000000,3/108=0140000000,'
          '3/119=0,3/120=0, 3/130=0, 3/131=0)')

def read_xyz(path):
    return pd.read_csv(path, sep='\t', header=None,
                       names=['atom', 'X', 'Y', 'Z']).round(6)

if __name__ == '__main__':
    mol, T, cutoff = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    run_dir = os.path.expanduser(f'~/MSci/runs/{mol}_{T}ps_{cutoff}A')
    out_dir = os.path.join(run_dir, 'CPJobs') + '/'      # write_CP_com needs the trailing /

    os.makedirs(run_dir, exist_ok=True)
    os.chdir(run_dir)                       # prepare_snapshot writes into the cwd
    if not os.path.exists('init.file'):
        gd.prepare_snapshot(mol, T, float(cutoff))
    #if the file doesnt already exist, generate the file


    init = np.loadtxt(os.path.join(run_dir, 'init.file'),       dtype=int, ndmin=1)
    fin  = np.loadtxt(os.path.join(run_dir, 'fin.file'),        dtype=int, ndmin=1)
    fl   = np.loadtxt(os.path.join(run_dir, 'finlattice.file'), dtype=int, ndmin=1)

    for i in range(len(init)):
        a = f"Mol{init[i]}_Lattice0"
        b = f"Mol{fin[i]}_Lattice{fl[i]}"
        geomA = read_xyz(os.path.join(run_dir, 'MoleculesPBC', a + '.xyz'))
        geomB = read_xyz(os.path.join(run_dir, 'MoleculesPBC', b + '.xyz'))
        cp_gjf.write_CP_com(geomA, geomB, f"{a}_{b}_{T}ps", out_dir,
                            nprocshared=16, mem='10GB', method=METHOD)
    print(f"{len(init)} pairs. Submit with:")
    print(f"qsub -J 1-{3*len(init)} -v MOL={mol},T={T},CUTOFF={cutoff} batch_cp.pbs")