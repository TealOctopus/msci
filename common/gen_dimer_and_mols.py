from prep import gro_to_atoms, truncate_and_cap, get_box_vectors, pair_finder, build_supercell
import os
import sys
import numpy as np
import pandas as pd



def write_monomer_xyz(df, mol, Latt):
    #df is the input, the atoms of one copy of one molecule picked out of bigmolecules
    #func writes a new file, with one line per atom
    fname = "MoleculesPBC/Mol{:.0f}_Lattice{:.0f}.xyz".format(mol, Latt)
    df[['Symbol', 'X', 'Y', 'Z']].to_csv(fname, sep='\t', header=False,
                                            index=False, float_format='%.7f')
    with open(fname, 'a') as f:
        f.write('\n')   # blank line for Gaussian


def write_dimer_xyz(mol1_df, mol2_df, mol1, mol2, fin_lattice):
    fname = f"MoleculesPBC/Mol{mol1:.0f}_Mol{mol2:.0f}_Lattice{fin_lattice:.0f}.xyz"
    dimer = pd.concat([mol1_df, mol2_df], ignore_index=True)
    dimer[['Symbol', 'X', 'Y', 'Z']].to_csv(fname, sep='\t', header=False,
                                        index=False, float_format='%.7f')
    with open(fname, 'a') as f:
        f.write('\n')   # blank line at end, required by Gaussian


def prepare_snapshot(mol_name, T, cutoff=5.0,save_val=False): #molecule eg 'Y6, snapshot time, cutoff distance
    
    datafile = os.path.expanduser(f"~/MSci/Structures/{mol_name}")
    fname = f"{mol_name}_300K_{T}ps.gro"
    file = os.path.join(datafile,fname)

    b_vectors = get_box_vectors(file)

    df = gro_to_atoms(file, mol_name) 
    #Reads a gromacs .gro file and returns a pd datatframe
    #df has one row per atom of the molecule
    #columns are what moleucles the atom belongs to, the element symbol, 
    #the atomic weight, and the x,y,z coordinates

    molecules = truncate_and_cap(df, mol_name)
    #strips the molecules down to the backbone and adds capping hydrogens
    #only contains the atoms of the original lattice
    bigmolecules = build_supercell(molecules, b_vectors)
    #bigmolecules dataframe - holds every atom of every capped molecule in all 27 boxes

    if save_val:
        unfilt_pairs = pair_finder(bigmolecules, limit=save_val)
        pairs = unfilt_pairs[unfilt_pairs.dist<=cutoff].reset_index(drop=True)
        unfilt_pairs.to_csv(f'pairs_{mol_name}_{T}ps_upto{save_val:g}A.csv', index=False, float_format='%.6f')
    else:
        pairs = pair_finder(bigmolecules, limit=cutoff)


    groups = {}
    for key, mol_df in bigmolecules.groupby(['Mol', 'Lattice']):
        #we have created groups, labelled by (ie has key value) (mol, lattice) - the key is a tuple
        #each group is a df containing the atoms
        groups[key] = mol_df
        # groupby doesnt actually create the groups, this func does

    needed = pd.concat([
        pd.DataFrame({'Mol': pairs.init, 'Lattice': 0}),
        pd.DataFrame({'Mol': pairs.fin,  'Lattice': pairs.fin_lattice}),
    ]).drop_duplicates()
    # needed is a df containing all the monomers the dimers will use

     # write the monomers
    os.makedirs('MoleculesPBC', exist_ok=True)
    for mol, latt in needed.itertuples(index=False):
        mol_df = groups[(mol, latt)]
        write_monomer_xyz(mol_df, mol, latt)

    # write the dimers
    for pair in pairs.itertuples(index=False):
        mol1_df = groups[(pair.init, 0)]
        mol2_df = groups[(pair.fin, pair.fin_lattice)]
        mol1 = pair.init
        mol2 = pair.fin
        fin_lattice = pair.fin_lattice
        write_dimer_xyz(mol1_df, mol2_df, mol1, mol2, fin_lattice)

    np.savetxt('init.file',       pairs['init'].values,        fmt='%d')
    np.savetxt('fin.file',        pairs['fin'].values,         fmt='%d')
    np.savetxt('finlattice.file', pairs['fin_lattice'].values, fmt='%d')
    np.savetxt('dist.txt',        pairs['dist'].values,        fmt='%.8f')



if __name__ == '__main__': 
    mol, T, cutoff = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
    save_val = float(sys.argv[4]) if len(sys.argv) > 4 else False
    prepare_snapshot(mol, T, cutoff, save_val)
                    

