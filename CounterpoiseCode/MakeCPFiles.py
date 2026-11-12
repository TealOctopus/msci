# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 09:48:48 2026

@author: Lucy
"""
import numpy as np
import pandas as pd
import os
import generate_CP_GJF as cp_gjf
import generate_crystal_rotations as gcrys

#%%

def save_molecules_to_xyz(df, save_path, PBC = 0, pairs = 0):
    for mol in df['Mol'].unique():
        if PBC == 0:
            df1 = df[df['Mol']==mol][['Symbol','X','Y','Z']]
            fname = "{}Mol{}.xyz".format(save_path, mol)
        elif PBC == 1:
            df1 = df[(df['Mol']==mol) & (df['Lattice']==0)][['Symbol','X','Y','Z']]
            fname = "{}Mol{}_Lattice0.xyz".format(save_path, mol)
        with open(fname, "w") as f:
            f.write("Symbol\tX\tY\tZ\n")
            for _, row in df1.iterrows():
                symbol, x, y, z = row
                f.write(f"{symbol}\t{x:.7f}\t{y:.7f}\t{z:.7f}\n")
    if PBC == 1:
        problem_children = pairs[pairs['fin_lattice'] != 0]
        for i in range(len(problem_children)):
            fin = problem_children['fin_lattice'].iloc[i]
            mol = problem_children['fin'].iloc[i]
            df1 = df[(df['Mol']==mol) & (df['Lattice']==fin)][['Symbol','X','Y','Z']]
            fname = "{}Mol{}_Lattice{}.xyz".format(save_path, mol, fin)
            with open(fname, "w") as f:
                f.write("Symbol\tX\tY\tZ\n")
                for _, row in df1.iterrows():
                    symbol, x, y, z = row
                    f.write(f"{symbol}\t{x:.7f}\t{y:.7f}\t{z:.7f}\n")
                


#%%
rng = np.random.default_rng()
t_ar = np.arange(20000,50001,6000)

savefile_xyz = 'GenCPFiles/xyz/'
savefile_CP = 'GenCPFiles/JobFiles/'
os.makedirs(savefile_xyz, exist_ok=True)
os.makedirs(savefile_CP, exist_ok=True)

makeCPfiles = 0
PBC = 1
save_csv = 0 

big_batch = 0
Nsamples = 33
max_dist = 6

backbone_stop = 81

which = 'L8BO'

#These will be different for L8BO and Y6
#NB, have subtracted 1 versus GROMACS numbering as python indexes from 0
if which == 'Y6':
    carbon_to_keep = [156,125,96,81]
    hydrogen_triples = [[75,76,156],[73,74,125],[71,72,96],[69,70,81]]
elif which == 'L8BO':
    carbon_to_keep = [83,112,82,81]
    hydrogen_triples = [[75,76,83],[73,74,112],[71,72,82],[69,70,81]]
    
atoms_pairs = [[21,34], [24,33], [14,32], [11,31]]
ref_hydrogen = [75,73,71,69]   
atoms_to_keep = np.arange(0, backbone_stop)
atoms_to_keep = np.append(atoms_to_keep,carbon_to_keep)

count = 0
for t in [t_ar[0]]:
    print(count)
    file = "FILLTHISIN".format(which,t)
    b_vectors = get_box_vectors(file)
    df = gro_to_atoms(file,which)
    for m,mol in enumerate(df['Mol'].unique()):
        df1 = df[df['Mol']==mol].iloc[atoms_to_keep].copy()
        df1.index = atoms_to_keep
        for j,c in enumerate(carbon_to_keep):
            df1.loc[c,'Symbol'] = "H"
            df1.loc[c,'weight'] = 1
            df1 = fix_hydrogens(df1,atoms_pairs[j],
                                ref_hydrogen[j],hydrogen_triples[j])       
        if m == 0:
            molecules = df1
        else:
            molecules = pd.concat([molecules,df1])
            
    if PBC:
        bigmolecules, pairs = gcrys.define_pairs_PBC(molecules, b_vectors)
        if save_csv:
            pairs.to_csv('Pairs/{}/pairs_total_PBC_{}ps.csv'.format(which,t), 
                         index=False)
            bigmolecules.to_csv('Atoms/{}/atoms_PBC_{}ps.csv'.format(which,t), 
                                index=False)
    else:
        pairs = gcrys.define_pairs(molecules)
        if save_csv:
            pairs.to_csv('Pairs/{}/pairs_total_{}ps.csv'.format(which,t), 
                         index=False)
            molecules.to_csv('Atoms/{}/atoms_{}ps.csv'.format(which,t), 
                             index=False)
    count += 1 
    
    if makeCPfiles == 1:
        neighbours = pairs[pairs['dist'] < max_dist]
        if PBC == 0:
            save_molecules_to_xyz(molecules, savefile_xyz)
        elif PBC == 1:
            save_molecules_to_xyz(bigmolecules, savefile_xyz, PBC=1, 
                                  pairs = neighbours)
        if big_batch == 1:
            sample_idx = rng.integers(low = 0, high = len(neighbours), 
                                      size = Nsamples)
            sample_idx = np.unique(sample_idx)
            while len(sample_idx) < Nsamples:
                new_values = rng.integers(low = 0, 
                                          high = len(neighbours), 
                                          size = Nsamples-len(sample_idx))
                sample_idx = np.concat([sample_idx, new_values])
            neighbours = neighbours.iloc[sample_idx]   
        cp_gjf.GaussInputFile_list(neighbours, savefile_CP,
                                   savefile_xyz,
                                   nprocshared=16,
                                   mem = '10GB', 
                                   method = 'LC-wHPBE/6-31G(d,p) IOp(3/107=0140000000,3/108=0140000000,3/119=0,3/120=0, 3/130=0, 3/131=0)',
                                   append_t = True,
                                   PBC = PBC)





