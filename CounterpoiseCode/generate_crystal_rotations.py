# -*- coding: utf-8 -*-
"""
Purpose: functions to generate the crystal structure and calculate interactions between pairs
"""
import numpy as np
import pandas as pd
         
def calculate_min_distance_between_C_atoms(molecule1, molecule2):
    C_atoms_mol1 = molecule1[molecule1['Symbol'] == 'C'][['X', 'Y', 'Z']].values
    C_atoms_mol2 = molecule2[molecule2['Symbol'] == 'C'][['X', 'Y', 'Z']].values
    pairwise_distances = np.linalg.norm(C_atoms_mol1[:, None] - C_atoms_mol2, axis=-1)
    min_distance = np.min(pairwise_distances)
    return min_distance

def find_centre(mol):
    COM = np.array(mol['weight']@mol[['X','Y','Z']]/sum(mol['weight']))
    return COM

def calculate_dist_centres(molecule1, molecule2):
    CoM_mol1 = find_centre(molecule1)
    CoM_mol2 = find_centre(molecule2)
    cen2cen_dis = abs(np.linalg.norm(CoM_mol1 - CoM_mol2))
    return cen2cen_dis
   
def define_pairs(df):
    molecules = df['Mol'].unique()
    init,fin,dist,dist_CC = [],[],[],[]
    for i in range(len(molecules)):
        mol1 = df[df['Mol'] == molecules[i]].copy()
        for j in range(i+1,len(molecules)):
            mol2 = df[df['Mol'] == molecules[j]].copy()
            init.append(molecules[i])
            fin.append(molecules[j])
            dis = calculate_min_distance_between_C_atoms(mol1,mol2)
            dis_CC = calculate_dist_centres(mol1,mol2)
            dist.append(np.round(dis,6))
            dist_CC.append(np.round(dis_CC,6))
    pairs = pd.DataFrame({'init' : init, 'fin' : fin,
                          'dist' : dist})
    return pairs

def define_pairs_PBC(input_df, latticevec):
    df = input_df.copy()
    df['Lattice'] = 0 
    molecules = df['Mol'].unique()  
    count = 1
    for i in [-1,0,1]:
        for j in [-1,0,1]:
            for k in [-1,0,1]:
                if i == 0 and j == 0 and k == 0:
                    continue
                temp = df[df['Lattice']==0].copy()
                temp['Lattice'] = count
                temp['X'] += i*latticevec[0,0] + j*latticevec[1,0] + k*latticevec[2,0]
                temp['Y'] += i*latticevec[0,1] + j*latticevec[1,1] + k*latticevec[2,1]
                temp['Z'] += i*latticevec[0,2] + j*latticevec[1,2] + k*latticevec[2,2]
                df = pd.concat([df,temp])
                count += 1
    init,fin,fin_lattice,dist,dist_CC = [],[],[],[],[]
    fin_X, fin_Y, fin_Z = [],[],[]
    for i in range(len(molecules)):
        mol1 = df[(df['Mol'] == molecules[i]) & (df['Lattice'] == 0)].copy()      
        for j in range(i+1,len(molecules)):
            seps = []
            temp = df[df['Mol'] == molecules[j]]
            for m in temp['Lattice'].unique():
                m2 = temp[temp['Lattice'] == m]
                sep = calculate_min_distance_between_C_atoms(mol1,m2)
                seps.append(sep)
            which_lattice = np.argmin(seps)
            mol2 = df[(df['Mol'] == molecules[j]) \
                      & (df['Lattice'] == which_lattice)].copy()
            init.append(molecules[i])
            fin.append(molecules[j])
            fin_lattice.append(which_lattice)
            fin_X.append(sum(mol2['X'].values)/len(mol2))
            fin_Y.append(sum(mol2['Y'].values)/len(mol2))
            fin_Z.append(sum(mol2['Z'].values)/len(mol2))
            dis = calculate_min_distance_between_C_atoms(mol1,mol2)
            dis_CC = calculate_dist_centres(mol1,mol2)
            dist.append(np.round(dis,6))
            dist_CC.append(np.round(dis_CC,6))
            #print(molecules[i], molecules[j])
    pairs = pd.DataFrame({'init' : init, 'fin' : fin,
                          'fin_lattice' : fin_lattice,
                          'fin_X': fin_X, 'fin_Y': fin_Y, 'fin_Z': fin_Z,
                          'dist' : dist})
    return df, pairs

def find_unique_pairs(pairs):
    df = pairs.copy()
    df1 = df.drop_duplicates(subset=['dist']).copy()
    return df1


    
    