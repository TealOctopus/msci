import numpy as np
import pandas as pd
import sys
import os

ele_to_num = {"H":1, "C":12,"N":14,"F":19,"S":32,"O":16}

def gro_to_atoms(gro_file):
    weight = []
    element = []
    res_ID = []
    x = []
    y = []
    z = []
    with open(gro_file, "r") as f:                       
         f.readline()
         f.readline()
         for line in f: 
             values = line.split()
             try:
                 res_ID.append(int(values[0].split("Y")[0])+100)
             except:
                continue
             weight.append(ele_to_num[values[1][0]])
             element.append(values[1][0])
             if len(values) == 6:
                 x.append(10*float(values[3])) 
                 y.append(10*float(values[4]))
                 z.append(10*float(values[5]))
             elif len(values) == 5:
                 x.append(10*float(values[2])) 
                 y.append(10*float(values[3]))
                 z.append(10*float(values[4]))  
    atoms = pd.DataFrame({'Mol':res_ID, 'Symbol':element, 'weight':weight,
                          'X':x, 'Y':y, 'Z':z})
    return atoms

def arb_rotation(axis, angle, vector):
    #Rodrigues’ rotation formula
    #Rotate vector by angle theta around axis (u,v,w)
    angle = angle*2*np.pi/360
    c = np.cos(angle)
    s = np.sin(angle)
    u = axis[0]
    v = axis[1]
    w = axis[2]
    rotation_mat = np.array([[u**2+(1-u**2)*c, u*v*(1-c)-w*s, u*w*(1-c)+v*s],
                             [u*v*(1-c)+w*s, v**2+(1-v**2)*c, v*w*(1-c)-u*s],
                             [u*w*(1-c)-v*s, v*w*(1-c)+u*s, w**2+(1-w**2)*c]])
    return rotation_mat@vector

def fix_hydrogens(df,norm_atoms,
                  ref_hydrogen,hydrogen_idx, 
                  theta=109.5,rCH=1.09): 
    theta = theta*(2*np.pi/360)
    loc_bb2 = np.asarray(df.loc[norm_atoms[1]][['X','Y','Z']], dtype = float)
    #Shift carbon atom to the origin
    df[['X','Y','Z']] -= loc_bb2    
    loc_bb1 = np.asarray(df.loc[norm_atoms[0]][['X','Y','Z']], dtype = float)
    rhat = (loc_bb1)/np.sqrt(sum([i**2 for i in loc_bb1]))
    loc_O = rCH*np.cos(theta)*rhat
    loc_H = np.asarray(df.loc[ref_hydrogen][['X','Y','Z']], dtype = float)
    coeff_mat = np.array([[2,0,0,rhat[0]],
                          [0,2,0,rhat[1]],
                          [0,0,2,rhat[2]],
                          [rhat[0],rhat[1],rhat[2],0]])
    constant_vec = np.array([2*loc_H[0],
                             2*loc_H[1],
                             2*loc_H[2],
                             np.dot(loc_O,rhat)])
    sol = np.linalg.solve(coeff_mat, constant_vec)
    X = sol[0:3]
    vecOX = X - loc_O
    vecOY = arb_rotation(-rhat, 120, vecOX)
    vecOZ = arb_rotation(-rhat, -120, vecOX)
    Y = loc_O + vecOY
    Z = loc_O + vecOZ
    df.loc[ref_hydrogen, ['X','Y','Z']] = X
    df.loc[hydrogen_idx[1], ['X','Y','Z']] = Y
    df.loc[hydrogen_idx[2], ['X','Y','Z']] = Z
    #shift coordinates back 
    df[['X','Y','Z']] += loc_bb2     
    return df

def write_xyz_file(df, mol, t):
    df1 = df[['Symbol','X','Y','Z']]
    fname = "Molecules/Mol{}.xyz".format(mol, t)
    with open(fname, "w") as f:
        for _, row in df1.iterrows():
            symbol, x, y, z = row
            f.write(f"{symbol}\t{x:.7f}\t{y:.7f}\t{z:.7f}\n")
        f.write('\n')
        
def calculate_min_distance_between_C_atoms(molecule1, molecule2):
    C_atoms_mol1 = molecule1[molecule1['Symbol'] == 'C'][['X', 'Y', 'Z']].values
    C_atoms_mol2 = molecule2[molecule2['Symbol'] == 'C'][['X', 'Y', 'Z']].values
    pairwise_distances = np.linalg.norm(C_atoms_mol1[:, None] - C_atoms_mol2, axis=-1)
    min_distance = np.min(pairwise_distances)
    return min_distance
            
if __name__ == '__main__':
    t = int(sys.argv[1])
    datafile = sys.argv[2]
    fname = "Y6_300K_{}ps.gro".format(t)
    file = os.path.join(datafile,fname)
    
    MOO_cutoff = 5
    
    backbone_stop = 81  
    #These will be different for L8BO and Y6
    #NB, have subtracted 1 versus GROMACS numbering as python indexes from 0
    carbon_to_keep = [156,125,96,81]
    atoms_pairs = [[21,34], [24,33], [14,32], [11,31]]
    ref_hydrogen = [75,73,71,69]
    hydrogen_triples = [[75,76,156],[73,74,125],[71,72,96],[69,70,81]]
    
    atoms_to_keep = np.arange(0, backbone_stop)
    atoms_to_keep = np.append(atoms_to_keep,carbon_to_keep)

    df = gro_to_atoms(file)
    for m in df['Mol'].unique():
        df1 = df[df['Mol']==m].iloc[atoms_to_keep].copy()
        df1.index = atoms_to_keep
        for j,c in enumerate(carbon_to_keep):
            df1.loc[c,'Symbol'] = "H"
            df1.loc[c,'weight'] = 1
            df1 = fix_hydrogens(df1,atoms_pairs[j],
                                    ref_hydrogen[j],hydrogen_triples[j])
        write_xyz_file(df1,m,t)
        
    init = []
    fin = []
    seps = []
    for m1 in df['Mol'].unique():
        print(m1)
        for m2 in df['Mol'].unique():
            if m1 < m2:
                molecule1 = df[df['Mol']==m1].iloc[atoms_to_keep].copy()
                molecule2 = df[df['Mol']==m2].iloc[atoms_to_keep].copy()
                min_sep = calculate_min_distance_between_C_atoms(molecule1, 
                                                                 molecule2)
                if min_sep <= MOO_cutoff: 
                    seps.append(min_sep)
                    init.append(m1)
                    fin.append(m2)
    
    np.savetxt('init.file', np.asarray(init), fmt='%.0f')
    np.savetxt('fin.file', np.asarray(fin), fmt='%.0f')
    np.savetxt('dist.txt', np.asarray(seps), fmt='%.8f')







