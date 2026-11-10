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

def write_xyz_file(df, mol, Latt):
    df1 = df[['Symbol','X','Y','Z']]
    fname = "MoleculesPBC/Mol{:.0f}_Lattice{:.0f}.xyz".format(mol, Latt)
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
    init,fin,fin_lattice,dist = [],[],[],[]
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
            dist.append(np.round(dis,6))
    pairs = pd.DataFrame({'init' : init, 'fin' : fin,
                          'fin_lattice' : fin_lattice,
                          'fin_X': fin_X, 'fin_Y': fin_Y, 'fin_Z': fin_Z,
                          'dist' : dist})
    return df, pairs
                          
def get_box_vectors(gro_file):
    with open(gro_file, "r") as f:
        lines = f.readlines()
        if lines:  # Ensure the file is not empty
            last_line = lines[-1]           
            values = last_line.split()
            box_vectors = np.zeros((3,3))
            box_vectors[0,0] = 10*float(values[0])
            box_vectors[0,1] = 10*float(values[3])
            box_vectors[0,2] = 10*float(values[4])
            box_vectors[1,0] = 10*float(values[5])
            box_vectors[1,1] = 10*float(values[1])
            box_vectors[1,2] = 10*float(values[6])
            box_vectors[2,0] = 10*float(values[7])
            box_vectors[2,1] = 10*float(values[8])
            box_vectors[2,2] = 10*float(values[2])
    return box_vectors
          
if __name__ == '__main__':
    t = int(sys.argv[1])
    datafile = sys.argv[2]
    fname = "Y6_300K_{}ps.gro".format(t)
    file = os.path.join(datafile,fname)
    b_vectors = get_box_vectors(file)
    
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
    count = 0
    for m in df['Mol'].unique():
        df1 = df[df['Mol']==m].iloc[atoms_to_keep].copy()
        df1.index = atoms_to_keep
        for j,c in enumerate(carbon_to_keep):
            df1.loc[c,'Symbol'] = "H"
            df1.loc[c,'weight'] = 1
            df1 = fix_hydrogens(df1,atoms_pairs[j],
                                    ref_hydrogen[j],hydrogen_triples[j])
        if count == 0:
            molecules = df1
        else:
            molecules = pd.concat([molecules,df1])
        count += 1
    bigmolecules, pairs = define_pairs_PBC(molecules, b_vectors)
            
    init = []
    fin = []
    seps = []
    fin_lattices = []
    for i, pair in pairs.iterrows():
        if pair.dist <= MOO_cutoff: 
            if pair.fin_lattice != 0:
                fin_lattice = pair.fin_lattice          
                molecule1 = bigmolecules[(bigmolecules['Mol'] == pair.init) & (bigmolecules['Lattice'] == 0)]
                write_xyz_file(molecule1,pair.init,0)
                molecule2 = bigmolecules[(bigmolecules['Mol'] == pair.fin) & (bigmolecules['Lattice'] == fin_lattice)]
                write_xyz_file(molecule2,pair.fin,fin_lattice)       
                fin_lattices.append(fin_lattice)
                seps.append(pair.dist)
                init.append(pair.init)
                fin.append(pair.fin)
    
    np.savetxt('init.file', np.asarray(init), fmt='%.0f')
    np.savetxt('fin.file', np.asarray(fin), fmt='%.0f')
    np.savetxt('dist.txt', np.asarray(seps), fmt='%.8f')
    np.savetxt('finlattice.file', np.asarray(fin_lattices), fmt='%.0f')







