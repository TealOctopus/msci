import numpy as np
import pandas as pd
import sys
import os
import re
from scipy.spatial import cKDTree

from molecule_dict import MOLECULES

_ele_to_num = {"H":1, "C":12,"N":14,"F":19,"S":32,"O":16}

GRO_RES = re.compile(r"(\d+)(\D\S*)")   # "23L8BO" -> ("23", "L8BO")
#regex gobbledegook. Apparantly splits into a number and a string. Used to get the residue number from the gro file.


def gro_to_atoms(gro_file, mol):
    #Reads a gromacs .gro file and returns a pd datatframe
    #df has one row per atom of the molecule
    #columns are what moleucles the atom belongs to, the element symbol, the atomic weight, and the x,y,z coordinates
    resname = MOLECULES[mol]["resname"]
    res_ID, element, weight, x, y, z = [], [], [], [], [], []
    # Added to GROMACS residue numbers to give molecule IDs
    RES_ID_OFFSET = 100

    with open(gro_file, "r") as f:                       
         f.readline() #title
         f.readline() #number of atoms
         for line in f: 
            values = line.split() # breaks the line into a list of strings, split by whitespace
            if values: #ie if line not empty
                m = GRO_RES.match(values[0])
                #tries to match the first value in the line to the regex pattern defined
                #If it matches, m.group holds the broken string; 
                #if it doesnt match, m will be None.
            else:
                m = None
            
            if m is None or m.group(2) != resname:
                continue   # skips the box line and any other residues

            elem_symbol = values[1][0] #values[1] = eg C1, H23, so values[1][0]= e.g. C or H

            res_ID.append(int(m.group(1)) + RES_ID_OFFSET) #ie 100+ molecule num
            #I'm not sure why this is done, maybe to keep everything to 3 figs?

            element.append(elem_symbol)
            weight.append(_ele_to_num[elem_symbol])
            # last three fields are x, y, z (assumes no velocities in the .gro)
            xi, yi, zi = (10 * float(v) for v in values[-3:])   # nm -> Å
            # the x y z coords are always the final 3 values
            x.append(xi); y.append(yi); z.append(zi)

    atoms = pd.DataFrame({"Mol": res_ID, "Symbol": element, "weight": weight,
                         "X": x, "Y": y, "Z": z})
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


def truncate_and_cap(df, mol):
    """
    Read a .gro snapshot, keep only the atoms in `atoms_to_keep` for each
    molecule, and turn the cut-point carbons into capping hydrogens.

    Returns a DataFrame of all truncated molecules (indexed by original
    atom number within each molecule; use the 'Mol' column to tell
    molecules apart).
    """
    from molecule_dict import get_trim_params

    p = get_trim_params(mol)
    backbone_stop    = p["backbone_stop"]
    carbon_to_keep   = p["carbon_to_keep"]
    atoms_pairs      = p["atoms_pairs"]
    ref_hydrogen     = p["ref_hydrogen"]
    hydrogen_triples = p["hydrogen_triples"]
    atoms_to_keep = np.arange(0, backbone_stop)
    atoms_to_keep = np.append(atoms_to_keep,carbon_to_keep)

    frames = [] #holds a list of dataframes
    for m in df['Mol'].unique():
        #ie for each moledule in the box
        # each row of df is one atom, and the 'Mol' column tells us which molecule it belongs to
        #.unique means we only get one copy of each molecule number
        is_this_mol = (df['Mol'] == m)           # True/False for every atom, True if atom is in mol m
        this_mol    = df[is_this_mol]           # creates a dataframe only keeping True entries, ie only keeping those 
        #atoms that belong to the current molecule
        kept        = this_mol.iloc[atoms_to_keep]   # remove the atoms we are not keeping
        df1         = kept.copy()   #and copies it to avoid modifying the original df
        
        df1.index = atoms_to_keep
        #before, atoms were indexed by their original atom number in the gro file, 
        #now they are indexed by their position in the list of atoms we are keeping

        for j, c in enumerate(carbon_to_keep):
            #truncation step, and capping
            df1.loc[c, 'Symbol'] = "H"
            df1.loc[c, 'weight'] = 1
            df1 = fix_hydrogens(df1, atoms_pairs[j],
                                ref_hydrogen[j], hydrogen_triples[j])
        frames.append(df1) #we no longer use whoe crystal index
    return pd.concat(frames) #joins all the dataframes in the list together into one big data frame

def build_supercell(input_df, latticevec):
    cell0 = input_df.copy()
    images = [cell0] #so we dont concat each time
    cell0['Lattice'] = 0
    count = 1
    for i in [-1,0,1]:
        for j in [-1,0,1]:
            for k in [-1,0,1]:
                if i == 0 and j == 0 and k == 0:
                    continue
                shift = i*latticevec[0] + j*latticevec[1] + k*latticevec[2]
                new_cell = cell0.copy()
                new_cell['Lattice'] = count
                new_cell['X'] += shift[0]
                new_cell['Y'] += shift[1]
                new_cell['Z'] += shift[2]
                images.append(new_cell)
                count += 1
    df = pd.concat(images, ignore_index=True)
    #our full df with all cells
    return df


def pair_finder(big, limit):
    C = big[big['Symbol'] == 'C'] # keeps only the carbon atoms
    C_lattice0 = C[C['Lattice'] == 0] #only keeps C_arr entries in lattice 0

    C_lattice0_xyz = C_lattice0[['X', 'Y', 'Z']].values
    all_xyz  = C[['X', 'Y', 'Z']].values

    lattice0_tree = cKDTree(C_lattice0_xyz)
    all_tree = cKDTree(all_xyz) 
    #cKDTree(points) builds a KD-tree. It takes all points and splits them into nested boxes
    #this makes it easier to judge distances, and eliminate points too far without calcs

    close_pairs = lattice0_tree.sparse_distance_matrix(all_tree, limit, output_type='coo_matrix')
    # every (lattice0 carbon, any carbon) pair closer than limit
    #coo_matrix: ie only store the pairs which are closer than limit

    i_lattice0 = close_pairs.row     # position of each pair's mol1 carbon in C_arr_lattice0
    i_all  = close_pairs.col     # position of the partner ie mol2 carbon in C_arr
    d  = close_pairs.data    # distance between them (Å)

    atom_pairs = pd.DataFrame({'init':        C_lattice0['Mol'].values[i_lattice0],
                               'fin':         C['Mol'].values[i_all],
                               'fin_lattice': C['Lattice'].values[i_all],
                               'dist':        d})
    #this pd contains the atom pairs- we want only relevant molecules
    
    # remove carbons paired with carbons in the same molecule
    same_molecule = (atom_pairs.init == atom_pairs.fin) & (atom_pairs.fin_lattice == 0)
    atom_pairs = atom_pairs[~same_molecule] #~column == NOT(column)

    dimer_groups = atom_pairs.groupby(['init', 'fin', 'fin_lattice'], as_index=False)
    #sorts the rows into groups, st each group member (atom pair) has same init mol and fin mol
    #ie each pair group represents atom pairs in the same dimer
    #each group carries a label (init, fin, fin_lattice)

    dimer_dists = dimer_groups['dist']
    # essentially makes a dimer_dists df structured like: group (101, 102, 0):   [4.21, 3.58, 5.02]
    #labels remain

    pairs = dimer_dists.min() #only want to keep min dist
    #the group label gets split back into columns, since implicitly as_index =False

    # each dimer appears as (A,B) and (B,A); we only want to keep one
    pairs = pairs[pairs.init < pairs.fin].reset_index(drop=True)
    return pairs