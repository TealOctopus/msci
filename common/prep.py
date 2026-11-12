import numpy as np
import pandas as pd
import sys
import os


ele_to_num = {"H":1, "C":12,"N":14,"F":19,"S":32,"O":16}

def gro_to_atoms(gro_file, mol = 'Y6'):
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
                if mol == 'Y6':
                     res_ID.append(int(values[0].split("Y")[0])+100)
                elif mol == 'L8BO':
                    res_ID.append(int(values[0].split("L")[0])+100)
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